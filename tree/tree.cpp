#include "tree.h"

#include <cassert>
#include <utils/matrix.h>

using namespace std;

Tree::Tree(Params&& n_params)
{
    params = std::move(n_params);
}

void Tree::setupFunctionPointers()
{
    switch (aln.num_states) {
        case 2:
            updatePartialLhPointer = &Tree::updatePartialLhTemplate<4>;
            calculateSamplePlacementCostPointer = &Tree::calculateSamplePlacementCostTemplate<2>;
            calculateSubTreePlacementCostPointer = &Tree::calculateSubTreePlacementCostTemplate<2>;
            break;
        case 4:
            updatePartialLhPointer = &Tree::updatePartialLhTemplate<4>;
            calculateSamplePlacementCostPointer = &Tree::calculateSamplePlacementCostTemplate<4>;
            calculateSubTreePlacementCostPointer = &Tree::calculateSubTreePlacementCostTemplate<4>;
            break;
        case 20:
            updatePartialLhPointer = &Tree::updatePartialLhTemplate<4>;
            calculateSamplePlacementCostPointer = &Tree::calculateSamplePlacementCostTemplate<20>;
            calculateSubTreePlacementCostPointer = &Tree::calculateSubTreePlacementCostTemplate<20>;
            break;
            
        default:
            outError("Sorry! currently we only support DNA data!");
            break;
    }
}

void Tree::setupBlengthThresh()
{
    default_blength = 1.0 / aln.ref_seq.size();
    min_blength = params->min_blength_factor * default_blength;
    max_blength = params->max_blength_factor * default_blength;
    min_blength_mid = params->min_blength_mid_factor * default_blength;
    min_blength_sensitivity = min_blength * 1e-5;
    half_min_blength_mid = min_blength_mid * 0.5;
    half_max_blength = max_blength * 0.5;
    double_min_blength = min_blength + min_blength;
}

void Tree::setup()
{
    setupFunctionPointers();
    setupBlengthThresh();
}

const string Tree::exportTreeString(const bool binary, const NumSeqsType node_vec_index) const
{
    const PhyloNode& node = nodes[node_vec_index];
    string output = "(";
    
    // if it's a leaf
    if (!node.isInternal())
        return node.exportString(binary, aln);
    // if it's an internal node
    else
    {
        /*bool add_comma = false;
         for (Index neighbor_index:node.getNeighborIndexes(TOP))
         {
         if (!add_comma)
         add_comma = true;
         else
         output += ",";
         output += exportTreeString(binary, neighbor_index.getVectorIndex());
         }*/
        output += exportTreeString(binary, node.getNeighborIndex(RIGHT).getVectorIndex());
        output += ",";
        output += exportTreeString(binary, node.getNeighborIndex(LEFT).getVectorIndex());
    }

    string length = node.getUpperLength() < 0 ? "0" : convertDoubleToString(node.getUpperLength());
    output += "):" + length;
    
    return output;
}

void Tree::updateMidBranchLh(const Index node_index, PhyloNode& node, const std::unique_ptr<SeqRegions>& parent_upper_regions, stack<Index> &node_stack, bool &update_blength)
{
    // update vector of regions at mid-branch point
    std::unique_ptr<SeqRegions> mid_branch_regions = nullptr;
    computeMidBranchRegions(node, mid_branch_regions, *parent_upper_regions);
    
    if (!mid_branch_regions)
        handleNullNewRegions(node_index, node, (node.getUpperLength() <= 1e-100), node_stack, update_blength, "inside updatePartialLh(), from parent: should not have happened since node->length > 0");
    // update likelihood at the mid-branch point
    else
        node.setMidBranchLh(std::move(mid_branch_regions));
}

std::unique_ptr<SeqRegions> Tree::computeUpperLeftRightRegions(const Index node_index, PhyloNode& node, const MiniIndex next_node_mini, const std::unique_ptr<SeqRegions>& parent_upper_regions, std::stack<Index> &node_stack, bool &update_blength)
{
    std::unique_ptr<SeqRegions> upper_left_right_regions = nullptr;
    const std::unique_ptr<SeqRegions>& lower_regions = getPartialLhAtNode(node.getNeighborIndex(next_node_mini));  // next_node->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob);
    // parent_upper_regions->mergeUpperLower(upper_left_right_regions, node.getUpperLength(), *lower_regions, next_node->length, aln, model, params->threshold_prob);
    parent_upper_regions->mergeUpperLower(upper_left_right_regions, node.getUpperLength(), *lower_regions, node.getCorrespondingLength(next_node_mini, nodes), aln, model, params->threshold_prob);
    
    // handle cases when new regions is null/empty
    if (!upper_left_right_regions || upper_left_right_regions->size() == 0)
        handleNullNewRegions(node_index, node, (node.getUpperLength() <= 0 && node.getCorrespondingLength(next_node_mini, nodes) <= 0), node_stack, update_blength, "Strange: None vector from non-zero distances in updatePartialLh() from parent direction.");
    
    return upper_left_right_regions;
}

void Tree::updateNewPartialIfDifferent(PhyloNode& node, const MiniIndex next_node_mini, std::unique_ptr<SeqRegions>& upper_left_right_regions, std::stack<Index> &node_stack, const PositionType seq_length)
{
    if (node.getPartialLh(next_node_mini)->areDiffFrom(*upper_left_right_regions, seq_length, aln.num_states, &params.value())) //(next_node->getPartialLhAtNode(aln, model, params->threshold_prob)->areDiffFrom(*upper_left_right_regions, seq_length, aln.num_states, &params.value()))
    {
        /*replacePartialLH(next_node->partial_lh, upper_left_right_regions);
        node_stack.push(next_node->neighbor);*/
        node.setPartialLh(next_node_mini, std::move(upper_left_right_regions));
        node_stack.push(node.getNeighborIndex(next_node_mini));
    }
}

void Tree::updatePartialLhFromParent(const Index index, PhyloNode& node, stack<Index> &node_stack, const std::unique_ptr<SeqRegions>& parent_upper_regions, const PositionType seq_length)
{
    bool update_blength = false;
    
    // if necessary, update the total probabilities at the mid node.
    if (node.getUpperLength() > 0)
    {
        // update vector of regions at mid-branch point
        updateMidBranchLh(index, node, parent_upper_regions, node_stack, update_blength);
        
        // if necessary, update the total probability vector.
        if (!update_blength)
        {
            // node->computeTotalLhAtNode(aln, model, params->threshold_prob, node == root);
            node.updateTotalLhAtNode(nodes[node.getNeighborIndex(TOP).getVectorIndex()], aln, model, params->threshold_prob, index.getVectorIndex() == root_vector_index);
            
            if (!node.getTotalLh() || node.getTotalLh()->size() == 0)
                outError("inside updatePartialLh(), from parent 2: should not have happened since node->length > 0");
        }
    }
    
    // at valid internal node, update upLeft and upRight, and if necessary add children to node_stack.
    if (node.isInternal() && !update_blength)//(node->next && !update_blength)
    {
        /*Node* next_node_1 = node->next;
        Node* next_node_2 = next_node_1->next;*/
        
        // compute new upper left/right for next_node_1
        std::unique_ptr<SeqRegions> upper_left_right_regions_1 = computeUpperLeftRightRegions(index, node, RIGHT, parent_upper_regions, node_stack, update_blength);//computeUpperLeftRightRegions(next_node_1, node, parent_upper_regions, node_stack, update_blength);
        std::unique_ptr<SeqRegions> upper_left_right_regions_2 = nullptr;
        
        // compute new upper left/right for next_node_1
        if (!update_blength)
            upper_left_right_regions_2 = computeUpperLeftRightRegions(index, node, LEFT, parent_upper_regions, node_stack, update_blength);
        
        if (!update_blength)
        {
            // update new partiallh for next_node_1
            updateNewPartialIfDifferent(node, RIGHT, upper_left_right_regions_2, node_stack, seq_length);
            
            // update new partiallh for next_node_2
            updateNewPartialIfDifferent(node, LEFT, upper_left_right_regions_1, node_stack, seq_length);
        }
    }
}

void Tree::updatePartialLhFromChildren(const Index index, PhyloNode& node, std::stack<Index> &node_stack, const std::unique_ptr<SeqRegions>& parent_upper_regions, const bool is_non_root, const PositionType seq_length)
{
    bool update_blength = false;
    /*Node* top_node = NULL;
    Node* other_next_node = NULL;
    Node* next_node = NULL;
    FOR_NEXT(node, next_node)
    {
        if (next_node->is_top)
            top_node = next_node;
        else
            other_next_node = next_node;
    }
    
    ASSERT(top_node && other_next_node);*/
    
    const NumSeqsType node_vec_index = index.getVectorIndex();
    const Index top_node_index = Index(node_vec_index, TOP);
    const MiniIndex node_mini = index.getMiniIndex();
    const MiniIndex other_next_node_mini = index.getFlipMiniIndex();
    
    const RealNumType this_node_distance = node.getCorrespondingLength(node_mini, nodes); //node->length;
    const RealNumType other_next_node_distance = node.getCorrespondingLength(other_next_node_mini, nodes); //other_next_node->length;
    const Index neighbor_index = node.getNeighborIndex(node_mini);
    PhyloNode& neighbor = nodes[neighbor_index.getVectorIndex()];
    const std::unique_ptr<SeqRegions>& this_node_lower_regions = neighbor.getPartialLh(neighbor_index.getMiniIndex()); // node->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob);
    
    // update lower likelihoods
    std::unique_ptr<SeqRegions> merged_two_lower_regions = nullptr;
    std::unique_ptr<SeqRegions> old_lower_regions = nullptr;
    // other_next_node->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob)->mergeTwoLowers(merged_two_lower_regions, other_next_node_distance, *this_node_lower_regions, this_node_distance, aln, model, params->threshold_prob);
    getPartialLhAtNode(node.getNeighborIndex(other_next_node_mini))->mergeTwoLowers(merged_two_lower_regions, other_next_node_distance, *this_node_lower_regions, this_node_distance, aln, model, params->threshold_prob);
    
    if (!merged_two_lower_regions || merged_two_lower_regions->size() == 0)
    {
        //handleNullNewRegions(node->neighbor, (this_node_distance <= 0 && other_next_node_distance <= 0), node_stack, update_blength, "Strange: None vector from non-zero distances in updatePartialLh() from child direction.");
        handleNullNewRegions(neighbor_index, neighbor, (this_node_distance <= 0 && other_next_node_distance <= 0), node_stack, update_blength, "Strange: None vector from non-zero distances in updatePartialLh() from child direction.");
    }
    else
    {
        /*replacePartialLH(old_lower_regions, top_node->partial_lh);
        top_node->partial_lh = merged_two_lower_regions;
        merged_two_lower_regions = NULL;*/
        old_lower_regions = std::move(node.getPartialLh(TOP));
        node.setPartialLh(TOP, std::move(merged_two_lower_regions));
    }

    // update total likelihood
    if (!update_blength)
    {
        if (node.getUpperLength() > 0 || node_vec_index == root_vector_index) //(top_node->length > 0 || top_node == root)
        {
            // SeqRegions* new_total_lh_regions = top_node->computeTotalLhAtNode(aln, model, params->threshold_prob, top_node == root, false);
            std::unique_ptr<SeqRegions> new_total_lh_regions = node.computeTotalLhAtNode(nodes[node.getNeighborIndex(TOP).getVectorIndex()], aln, model, params->threshold_prob, node_vec_index == root_vector_index);
            
            if (!new_total_lh_regions)
            {
                // handleNullNewRegions(top_node, (node.getUpperLength() <= 0), node_stack, update_blength, "Strange: None vector from non-zero distances in updatePartialLh() from child direction while doing overall likelihood.");
                handleNullNewRegions(top_node_index, node, (node.getUpperLength() <= 0), node_stack, update_blength, "Strange: None vector from non-zero distances in updatePartialLh() from child direction while doing overall likelihood.");
            }
            else
            {
                // replacePartialLH(top_node->total_lh, new_total_lh_regions);
                node.setTotalLh(std::move(new_total_lh_regions));
            }
        }
    }
    
    // update total mid-branches likelihood
    if (!update_blength && node.getUpperLength() > 0 && is_non_root) //(!update_blength && top_node->length > 0 && is_non_root)
        //updateMidBranchLh(top_node, parent_upper_regions, node_stack, update_blength);
        updateMidBranchLh(top_node_index, node, parent_upper_regions, node_stack, update_blength);
    
    if (!update_blength)
    {
        // update likelihoods at parent node
        if (node.getPartialLh(TOP)->areDiffFrom(*old_lower_regions, seq_length, aln.num_states, &params.value()) && node_vec_index != root_vector_index) //(top_node->getPartialLhAtNode(aln, model, params->threshold_prob)->areDiffFrom(*old_lower_regions, seq_length, aln.num_states, &params.value()) && root != top_node)
            //node_stack.push(top_node->neighbor);
            node_stack.push(node.getNeighborIndex(TOP));

        // update likelihoods at sibling node
        std::unique_ptr<SeqRegions> new_upper_regions = nullptr;
        if (is_non_root)
            parent_upper_regions->mergeUpperLower(new_upper_regions, node.getUpperLength(), *this_node_lower_regions, this_node_distance, aln, model, params->threshold_prob);
        else
            //new_upper_regions = node->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob)->computeTotalLhAtRoot(aln.num_states, model, this_node_distance);
            new_upper_regions = getPartialLhAtNode(neighbor_index)->computeTotalLhAtRoot(aln.num_states, model, this_node_distance);
        
        if (!new_upper_regions || new_upper_regions->size() == 0)
        {
            //handleNullNewRegions(top_node, (top_node->length <= 0 && this_node_distance <= 0), node_stack, update_blength, "Strange: None vector from non-zero distances in updatePartialLh() from child direction, new_upper_regions.");
            handleNullNewRegions(top_node_index, node, (node.getUpperLength() <= 0 && this_node_distance <= 0), node_stack, update_blength, "Strange: None vector from non-zero distances in updatePartialLh() from child direction, new_upper_regions.");
        }
        // update partiallh of other_next_node if the new one is different from the current one
        else
            // updateNewPartialIfDifferent(other_next_node, new_upper_regions, node_stack, seq_length);
            updateNewPartialIfDifferent(node, other_next_node_mini, new_upper_regions, node_stack, seq_length);
            
        // delete new_upper_regions
       // if (new_upper_regions) delete new_upper_regions;
    }
            
    // delete old_lower_regions
    // if (old_lower_regions) delete old_lower_regions;
}

void Tree::updatePartialLh(stack<Index> &node_stack)
{
    (this->*updatePartialLhPointer)(node_stack);
}

template <const StateType num_states>
void Tree::updatePartialLhTemplate(stack<Index> &node_stack)
{
    const PositionType seq_length = aln.ref_seq.size();
    
    while (!node_stack.empty())
    {
        Index node_index = node_stack.top();
        node_stack.pop();
        PhyloNode& node = nodes[node_index.getVectorIndex()];
        
        // NHANLT: debug
        //if (node->next && (node->neighbor->seq_name == "25" || (node->next->neighbor && node->next->neighbor->seq_name == "25") || //(node->next->next->neighbor && node->next->next->neighbor->seq_name == "25")))
        //if (node->seq_name == "25")
         //   cout << "dsdas";
        
        node.setOutdated(true);
        
        std::unique_ptr<SeqRegions> null_seqregions_ptr = nullptr;
        bool is_non_root = root_vector_index != node_index.getVectorIndex();
        /*if (is_non_root)
        {
            //parent_upper_regions = node->getTopNode()->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob);
            SeqRegions parent_upper_regions_clone = SeqRegions(getPartialLhAtNode(node.getNeighborIndex(TOP)));
            parent_upper_regions = std::make_unique<SeqRegions>(std::move(parent_upper_regions_clone));
        }*/
        const std::unique_ptr<SeqRegions>& parent_upper_regions = is_non_root ? getPartialLhAtNode(node.getNeighborIndex(TOP)) :  null_seqregions_ptr;
            
        // change in likelihoods is coming from parent node
        if (node_index.getMiniIndex() == TOP)
        {
            ASSERT(is_non_root);
            updatePartialLhFromParent(node_index, node, node_stack, parent_upper_regions, seq_length);
        }
        // otherwise, change in likelihoods is coming from a child.
        else
        {
            updatePartialLhFromChildren(node_index, node, node_stack, parent_upper_regions, is_non_root, seq_length);
        }
    }
}

void Tree::examineSamplePlacementMidBranch(Index& selected_node_index, const std::unique_ptr<SeqRegions>& mid_branch_lh, RealNumType &best_lh_diff, bool& is_mid_branch, RealNumType& lh_diff_mid_branch, TraversingNode& current_extended_node, const std::unique_ptr<SeqRegions>& sample_regions)
{
    // compute the placement cost
    lh_diff_mid_branch = calculateSamplePlacementCost(nodes[current_extended_node.getIndex().getVectorIndex()].getMidBranchLh(), sample_regions, default_blength);
    
    // record the best_lh_diff if lh_diff_mid_branch is greater than the best_lh_diff ever
    if (lh_diff_mid_branch > best_lh_diff)
    {
        best_lh_diff = lh_diff_mid_branch;
        selected_node_index = current_extended_node.getIndex();
        current_extended_node.setFailureCount(0);
        is_mid_branch = true;
    }
}

void Tree::examineSamplePlacementAtNode(Index& selected_node_index, const std::unique_ptr<SeqRegions>& total_lh, RealNumType &best_lh_diff, bool& is_mid_branch, RealNumType& lh_diff_at_node, RealNumType& lh_diff_mid_branch, RealNumType &best_up_lh_diff, RealNumType &best_down_lh_diff, Index& best_child_index, TraversingNode& current_extended_node, const std::unique_ptr<SeqRegions>& sample_regions)
{
    // compute the placement cost
    lh_diff_at_node = calculateSamplePlacementCost(total_lh, sample_regions, default_blength);
    
    // record the best_lh_diff if lh_diff_at_node is greater than the best_lh_diff ever
    if (lh_diff_at_node > best_lh_diff)
    {
        best_lh_diff = lh_diff_at_node;
        selected_node_index = current_extended_node.getIndex();
        current_extended_node.setFailureCount(0);
        is_mid_branch = false;
        best_up_lh_diff = lh_diff_mid_branch;
    }
    else if (lh_diff_mid_branch >= (best_lh_diff - params->threshold_prob))
    {
        best_up_lh_diff = current_extended_node.getLhDiff();
        best_down_lh_diff = lh_diff_at_node;
        best_child_index = current_extended_node.getIndex();
    }
    // placement at current node is considered failed if placement likelihood is not improved by a certain margin compared to best placement so far for the nodes above it.
    else if (lh_diff_at_node < (current_extended_node.getLhDiff() - params->thresh_log_lh_failure))
        current_extended_node.increaseFailureCount();
}

void Tree::finetuneSamplePlacementAtNode(const PhyloNode& selected_node, RealNumType &best_down_lh_diff, Index& best_child_index, const std::unique_ptr<SeqRegions>& sample_regions)
{
    // current node might be part of a polytomy (represented by 0 branch lengths) so we want to explore all the children of the current node to find out if the best placement is actually in any of the branches below the current node.
    // Node* neighbor_node;
    stack<Index> node_stack;
    /*for (Index neighbor_index:nodes[selected_node_index.getVectorIndex()].getNeighborIndexes(selected_node_index.getMiniIndex()))
        node_stack.push(neighbor_index);*/
    // ASSERT(selected_node_index.getMiniIndex() == TOP);
    if (selected_node.isInternal())
    {
        node_stack.push(selected_node.getNeighborIndex(RIGHT));
        node_stack.push(selected_node.getNeighborIndex(LEFT));
    }
    
    
    while (!node_stack.empty())
    {
        Index node_index = node_stack.top();
        // MiniIndex node_mini_index = node_index.getMiniIndex();
        node_stack.pop();
        ASSERT(node_index.getMiniIndex() == TOP);
        PhyloNode& node = nodes[node_index.getVectorIndex()];
        // const RealNumType current_blength = node.getCorrespondingLength(node_mini_index, nodes);
        const RealNumType current_blength = node.getUpperLength();

        if (current_blength <= 0)
        {
            /*for (Index neighbor_index:node.getNeighborIndexes(node_index.getMiniIndex()))
                node_stack.push(neighbor_index);*/
            if (node.isInternal())
            {
                node_stack.push(node.getNeighborIndex(RIGHT));
                node_stack.push(node.getNeighborIndex(LEFT));
            }
        }
        else
        {
            // now try to place on the current branch below the best node, at an height above the mid-branch.
            RealNumType new_blength = current_blength * 0.5;
            RealNumType new_best_lh_mid_branch = MIN_NEGATIVE;
            // node->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob);
            // const std::unique_ptr<SeqRegions>& upper_lr_regions = getPartialLhAtNode(node.getNeighborIndex(node_mini_index));
            const std::unique_ptr<SeqRegions>& upper_lr_regions = getPartialLhAtNode(node.getNeighborIndex(TOP));
            //SeqRegions* lower_regions = node->getPartialLhAtNode(aln, model, params->threshold_prob);
            // const std::unique_ptr<SeqRegions>& lower_regions = node.getPartialLh(node_mini_index);
            const std::unique_ptr<SeqRegions>& lower_regions = node.getPartialLh(TOP);
            RealNumType new_lh_mid_branch = calculateSamplePlacementCost(node.getMidBranchLh(), sample_regions, default_blength);
            std::unique_ptr<SeqRegions> mid_branch_regions = nullptr;

            // try to place new sample along the upper half of the current branch
            while (true)
            {
                // record new_best_lh_mid_branch
                if (new_lh_mid_branch > new_best_lh_mid_branch)
                    new_best_lh_mid_branch = new_lh_mid_branch;
                // otherwise, stop trying along the current branch
                else
                    break;
                
                // stop trying if reaching the minimum branch length
                if (new_blength <= min_blength_mid)
                    break;
             
                // try at different position along the current branch
                new_blength *= 0.5;

                // get new mid_branch_regions based on the new_blength
                upper_lr_regions->mergeUpperLower(mid_branch_regions, new_blength, *lower_regions, current_blength - new_blength, aln, model, params->threshold_prob);
                
                // compute the placement cost
                new_lh_mid_branch = calculateSamplePlacementCost(mid_branch_regions, sample_regions, default_blength);
            }
            
            // record new best_down_lh_diff
            if (new_best_lh_mid_branch > best_down_lh_diff)
            {
                best_down_lh_diff = new_best_lh_mid_branch;
                best_child_index = node_index;
            }
        }
    }
}

void Tree::seekSamplePlacement(const Index start_node_index, const NumSeqsType seq_name_index, const std::unique_ptr<SeqRegions>& sample_regions, Index& selected_node_index, RealNumType &best_lh_diff, bool &is_mid_branch, RealNumType &best_up_lh_diff, RealNumType &best_down_lh_diff, Index& best_child_index)
{
    // init variables
    // output variables
    selected_node_index = start_node_index;
    // dummy variables
    RealNumType lh_diff_mid_branch = 0;
    RealNumType lh_diff_at_node = 0;
    // stack of nodes to examine positions
    stack<TraversingNode> extended_node_stack;
    extended_node_stack.push(TraversingNode(start_node_index, 0, MIN_NEGATIVE));
    
    // recursively examine positions for placing the new sample
    while (!extended_node_stack.empty())
    {
        TraversingNode current_extended_node = std::move(extended_node_stack.top());
        extended_node_stack.pop();
        const NumSeqsType current_node_vec = current_extended_node.getIndex().getVectorIndex();
        PhyloNode& current_node = nodes[current_node_vec];
        const bool& is_internal = current_node.isInternal();
        
        // NHANLT: debug
        //if (current_node->next && ((current_node->next->neighbor && current_node->next->neighbor->seq_name == "25")
          //                        || (current_node->next->next->neighbor && current_node->next->next->neighbor->seq_name == "25")))
            //cout << "fdsfsd";
    
        // if the current node is a leaf AND the new sample/sequence is strictly less informative than the current node
        // -> add the new sequence into the list of minor sequences of the current node + stop seeking the placement
        if ((!is_internal) && (current_node.getPartialLh(TOP)->compareWithSample(*sample_regions, aln.ref_seq.size(), aln.num_states) == 1))
        {
            current_node.addLessInfoSeqs(seq_name_index);
            selected_node_index = Index();
            return;
        }
        
        const RealNumType current_node_blength = current_node.getUpperLength();
        
        // 1. try first placing as a descendant of the mid-branch point of the branch above the current node
        if (current_node_vec != root_vector_index && current_node_blength > 0)
        {
            examineSamplePlacementMidBranch(selected_node_index, current_node.getMidBranchLh(), best_lh_diff, is_mid_branch, lh_diff_mid_branch, current_extended_node, sample_regions);
        }
        // otherwise, don't consider mid-branch point
        else
            lh_diff_mid_branch = MIN_NEGATIVE;

        // 2. try to place as descendant of the current node (this is skipped if the node has top branch length 0 and so is part of a polytomy).
        if (current_node_vec == root_vector_index || current_node_blength > 0)
        {
            examineSamplePlacementAtNode(selected_node_index, current_node.getTotalLh(), best_lh_diff, is_mid_branch, lh_diff_at_node, lh_diff_mid_branch, best_up_lh_diff, best_down_lh_diff, best_child_index, current_extended_node, sample_regions);
        }
        else
            lh_diff_at_node = current_extended_node.getLhDiff();
        
        // keep trying to place at children nodes, unless the number of attempts has reaches the failure limit
        const short int failure_count = current_extended_node.getFailureCount();
        if ((params->strict_stop_seeking_placement_sample
             && failure_count < params->failure_limit_sample
             && lh_diff_at_node > (best_lh_diff - params->thresh_log_lh_sample))
            || (!params->strict_stop_seeking_placement_sample
                && (failure_count < params->failure_limit_sample
                    || lh_diff_at_node > (best_lh_diff - params->thresh_log_lh_sample))))
        {
            /*for (Index neighbor_index:current_node.getNeighborIndexes(TOP))
                extended_node_stack.push(TraversingNode(neighbor_index, current_extended_node.getFailureCount(), lh_diff_at_node));*/
            if (is_internal)
            {
                extended_node_stack.push(TraversingNode(current_node.getNeighborIndex(RIGHT), failure_count, lh_diff_at_node));
                extended_node_stack.push(TraversingNode(current_node.getNeighborIndex(LEFT), failure_count, lh_diff_at_node));
            }
        }
    }

    // exploration of the tree is finished, and we are left with the node found so far with the best appending likelihood cost. Now we explore placement just below this node for more fine-grained placement within its descendant branches.
    best_down_lh_diff = MIN_NEGATIVE;
    best_child_index = Index();;
    
    // if best position so far is the descendant of a node -> explore further at its children
    if (!is_mid_branch)
    {
        finetuneSamplePlacementAtNode(nodes[selected_node_index.getVectorIndex()], best_down_lh_diff, best_child_index, sample_regions);
    }
}

void Tree::addStartingNodes(const Index& node_index, PhyloNode& node, const Index& other_child_node_index, const RealNumType best_lh_diff, std::stack<std::unique_ptr<UpdatingNode>> &node_stack)
{
    // dummy variables
    const NumSeqsType vec_index = node_index.getVectorIndex();
    PhyloNode& other_child_node = nodes[other_child_node_index.getVectorIndex()];
    
    // node is not the root
    if (vec_index != root_vector_index)
    {
        std::unique_ptr<SeqRegions>& parent_upper_lr_regions = getPartialLhAtNode(node.getNeighborIndex(TOP)); // node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
        std::unique_ptr<SeqRegions>& other_child_node_regions = other_child_node.getPartialLh(TOP); // other_child_node->getPartialLhAtNode(aln, model, threshold_prob);
       
        // add nodes (sibling and parent of the current node) into node_stack which we will need to traverse to update their regions due to the removal of the sub tree
        RealNumType branch_length = other_child_node.getUpperLength(); // other_child_node->length;
        if (node.getUpperLength() > 0)
            branch_length = branch_length > 0 ? branch_length + node.getUpperLength() : node.getUpperLength();
        
        /*node_stack.push(new UpdatingNode(node->neighbor, other_child_node_regions, branch_length, true, best_lh_diff, 0, false));
        node_stack.push(new UpdatingNode(other_child_node, parent_upper_lr_regions, branch_length, true, best_lh_diff, 0, false));*/
        std::unique_ptr<SeqRegions> null_seqregions_ptr1 = nullptr;
        node_stack.push(std::make_unique<UpdatingNode>(std::move(UpdatingNode(node.getNeighborIndex(TOP), std::move(null_seqregions_ptr1), other_child_node_regions, branch_length, true, best_lh_diff, 0))));
        std::unique_ptr<SeqRegions> null_seqregions_ptr2 = nullptr;
        node_stack.push(std::make_unique<UpdatingNode>(std::move(UpdatingNode(other_child_node_index, std::move(null_seqregions_ptr2),  parent_upper_lr_regions, branch_length, true, best_lh_diff, 0))));
    }
    // node is root
    else
    {
        // there is only one sample outside of the subtree doesn't need to be considered
        if (other_child_node.isInternal()) // other_child_node->next)
        {
            // add nodes (children of the sibling of the current node) into node_stack which we will need to traverse to update their regions due to the removal of the sub tree
            const Index grand_child_1_index = other_child_node.getNeighborIndex(RIGHT);
            const Index grand_child_2_index = other_child_node.getNeighborIndex(LEFT);
            PhyloNode& grand_child_1 = nodes[grand_child_1_index.getVectorIndex()]; // other_child_node->next->neighbor;
            PhyloNode& grand_child_2 = nodes[grand_child_2_index.getVectorIndex()]; // other_child_node->next->next->neighbor;
            
            // always unique_ptr<SeqRegions>& -> always automatically delete
            // SeqRegions* up_lr_regions_1 = grand_child_2->computeTotalLhAtNode(aln, model, threshold_prob, true, false, grand_child_2->length);
            // std::unique_ptr<SeqRegions> up_lr_regions_1 = grand_child_2.computeTotalLhAtNode(other_child_node, aln, model, threshold_prob, true, grand_child_2.getUpperLength());
            std::unique_ptr<SeqRegions> up_lr_regions_1 = std::move(grand_child_2.getPartialLh(TOP)->computeTotalLhAtRoot(aln.num_states, model, grand_child_2.getUpperLength()));
            
            // node_stack.push(new UpdatingNode(grand_child_1, up_lr_regions_1, grand_child_1->length, true, best_lh_diff, 0, true));
            std::unique_ptr<SeqRegions> null_seqregions_ptr1 = nullptr;
            node_stack.push(std::make_unique<UpdatingNode>(std::move(UpdatingNode(grand_child_1_index, std::move(up_lr_regions_1), null_seqregions_ptr1, grand_child_1.getUpperLength(), true, best_lh_diff, 0))));
            
            
            // SeqRegions* up_lr_regions_2 = grand_child_1->computeTotalLhAtNode(aln, model, threshold_prob, true, false, grand_child_1->length);
            std::unique_ptr<SeqRegions> up_lr_regions_2 = std::move(grand_child_1.getPartialLh(TOP)->computeTotalLhAtRoot(aln.num_states, model, grand_child_1.getUpperLength()));
            
            // node_stack.push(new UpdatingNode(grand_child_2, up_lr_regions_2, grand_child_2->length, true, best_lh_diff, 0, true));
            std::unique_ptr<SeqRegions> null_seqregions_ptr2 = nullptr;
            node_stack.push(std::make_unique<UpdatingNode>(std::move(UpdatingNode(grand_child_2_index, std::move(up_lr_regions_2), null_seqregions_ptr2, grand_child_2.getUpperLength(), true, best_lh_diff, 0))));
        }
    }
}

// NOTE: top_node != null <=> case when crawling up from child to parent
// otherwise, top_node == null <=> case we are moving from a parent to a child
bool Tree::examineSubtreePlacementMidBranch(Index& best_node_index, PhyloNode& current_node, RealNumType& best_lh_diff, bool& is_mid_branch, RealNumType& lh_diff_at_node, RealNumType& lh_diff_mid_branch, RealNumType& best_up_lh_diff, RealNumType& best_down_lh_diff, std::unique_ptr<UpdatingNode>& updating_node, const std::unique_ptr<SeqRegions>& subtree_regions, const RealNumType threshold_prob, const RealNumType removed_blength, const Index top_node_index, std::unique_ptr<SeqRegions>& bottom_regions)
{
    const bool top_node_exists = (top_node_index.getMiniIndex() != UNDEFINED);
    const Index updating_node_index = updating_node->getIndex();
    const Index at_node_index = top_node_exists ? top_node_index : updating_node_index;
    const NumSeqsType at_node_vec = at_node_index.getVectorIndex();
    PhyloNode& at_node = top_node_exists ? nodes[at_node_vec] : current_node;
    
    std::unique_ptr<SeqRegions> new_mid_branch_regions = nullptr;
    // get or recompute the lh regions at the mid-branch position
    if (updating_node->needUpdate())
    {
        // recompute mid_branch_regions in case when crawling up from child to parent
        if (top_node_exists)
        {
            /*Node* other_child = updating_node->node->getOtherNextNode()->neighbor;
            SeqRegions* other_child_lower_regions = other_child->getPartialLhAtNode(aln, model, threshold_prob);
            other_child_lower_regions->mergeTwoLowers(bottom_regions, other_child->length, *updating_node->incoming_regions, updating_node->branch_length, aln, model, threshold_prob);*/
            const Index other_child_index = current_node.getNeighborIndex(updating_node_index.getFlipMiniIndex());
            PhyloNode& other_child = nodes[other_child_index.getVectorIndex()]; // updating_node->node->getOtherNextNode()->neighbor;
            const std::unique_ptr<SeqRegions>& other_child_lower_regions = other_child.getPartialLh(TOP); // other_child->getPartialLhAtNode(aln, model, threshold_prob);
            // other_child_lower_regions->mergeTwoLowers(bottom_regions, other_child->length, *updating_node->incoming_regions, updating_node->branch_length, aln, model, threshold_prob);
            other_child_lower_regions->mergeTwoLowers(bottom_regions, other_child.getUpperLength(), *updating_node->getIncomingRegions(), updating_node->getBranchLength(), aln, model, threshold_prob);
                                
            // skip if bottom_regions is null (inconsistent)
            if (!bottom_regions)
            {
                // delete updating_node
                // delete updating_node;

                return false; // continue;
            }
           
            // compute new mid-branch regions
            /*SeqRegions* upper_lr_regions = top_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
            RealNumType mid_branch_length = top_node->length * 0.5;
            upper_lr_regions->mergeUpperLower(mid_branch_regions, mid_branch_length, *bottom_regions, mid_branch_length, aln, model, threshold_prob);*/
            
            const std::unique_ptr<SeqRegions>& upper_lr_regions = getPartialLhAtNode(at_node.getNeighborIndex(TOP)); // top_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
            const RealNumType mid_branch_length = at_node.getUpperLength() * 0.5;
            upper_lr_regions->mergeUpperLower(new_mid_branch_regions, mid_branch_length, *bottom_regions, mid_branch_length, aln, model, threshold_prob);
        }
        // recompute mid_branch_regions in case we are moving from a parent to a child
        else
        {
            /*SeqRegions* lower_regions = updating_node->node->getPartialLhAtNode(aln, model, threshold_prob);
            RealNumType mid_branch_length = updating_node->branch_length * 0.5;
            updating_node->incoming_regions->mergeUpperLower(new_mid_branch_regions, mid_branch_length, *lower_regions, mid_branch_length, aln, model, threshold_prob);*/
            const std::unique_ptr<SeqRegions>& lower_regions = current_node.getPartialLh(TOP); // getPartialLhAtNode(updating_node_index);
            const RealNumType mid_branch_length = updating_node->getBranchLength() * 0.5;
            updating_node->getIncomingRegions()->mergeUpperLower(new_mid_branch_regions, mid_branch_length, *lower_regions, mid_branch_length, aln, model, threshold_prob);
        }
    }
    
    std::unique_ptr<SeqRegions>& mid_branch_regions = updating_node->needUpdate() ? new_mid_branch_regions : at_node.getMidBranchLh();
    
    // skip if mid_branch_regions is null (branch length == 0)
    if (!mid_branch_regions)
    {
        // delete bottom_regions if it's existed
        // if (bottom_regions) delete bottom_regions;
                            
        // delete updating_node
        // delete updating_node;
        
        return false; // continue;
    }
    
    // compute the placement cost
    // if (search_subtree_placement)
    lh_diff_mid_branch = calculateSubTreePlacementCost(mid_branch_regions, subtree_regions, removed_blength);
    //else
       // lh_diff_mid_branch = calculateSamplePlacementCost( mid_branch_regions, subtree_regions, removed_blength);
    
    if (top_node_exists && best_node_index.getVectorIndex() == at_node_vec)// top_node && best_node == top_node) // only update in case when crawling up from child to parent
        best_up_lh_diff = lh_diff_mid_branch;
    
    // if this position is better than the best position found so far -> record it
    if (lh_diff_mid_branch > best_lh_diff)
    {
        best_node_index = at_node_index;
        best_lh_diff = lh_diff_mid_branch;
        is_mid_branch = true;
        updating_node->setFailureCount(0);
        if (top_node_exists) best_down_lh_diff = lh_diff_at_node; // only update in case when crawling up from child to parent
    }
    else if (top_node_exists && lh_diff_at_node >= (best_lh_diff - threshold_prob))
        best_up_lh_diff = lh_diff_mid_branch;
        
    // delete mid_branch_regions
    // if (updating_node->need_updating) delete mid_branch_regions;
    
    // no error
    return true;
}

// NOTE: top_node != null <=> case when crawling up from child to parent
// otherwise, top_node == null <=> case we are moving from a parent to a child
bool Tree::examineSubTreePlacementAtNode(Index& best_node_index, PhyloNode& current_node, RealNumType &best_lh_diff, bool& is_mid_branch, RealNumType& lh_diff_at_node, RealNumType& lh_diff_mid_branch, RealNumType &best_up_lh_diff, RealNumType &best_down_lh_diff, std::unique_ptr<UpdatingNode>& updating_node, const std::unique_ptr<SeqRegions>& subtree_regions, const RealNumType threshold_prob, const RealNumType removed_blength, const Index top_node_index)
{
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    
    // Node* at_node = top_node ? top_node: updating_node->node;
    const bool top_node_exits = top_node_index.getMiniIndex() != UNDEFINED;
    const Index updating_node_index = updating_node->getIndex();
    const Index at_node_index = top_node_exits ? top_node_index : updating_node_index;
    PhyloNode& at_node = top_node_exits ? nodes[at_node_index.getVectorIndex()] : current_node;
    
    std::unique_ptr<SeqRegions> new_at_node_regions = nullptr;
    const bool need_updating = updating_node->needUpdate();
    if (updating_node->needUpdate())
    {
        // get or recompute the lh regions at the current node position
        const std::unique_ptr<SeqRegions>& updating_node_partial = getPartialLhAtNode(updating_node_index); // updating_node->node->getPartialLhAtNode(aln, model, threshold_prob);
        if (top_node_exits)
        {
            updating_node_partial->mergeUpperLower(new_at_node_regions, -1, *updating_node->getIncomingRegions(), updating_node->getBranchLength(), aln, model, threshold_prob);
        }
        else
        {
            updating_node->getIncomingRegions()->mergeUpperLower(new_at_node_regions, updating_node->getBranchLength(), *updating_node_partial, -1, aln, model, threshold_prob);
        }
        
        // skip if at_node_regions is null (branch length == 0)
        if (!new_at_node_regions)
        {
            // delete updating_node
            // delete updating_node;
            
            // continue;
            return false;
        }
        
        // stop updating if the difference between the new and old regions is insignificant
        if  (!new_at_node_regions->areDiffFrom(*at_node.getTotalLh(), seq_length, num_states, &params.value()))
            updating_node->setUpdate(false);
    }
    // else
    const std::unique_ptr<SeqRegions>& at_node_regions = need_updating ? new_at_node_regions : at_node.getTotalLh();
    
    //if (search_subtree_placement)
    lh_diff_at_node = calculateSubTreePlacementCost(at_node_regions, subtree_regions, removed_blength);
    //else
        //lh_diff_at_node = calculateSamplePlacementCost(at_node_regions, subtree_regions, removed_blength);
    
    // if this position is better than the best position found so far -> record it
    if (lh_diff_at_node > best_lh_diff)
    {
        best_node_index = at_node_index;
        best_lh_diff = lh_diff_at_node;
        is_mid_branch = false;
        updating_node->setFailureCount(0);
        if (!top_node_exits) best_up_lh_diff = lh_diff_mid_branch; // only update in case we are moving from a parent to a child
    }
    else if (!top_node_exits && lh_diff_mid_branch >= (best_lh_diff - threshold_prob)) // only update in case we are moving from a parent to a child
    {
        best_up_lh_diff = updating_node->getLhDiff();
        best_down_lh_diff = lh_diff_at_node;
    }
    // placement at current node is considered failed if placement likelihood is not improved by a certain margin compared to best placement so far for the nodes above it.
    else if (lh_diff_at_node < (updating_node->getLhDiff() - params->thresh_log_lh_failure))
        updating_node->increaseFailureCount();
        
    // delete at_node_regions
    // if (delete_at_node_regions) delete at_node_regions;
    
    // no error
    return true;
}

bool keepTraversing(const RealNumType& best_lh_diff, const RealNumType& lh_diff_at_node, const bool& strict_stop_seeking_placement_subtree, const std::unique_ptr<UpdatingNode>& updating_node, const int& failure_limit_subtree, const RealNumType& thresh_log_lh_subtree, const bool able2traverse = true)
{
    //if (search_subtree_placement)
    //{
    if (strict_stop_seeking_placement_subtree)
    {
        if (updating_node->getFailureCount() <= failure_limit_subtree && lh_diff_at_node > (best_lh_diff - thresh_log_lh_subtree) && able2traverse)
            return true;
    }
    else
    {
        if ((updating_node->getFailureCount() <= failure_limit_subtree || lh_diff_at_node > (best_lh_diff - thresh_log_lh_subtree))
            && able2traverse)
            return true;
    }
    //}
    //else
    //{
      //  if (params->strict_stop_seeking_placement_sample)
        //{
          //  if (updating_node->failure_count <= params->failure_limit_sample && lh_diff_at_node > (best_lh_diff - params->thresh_log_lh_sample)
            //    && updating_node->node->next)
              //      return true;
       // }
       // else
        //{
          //  if ((updating_node->failure_count <= params->failure_limit_sample || lh_diff_at_node > (best_lh_diff - params->thresh_log_lh_sample))
            //    && updating_node->node->next)
              //      return true;
        //}
    //}
    
    // default
    return false;
}

void Tree::addChildSeekSubtreePlacement(const Index child_1_index, const Index child_2_index, PhyloNode& child_1, PhyloNode& child_2, const RealNumType& lh_diff_at_node, const std::unique_ptr<UpdatingNode>& updating_node, std::stack<std::unique_ptr<UpdatingNode>>& node_stack, const RealNumType threshold_prob)
{
    // get or recompute the upper left/right regions of the children node
    if (updating_node->needUpdate())
    {
        std::unique_ptr<SeqRegions> upper_lr_regions = nullptr;
        const std::unique_ptr<SeqRegions>& lower_regions = child_2.getPartialLh(TOP); // ->getPartialLhAtNode(aln, model, threshold_prob);
        
        updating_node->getIncomingRegions()->mergeUpperLower(upper_lr_regions, updating_node->getBranchLength(), *lower_regions, child_2.getUpperLength(), aln, model, threshold_prob);
        
        // traverse to this child's subtree
        if (upper_lr_regions)
        {
            // node_stack.push(new UpdatingNode(child_1, upper_lr_regions, child_1->length, updating_node->need_updating, lh_diff_at_node, updating_node->failure_count, updating_node->need_updating));
            
            std::unique_ptr<SeqRegions> null_seqregions_ptr = nullptr;
            node_stack.push(std::make_unique<UpdatingNode>(std::move(UpdatingNode(child_1_index, std::move(upper_lr_regions), null_seqregions_ptr, child_1.getUpperLength(), updating_node->needUpdate(), lh_diff_at_node, updating_node->getFailureCount()))));
        }
    }
    else
    {
        std::unique_ptr<SeqRegions>& upper_lr_regions = getPartialLhAtNode(child_1.getNeighborIndex(TOP));
        const std::unique_ptr<SeqRegions> null_seqregions_ptr = nullptr;
        /*if (child_1->neighbor->partial_lh)
            upper_lr_regions = child_1->neighbor->getPartialLhAtNode(aln, model, threshold_prob);*/
        
        // traverse to this child's subtree
        if (upper_lr_regions)
        {
            // node_stack.push(new UpdatingNode(child_1, upper_lr_regions, child_1->length, updating_node->need_updating, lh_diff_at_node, updating_node->failure_count, updating_node->need_updating));
            std::unique_ptr<SeqRegions> null_seqregions_ptr = nullptr;
            node_stack.push(std::make_unique<UpdatingNode>(std::move(UpdatingNode(child_1_index, std::move(null_seqregions_ptr), upper_lr_regions, child_1.getUpperLength(), updating_node->needUpdate(), lh_diff_at_node, updating_node->getFailureCount()))));
        }
    }
}

bool Tree::addNeighborsSeekSubtreePlacement(PhyloNode& current_node, const Index other_child_index, std::unique_ptr<SeqRegions>&& bottom_regions, const RealNumType& lh_diff_at_node, const std::unique_ptr<UpdatingNode>& updating_node, std::stack<std::unique_ptr<UpdatingNode>>& node_stack, const RealNumType threshold_prob)
{
    ASSERT(other_child_index.getMiniIndex() == TOP);
    PhyloNode& other_child = nodes[other_child_index.getVectorIndex()];
    const Index updating_node_index = updating_node->getIndex();
    const MiniIndex updating_node_mini = updating_node_index.getMiniIndex();
    
    // keep crawling up into parent and sibling node
    // case the node is not the root
    if (updating_node_index.getVectorIndex() != root_vector_index)
    {
        // first pass the crawling down the other child (sibling)
        
        // get or recompute the upper left/right regions of the sibling node
        if (updating_node->needUpdate())
        {
            const std::unique_ptr<SeqRegions>& parent_upper_lr_regions = getPartialLhAtNode(current_node.getNeighborIndex(TOP)); // top_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
            std::unique_ptr<SeqRegions> upper_lr_regions = nullptr;
            
            parent_upper_lr_regions->mergeUpperLower(upper_lr_regions, current_node.getUpperLength(), *updating_node->getIncomingRegions(), updating_node->getBranchLength(), aln, model, threshold_prob);
            
            if (!upper_lr_regions)
            {
                // delete bottom_regions if it's existed
                // if (bottom_regions) delete bottom_regions;
                
                // delete updating_node
                // delete updating_node;
                
                return false; // continue;
            }
            else
            {
                std::unique_ptr<SeqRegions> null_seqregions_ptr = nullptr;
                //node_stack.push(std::make_unique<UpdatingNode>(std::move(UpdatingNode(other_child_index, std::move(upper_lr_regions), null_seqregions_ptr, other_child->length, updating_node->need_updating_, lh_diff_at_node, updating_node->failure_count_))));
                node_stack.push(std::make_unique<UpdatingNode>(std::move(UpdatingNode(other_child_index, std::move(upper_lr_regions), null_seqregions_ptr, other_child.getUpperLength(), updating_node->needUpdate(), lh_diff_at_node, updating_node->getFailureCount()))));
            }
        }
        else
        {
            std::unique_ptr<SeqRegions>& upper_lr_regions = current_node.getPartialLh(updating_node_mini);
            
            if (!upper_lr_regions)//updating_node->node->partial_lh)
            {
                // delete bottom_regions if it's existed
                // if (bottom_regions) delete bottom_regions;
                
                // delete updating_node
                // delete updating_node;
                
                return false; // continue;
            }
            else
            {
                /*upper_lr_regions = updating_node->node->getPartialLhAtNode(aln, model, threshold_prob);
                node_stack.push(new UpdatingNode(other_child, upper_lr_regions, other_child->length, updating_node->need_updating, lh_diff_at_node, updating_node->failure_count, updating_node->need_updating));*/
                std::unique_ptr<SeqRegions> null_seqregions_ptr = nullptr;
                node_stack.push(std::make_unique<UpdatingNode>(std::move(UpdatingNode(other_child_index, std::move(null_seqregions_ptr), upper_lr_regions, other_child.getUpperLength(), updating_node->needUpdate(), lh_diff_at_node, updating_node->getFailureCount()))));
            }
        }

        // add sibling node to node_stack for traversing later; skip if upper_lr_regions is null (inconsistent)
        /*if (!upper_lr_regions)
        {
            // delete bottom_regions if it's existed
            // if (bottom_regions) delete bottom_regions;
            
            // delete updating_node
            // delete updating_node;
            
            return false; // continue;
        }
        else
            node_stack.push(new UpdatingNode(other_child, upper_lr_regions, other_child->length, updating_node->need_updating, lh_diff_at_node, updating_node->failure_count, updating_node->need_updating));*/
        
        // now pass the crawling up to the parent node
        // get or recompute the bottom regions (comming from 2 children) of the parent node
        if (updating_node->needUpdate())
        {
            if (!bottom_regions)
            {
                const std::unique_ptr<SeqRegions>& other_child_lower_regions = other_child.getPartialLh(TOP); // other_child->getPartialLhAtNode(aln, model, threshold_prob);
                // other_child_lower_regions->mergeTwoLowers(bottom_regions, other_child->length, *updating_node->incoming_regions, updating_node->branch_length_, aln, model, threshold_prob);
                other_child_lower_regions->mergeTwoLowers(bottom_regions, other_child.getUpperLength(), *updating_node->getIncomingRegions(), updating_node->getBranchLength(), aln, model, threshold_prob);
                
                // skip if bottom_regions is null (inconsistent)
                if (!bottom_regions)
                {
                    // delete updating_node
                    // delete updating_node;
                    
                    return false; // continue;
                }
            }
            
            // node_stack.push(std::make_unique<UpdatingNode>(std::move(UpdatingNode(top_node->neighbor, bottom_regions, top_node->length, updating_node->need_updating, lh_diff_at_node, updating_node->failure_count, updating_node->need_updating))));
            std::unique_ptr<SeqRegions> null_seqregions_ptr = nullptr;
            node_stack.push(std::make_unique<UpdatingNode>(std::move(UpdatingNode(current_node.getNeighborIndex(TOP), std::move(bottom_regions), null_seqregions_ptr, current_node.getUpperLength(), updating_node->needUpdate(), lh_diff_at_node, updating_node->getFailureCount()))));
        }
        else
        {
            // if (bottom_regions) delete bottom_regions;
            bottom_regions = nullptr;
            std::unique_ptr<SeqRegions>& bottom_regions_ref = current_node.getPartialLh(TOP); // top_node->getPartialLhAtNode(aln, model, threshold_prob);
            
            std::unique_ptr<SeqRegions> null_seqregions_ptr = nullptr;
            // node_stack.push(std::make_unique<UpdatingNode>(std::move(UpdatingNode(top_node->neighbor, bottom_regions, top_node->length, updating_node->need_updating, lh_diff_at_node, updating_node->failure_count, updating_node->need_updating))));
            node_stack.push(std::make_unique<UpdatingNode>(std::move(UpdatingNode(current_node.getNeighborIndex(TOP), std::move(null_seqregions_ptr), bottom_regions_ref, current_node.getUpperLength(), updating_node->needUpdate(), lh_diff_at_node, updating_node->getFailureCount()))));
        }
        
        // add the parent node to node_stack for traversing later
        /* node_stack.push(new UpdatingNode(top_node->neighbor, bottom_regions, top_node->length, updating_node->need_updating, lh_diff_at_node, updating_node->failure_count, updating_node->need_updating));*/
    }
    // now consider case of root node -> only need to care about the sibling node
    else
    {
        // get or recompute the upper left/right regions of the sibling node
        if (updating_node->needUpdate())
        {
            std::unique_ptr<SeqRegions> upper_lr_regions = updating_node->getIncomingRegions()->computeTotalLhAtRoot(aln.num_states, model, updating_node->getBranchLength());
            
            std::unique_ptr<SeqRegions> null_seqregions_ptr = nullptr;
            node_stack.push(std::make_unique<UpdatingNode>(std::move(UpdatingNode(other_child_index, std::move(upper_lr_regions), null_seqregions_ptr, other_child.getUpperLength(), updating_node->needUpdate(), lh_diff_at_node, updating_node->getFailureCount()))));
        }
        else
        {
            std::unique_ptr<SeqRegions>& upper_lr_regions = current_node.getPartialLh(updating_node_mini); // updating_node->node->getPartialLhAtNode(aln, model, threshold_prob);
            
            std::unique_ptr<SeqRegions> null_seqregions_ptr = nullptr;
            node_stack.push(std::make_unique<UpdatingNode>(std::move(UpdatingNode(other_child_index, std::move(null_seqregions_ptr), upper_lr_regions, other_child.getUpperLength(), updating_node->needUpdate(), lh_diff_at_node, updating_node->getFailureCount()))));
        }
        
        // add the sibling node to node_stack for traversing later
        // node_stack.push(new UpdatingNode(other_child, upper_lr_regions, other_child->length, updating_node->need_updating, lh_diff_at_node, updating_node->failure_count, updating_node->need_updating));
        
        // delete bottom_regions
        // if (bottom_regions) delete bottom_regions;
    }
    
    return true;
}

void Tree::seekSubTreePlacement(Index& best_node_index, RealNumType &best_lh_diff, bool &is_mid_branch, RealNumType &best_up_lh_diff, RealNumType &best_down_lh_diff, Index& best_child_index, const bool short_range_search, const Index child_node_index, RealNumType &removed_blength) //, bool search_subtree_placement, SeqRegions* sample_regions)
{
    // init variables
    PhyloNode& child_node = nodes[child_node_index.getVectorIndex()];
    const Index node_index = child_node.getNeighborIndex(TOP);
    const NumSeqsType vec_index = node_index.getVectorIndex();
    PhyloNode& node = nodes[vec_index]; // child_node->neighbor->getTopNode();
    const Index other_child_node_index = node.getNeighborIndex(node_index.getFlipMiniIndex()); // child_node->neighbor->getOtherNextNode()->neighbor;
    best_node_index = node_index;
    const std::unique_ptr<SeqRegions>& subtree_regions = child_node.getPartialLh(TOP); // child_node->getPartialLhAtNode(aln, model, threshold_prob); nullptr;
    // stack of nodes to examine positions
    stack<std::unique_ptr<UpdatingNode>> node_stack;
    // dummy variables
    const RealNumType threshold_prob = params->threshold_prob;
    RealNumType lh_diff_mid_branch = 0;
    RealNumType lh_diff_at_node = 0;
    // const std::unique_ptr<SeqRegions> null_seqregions_ptr = nullptr;
    // const std::unique_ptr<SeqRegions>& parent_upper_lr_regions = vec_index == root_vector_index ? null_seqregions_ptr : getPartialLhAtNode(node.getNeighborIndex(TOP));
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    
    // get/init approximation params
    bool strict_stop_seeking_placement_subtree = params->strict_stop_seeking_placement_subtree;
    int failure_limit_subtree = params->failure_limit_subtree;
    RealNumType thresh_log_lh_subtree = params->thresh_log_lh_subtree;
    
    // for short range topology search
    if (short_range_search)
    {
        strict_stop_seeking_placement_subtree = params->strict_stop_seeking_placement_subtree_short_search;
        failure_limit_subtree = params->failure_limit_subtree_short_search;
        thresh_log_lh_subtree = params->thresh_log_lh_subtree_short_search;
    }

    // search a placement for a subtree
    //if (search_subtree_placement)
    //{
    // get the lower regions of the child node
    // subtree_regions = child_node->getPartialLhAtNode(aln, model, threshold_prob);
    
    // add starting nodes to start seek placement for the subtree
    addStartingNodes(node_index, node, other_child_node_index, best_lh_diff, node_stack);
    
    //}
     // search a placement for a new sample
    //else
    //{
        // get the regions of the input sample
      //  subtree_regions = sample_regions;
        //RealNumType down_lh = is_mid_branch ? best_down_lh_diff : best_lh_diff;
        
        // node is not the root
  //      if (node != root)
    //    {
      //      SeqRegions* lower_regions = new SeqRegions(node->getPartialLhAtNode(aln, model, threshold_prob), num_states);
        //    parent_upper_lr_regions = node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
            
            // add the parent node of the current node into node_stack for traversing to seek the placement for the new sample
    //        node_stack.push(new UpdatingNode(node->neighbor, lower_regions, node->length, false, down_lh, 0));
      //  }
        
        // add the children nodes of the current node into node_stack for traversing to seek the placement for the new sample
    //    Node* neighbor_node;
      //  FOR_NEIGHBOR(node, neighbor_node)
        //{
         //   SeqRegions* upper_lr_regions = new SeqRegions(neighbor_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob), num_states);
          //  node_stack.push(new UpdatingNode(neighbor_node, upper_lr_regions, neighbor_node->length, false, down_lh, 0));
        //}
     //}
    
    // examine each node in the node stack to seek the "best" placement
    while (!node_stack.empty())
    {
        // extract updating_node from stack
        std::unique_ptr<UpdatingNode> updating_node = std::move(node_stack.top());
        node_stack.pop();
        const Index current_node_index = updating_node->getIndex();
        const NumSeqsType current_node_vec = current_node_index.getVectorIndex();
        PhyloNode& current_node = nodes[current_node_vec];
        
        // consider the case we are moving from a parent to a child
        if (current_node_index.getMiniIndex() == TOP)
        {
            if (current_node.getUpperLength() > 0) // updating_node->node->length > 0)
            {
                //  try to append mid-branch
                // avoid adding to the old position where the subtree was just been removed from
                if (current_node_vec != root_vector_index && current_node.getNeighborIndex(TOP).getVectorIndex() != vec_index) // updating_node->node != root && updating_node->node->neighbor->getTopNode() != node)
                {
                    std::unique_ptr<SeqRegions> bottom_regions = nullptr;
                    if (!examineSubtreePlacementMidBranch(best_node_index, current_node, best_lh_diff, is_mid_branch, lh_diff_at_node, lh_diff_mid_branch, best_up_lh_diff, best_down_lh_diff, updating_node, subtree_regions, threshold_prob, removed_blength, Index(), bottom_regions)) continue;
                }
                // set the placement cost at the mid-branch position the most negative value if branch length is zero -> we can't place the subtree on that branch
                else
                    lh_diff_mid_branch = MIN_NEGATIVE;
                    
                // now try appending exactly at node
                if(!examineSubTreePlacementAtNode(best_node_index, current_node, best_lh_diff, is_mid_branch, lh_diff_at_node, lh_diff_mid_branch, best_up_lh_diff, best_down_lh_diff, updating_node, subtree_regions, threshold_prob, removed_blength, Index())) continue;
            }
            // set the placement cost at the current node position at the most negative value if branch length is zero -> we can't place the subtree on that branch
            else
                lh_diff_at_node = updating_node->getLhDiff();
            
            // keep crawling down into children nodes unless the stop criteria for the traversal are satisfied.
            // check the stop criteria
            // keep traversing further down to the children
            if (keepTraversing(best_lh_diff, lh_diff_at_node, strict_stop_seeking_placement_subtree, updating_node, failure_limit_subtree, thresh_log_lh_subtree, current_node.isInternal()))// updating_node->node->next))
            {
                /*Node* child_1 = updating_node->node->getOtherNextNode()->neighbor;
                Node* child_2 = child_1->neighbor->getOtherNextNode()->neighbor;*/
                const Index child_1_index = current_node.getNeighborIndex(RIGHT); // updating_node->node->getOtherNextNode()->neighbor;
                const Index child_2_index = current_node.getNeighborIndex(LEFT); // child_1->neighbor->getOtherNextNode()->neighbor;
                PhyloNode& child_1 = nodes[child_1_index.getVectorIndex()];
                PhyloNode& child_2 = nodes[child_2_index.getVectorIndex()];
                
                // add child_1 to node_stack
                addChildSeekSubtreePlacement(child_1_index, child_2_index, child_1, child_2, lh_diff_at_node, updating_node, node_stack, threshold_prob);
                
                // add child_2 to node_stack
                addChildSeekSubtreePlacement(child_2_index, child_1_index, child_2, child_1, lh_diff_at_node, updating_node, node_stack, threshold_prob);
            }
        }
        // case when crawling up from child to parent
        else
        {
            // Node* top_node = updating_node->node->getTopNode();
            const Index top_node_index = Index(current_node_vec, TOP);
            
            // append directly at the node
            if (current_node.getUpperLength() > 0 || current_node_vec == root_vector_index) // top_node->length > 0 || top_node == root)
            {
                if (!examineSubTreePlacementAtNode(best_node_index, current_node, best_lh_diff, is_mid_branch, lh_diff_at_node, lh_diff_mid_branch, best_up_lh_diff, best_down_lh_diff, updating_node, subtree_regions, threshold_prob, removed_blength, top_node_index)) continue;
            }
            // if placement cost at new position gets worse -> restore to the old one
            else
                lh_diff_at_node = updating_node->getLhDiff();

            // try appending mid-branch
            const Index other_child_index = current_node.getNeighborIndex(current_node_index.getFlipMiniIndex()); // updating_node->node->getOtherNextNode()->neighbor;
            std::unique_ptr<SeqRegions> bottom_regions = nullptr;
            if (current_node.getUpperLength() > 0 && current_node_vec != root_vector_index) // top_node->length > 0 && top_node != root)
            {
                if (!examineSubtreePlacementMidBranch(best_node_index, current_node, best_lh_diff, is_mid_branch, lh_diff_at_node, lh_diff_mid_branch, best_up_lh_diff, best_down_lh_diff, updating_node, subtree_regions, threshold_prob, removed_blength, top_node_index, bottom_regions)) continue;
            }
            // set the placement cost at the mid-branch position at the most negative value if branch length is zero -> we can't place the subtree on that branch
            // NHANLT: we actually don't need to do that since lh_diff_mid_branch will never be read
            // else
                // lh_diff_mid_branch = MIN_NEGATIVE;
            
            // check stop rule of the traversal process
            // keep traversing upwards
            if (keepTraversing(best_lh_diff, lh_diff_at_node, strict_stop_seeking_placement_subtree, updating_node, failure_limit_subtree, thresh_log_lh_subtree))
            {
                // if(!addNeighborsSeekSubtreePlacement(top_node_index, other_child_index, parent_upper_lr_regions, bottom_regions, lh_diff_at_node, updating_node, node_stack, threshold_prob)) continue;
                if(!addNeighborsSeekSubtreePlacement(current_node, other_child_index, std::move(bottom_regions), lh_diff_at_node, updating_node, node_stack, threshold_prob)) continue;
            }
            /*else
            {
                // delete bottom_regions if it's existed
                if (bottom_regions) delete bottom_regions;
            }*/
        }
        
        // delete updating_node
        // delete updating_node;
    }
    
    // ############ KEEP this section DISABLE/COMMENTED OUT ############
    // exploration of the tree is finished, and we are left with the node found so far with the best appending likelihood cost.
    // Now we explore placement just below this node for more fine-grained placement within its descendant branches.
    /*if (!search_subtree_placement)
    {
        best_down_lh_diff = MIN_NEGATIVE;
        best_child = NULL;
        
        if (is_mid_branch)
        {
            // go upward until we reach the parent node of a polytomy
            Node* parent_node = best_node->neighbor->getTopNode();
            while (parent_node->length <= 0 && parent_node != root)
                parent_node = parent_node->neighbor->getTopNode();
            
            best_up_lh_diff = calculateSamplePlacementCost(parent_node->total_lh, subtree_regions, removed_blength);
            child_node = best_node;
        }
        else
        {
            // current node might be part of a polytomy (represented by 0 branch lengths) so we want to explore all the children of the current node to find out if the best placement is actually in any of the branches below the current node.
            Node* neighbor_node;
            stack<Node*> new_node_stack;
            FOR_NEIGHBOR(best_node, neighbor_node)
                new_node_stack.push(neighbor_node);
            
            while (!new_node_stack.empty())
            {
                Node* node = new_node_stack.top();
                new_node_stack.pop();
                
                if (node->length <= 0)
                {
                    FOR_NEIGHBOR(node, neighbor_node)
                        new_node_stack.push(neighbor_node);
                }
                // now try to place on the current branch below the best node, at an height above the mid-branch.
                else
                {
                    // now try to place on the current branch below the best node, at an height above the mid-branch.
                    RealNumType new_blength = node->length * 0.5;
                    RealNumType new_best_lh_mid_branch = MIN_NEGATIVE;
                    SeqRegions* upper_lr_regions = node->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob);
                    SeqRegions* lower_regions = node->getPartialLhAtNode(aln, model, params->threshold_prob);
                    SeqRegions* mid_branch_regions = new SeqRegions(node->mid_branch_lh, aln.num_states);

                    // try to place new sample along the upper half of the current branch
                    while (true)
                    {
                        // compute the placement cost
                        RealNumType new_lh_mid_branch = calculateSamplePlacementCost(mid_branch_regions, subtree_regions, removed_blength);
                        
                        // record new_best_lh_mid_branch
                        if (new_lh_mid_branch > new_best_lh_mid_branch)
                            new_best_lh_mid_branch = new_lh_mid_branch;
                        // otherwise, stop trying along the current branch
                        else
                            break;
                        
                        // stop trying if reaching the minimum branch length
                        if (new_blength <= min_blength_mid)
                            break;
     
                         // try at different position along the current branch
                         new_blength *= 0.5;
                     
                        // get new mid_branch_regions based on the new_blength
                        upper_lr_regions->mergeUpperLower(mid_branch_regions, new_blength, lower_regions, node->length - new_blength, aln, model, params->threshold_prob);
                    }
                    
                    //RealNumType new_best_lh_mid_branch = calculateSamplePlacementCost(node->mid_branch_lh, sample_regions, default_blength);
                    
                    // record new best_down_lh_diff
                    if (new_best_lh_mid_branch > best_down_lh_diff)
                    {
                        best_down_lh_diff = new_best_lh_mid_branch;
                        best_child = node;
                    }
                }
            }
        }
    }*/
    // ############ KEEP this section DISABLE/COMMENTED OUT ############
}

void Tree::applySPR(const Index subtree_index, PhyloNode& subtree, const Index best_node_index, const bool is_mid_branch, const RealNumType branch_length, const RealNumType best_lh_diff)
{
    // record the SPR applied at this subtree
    subtree.setSPRApplied(true);

    // remove subtree from the tree
    const Index parent_index = subtree.getNeighborIndex(TOP);
    PhyloNode& parent_subtree = nodes[parent_index.getVectorIndex()]; // subtree->neighbor->getTopNode();
    const Index sibling_index = parent_subtree.getNeighborIndex(parent_index.getFlipMiniIndex());
    PhyloNode& sibling_subtree = nodes[sibling_index.getVectorIndex()]; // subtree->neighbor->getOtherNextNode()->neighbor;
    const RealNumType threshold_prob = params->threshold_prob;
    const StateType num_states = aln.num_states;
    
    // connect grandparent to sibling
    const Index grandparent_index = parent_subtree.getNeighborIndex(TOP);
    if (parent_index.getVectorIndex() != root_vector_index)
    {
        // parent_subtree->neighbor->neighbor = sibling_subtree;
        nodes[grandparent_index.getVectorIndex()].setNeighborIndex(grandparent_index.getMiniIndex(), sibling_index);
    }
    // sibling_subtree->neighbor = parent_subtree->neighbor;
    sibling_subtree.setNeighborIndex(TOP, grandparent_index);
    
    // update the length of the branch connecting grandparent to sibling
    if (sibling_subtree.getUpperLength() > 0)// sibling_subtree->length > 0)
    {
        if (parent_subtree.getUpperLength() > 0)
            sibling_subtree.setUpperLength(sibling_subtree.getUpperLength() + parent_subtree.getUpperLength()); // sibling_subtree->length += parent_subtree->length;
    }
    else
        sibling_subtree.setUpperLength(parent_subtree.getUpperLength()); // sibling_subtree->length = parent_subtree->length;
    
    // update likelihood lists after subtree removal
    // case when the sibling_subtree becomes the new root
    if (parent_index.getVectorIndex() == root_vector_index) //!sibling_subtree->neighbor)
    {
        // update root
        root_vector_index = sibling_index.getVectorIndex();
        
        // delete mid_branch_lh
        /*if (root->mid_branch_lh)
        {
            delete sibling_subtree->mid_branch_lh;
            sibling_subtree->mid_branch_lh = NULL;
        }*/
        sibling_subtree.setMidBranchLh(nullptr);
        
        // reset branch length (to 0) if sibling_subtree is root
        sibling_subtree.setUpperLength(0);

        // recompute the total lh regions at sibling
        //sibling_subtree->computeTotalLhAtNode(aln, model, threshold_prob, true);
        sibling_subtree.setTotalLh(std::move(sibling_subtree.getPartialLh(TOP)->computeTotalLhAtRoot(aln.num_states, model)));
        
        // traverse downwards (to childrens of the sibling) to update their lh regions
        if (sibling_subtree.isInternal()) //->next)
        {
            // update upper left/right regions
            /*Node* next_node_1 = sibling_subtree->next; // right
            Node* next_node_2 = next_node_1->next;*/ // left
            const Index neighbor_1_index = sibling_subtree.getNeighborIndex(RIGHT);
            PhyloNode& neighbor_1 = nodes[neighbor_1_index.getVectorIndex()];
            const Index neighbor_2_index = sibling_subtree.getNeighborIndex(LEFT);
            PhyloNode& neighbor_2 = nodes[neighbor_2_index.getVectorIndex()];
            
            // if (next_node_1->partial_lh) delete next_node_1->partial_lh;
            const std::unique_ptr<SeqRegions>& lower_regions_2 = neighbor_2.getPartialLh(TOP); // next_node_2->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
            // next_node_1->partial_lh = lower_reions->computeTotalLhAtRoot(num_states, model, next_node_2->length);
            sibling_subtree.setPartialLh(RIGHT, lower_regions_2->computeTotalLhAtRoot(num_states, model, neighbor_2.getUpperLength()));
            
            // if (next_node_2->partial_lh) delete next_node_2->partial_lh;
            /*lower_reions = next_node_1->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
            next_node_2->partial_lh = lower_reions->computeTotalLhAtRoot(num_states, model, next_node_1->length);*/
            const std::unique_ptr<SeqRegions>& lower_regions_1 = neighbor_1.getPartialLh(TOP);
            sibling_subtree.setPartialLh(LEFT, lower_regions_1->computeTotalLhAtRoot(num_states, model, neighbor_1.getUpperLength()));
            
            // add children to node_stack for further traversing and updating likelihood regions
            stack<Index> node_stack;
            node_stack.push(neighbor_1_index);
            node_stack.push(neighbor_2_index);
            updatePartialLh(node_stack);
        }
    }
    // case when the sibling_subtree is non-root node
    else
    {
        // update branch length from the grandparent side
        // sibling_subtree->neighbor->length = sibling_subtree->length;
        
        // traverse to parent and sibling node to update their likelihood regions due to subtree remova
        stack<Index> node_stack;
        node_stack.push(sibling_index);
        node_stack.push(grandparent_index); // sibling_subtree->neighbor);
        updatePartialLh(node_stack);
    }
    
    // replace the node and re-update the vector lists
    const std::unique_ptr<SeqRegions>& subtree_lower_regions = subtree.getPartialLh(TOP); // ->getPartialLhAtNode(aln, model, threshold_prob);
    // try to place the new sample as a descendant of a mid-branch point
    if (is_mid_branch && best_node_index.getVectorIndex() != root_vector_index)
        placeSubTreeMidBranch(best_node_index, subtree_index, subtree, subtree_lower_regions, branch_length, best_lh_diff);
    // otherwise, best lk so far is for appending directly to existing node
    else
        placeSubTreeAtNode(best_node_index, subtree_index, subtree, subtree_lower_regions, branch_length, best_lh_diff);
}

void Tree::updateRegionsPlaceSubTree(PhyloNode& subtree, PhyloNode& sibling_node, PhyloNode& internal, std::unique_ptr<SeqRegions>&& best_child_regions, const std::unique_ptr<SeqRegions>& subtree_regions, const std::unique_ptr<SeqRegions>& upper_left_right_regions, const std::unique_ptr<SeqRegions>& lower_regions, RealNumType& best_blength)
{
    // update next_node_1->partial_lh
    // replacePartialLH(next_node_1->partial_lh, best_child_regions);
    internal.setPartialLh(LEFT, std::move(best_child_regions));

    // update new_internal_node->partial_lh
    // sibling_node.getPartialLhAtNode(aln, model, params->threshold_prob)->mergeTwoLowers(new_internal_node->partial_lh, sibling_node->length, *subtree_regions, best_blength, aln, model, params->threshold_prob);
    sibling_node.getPartialLh(TOP)->mergeTwoLowers(internal.getPartialLh(TOP), sibling_node.getUpperLength(), *subtree_regions, best_blength, aln, model, params->threshold_prob);
}

void Tree::updateRegionsPlaceSubTreeAbove(PhyloNode& subtree, PhyloNode& sibling_node, PhyloNode& internal, std::unique_ptr<SeqRegions>&& best_child_regions, const std::unique_ptr<SeqRegions>& subtree_regions, const std::unique_ptr<SeqRegions>& upper_left_right_regions, const std::unique_ptr<SeqRegions>& lower_regions, RealNumType& best_length)
{
    // sibling_node->getPartialLhAtNode(aln, model, params->threshold_prob)->mergeTwoLowers(new_internal_node->partial_lh, sibling_node->length, *subtree_regions, best_length, aln, model, params->threshold_prob);
    sibling_node.getPartialLh(TOP)->mergeTwoLowers(internal.getPartialLh(TOP), sibling_node.getUpperLength(), *subtree_regions, best_length, aln, model, params->threshold_prob);
    
    if (!internal.getPartialLh(TOP)) // new_internal_node->partial_lh)
    {
        outWarning("Problem, non lower likelihood while placing subtree -> set best branch length to min length");
        best_length = min_blength;
        /*subtree->length = best_length;
        subtree->neighbor->length = best_length;*/
        subtree.setUpperLength(best_length);
        // lower_regions->mergeTwoLowers(new_internal_node->partial_lh, sibling_node->length, *subtree_regions, best_length, aln, model, params->threshold_prob);
        lower_regions->mergeTwoLowers(internal.getPartialLh(TOP), sibling_node.getUpperLength(), *subtree_regions, best_length, aln, model, params->threshold_prob);
    }
    // upper_left_right_regions->mergeUpperLower(next_node_1->partial_lh, new_internal_node->length, *lower_regions, sibling_node->length, aln, model, params->threshold_prob);
    upper_left_right_regions->mergeUpperLower(internal.getPartialLh(LEFT), internal.getUpperLength(), *lower_regions, sibling_node.getUpperLength(), aln, model, params->threshold_prob);
}

template<void (Tree::*updateRegionsSubTree)(PhyloNode&, PhyloNode&, PhyloNode&, std::unique_ptr<SeqRegions>&&, const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, RealNumType&)>
void Tree::connectSubTree2Branch(const std::unique_ptr<SeqRegions>& subtree_regions, const std::unique_ptr<SeqRegions>& lower_regions, const Index subtree_index, PhyloNode& subtree, const Index sibling_node_index, PhyloNode& sibling_node, const RealNumType top_distance, const RealNumType down_distance, RealNumType &best_blength, std::unique_ptr<SeqRegions>&& best_child_regions, const std::unique_ptr<SeqRegions>& upper_left_right_regions)
{
    const RealNumType threshold_prob = params->threshold_prob;
    ASSERT(sibling_node_index.getMiniIndex() == TOP);
    const NumSeqsType internal_vec = subtree.getNeighborIndex(TOP).getVectorIndex();
    PhyloNode& internal = nodes[internal_vec];
        
    // re-use internal nodes
    /*Node* next_node_1 = subtree->neighbor;
    Node* new_internal_node = next_node_1->getTopNode();
    Node* next_node_2 = next_node_1->getOtherNextNode();*/
    
    // NHANLT NOTES: UNNECESSARY
    // re-order next circle (not neccessary, just to make it consistent with Python code)
    /*new_internal_node->next = next_node_2;
    next_node_2->next = next_node_1;
    next_node_1->next = new_internal_node;*/
    
    // connect new_internal_node to the parent of the selected node
    /*new_internal_node->outdated = true;
    sibling_node->neighbor->neighbor = new_internal_node;
    new_internal_node->neighbor = sibling_node->neighbor;
    new_internal_node->length = top_distance;
    new_internal_node->neighbor->length = top_distance;*/
    
    internal.setOutdated(true);
    const Index parent_index = sibling_node.getNeighborIndex(TOP);
    const NumSeqsType parent_vec = parent_index.getVectorIndex();
    nodes[parent_vec].setNeighborIndex(parent_index.getMiniIndex(), Index(internal_vec, TOP));
    internal.setNeighborIndex(TOP, parent_index);
    internal.setUpperLength(top_distance);
    
    // connect the selected_node to new_internal_node (via next_node_2)
    /*sibling_node->neighbor = next_node_2;
    next_node_2->neighbor = sibling_node;
    sibling_node->length = down_distance;
    sibling_node->neighbor->length = down_distance;*/
    
    sibling_node.setNeighborIndex(TOP, Index(internal_vec, RIGHT));
    internal.setNeighborIndex(RIGHT, sibling_node_index);
    sibling_node.setUpperLength(down_distance);
    
    // subtree already connected to new_internal_node (via next_node_1)
    /*subtree->length = best_blength;
    subtree->neighbor->length = best_blength;*/
    
    subtree.setNeighborIndex(TOP, Index(internal_vec, LEFT));
    internal.setNeighborIndex(LEFT, subtree_index);
    subtree.setUpperLength(best_blength);
            
    // update all likelihood regions
    (this->*updateRegionsSubTree)(subtree, sibling_node, internal, std::move(best_child_regions), subtree_regions, upper_left_right_regions, lower_regions, best_blength);
    
    // upper_left_right_regions->mergeUpperLower(next_node_2->partial_lh, new_internal_node->length, *subtree_regions, best_blength, aln, model, threshold_prob);
    upper_left_right_regions->mergeUpperLower(internal.getPartialLh(RIGHT), internal.getUpperLength(), *subtree_regions, best_blength, aln, model, threshold_prob);
    
    RealNumType mid_branch_length = internal.getUpperLength() * 0.5; // new_internal_node->length * 0.5;
    // upper_left_right_regions->mergeUpperLower(new_internal_node->mid_branch_lh, mid_branch_length, *new_internal_node->partial_lh, mid_branch_length, aln, model, threshold_prob);
    upper_left_right_regions->mergeUpperLower(internal.getMidBranchLh(), mid_branch_length, *internal.getPartialLh(TOP), mid_branch_length, aln, model, threshold_prob);
    
    // new_internal_node->computeTotalLhAtNode(aln, model, threshold_prob, new_internal_node == root);
    internal.updateTotalLhAtNode(nodes[parent_vec], aln, model, threshold_prob, internal_vec == root_vector_index);
    
    if (!internal.getTotalLh()) //->total_lh || new_internal_node->total_lh->size() == 0)
        outError("Problem, None vector when re-placing sample, placing subtree at mid-branch point");
    
    //if distTop>=2*min_blengthForMidNode:
    // createFurtherMidNodes(newInternalNode,upper_left_right_regions)
    
    // NHANLT: LOGS FOR DEBUGGING
    /*if (params->debug)
    {
        cout << "2Branch " << (best_blength > 0 ? subtree.getTotalLh()->size():0)<< " " << (best_blength > 0 ? subtree.getMidBranchLh()->size():0)<< " " << subtree.getPartialLh(TOP)->size() << " " << internal.getTotalLh()->size() << " " << internal.getMidBranchLh()->size()<< " " << internal.getPartialLh(TOP)->size() << " " << internal.getPartialLh(LEFT)->size() << " " << internal.getPartialLh(RIGHT)->size() << std::endl;
         cout << std::setprecision(20) << internal.getUpperLength() << " " << internal.getCorrespondingLength(RIGHT, nodes) << " " << internal.getCorrespondingLength(LEFT, nodes) << std::endl;
    }*/

    // iteratively traverse the tree to update partials from the current node
    stack<Index> node_stack;
    node_stack.push(sibling_node_index);
    node_stack.push(subtree_index);
    node_stack.push(parent_index); // new_internal_node->neighbor);
    updatePartialLh(node_stack);
}

void Tree::placeSubTreeMidBranch(const Index selected_node_index, const Index subtree_index, PhyloNode& subtree, const std::unique_ptr<SeqRegions>& subtree_regions, const RealNumType new_branch_length, const RealNumType new_lh)
{
    PhyloNode& selected_node = nodes[selected_node_index.getVectorIndex()];
    const RealNumType threshold_prob = params->threshold_prob;
    const std::unique_ptr<SeqRegions>& upper_left_right_regions = getPartialLhAtNode(selected_node.getNeighborIndex(TOP)); // selected_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
    // RealNumType best_split = 0.5;
    RealNumType best_blength_split = selected_node.getUpperLength() * 0.5;
    RealNumType best_split_lh = new_lh;
    // RealNumType new_split = 0.25;
    std::unique_ptr<SeqRegions> best_child_regions = nullptr; // std::make_unique<SeqRegions>(std::move(SeqRegions(selected_node.getMidBranchLh())));
    const std::unique_ptr<SeqRegions>& lower_regions = selected_node.getPartialLh(TOP);

    // try different positions on the existing branch
    bool found_new_split = tryShorterBranch<&Tree::calculateSubTreePlacementCost>(selected_node.getUpperLength(), best_child_regions, subtree_regions, upper_left_right_regions, lower_regions, best_split_lh, best_blength_split, new_branch_length, true);
    
    if (!found_new_split)
    {
        found_new_split = tryShorterBranch<&Tree::calculateSubTreePlacementCost>(selected_node.getUpperLength(), best_child_regions, subtree_regions, upper_left_right_regions, lower_regions, best_split_lh, best_blength_split, new_branch_length, false);
        
        if (found_new_split)
            best_blength_split = selected_node.getUpperLength() - best_blength_split;
    }
    
    // Delay cloning SeqRegions
    if (!best_child_regions)
        best_child_regions = std::make_unique<SeqRegions>(std::move(SeqRegions(selected_node.getMidBranchLh())));
    
    // now try different lengths for the new branch
    RealNumType best_blength = new_branch_length;
    estimateLengthNewBranch<&Tree::calculateSubTreePlacementCost>(best_split_lh, best_child_regions, subtree_regions, best_blength, max_blength, double_min_blength, (new_branch_length <= 0));
    
    // attach subtree to the branch above the selected node
    connectSubTree2Branch<&Tree::updateRegionsPlaceSubTree>(subtree_regions, nullptr, subtree_index, subtree, selected_node_index, selected_node, best_blength_split, selected_node.getUpperLength() - best_blength_split, best_blength, std::move(best_child_regions), upper_left_right_regions);
    
    // delete best_child_regions
    /*if (best_child_regions)
        delete best_child_regions;*/
}

void Tree::connectSubTree2Root(const Index subtree_index, PhyloNode& subtree, const std::unique_ptr<SeqRegions>& subtree_regions, const std::unique_ptr<SeqRegions>& lower_regions, const Index sibling_node_index, PhyloNode& sibling_node, const RealNumType best_root_blength, const RealNumType best_length2, std::unique_ptr<SeqRegions>&& best_parent_regions)
{
    ASSERT(sibling_node_index.getMiniIndex() == TOP);
    const NumSeqsType new_root_vec = subtree.getNeighborIndex(TOP).getVectorIndex();
    PhyloNode& new_root = nodes[new_root_vec];
    
    // re-use internal nodes
    /*Node* next_node_1 = subtree->neighbor;
    Node* new_root = next_node_1->getTopNode();
    Node* next_node_2 = next_node_1->getOtherNextNode();
    
    // NHANLT NOTES: UNNECESSARY
    // re-order next circle (not neccessary, just to make it consistent with Python code)
    new_root->next = next_node_2;
    next_node_2->next = next_node_1;
    next_node_1->next = new_root;*/
    
    // connect new_internal_node to the parent of the selected node
    /*new_root->outdated = true;
    new_root->neighbor = sibling_node->neighbor; // actually NULL since selected_node is root
    new_root->length = 0;*/
    
    new_root.setOutdated(true);
    new_root.setNeighborIndex(TOP, Index());
    new_root.setUpperLength(0);
    
    // connect the selected_node to new_internal_node (via next_node_2)
    /*sibling_node->neighbor = next_node_2;
    next_node_2->neighbor = sibling_node;
    sibling_node->length = best_root_blength;
    sibling_node->neighbor->length = best_root_blength;*/
    
    sibling_node.setNeighborIndex(TOP, Index(new_root_vec, RIGHT));
    new_root.setNeighborIndex(RIGHT, sibling_node_index);
    sibling_node.setUpperLength(best_root_blength);
    
    if (best_root_blength <= 0)
    {
        /*delete sibling_node->total_lh;
        sibling_node->total_lh = NULL;
        
        if (sibling_node->mid_branch_lh) delete sibling_node->mid_branch_lh;
        sibling_node->mid_branch_lh = NULL;*/
        
        sibling_node.setTotalLh(nullptr);
        sibling_node.setMidBranchLh(nullptr);
        
        //selected_node.furtherMidNodes=None
    }
    
    // subtree already connected to new_internal_node (via next_node_1)
    /*subtree->length = best_length2;
    subtree->neighbor->length = best_length2;*/
    
    subtree.setNeighborIndex(TOP, Index(new_root_vec, LEFT));
    new_root.setNeighborIndex(LEFT, subtree_index);
    subtree.setUpperLength(best_length2);
            
    // update all likelihood regions
    /*replacePartialLH(new_root->partial_lh, best_parent_regions);
    if (new_root->mid_branch_lh) delete new_root->mid_branch_lh;
    new_root->mid_branch_lh = NULL;*/
    
    new_root.setPartialLh(TOP, std::move(best_parent_regions));
    new_root.setMidBranchLh(nullptr);
    
    //new_root->computeTotalLhAtNode(aln, model, params->threshold_prob, true);
    new_root.setTotalLh(std::move(new_root.getPartialLh(TOP)->computeTotalLhAtRoot(aln.num_states, model)));

    /*if (next_node_1->partial_lh) delete next_node_1->partial_lh;
    next_node_1->partial_lh = lower_regions->computeTotalLhAtRoot(aln.num_states, model, best_root_blength);*/
    new_root.setPartialLh(LEFT, lower_regions->computeTotalLhAtRoot(aln.num_states, model, best_root_blength));
    
    /*if (next_node_2->partial_lh) delete next_node_2->partial_lh;
    next_node_2->partial_lh = subtree_regions->computeTotalLhAtRoot(aln.num_states, model, best_length2);*/
    new_root.setPartialLh(RIGHT, subtree_regions->computeTotalLhAtRoot(aln.num_states, model, best_length2));
    
    if (!new_root.getTotalLh()) // ->total_lh || new_root->total_lh->size() == 0)
        outWarning("Problem, None vector when re-placing sample, position root");
    
    // NHANLT: LOGS FOR DEBUGGING
    /*if (params->debug)
    {
        cout << "2Root " << (best_length2 > 0 ? subtree.getTotalLh()->size():0)<< " " << (best_length2 > 0 ? subtree.getMidBranchLh()->size():0)<< " " << subtree.getPartialLh(TOP)->size()<< " " << new_root.getTotalLh()->size()<< " " << new_root.getPartialLh(TOP)->size() << " " << new_root.getPartialLh(LEFT)->size() << " " << new_root.getPartialLh(RIGHT)->size() << std::endl;
        cout << std::setprecision(20) << new_root.getCorrespondingLength(RIGHT, nodes) << " " << new_root.getCorrespondingLength(LEFT, nodes) << std::endl;
    }*/
    
    // update tree->root;
    root_vector_index = new_root_vec;
    
    // iteratively traverse the tree to update partials from the current node
    stack<Index> node_stack;
    node_stack.push(sibling_node_index);
    node_stack.push(subtree_index);
    updatePartialLh(node_stack);
}

void Tree::handlePolytomyPlaceSubTree(const Index selected_node_index, PhyloNode& selected_node, const std::unique_ptr<SeqRegions>& subtree_regions, const RealNumType new_branch_length, RealNumType& best_down_lh_diff, Index& best_child_index, RealNumType& best_child_blength_split, std::unique_ptr<SeqRegions>& best_child_regions)
{
    // current node might be part of a polytomy (represented by 0 branch lengths) so we want to explore all the children of the current node to find out if the best placement is actually in any of the branches below the current node.
    const RealNumType threshold_prob = params->threshold_prob;
    stack<Index> new_node_stack;
    /*for (Index neighbor_index:nodes[selected_node_index.getVectorIndex()].getNeighborIndexes(TOP))
        new_node_stack.push(neighbor_index);*/
    if (selected_node.isInternal())
    {
        new_node_stack.push(selected_node.getNeighborIndex(RIGHT));
        new_node_stack.push(selected_node.getNeighborIndex(LEFT));
    }
    
    while (!new_node_stack.empty())
    {
        const Index node_index = new_node_stack.top();
        new_node_stack.pop();
        PhyloNode& node = nodes[node_index.getVectorIndex()];


        // add all nodes in polytomy
        if (node.getUpperLength() <= 0) //node->length <= 0)
        {
            /*FOR_NEIGHBOR(node, neighbor_node)
                new_node_stack.push(neighbor_node);*/
            /*for (Index neighbor_index:node.getNeighborIndexes(TOP))
                new_node_stack.push(neighbor_index);*/
            if (node.isInternal())
            {
                new_node_stack.push(node.getNeighborIndex(RIGHT));
                new_node_stack.push(node.getNeighborIndex(LEFT));
            }
        }
        else
        {
            // now try to place on the current branch below the best node, at an height above or equal to the mid-branch.
            std::unique_ptr<SeqRegions> mid_branch_regions = nullptr; // std::make_unique<SeqRegions>(std::move(SeqRegions(node.getMidBranchLh()))); // new SeqRegions(node->mid_branch_lh);
            const std::unique_ptr<SeqRegions>& parent_upper_lr_regions = getPartialLhAtNode(node.getNeighborIndex(TOP)); // node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
            const std::unique_ptr<SeqRegions>& lower_regions = node.getPartialLh(TOP); // node->getPartialLhAtNode(aln, model, threshold_prob);
            RealNumType new_branch_length_split = 0.5 * node.getUpperLength(); // node->length;
            RealNumType tmp_lh_diff = calculateSubTreePlacementCost(node.getMidBranchLh(), subtree_regions, new_branch_length);;
            
            while (true)
            {
                // if better placement found -> record it
                if (tmp_lh_diff > best_down_lh_diff)
                {
                    best_down_lh_diff = tmp_lh_diff;
                    best_child_index = node_index;
                    best_child_blength_split = new_branch_length_split;
                    new_branch_length_split *= 0.5;
                    
                    // replacePartialLH(best_child_regions, mid_branch_regions);
                    // Delay cloning SeqRegions
                    best_child_regions = mid_branch_regions ? std::move(mid_branch_regions) : std::make_unique<SeqRegions>(std::move(SeqRegions(node.getMidBranchLh())));
                    
                    if (new_branch_length_split <= half_min_blength_mid)
                        break;
                    
                    // compute mid_branch_regions
                    parent_upper_lr_regions->mergeUpperLower(mid_branch_regions, new_branch_length_split, *lower_regions, node.getUpperLength() - new_branch_length_split, aln, model, threshold_prob);
                    
                    tmp_lh_diff = calculateSubTreePlacementCost(mid_branch_regions, subtree_regions, new_branch_length);
                }
                else
                    break;
            }
            
            // delete mid_branch_regions
            // if (mid_branch_regions) delete mid_branch_regions;
        }
    }
}

void Tree::placeSubTreeAtNode(const Index selected_node_index, const Index subtree_index, PhyloNode& subtree, const std::unique_ptr<SeqRegions>& subtree_regions, const RealNumType new_branch_length, const RealNumType new_lh)
{
    // dummy variables
    const StateType num_states = aln.num_states;
    const RealNumType threshold_prob = params->threshold_prob;
    RealNumType best_child_lh;
    RealNumType best_child_blength_split = -1;
    RealNumType best_parent_lh;
    RealNumType best_parent_blength_split = 0;
    std::unique_ptr<SeqRegions> best_parent_regions = nullptr;
    RealNumType best_root_blength = -1;
    std::unique_ptr<SeqRegions> best_child_regions = nullptr;
    RealNumType best_down_lh_diff = MIN_NEGATIVE;
    Index best_child_index;
    const NumSeqsType selected_node_vec = selected_node_index.getVectorIndex();
    PhyloNode& selected_node = nodes[selected_node_vec];
    
    // We first explore placement just below the best placement node for more fine-grained placement within its descendant branches (accounting for polytomies).
    handlePolytomyPlaceSubTree(selected_node_index, selected_node, subtree_regions, new_branch_length, best_down_lh_diff, best_child_index, best_child_blength_split, best_child_regions);
    
    // place the new sample as a descendant of an existing node
    if (best_child_index.getMiniIndex() != UNDEFINED)
    {
        PhyloNode& best_child = nodes[best_child_index.getVectorIndex()];
        const std::unique_ptr<SeqRegions>& upper_left_right_regions = getPartialLhAtNode(best_child.getNeighborIndex(TOP)); // best_child->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
        const std::unique_ptr<SeqRegions>& lower_regions = best_child.getPartialLh(TOP); // ->getPartialLhAtNode(aln, model, threshold_prob);
        best_child_lh = best_down_lh_diff;
        best_child_blength_split = (best_child_blength_split == -1) ? (0.5 * best_child.getUpperLength()) : best_child_blength_split;
        
        // try with a shorter branch length
        tryShorterBranch<&Tree::calculateSubTreePlacementCost>(best_child.getUpperLength(), best_child_regions, subtree_regions, upper_left_right_regions, lower_regions, best_child_lh, best_child_blength_split, new_branch_length, true);
    }
    else
        best_child_lh = MIN_NEGATIVE;
    
    // if node is root, try to place as sibling of the current root.
    RealNumType old_root_lh = MIN_NEGATIVE;
    if (root_vector_index == selected_node_vec)
    {
        const std::unique_ptr<SeqRegions>& lower_regions = selected_node.getPartialLh(TOP); // ->getPartialLhAtNode(aln, model, threshold_prob);
        old_root_lh = lower_regions->computeAbsoluteLhAtRoot(num_states, model);
        
        // merge 2 lower vector into one
        best_parent_lh = lower_regions->mergeTwoLowers(best_parent_regions, default_blength, *subtree_regions, new_branch_length, aln, model, threshold_prob, true);
        
        best_parent_lh += best_parent_regions->computeAbsoluteLhAtRoot(num_states, model);
        
        // Try shorter branch lengths at root
        best_root_blength = default_blength;
        tryShorterBranchAtRoot(subtree_regions, lower_regions, best_parent_regions, best_root_blength, best_parent_lh, new_branch_length);
        
        // update best_parent_lh (taking into account old_root_lh)
        best_parent_lh -= old_root_lh;
    }
    // selected_node is not root
    // try to append just above node
    else
    {
        const std::unique_ptr<SeqRegions>& upper_left_right_regions = getPartialLhAtNode(selected_node.getNeighborIndex(TOP)); // selected_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
        const std::unique_ptr<SeqRegions>& lower_regions = selected_node.getPartialLh(TOP); // selected_node->getPartialLhAtNode(aln, model, threshold_prob);
        best_parent_regions = nullptr; // std::make_unique<SeqRegions>(std::move(SeqRegions(selected_node.getMidBranchLh())));
        best_parent_lh = calculateSubTreePlacementCost(selected_node.getMidBranchLh(), subtree_regions, new_branch_length);
        best_parent_blength_split = 0.5 * selected_node.getUpperLength();
        
        // try with a shorter split
        tryShorterBranch<&Tree::calculateSubTreePlacementCost>(selected_node.getUpperLength(), best_parent_regions, subtree_regions, upper_left_right_regions, lower_regions, best_parent_lh, best_parent_blength_split, new_branch_length, false);
        
        // Delay cloning SeqRegions
        if (!best_parent_regions)
            best_parent_regions = std::make_unique<SeqRegions>(std::move(SeqRegions(selected_node.getMidBranchLh())));
    }
    
    // if the best placement is below the selected_node => add an internal node below the selected_node
    // now we have three likelihood costs,
    // best_child_lh is the likelihood score of appending below node;
    // best_parent_lh is the likelihood score of appending above node;
    // new_lh is the likelihood cost of appending exactly at node.
    if (best_child_lh >= best_parent_lh && best_child_lh >= new_lh)
    {
        ASSERT(best_child_index.getMiniIndex() != UNDEFINED);
        PhyloNode& best_child = nodes[best_child_index.getVectorIndex()];
        const std::unique_ptr<SeqRegions>& upper_left_right_regions = getPartialLhAtNode(best_child.getNeighborIndex(TOP)); // best_child->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
        
        // now try different lengths for the new branch
        RealNumType best_length = new_branch_length;
        estimateLengthNewBranch<&Tree::calculateSubTreePlacementCost>(best_child_lh, best_child_regions, subtree_regions, best_length, max_blength, double_min_blength, (new_branch_length <= 0));
        
        // attach subtree to the phylogenetic tree (below the selected_node ~ above the child node)
        connectSubTree2Branch<&Tree::updateRegionsPlaceSubTree>(subtree_regions, nullptr, subtree_index, subtree, best_child_index, best_child, best_child_blength_split, best_child.getUpperLength() - best_child_blength_split, best_length, std::move(best_child_regions), upper_left_right_regions);

    }
    // otherwise, add new parent to the selected_node
    else
    {
        const std::unique_ptr<SeqRegions>& lower_regions = selected_node.getPartialLh(TOP); // ->getPartialLhAtNode(aln, model, threshold_prob);
        
        // new parent is actually part of a polytomy since best placement is exactly at the node
        if (new_lh >= best_parent_lh)
        {
            best_root_blength = -1;
            best_parent_blength_split = -1;
            best_parent_lh = new_lh;
            /*if (best_parent_regions) delete best_parent_regions;
            best_parent_regions = NULL;*/
            best_parent_regions = nullptr;
            
            if (selected_node_vec == root_vector_index)
                lower_regions->mergeTwoLowers(best_parent_regions, -1, *subtree_regions, new_branch_length, aln, model, threshold_prob);
            else
                best_parent_regions = std::make_unique<SeqRegions>(std::move(selected_node.getTotalLh()));
        }

        // add parent to the root
        if (selected_node_vec == root_vector_index)
        {
            // remove old_root_lh from best_parent_lh before estimating the new branch length
            best_parent_lh += old_root_lh;
            
            // estimate the new branch length
            RealNumType best_length2 = new_branch_length;
            estimateLengthNewBranchAtRoot(subtree_regions, lower_regions, best_parent_regions, best_length2, best_parent_lh, best_root_blength, double_min_blength, new_branch_length <= 0);
            
            // update best_parent_lh (taking into account old_root_lh)
            best_parent_lh -= old_root_lh;
            
            // attach subtree to the phylogenetic tree (exactly at the seleted root node)
            connectSubTree2Root(subtree_index, subtree, subtree_regions, lower_regions, selected_node_index, selected_node, best_root_blength, best_length2, std::move(best_parent_regions));
        }
        //add parent to non-root node (place subtree exactly at the selected non-root node)
        else
        {
            const std::unique_ptr<SeqRegions>& upper_left_right_regions = getPartialLhAtNode(selected_node.getNeighborIndex(TOP)); // selected_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
            
            // estimate the length for the new branch
            RealNumType best_length = new_branch_length;
            estimateLengthNewBranch<&Tree::calculateSubTreePlacementCost>(best_parent_lh, best_parent_regions, subtree_regions, best_length, new_branch_length * 10, double_min_blength, (new_branch_length <= 0));
            
            // attach subtree to the phylogenetic tree (exactly at the selected non-root node)
            RealNumType down_distance = best_parent_blength_split;
            RealNumType top_distance = selected_node.getUpperLength() - down_distance; // selected_node->length - down_distance;
            if (best_parent_blength_split <= 0)
            {
                down_distance = -1;
                top_distance = selected_node.getUpperLength(); // ->length;
                
                /*if (selected_node->total_lh) delete selected_node->total_lh;
                selected_node->total_lh = NULL;*/
                selected_node.setTotalLh(nullptr);
                
                /*if (selected_node->mid_branch_lh) delete selected_node->mid_branch_lh;
                selected_node->mid_branch_lh = NULL;*/
                selected_node.setMidBranchLh(nullptr);
            }
            
            connectSubTree2Branch<&Tree::updateRegionsPlaceSubTreeAbove>(subtree_regions, lower_regions, subtree_index, subtree, selected_node_index, selected_node, top_distance, down_distance, best_length, std::move(best_child_regions), upper_left_right_regions);
        }
    }
    
    // delete best_parent_regions and best_child_regions
    /*if (best_parent_regions)
        delete best_parent_regions;
    if (best_child_regions)
        delete best_child_regions;*/
}

template <RealNumType(Tree::*calculatePlacementCost)(const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, const RealNumType)>
bool Tree::tryShorterBranch(const RealNumType current_blength, std::unique_ptr<SeqRegions>& best_child_regions, const std::unique_ptr<SeqRegions>& sample, const std::unique_ptr<SeqRegions>& upper_left_right_regions, const std::unique_ptr<SeqRegions>& lower_regions, RealNumType &best_split_lh, RealNumType &best_branch_length_split, const RealNumType new_branch_length, const bool try_first_branch)
{
    std::unique_ptr<SeqRegions> new_parent_regions = nullptr;
    bool found_new_split = false;
    RealNumType new_branch_length_split = 0.5 * best_branch_length_split;
    
    while (new_branch_length_split > min_blength)
    {
        // try on the first or second branch
        if (try_first_branch)
            upper_left_right_regions->mergeUpperLower(new_parent_regions, new_branch_length_split, *lower_regions,  current_blength - new_branch_length_split, aln, model, params->threshold_prob);
        else
            upper_left_right_regions->mergeUpperLower(new_parent_regions, current_blength - new_branch_length_split, *lower_regions, new_branch_length_split, aln, model, params->threshold_prob);
        
        // calculate placement_cost
        RealNumType placement_cost = (this->*calculatePlacementCost)(new_parent_regions, sample, new_branch_length);
        
        if (placement_cost > best_split_lh)
        {
            best_split_lh = placement_cost;
            best_branch_length_split = new_branch_length_split;
            new_branch_length_split *= 0.5;
            found_new_split = true;
            
            best_child_regions = std::move(new_parent_regions);
        }
        else
            break;
    }
    
    return found_new_split;
}

template <RealNumType(Tree::*calculatePlacementCost)(const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, const RealNumType)>
bool Tree::tryShorterNewBranch(const std::unique_ptr<SeqRegions>& best_child_regions, const std::unique_ptr<SeqRegions>& sample, RealNumType &best_blength, RealNumType &new_branch_lh, const RealNumType short_blength_thresh)
{
    bool found_new_split = false;
    RealNumType new_blength = best_blength;
    RealNumType placement_cost;
    
    while (best_blength > short_blength_thresh)
    {
        new_blength *= 0.5;
        placement_cost = (this->*calculatePlacementCost)(best_child_regions, sample, new_blength);
        
        if (placement_cost > new_branch_lh)
        {
            new_branch_lh = placement_cost;
            best_blength = new_blength;
            found_new_split = true;
        }
        else
            break;
    }
    
    return found_new_split;
}

template <RealNumType(Tree::*calculatePlacementCost)(const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, const RealNumType)>
void Tree::tryLongerNewBranch(const std::unique_ptr<SeqRegions>& best_child_regions, const std::unique_ptr<SeqRegions>& sample, RealNumType &best_blength, RealNumType &new_branch_lh, const RealNumType long_blength_thresh)
{
    RealNumType new_blength = best_blength;
    RealNumType placement_cost;
    
    while (best_blength < long_blength_thresh)
    {
        new_blength += new_blength;
        placement_cost = (this->*calculatePlacementCost)(best_child_regions, sample, new_blength);
        if (placement_cost > new_branch_lh)
        {
            new_branch_lh = placement_cost;
            best_blength = new_blength;
        }
        else
            break;
    }
}

template <RealNumType(Tree::*calculatePlacementCost)(const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, const RealNumType)>
void Tree::estimateLengthNewBranch(const RealNumType best_split_lh, const std::unique_ptr<SeqRegions>& best_child_regions, const std::unique_ptr<SeqRegions>& sample, RealNumType &best_blength, const RealNumType long_blength_thresh, const RealNumType short_blength_thresh, const bool optional_check)
{
    RealNumType new_branch_lh = best_split_lh;
    
    // change zero branch length to min branch length
    if (optional_check)
    {
        best_blength = min_blength;
        new_branch_lh = (this->*calculatePlacementCost)(best_child_regions, sample, best_blength);
    }
    
    // try shorter lengths for the new branch
    bool found_new_blength = tryShorterNewBranch<calculatePlacementCost>(best_child_regions, sample, best_blength, new_branch_lh, min_blength);
    
    // try longer lengths for the new branch
    if (optional_check || !found_new_blength) // (best_blength > 0.7 * default_blength)
        tryLongerNewBranch<calculatePlacementCost>(best_child_regions, sample, best_blength, new_branch_lh, long_blength_thresh);
    
    // try zero-length for the new branch
    if (best_blength < short_blength_thresh)
    {
        RealNumType zero_branch_lh = (this->*calculatePlacementCost)(best_child_regions, sample, -1);
        if (zero_branch_lh > new_branch_lh)
            best_blength = -1;
    }
}

void Tree::connectNewSample2Branch(std::unique_ptr<SeqRegions>& sample, const NumSeqsType seq_name_index, const Index sibling_node_index, PhyloNode& sibling_node, const RealNumType top_distance, const RealNumType down_distance, const RealNumType best_blength, std::unique_ptr<SeqRegions>& best_child_regions, const std::unique_ptr<SeqRegions>& upper_left_right_regions)
{
    const RealNumType threshold_prob = params->threshold_prob;
    
    // create new internal node and append child to it
    /*Node* new_internal_node = new Node(true);
    Node* next_node_1 = new Node();
    Node* next_node_2 = new Node();
    Node* new_sample_node = new Node(seq_name);
    
    new_internal_node->next = next_node_2;
    next_node_2->next = next_node_1;
    next_node_1->next = new_internal_node;*/
    createAnInternalNode();
    createALeafNode(seq_name_index);
    NumSeqsType leaf_vec_index = nodes.size() - 1;
    PhyloNode& leaf = nodes[leaf_vec_index];
    NumSeqsType internal_vec_index = leaf_vec_index - 1;
    PhyloNode& internal = nodes[internal_vec_index];
    
    /*new_internal_node->neighbor = sibling_node->neighbor;
     sibling_node->neighbor->neighbor = new_internal_node;
     new_internal_node->length = top_distance;
     new_internal_node->neighbor->length = top_distance;*/
    
    Index parent_index = sibling_node.getNeighborIndex(TOP);
    PhyloNode& parent_node = nodes[parent_index.getVectorIndex()];
    internal.setNeighborIndex(TOP, parent_index);
    parent_node.setNeighborIndex(parent_index.getMiniIndex(), Index(internal_vec_index, TOP));
    internal.setUpperLength(top_distance);
    
    /*sibling_node->neighbor = next_node_2;
     next_node_2->neighbor = sibling_node;
     sibling_node->length = down_distance;
     sibling_node->neighbor->length = down_distance;*/
    
    internal.setNeighborIndex(RIGHT, sibling_node_index);
    sibling_node.setNeighborIndex(TOP, Index(internal_vec_index, RIGHT));
    sibling_node.setUpperLength(down_distance);
    
    /*new_sample_node->neighbor = next_node_1;
     next_node_1->neighbor = new_sample_node;
     new_sample_node->length = best_blength;
     new_sample_node->neighbor->length = best_blength;*/
    
    internal.setNeighborIndex(LEFT, Index(leaf_vec_index, TOP));
    leaf.setNeighborIndex(TOP, Index(internal_vec_index, LEFT));
    leaf.setUpperLength(best_blength);
    
    /*new_sample_node->partial_lh = sample;
    next_node_1->partial_lh = best_child_regions;
    best_child_regions = NULL;
    upper_left_right_regions->mergeUpperLower(next_node_2->partial_lh, new_internal_node->length, *sample, best_blength, aln, model, threshold_prob);
    sibling_node->getPartialLhAtNode(aln, model, threshold_prob)->mergeTwoLowers(new_internal_node->partial_lh, sibling_node->length, *sample, best_blength, aln, model, threshold_prob);
    RealNumType half_branch_length = new_internal_node->length * 0.5;
    upper_left_right_regions->mergeUpperLower(new_internal_node->mid_branch_lh, half_branch_length, *new_internal_node->partial_lh, half_branch_length, aln, model, threshold_prob);*/
    
    leaf.setPartialLh(TOP, std::move(sample));
    SeqRegions& leaf_lower_regions = *leaf.getPartialLh(TOP);
    internal.setPartialLh(LEFT, std::move(best_child_regions));
    upper_left_right_regions->mergeUpperLower(internal.getPartialLh(RIGHT), internal.getUpperLength(), leaf_lower_regions, best_blength, aln, model, threshold_prob);
    sibling_node.getPartialLh(TOP)->mergeTwoLowers(internal.getPartialLh(TOP), sibling_node.getUpperLength(), leaf_lower_regions, best_blength, aln, model, threshold_prob);
    RealNumType half_branch_length = internal.getUpperLength() * 0.5;
    upper_left_right_regions->mergeUpperLower(internal.getMidBranchLh(), half_branch_length, *(internal.getPartialLh(TOP)), half_branch_length, aln, model, threshold_prob);
    
    //new_internal_node->computeTotalLhAtNode(aln, model, threshold_prob, new_internal_node == root);
    internal.updateTotalLhAtNode(parent_node, aln, model, threshold_prob, internal_vec_index == root_vector_index);
    
    //if (!internal.getTotalLh() || internal.getTotalLh()->empty())
    if (!internal.getTotalLh())
        outError("Problem, None vector when placing sample, below node");
    
    if (best_blength > 0)
    {
        // new_sample_node->computeTotalLhAtNode(aln, model, threshold_prob, new_sample_node == root);
        leaf.updateTotalLhAtNode(internal, aln, model, threshold_prob, leaf_vec_index == root_vector_index);
        
        /*RealNumType half_branch_length = new_sample_node->length * 0.5;
        next_node_1->getPartialLhAtNode(aln, model, threshold_prob)->mergeUpperLower(new_sample_node->mid_branch_lh, half_branch_length, *sample, half_branch_length, aln, model, threshold_prob);*/
        RealNumType half_branch_length = leaf.getUpperLength() * 0.5;
        internal.getPartialLh(LEFT)->mergeUpperLower(leaf.getMidBranchLh(), half_branch_length, leaf_lower_regions, half_branch_length, aln, model, threshold_prob);
    }
    
    // NHANLT: LOGS FOR DEBUGGING
    /*if (params->debug)
    {
        cout << "2Branch " << aln.data[seq_name_index].seq_name << " " << (best_blength > 0 ? leaf.getTotalLh()->size():0)<< " " << (best_blength > 0 ? leaf.getMidBranchLh()->size():0)<< " " << leaf.getPartialLh(TOP)->size() << " " << internal.getTotalLh()->size() << " " << internal.getMidBranchLh()->size()<< " " << internal.getPartialLh(TOP)->size() << " " << internal.getPartialLh(LEFT)->size() << " " << internal.getPartialLh(RIGHT)->size() << std::endl;
        cout << std::setprecision(20) << internal.getUpperLength() << " " << internal.getCorrespondingLength(RIGHT, nodes) << " " << internal.getCorrespondingLength(LEFT, nodes) << std::endl;
    }*/
    
    // update pseudo_count
    model.updatePesudoCount(aln, *internal.getPartialLh(LEFT), *leaf.getPartialLh(TOP));

    // iteratively traverse the tree to update partials from the current node
    stack<Index> node_stack;
    node_stack.push(sibling_node_index);
    node_stack.push(parent_index);
    updatePartialLh(node_stack);
}

void Tree::placeNewSampleMidBranch(const Index& selected_node_index, std::unique_ptr<SeqRegions>& sample, const NumSeqsType seq_name_index, const RealNumType best_lh_diff)
{
    // dummy variables
    // const RealNumType threshold_prob = params->threshold_prob;
    std::unique_ptr<SeqRegions> best_child_regions = nullptr;
    
    // selected_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
    // const MiniIndex seleted_node_mini_index = selected_node_index.getMiniIndex();
    ASSERT(selected_node_index.getMiniIndex() == TOP);
    PhyloNode& selected_node = nodes[selected_node_index.getVectorIndex()];
    const std::unique_ptr<SeqRegions>& upper_left_right_regions = getPartialLhAtNode(selected_node.getNeighborIndex(TOP));
    RealNumType best_split_lh = best_lh_diff;
    const RealNumType selected_node_blength = selected_node.getUpperLength();
    RealNumType best_branch_length_split = 0.5 * selected_node_blength;
    best_child_regions = nullptr; // std::make_unique<SeqRegions>(std::move(SeqRegions(selected_node.getMidBranchLh())));
    //selected_node->getPartialLhAtNode(aln, model, threshold_prob);
    const std::unique_ptr<SeqRegions>& lower_regions = selected_node.getPartialLh(TOP);
    
    // try different positions on the existing branch
    bool found_new_split = tryShorterBranch<&Tree::calculateSamplePlacementCost>(selected_node_blength, best_child_regions, sample, upper_left_right_regions, lower_regions, best_split_lh, best_branch_length_split, default_blength, true);
    
    if (!found_new_split)
    {
        // try on the second branch
        found_new_split = tryShorterBranch<&Tree::calculateSamplePlacementCost>(selected_node_blength, best_child_regions, sample, upper_left_right_regions, lower_regions, best_split_lh, best_branch_length_split, default_blength, false);
        
        if (found_new_split)
            best_branch_length_split = selected_node_blength - best_branch_length_split;
    }
    
    // Delay cloning SeqRegions
    if (!best_child_regions)
        best_child_regions = std::make_unique<SeqRegions>(std::move(SeqRegions(selected_node.getMidBranchLh())));
    
    // now try different lengths for the new branch
    RealNumType best_blength = default_blength;
    estimateLengthNewBranch<&Tree::calculateSamplePlacementCost>(best_split_lh, best_child_regions, sample, best_blength, max_blength, min_blength, false);
    
    // create new internal node and append child to it
    connectNewSample2Branch(sample, seq_name_index, selected_node_index, selected_node, best_branch_length_split, selected_node_blength - best_branch_length_split, best_blength, best_child_regions, upper_left_right_regions);
}

void Tree::tryShorterBranchAtRoot(const std::unique_ptr<SeqRegions>& sample, const std::unique_ptr<SeqRegions>& lower_regions, std::unique_ptr<SeqRegions>& best_parent_regions, RealNumType &best_root_blength, RealNumType &best_parent_lh, const RealNumType fixed_blength)
{
    std::unique_ptr<SeqRegions> merged_root_sample_regions = nullptr;
    RealNumType new_blength = 0.5 * best_root_blength;
    RealNumType new_root_lh;
    
    while (new_blength > min_blength)
    {
        // merge 2 lower vector into one
        new_root_lh = lower_regions->mergeTwoLowers(merged_root_sample_regions, new_blength, *sample, fixed_blength, aln, model, params->threshold_prob, true);
        new_root_lh += merged_root_sample_regions->computeAbsoluteLhAtRoot(aln.num_states, model);
        
        if (new_root_lh > best_parent_lh)
        {
            best_parent_lh = new_root_lh;
            best_root_blength = new_blength;
            new_blength *= 0.5;
            
            // replacePartialLH(best_parent_regions, merged_root_sample_regions);
            best_parent_regions = std::move(merged_root_sample_regions);
        }
        else
            break;
    }
    
    // delete merged_root_sample_regions
    // if (merged_root_sample_regions) delete merged_root_sample_regions;
}

bool Tree::tryShorterNewBranchAtRoot(const std::unique_ptr<SeqRegions>& sample, const std::unique_ptr<SeqRegions>& lower_regions, std::unique_ptr<SeqRegions>&best_parent_regions, RealNumType &best_length, RealNumType &best_parent_lh, const RealNumType fixed_blength)
{
    std::unique_ptr<SeqRegions> new_root_lower_regions = nullptr;
    bool found_new_split = false;
    RealNumType new_blength = best_length;
    RealNumType new_root_lh = 0;
    
    while (best_length > min_blength)
    {
        new_blength *= 0.5;
        
        new_root_lh = lower_regions->mergeTwoLowers(new_root_lower_regions, fixed_blength, *sample, new_blength, aln, model, params->threshold_prob, true);
        new_root_lh += new_root_lower_regions->computeAbsoluteLhAtRoot(aln.num_states, model);
        
        if (new_root_lh > best_parent_lh)
        {
            best_parent_lh = new_root_lh;
            best_length = new_blength;
            found_new_split = true;
            
            // replacePartialLH(best_parent_regions, new_root_lower_regions);
            best_parent_regions = std::move(new_root_lower_regions);
        }
        else
            break;
    }
    
    // delete new_root_lower_regions
    // if (new_root_lower_regions) delete new_root_lower_regions;
    
    return found_new_split;
}

bool Tree::tryLongerNewBranchAtRoot(const std::unique_ptr<SeqRegions>& sample, const std::unique_ptr<SeqRegions>& lower_regions, std::unique_ptr<SeqRegions>& best_parent_regions, RealNumType &best_length, RealNumType &best_parent_lh, const RealNumType fixed_blength)
{
    std::unique_ptr<SeqRegions> new_root_lower_regions = nullptr;
    bool found_new_split = false;
    RealNumType new_blength = best_length;
    RealNumType new_root_lh = 0;
    
    while (best_length < max_blength)
    {
        new_blength += new_blength;
        
        new_root_lh = lower_regions->mergeTwoLowers(new_root_lower_regions, fixed_blength, *sample, new_blength, aln, model, params->threshold_prob, true);
        new_root_lh += new_root_lower_regions->computeAbsoluteLhAtRoot(aln.num_states, model);
        
        if (new_root_lh > best_parent_lh)
        {
            best_parent_lh = new_root_lh;
            best_length = new_blength;
            found_new_split = true;
            
            // replacePartialLH(best_parent_regions, new_root_lower_regions);
            best_parent_regions = std::move(new_root_lower_regions);
        }
        else
            break;
    }
    
    // delete new_root_lower_regions
    // if (new_root_lower_regions) delete new_root_lower_regions;
    
    return found_new_split;
}

void Tree::estimateLengthNewBranchAtRoot(const std::unique_ptr<SeqRegions>& sample, const std::unique_ptr<SeqRegions>& lower_regions, std::unique_ptr<SeqRegions>& best_parent_regions, RealNumType &best_length, RealNumType &best_parent_lh, const RealNumType fixed_blength, const RealNumType short_blength_thresh, const bool optional_check)
{
    if (optional_check)
    {
        std::unique_ptr<SeqRegions> new_root_lower_regions = nullptr;
        best_length = min_blength;
        best_parent_lh = lower_regions->mergeTwoLowers(new_root_lower_regions, fixed_blength, *sample, best_length, aln, model, params->threshold_prob, true);
        
        best_parent_lh += new_root_lower_regions->computeAbsoluteLhAtRoot(aln.num_states, model);
        
        // replacePartialLH(best_parent_regions, new_root_lower_regions);
        best_parent_regions = std::move(new_root_lower_regions);
    }
    
    // try shorter lengths
    bool found_new_split = tryShorterNewBranchAtRoot(sample, lower_regions, best_parent_regions, best_length, best_parent_lh, fixed_blength);
    
    // try longer lengths
    if (optional_check || !found_new_split)
        tryLongerNewBranchAtRoot(sample, lower_regions, best_parent_regions, best_length, best_parent_lh, fixed_blength);
    
    
    // try with length zero
    if (best_length < short_blength_thresh)
    {
        std::unique_ptr<SeqRegions> new_root_lower_regions = nullptr;
        
        RealNumType new_root_lh = lower_regions->mergeTwoLowers(new_root_lower_regions, fixed_blength, *sample, -1, aln, model, params->threshold_prob, true);
        new_root_lh += new_root_lower_regions->computeAbsoluteLhAtRoot(aln.num_states, model);
        
        if (new_root_lh > best_parent_lh)
        {
            best_length = -1;
            // replacePartialLH(best_parent_regions, new_root_lower_regions);
            best_parent_regions = std::move(new_root_lower_regions);
        }
        
        // delete new_root_lower_regions
        // if (new_root_lower_regions) delete new_root_lower_regions;
    }
}

void Tree::connectNewSample2Root(std::unique_ptr<SeqRegions>& sample, const NumSeqsType seq_name_index, const Index sibling_node_index, PhyloNode& sibling_node, const RealNumType best_root_blength, const RealNumType best_length2, std::unique_ptr<SeqRegions>& best_parent_regions)
{
    const RealNumType threshold_prob = params->threshold_prob;
    // const MiniIndex sibling_node_mini_index = sibling_node_index.getMiniIndex();
    
    /*Node* new_root = new Node(true);
    Node* next_node_1 = new Node();
    Node* next_node_2 = new Node();
    Node* new_sample_node = new Node(seq_name);
    
    new_root->next = next_node_2;
    next_node_2->next = next_node_1;
    next_node_1->next = new_root;*/
    
    createAnInternalNode();
    createALeafNode(seq_name_index);
    NumSeqsType leaf_vec_index = nodes.size() - 1;
    PhyloNode& leaf = nodes[leaf_vec_index];
    NumSeqsType new_root_vec_index = leaf_vec_index - 1;
    PhyloNode& new_root = nodes[new_root_vec_index];
    
    new_root.setNeighborIndex(TOP, Index());
    
    // attach the left child
    /*sibling_node->neighbor = next_node_2;
    next_node_2->neighbor = sibling_node;
    sibling_node->length = best_root_blength;
    sibling_node->neighbor->length = best_root_blength;*/
    
    new_root.setNeighborIndex(RIGHT, sibling_node_index);
    sibling_node.setNeighborIndex(TOP, Index(new_root_vec_index, RIGHT));
    sibling_node.setUpperLength(best_root_blength);
    
    if (best_root_blength <= 0)
    {
        /*if (sibling_node->total_lh) delete sibling_node->total_lh;
        sibling_node->total_lh = NULL;
        
        if (sibling_node->mid_branch_lh) delete sibling_node->mid_branch_lh;
        sibling_node->mid_branch_lh = NULL;*/
        sibling_node.setTotalLh(nullptr);
        sibling_node.setMidBranchLh(nullptr);
        //selected_node.furtherMidNodes=None
    }
    
    // attach the right child
    /*new_sample_node->neighbor = next_node_1;
    next_node_1->neighbor = new_sample_node;
    new_sample_node->length = best_length2;
    new_sample_node->neighbor->length = best_length2;*/
    
    new_root.setNeighborIndex(LEFT, Index(leaf_vec_index, TOP));
    leaf.setNeighborIndex(TOP, Index(new_root_vec_index, LEFT));
    leaf.setUpperLength(best_length2);
    
    /*new_root->partial_lh = best_parent_regions;
    best_parent_regions = NULL;*/
    new_root.setPartialLh(TOP, std::move(best_parent_regions));
    
    // new_root->total_lh = new_root->computeTotalLhAtNode(aln, model, threshold_prob, true);
    new_root.setTotalLh(std::move(new_root.getPartialLh(TOP)->computeTotalLhAtRoot(aln.num_states, model)));

    /*next_node_1->partial_lh = sibling_node->getPartialLhAtNode(aln, model, threshold_prob)->computeTotalLhAtRoot(aln.num_states, model, best_root_blength);
    next_node_2->partial_lh = sample->computeTotalLhAtRoot(aln.num_states, model, best_length2);*/
    new_root.setPartialLh(LEFT, sibling_node.getPartialLh(TOP)->computeTotalLhAtRoot(aln.num_states, model, best_root_blength));
    new_root.setPartialLh(RIGHT, sample->computeTotalLhAtRoot(aln.num_states, model, best_length2));
    
    // new_sample_node->partial_lh = sample;
    leaf.setPartialLh(TOP, std::move(sample));
    
    if (!new_root.getTotalLh()) //(!new_root.getTotalLh() || new_root->total_lh->size() == 0)
    {
        outError("Problem, None vector when placing sample, new root");
        //print(merged_root_sample_regions)
        //print(node.probVect)
        //print(sample)
        //print(best_length2)
        //print(best_root_blength)
    }
    
    if (best_length2 > 0)
    {
        //new_sample_node->computeTotalLhAtNode(aln, model, threshold_prob, new_sample_node == root);
        leaf.updateTotalLhAtNode(new_root, aln, model, threshold_prob, leaf_vec_index = root_vector_index);
        
        RealNumType half_branch_length = leaf.getUpperLength() * 0.5; // new_sample_node->length * 0.5;
        // next_node_1->getPartialLhAtNode(aln, model, threshold_prob)->mergeUpperLower(new_sample_node->mid_branch_lh, half_branch_length, *sample, half_branch_length, aln, model, threshold_prob);
        new_root.getPartialLh(LEFT)->mergeUpperLower(leaf.getMidBranchLh(), half_branch_length, *leaf.getPartialLh(TOP), half_branch_length, aln, model, threshold_prob);
        
        //if best_length2>=2*min_blengthForMidNode:
          //  createFurtherMidNodes(new_root.children[1],new_root.probVectUpLeft)
    }
    
    // NHANLT: LOGS FOR DEBUGGING
    /*if (params->debug)
    {
        cout << "2Root " << aln.data[seq_name_index].seq_name << " " << (best_length2 > 0 ? leaf.getTotalLh()->size():0)<< " " << (best_length2 > 0 ? leaf.getMidBranchLh()->size():0)<< " " << leaf.getPartialLh(TOP)->size()<< " " << new_root.getTotalLh()->size()<< " " << new_root.getPartialLh(TOP)->size() << " " << new_root.getPartialLh(LEFT)->size() << " " << new_root.getPartialLh(RIGHT)->size() << std::endl;
        cout << std::setprecision(20) << new_root.getCorrespondingLength(RIGHT, nodes) << " " << new_root.getCorrespondingLength(LEFT, nodes) << std::endl;
    }*/
    
    // update tree->root;
    root_vector_index = new_root_vec_index;
    
    // iteratively traverse the tree to update partials from the current node
    stack<Index> node_stack;
    node_stack.push(sibling_node_index);
    updatePartialLh(node_stack);
}

void Tree::placeNewSampleAtNode(const Index selected_node_index, std::unique_ptr<SeqRegions>& sample, const NumSeqsType seq_name_index, const RealNumType best_lh_diff, const RealNumType best_up_lh_diff, const RealNumType best_down_lh_diff, const Index best_child_index)
{
    // dummy variables
    RealNumType best_child_lh = MIN_NEGATIVE;
    RealNumType best_child_blength_split = 0;
    RealNumType best_parent_lh;
    RealNumType best_parent_blength_split = 0;
    RealNumType best_root_blength = -1;
    const StateType num_states = aln.num_states;
    const RealNumType threshold_prob = params->threshold_prob;
    std::unique_ptr<SeqRegions> best_parent_regions = nullptr;
    std::unique_ptr<SeqRegions> best_child_regions = nullptr;
    
    ASSERT(selected_node_index.getMiniIndex() == TOP);
    const NumSeqsType selected_node_vec_index = selected_node_index.getVectorIndex();
    PhyloNode& selected_node = nodes[selected_node_vec_index];
    
    // place the new sample as a descendant of an existing node
    if (best_child_index.getMiniIndex() != UNDEFINED)
    {
        PhyloNode& best_child = nodes[best_child_index.getVectorIndex()];
        best_child_lh = best_down_lh_diff;
        ASSERT(best_child_index.getMiniIndex() == TOP);
        best_child_blength_split = 0.5 * best_child.getUpperLength(); // best_child->length;
        const std::unique_ptr<SeqRegions>& upper_left_right_regions = getPartialLhAtNode(best_child.getNeighborIndex(TOP)); // best_child->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
        const std::unique_ptr<SeqRegions>& lower_regions = best_child.getPartialLh(TOP); // ->getPartialLhAtNode(aln, model, threshold_prob);
        // best_child_regions = new SeqRegions(best_child->mid_branch_lh);
        // SeqRegions best_child_mid_clone = SeqRegions(best_child.getMidBranchLh());
        best_child_regions = nullptr; // std::make_unique<SeqRegions>(std::move(best_child_mid_clone));
        
        // try a shorter split
        // tryShorterBranch<&Tree::calculateSamplePlacementCost>(best_child->length, best_child_regions, sample, upper_left_right_regions, lower_regions, best_child_lh, best_child_blength_split, default_blength, true);
        tryShorterBranch<&Tree::calculateSamplePlacementCost>(best_child.getUpperLength(), best_child_regions, sample, upper_left_right_regions, lower_regions, best_child_lh, best_child_blength_split, default_blength, true);
        
        // Delay cloning SeqRegions
        if (!best_child_regions)
            best_child_regions = std::make_unique<SeqRegions>(std::move(SeqRegions(best_child.getMidBranchLh())));
    }
    
    // if node is root, try to place as sibling of the current root.
    RealNumType old_root_lh = MIN_NEGATIVE;
    if (root_vector_index == selected_node_vec_index)
    {
        /*old_root_lh = selected_node->getPartialLhAtNode(aln, model, threshold_prob)->computeAbsoluteLhAtRoot(num_states, model);
        SeqRegions* lower_regions = selected_node->getPartialLhAtNode(aln, model, threshold_prob);*/
        const std::unique_ptr<SeqRegions>& lower_regions = selected_node.getPartialLh(TOP);
        old_root_lh = lower_regions->computeAbsoluteLhAtRoot(num_states, model);
        
        // merge 2 lower vector into one
        RealNumType new_root_lh = lower_regions->mergeTwoLowers(best_parent_regions, default_blength, *sample, default_blength, aln, model, threshold_prob, true);
        
        new_root_lh += best_parent_regions->computeAbsoluteLhAtRoot(num_states, model);
        best_parent_lh = new_root_lh;
        
        // try shorter branch lengths
        best_root_blength = default_blength;
        tryShorterBranchAtRoot(sample, lower_regions, best_parent_regions, best_root_blength, best_parent_lh, default_blength);
        
        // update best_parent_lh (taking into account old_root_lh)
        best_parent_lh -= old_root_lh;
    }
    // selected_node is not root
    else
    {
        best_parent_lh = best_up_lh_diff;
        best_parent_blength_split = 0.5 * selected_node.getUpperLength(); // selected_node->length;
        /*SeqRegions* upper_left_right_regions = selected_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
        SeqRegions* lower_regions = selected_node->getPartialLhAtNode(aln, model, threshold_prob);
        best_parent_regions = new SeqRegions(selected_node->mid_branch_lh);*/
        const std::unique_ptr<SeqRegions>& upper_left_right_regions = getPartialLhAtNode(selected_node.getNeighborIndex(TOP));
        const std::unique_ptr<SeqRegions>& lower_regions = selected_node.getPartialLh(TOP);
        // SeqRegions seq_regions_clone = SeqRegions(selected_node.getMidBranchLh());
        best_parent_regions = nullptr; // std::make_unique<SeqRegions>(std::move(seq_regions_clone));
        
        // try a shorter split
        // tryShorterBranch<&Tree::calculateSamplePlacementCost>(selected_node->length, best_parent_regions, sample, upper_left_right_regions, lower_regions, best_parent_lh, best_parent_blength_split, default_blength, false);
        tryShorterBranch<&Tree::calculateSamplePlacementCost>(selected_node.getUpperLength(), best_parent_regions, sample, upper_left_right_regions, lower_regions, best_parent_lh, best_parent_blength_split, default_blength, false);
        
        // Delay cloning SeqRegions
        if (!best_parent_regions)
            best_parent_regions = std::make_unique<SeqRegions>(std::move(SeqRegions(selected_node.getMidBranchLh())));
    }
    
    // if the best placement is below the selected_node => add an internal node below the selected_node
    if (best_child_lh >= best_parent_lh && best_child_lh >= best_lh_diff)
    {
        ASSERT(best_child_index.getMiniIndex() == TOP);
        PhyloNode& best_child = nodes[best_child_index.getVectorIndex()];
        const std::unique_ptr<SeqRegions>& upper_left_right_regions = getPartialLhAtNode(best_child.getNeighborIndex(TOP)); // best_child->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
        
        // Estimate the length for the new branch
        RealNumType best_length = default_blength;
        estimateLengthNewBranch<&Tree::calculateSamplePlacementCost>(best_child_lh, best_child_regions, sample, best_length, max_blength, min_blength, false);
        
        // create new internal node and append child to it
        // connectNewSample2Branch(sample, seq_name_index, best_child, best_child_blength_split, best_child->length - best_child_blength_split, best_length, best_child_regions, upper_left_right_regions);
        connectNewSample2Branch(sample, seq_name_index, best_child_index, best_child, best_child_blength_split, best_child.getUpperLength() - best_child_blength_split, best_length, best_child_regions, upper_left_right_regions);
    }
    // otherwise, add new parent to the selected_node
    else
    {
        // new parent is actually part of a polytomy since best placement is exactly at the node
        if (best_lh_diff >= best_parent_lh)
        {
            best_root_blength = -1;
            best_parent_blength_split = -1;
            best_parent_lh = best_lh_diff;
            /*if (best_parent_regions) delete best_parent_regions;
            best_parent_regions = NULL;*/
            best_parent_regions = nullptr;
            
            if (selected_node_vec_index == root_vector_index)
                // selected_node->getPartialLhAtNode(aln, model, threshold_prob)->mergeTwoLowers(best_parent_regions, -1, *sample, default_blength, aln, model, threshold_prob);
                selected_node.getPartialLh(TOP)->mergeTwoLowers(best_parent_regions, -1, *sample, default_blength, aln, model, threshold_prob);
            else
            {
                // best_parent_regions = new SeqRegions(selected_node->total_lh);
                SeqRegions seq_regions_clone = SeqRegions(selected_node.getTotalLh());
                best_parent_regions = std::make_unique<SeqRegions>(std::move(seq_regions_clone));
            }
        }

        // add parent to the root
        if (selected_node_vec_index == root_vector_index)
        {
            // now try different lengths for right branch
            best_parent_lh += old_root_lh;
            RealNumType best_length2 = default_blength;
            const std::unique_ptr<SeqRegions>& lower_regions = selected_node.getPartialLh(TOP); // ->getPartialLhAtNode(aln, model, threshold_prob);
            
            estimateLengthNewBranchAtRoot(sample, lower_regions, best_parent_regions, best_length2, best_parent_lh, best_root_blength, min_blength, false);
            
            // update best_parent_lh (taking into account old_root_lh)
            best_parent_lh -= old_root_lh;
            
            // add new sample to a new root
            connectNewSample2Root(sample, seq_name_index, selected_node_index, selected_node, best_root_blength, best_length2, best_parent_regions);
        }
        //add parent to non-root node
        else
        {
            const std::unique_ptr<SeqRegions>& upper_left_right_regions = getPartialLhAtNode(selected_node.getNeighborIndex(TOP)); // selected_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
            
            // now try different lengths for the new branch
            RealNumType new_branch_length_lh = best_parent_lh;
            RealNumType best_length = default_blength;
            estimateLengthNewBranch<&Tree::calculateSamplePlacementCost>(best_parent_lh, best_parent_regions, sample, best_length, max_blength, min_blength, false);
            
            // now create new internal node and append child to it
            RealNumType down_distance = best_parent_blength_split;
            RealNumType top_distance = selected_node.getUpperLength() - down_distance; // selected_node->length - down_distance;
            if (best_parent_blength_split < 0)
            {
                down_distance = -1;
                top_distance = selected_node.getUpperLength(); // selected_node->length;
                
                /*if (selected_node->total_lh) delete selected_node->total_lh;
                selected_node->total_lh = NULL;
                
                if (selected_node->mid_branch_lh) delete selected_node->mid_branch_lh;
                selected_node->mid_branch_lh = NULL;*/
                selected_node.setTotalLh(nullptr);
                selected_node.setMidBranchLh(nullptr);
                
                // node.furtherMidNodes=None
            }
            connectNewSample2Branch(sample, seq_name_index, selected_node_index, selected_node, top_distance, down_distance, best_length, best_parent_regions, upper_left_right_regions);
        }
    }
    
    // delete best_parent_regions and best_child_regions
    /*if (best_parent_regions)
        delete best_parent_regions;
    if (best_child_regions)
        delete best_child_regions;*/
}

void Tree::refreshAllLhs()
{
    // 1. update all the lower lhs along the tree
    performDFS<&Tree::updateLowerLh>();
    
    // 2. update all the non-lower lhs along the tree
    refreshAllNonLowerLhs();
}

void Tree::refreshUpperLR(const Index node_index, PhyloNode& node, const Index neighbor_index, std::unique_ptr<SeqRegions>& replaced_regions, const SeqRegions& parent_upper_lr_lh)
{
    // recalculate the upper left/right lh of the current node
    std::unique_ptr<SeqRegions> new_upper_lr_lh = nullptr;
    PhyloNode& neighbor = nodes[neighbor_index.getVectorIndex()];
    const std::unique_ptr<SeqRegions>& lower_lh = neighbor.getPartialLh(TOP); // next_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
    // parent_upper_lr_lh.mergeUpperLower(new_upper_lr_lh, node->length, *lower_lh, next_node->length, aln, model, threshold_prob);
    parent_upper_lr_lh.mergeUpperLower(new_upper_lr_lh, node.getUpperLength(), *lower_lh, neighbor.getUpperLength(), aln, model, params->threshold_prob);
    
    // if the upper left/right lh is null -> try to increase the branch length
    if (!new_upper_lr_lh)
    {
        if (neighbor.getUpperLength() <= 0) // next_node->length <= 0)
        {
            stack<Index> node_stack;
            // updateZeroBlength(next_node->neighbor, node_stack, threshold_prob);
            updateZeroBlength(neighbor_index, neighbor, node_stack);
            updatePartialLh(node_stack);
        }
        else if (node.getUpperLength() <= 0) // node->length <= 0)
        {
            stack<Index> node_stack;
            // updateZeroBlength(node, node_stack, threshold_prob);
            updateZeroBlength(node_index, node, node_stack);
            updatePartialLh(node_stack);
        }
        else
            outError("Strange, inconsistent upper left/right lh creation in refreshAllNonLowerLhs()");
    }
    // otherwise, everything is good -> update upper left/right lh of the current node
    else
        // replacePartialLH(replaced_regions, new_upper_lr_lh);
        replaced_regions = std::move(new_upper_lr_lh);
    
    // delete new_upper_lr_lh
    // if (new_upper_lr_lh) delete new_upper_lr_lh;
}

void Tree::refreshNonLowerLhsFromParent(Index& node_index, Index& last_node_index)
{
    PhyloNode& node = nodes[node_index.getVectorIndex()];
    const RealNumType threshold_prob = params->threshold_prob;
    const Index parent_index = node.getNeighborIndex(TOP);
    PhyloNode& parent_node = nodes[parent_index.getVectorIndex()];
    const std::unique_ptr<SeqRegions>& parent_upper_lr_lh = parent_node.getPartialLh(parent_index.getMiniIndex()); // node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
    
    // update the total lh, total lh at the mid-branch point of the current node
    if (node.getUpperLength() > 0) // node->length > 0)
    {
        // update the total lh
        // node->computeTotalLhAtNode(aln, model, threshold_prob, node == root);
        node.updateTotalLhAtNode(parent_node, aln, model, threshold_prob, node_index.getVectorIndex() == root_vector_index);
        
        if (!node.getTotalLh())
            outError("Strange, inconsistent total lh creation in refreshAllNonLowerLhs()");

        // update mid_branch_lh
        computeMidBranchRegions(node, node.getMidBranchLh(), *parent_upper_lr_lh);
        
        // NHANLT: LOGS FOR DEBUGGING
        /*if (params->debug)
        {
            std::cout << "TotalLh " << node.getTotalLh()->size() << std::endl;
            std::cout << "MidbranchLh " << node.getMidBranchLh()->size() << std::endl;
        }*/
    }
    
    // if the current node is an internal node (~having children) -> update its upper left/right lh then traverse downward to update non-lower lhs of other nodes
    if (node.isInternal()) // !node->isLeave())
    {
        /*Node* next_node_1 = node->next;
        Node* next_node_2 = next_node_1->next;*/
        const Index neighbor_1_index = node.getNeighborIndex(RIGHT);
        const Index neighbor_2_index = node.getNeighborIndex(LEFT);
        
        // recalculate the FIRST upper left/right lh of the current node
        // refreshUpperLR(node, next_node_2, next_node_1->partial_lh, *parent_upper_lr_lh);
        refreshUpperLR(node_index, node, neighbor_2_index, node.getPartialLh(RIGHT), *parent_upper_lr_lh);
        
        // recalculate the SECOND upper left/right lh of the current node
        // refreshUpperLR(node, next_node_1, next_node_2->partial_lh, *parent_upper_lr_lh);
        refreshUpperLR(node_index, node, neighbor_1_index, node.getPartialLh(LEFT), *parent_upper_lr_lh);
        
        // NHANLT: LOGS FOR DEBUGGING
        /*if (params->debug)
        {
            std::cout << "next_node_1->partial_lh " << node.getPartialLh(RIGHT)->size() << std::endl;
            std::cout << "next_node_2->partial_lh " << node.getPartialLh(LEFT)->size() << std::endl;
        }*/
        
        // keep traversing downward to its firt child
        node_index = neighbor_1_index;
    }
    // if the current node is a leaf -> traverse upward to its parent
    else
    {
        last_node_index = node_index;
        node_index = node.getNeighborIndex(TOP); // node->neighbor;
    }
}

void Tree::refreshAllNonLowerLhs()
{
    // dummy variables
    const StateType num_states = aln.num_states;
    
    // start from the root
    // update the total lh at root
    //node->computeTotalLhAtNode(aln, model, params->threshold_prob, true);
    PhyloNode& root = nodes[root_vector_index];
    root.setTotalLh(std::move(root.getPartialLh(TOP)->computeTotalLhAtRoot(aln.num_states, model)));
    
    // if the root has children -> update its upper left/right lh then traverse downward to update non-lower lhs of other nodes
    if (root.isInternal()) // !node->isLeave())
    {
        // update upper left/right lh of the root
        /*Node* next_node_1 = node->next;
        Node* next_node_2 = next_node_1->next;*/
        Index neighbor_1_index = root.getNeighborIndex(RIGHT);
        PhyloNode& neighbor_1 = nodes[neighbor_1_index.getVectorIndex()];
        PhyloNode& neighbor_2 = nodes[root.getNeighborIndex(LEFT).getVectorIndex()];
        
        /*delete next_node_1->partial_lh;
        next_node_1->partial_lh = next_node_2->neighbor->getPartialLhAtNode(aln, model, threshold_prob)->computeTotalLhAtRoot(num_states, model, next_node_2->length);*/
        root.setPartialLh(RIGHT, neighbor_2.getPartialLh(TOP)->computeTotalLhAtRoot(num_states, model, neighbor_2.getUpperLength()));
        
        /*delete next_node_2->partial_lh;
        next_node_2->partial_lh = next_node_1->neighbor->getPartialLhAtNode(aln, model, threshold_prob)->computeTotalLhAtRoot(num_states, model, next_node_1->length);*/
        root.setPartialLh(LEFT, neighbor_1.getPartialLh(TOP)->computeTotalLhAtRoot(num_states, model, neighbor_1.getUpperLength()));
        
        // NHANLT: LOGS FOR DEBUGGING
        /*if (params->debug)
        {
            std::cout << "root.setPartialLh " << root.getPartialLh(RIGHT)->size() << std::endl;
            std::cout << "root.setPartialLh " << root.getPartialLh(LEFT)->size() << std::endl;
        }*/
        
        // traverse the tree downward and update the non-lower genome lists for all other nodes of the tree.
        /*Node* last_node = NULL;
        node = next_node_1->neighbor;*/
        Index last_node_index;
        Index node_index = neighbor_1_index;
        while (node_index.getMiniIndex() != UNDEFINED)
        {
            // we reach a top node by a downward traversing
            if (node_index.getMiniIndex() == TOP) // node->is_top)
                refreshNonLowerLhsFromParent(node_index, last_node_index);
            // we reach the current node by an upward traversing from its children
            else
            {
                /*Node* top_node = node->getTopNode();
                Node* next_node_1 = top_node->next;
                Node* next_node_2 = next_node_1->next;*/
                NumSeqsType node_vec = node_index.getVectorIndex();
                PhyloNode& node = nodes[node_vec];
                Index neighbor_1_index = node.getNeighborIndex(RIGHT);
                Index neighbor_2_index = node.getNeighborIndex(LEFT);
                
                // if we reach the current node by an upward traversing from its first children -> traversing downward to its second children
                if (last_node_index == neighbor_1_index) // next_node_1->neighbor)
                    node_index = neighbor_2_index; // node = next_node_2->neighbor;
                // otherwise, all children of the current node are updated -> update the lower lh of the current node
                else
                {
                    /*last_node_index = top_node;
                    node = top_node->neighbor;*/
                    last_node_index = Index(node_vec, TOP);
                    node_index = node.getNeighborIndex(TOP);
                }
            }
        }
    }
}

void Tree::resetSPRFlags(const bool reset_outdated)
{
    // start from the root
    stack<NumSeqsType> node_stack;
    node_stack.push(root_vector_index);
    
    // traverse downward to set all descentdant outdated
    while (!node_stack.empty())
    {
        // pick the top node from the stack
        PhyloNode& node = nodes[node_stack.top()];
        node_stack.pop();
        
        // set the current node outdated
        if (reset_outdated)
            node.setOutdated(true);
        node.setSPRApplied(false);
        
        // traverse downward
        /*for (Index neighbor_index:node.getNeighborIndexes(TOP))
            node_stack.push(neighbor_index.getVectorIndex());*/
        if (node.isInternal())
        {
            node_stack.push(node.getNeighborIndex(RIGHT).getVectorIndex());
            node_stack.push(node.getNeighborIndex(LEFT).getVectorIndex());
        }
    }
}

RealNumType Tree::improveEntireTree(bool short_range_search)
{
    // start from the root
    stack<Index> node_stack;
    node_stack.push(Index(root_vector_index, TOP));
    
    // dummy variables
    RealNumType total_improvement = 0;
    PositionType num_nodes = 0;
    
    // traverse downward the tree
    while (!node_stack.empty())
    {
        // pick the top node from the stack
        Index index = node_stack.top();
        node_stack.pop();
        PhyloNode& node = nodes[index.getVectorIndex()];
        // MiniIndex mini_index = index.getMiniIndex();

        // add all children of the current nodes to the stack for further traversing later
        /*for (Index neighbor_index:node.getNeighborIndexes(mini_index))
            node_stack.push(neighbor_index);*/
        ASSERT(index.getMiniIndex() == TOP);
        if (node.isInternal())
        {
            node_stack.push(node.getNeighborIndex(RIGHT));
            node_stack.push(node.getNeighborIndex(LEFT));
        }
        
        // only process outdated node to avoid traversing the same part of the tree multiple times
        if (node.isOutdated() && !node.isSPRApplied())
        {
            node.setOutdated(false);
            
            //if checkEachSPR:
              //  root=node
               // while root.up!=None:
               //     root=root.up
               // #print("Pre-SPR tree: "+createBinaryNewick(root))
               // oldTreeLK=calculateTreeLikelihood(root,mutMatrix,checkCorrectness=True)
               // #print("Pre-SPR tree likelihood: "+str(oldTreeLK))
                //reCalculateAllGenomeLists(root,mutMatrix, checkExistingAreCorrect=True)
            
            // do SPR moves to improve the tree
            RealNumType improvement = improveSubTree(index, node, short_range_search);
            
            //if checkEachSPR:
      //          #print(" apparent improvement "+str(improvement))
        //        root=node
         //       while root.up!=None:
           //         root=root.up
     //           #print("Post-SPR tree: "+createBinaryNewick(root))
       //         newTreeLK=calculateTreeLikelihood(root,mutMatrix)
         //       reCalculateAllGenomeLists(root,mutMatrix, checkExistingAreCorrect=True)
           //     #print("Post-SPR tree likelihood: "+str(newTreeLK))
             //   if newTreeLK-oldTreeLK < improvement-1.0:
      //              print("In startTopologyUpdates, LK score of improvement "+str(newTreeLK)+" - "+str(oldTreeLK)+" = "+str(newTreeLK-oldTreeLK)+" is less than //what is supposd to be "+str(improvement))
            //        exit()
             //
            
            // update total_improvement
            total_improvement += improvement;
            
            // NHANLT: LOGS FOR DEBUGGING
            /*if (params->debug && improvement > 0)
                std::cout << num_nodes << ": " << std::setprecision(20) << total_improvement << std::endl;*/
            
            // Show log every 1000 nodes
            num_nodes += 1;
            if (num_nodes % 1000 == 0)
                cout << "Processed topology for " << convertIntToString(num_nodes) << " nodes." << endl;
        }
    }
    
    return total_improvement;
}

PositionType Tree::optimizeBranchLengths()
{
    // start from the root's children
    stack<Index> node_stack;
    PhyloNode& root = nodes[root_vector_index];
    if (!root.isInternal()) // !root || !root->next)
        return 0;
    /*Node* neighbor_node = NULL;
    FOR_NEIGHBOR(root, neighbor_node)
        node_stack.push(neighbor_node);*/
    node_stack.push(root.getNeighborIndex(RIGHT));
    node_stack.push(root.getNeighborIndex(LEFT));
    
    // dummy variables
    PositionType num_improvement = 0;
    const RealNumType threshold_prob = params->threshold_prob;
    
    // traverse downward the tree
    while (!node_stack.empty())
    {
        // pick the top node from the stack
        Index node_index = node_stack.top();
        node_stack.pop();
        PhyloNode& node = nodes[node_index.getVectorIndex()];
        
        const Index parent_index = node.getNeighborIndex(TOP);
        const std::unique_ptr<SeqRegions>& upper_lr_regions = getPartialLhAtNode(parent_index); // node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
        const std::unique_ptr<SeqRegions>& lower_regions = node.getPartialLh(TOP); // node->getPartialLhAtNode(aln, model, threshold_prob);
        
        // add all children of the current nodes to the stack for further traversing later
        /*neighbor_node = NULL;
        FOR_NEIGHBOR(node, neighbor_node)
            node_stack.push(neighbor_node);*/
        /*for (Index neighbor_index:node.getNeighborIndexes(TOP))
            node_stack.push(neighbor_index);*/
        if (node.isInternal())
        {
            node_stack.push(node.getNeighborIndex(RIGHT));
            node_stack.push(node.getNeighborIndex(LEFT));
        }
        
        // only process outdated node to avoid traversing the same part of the tree multiple times
        if (node.isOutdated())
        {
            // estimate the branch length
            RealNumType best_length = estimateBranchLength(upper_lr_regions, lower_regions);
            
            if (best_length > 0 || node.getUpperLength() > 0)
            {
                RealNumType diff_thresh = 0.01 * best_length;
                if (best_length <= 0 || node.getUpperLength() <= 0 || (node.getUpperLength() > (best_length + diff_thresh)) || (node.getUpperLength() < (best_length - diff_thresh)))
                {
                    node.setUpperLength(best_length);
                    // node->neighbor->length = node->length;
                    ++num_improvement;
                    
                    // update partial likelihood regions
                    stack<Index> new_node_stack;
                    new_node_stack.push(node_index);
                    new_node_stack.push(parent_index);
                    updatePartialLh(new_node_stack);
                }
            }
        }
    }
    
    return num_improvement;
}

void Tree::estimateBlength_R_O(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType total_blength, const PositionType end_pos, RealNumType &coefficient, std::vector<RealNumType> &coefficient_vec)
{
    const StateType seq1_state = aln.ref_seq[end_pos];
    RealNumType* mutation_mat_row = model.mutation_mat + model.row_index[seq1_state];
    RealNumType coeff0 = seq2_region.getLH(seq1_state);
    RealNumType coeff1 = 0;

    if (seq1_region.plength_observation2root >= 0)
    {
      coeff0 *= model.root_freqs[seq1_state];

      RealNumType* transposed_mut_mat_row = model.transposed_mut_mat + model.row_index[seq1_state];

      assert(aln.num_states == 4);
      updateCoeffs<4>(model.root_freqs, transposed_mut_mat_row, &(*seq2_region.likelihood)[0], mutation_mat_row, seq1_region.plength_observation2node, coeff0, coeff1);

      coeff1 *= model.root_freqs[seq1_state];
    }
    else
    {
      // NHANLT NOTES:
      // x = seq1_state
      // l = log(1 + q_xx * t + sum(q_xy * t)
      // l' = [q_xx + sum(q_xy)]/[1 + q_xx * t + sum(q_xy * t)]
      // coeff1 = numerator = q_xx + sum(q_xy)
        assert(aln.num_states == 4);
        coeff1 += dotProduct<4>(&(*seq2_region.likelihood)[0], mutation_mat_row);
    }

    // NHANLT NOTES:
    // l = log(1 + q_xx * t + sum(q_xy * t)
    // l' = [q_xx + sum(q_xy)]/[1 + q_xx * t + sum(q_xy * t)]
    // coeff0 = denominator = 1 + q_xx * t + sum(q_xy * t)
    if (total_blength > 0)
      coeff0 += coeff1 * total_blength;

    // NHANLT NOTES:
    // l' = [q_xx + sum(q_xy)]/[1 + q_xx * t + sum(q_xy * t)] = coeff1 / coeff0
    if (coeff1 < 0)
      coefficient += coeff1 / coeff0;
    else
      coefficient_vec.push_back(coeff0 / coeff1);
}

void Tree::estimateBlength_R_ACGT(const SeqRegion& seq1_region, const StateType seq2_state, const RealNumType total_blength, const PositionType end_pos, std::vector<RealNumType> &coefficient_vec)
{
    if (seq1_region.plength_observation2root >= 0)
    {
        StateType seq1_state = aln.ref_seq[end_pos];
        
        RealNumType coeff1 = model.root_freqs[seq1_state] * model.mutation_mat[model.row_index[seq1_state] + seq2_state];
        RealNumType coeff0 = model.root_freqs[seq2_state] * model.mutation_mat[model.row_index[seq2_state] + seq1_state] * seq1_region.plength_observation2node;
        
        if (total_blength > 0)
            coeff0 += coeff1 * total_blength;
        
        coefficient_vec.push_back(coeff0 / coeff1);
    }
    // NHANLT: add else here, otherwise, coefficient_vec.push_back(total_blength > 0 ? total_blength : 0) is called even when (seq1_region->plength_observation2root >= 0)
    else
        // NHANLT NOTES:
        // l = log(q_xy * t)
        // l' = q_xy / (q_xy * t) = 1 / t
        coefficient_vec.push_back(total_blength > 0 ? total_blength : 0);
}

void Tree::estimateBlength_O_X(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType total_blength, const PositionType end_pos, RealNumType &coefficient, std::vector<RealNumType> &coefficient_vec)
{
    const StateType num_states = aln.num_states;
    RealNumType coeff0 = 0;
    RealNumType coeff1 = 0;
    
    // 3.1. e1.type = O and e2.type = O
    if (seq2_region.type == TYPE_O)
    {
        RealNumType* mutation_mat_row = model.mutation_mat;
        
        // NHANLT NOTES:
        // l = log(sum_x(1 + q_xx * t + sum_y(q_xy * t)))
        // l' = [sum_x(q_xx + sum_y(q_xy))]/[sum_x(1 + q_xx * t + sum_y(q_xy * t))]
        // coeff1 = numerator = sum_x(q_xx + sum_y(q_xy))
        // coeff0 = denominator = sum_x(1 + q_xx * t + sum_y(q_xy * t))
        for (StateType i = 0; i < num_states; ++i, mutation_mat_row += num_states)
        {
            RealNumType seq1_lh_i = seq1_region.getLH(i);
            coeff0 += seq1_lh_i * seq2_region.getLH(i);
            
            for (StateType j = 0; j < num_states; ++j)
                coeff1 += seq1_lh_i * seq2_region.getLH(j) * mutation_mat_row[j];
        }
    }
    // 3.2. e1.type = O and e2.type = R or A/C/G/T
    else
    {
        StateType seq2_state = seq2_region.type;
        if (seq2_state == TYPE_R)
            seq2_state = aln.ref_seq[end_pos];
        
        coeff0 = seq1_region.getLH(seq2_state);

        // NHANLT NOTES:
        // y = seq2_state
        // l = log(1 + q_yy * t + sum_x(q_xy * t)
        // l' = [q_yy + sum_x(q_xy))]/[1 + q_xx * t + sum_y(q_xy * t)]
        // coeff1 = numerator = q_yy + sum_x(q_xy))
        // coeff0 = denominator = 1 + q_xx * t + sum_y(q_xy * t)
        RealNumType* transposed_mut_mat_row = model.transposed_mut_mat + model.row_index[seq2_state];
        assert(num_states == 4);
        coeff1 += dotProduct<4>(&(*seq1_region.likelihood)[0], transposed_mut_mat_row);
    }
    
    if (total_blength > 0)
        coeff0 += coeff1 * total_blength;
    
    // NHANLT NOTES:
    // l' = coeff1 / coeff0
    if (coeff1 < 0)
        coefficient += coeff1 / coeff0;
    else
        coefficient_vec.push_back(coeff0 / coeff1);
}

void Tree::estimateBlength_ACGT_O(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType total_blength, RealNumType &coefficient, std::vector<RealNumType> &coefficient_vec)
{
    StateType seq1_state = seq1_region.type;
    RealNumType coeff0 = seq2_region.getLH(seq1_state);
    RealNumType coeff1 = 0;
    
    RealNumType* mutation_mat_row = model.mutation_mat + model.row_index[seq1_state];
    
    if (seq1_region.plength_observation2root >= 0)
    {
        coeff0 *= model.root_freqs[seq1_state];

        RealNumType* transposed_mut_mat_row = model.transposed_mut_mat + model.row_index[seq1_state];
                            
        assert(aln.num_states == 4);
        updateCoeffs<4>(model.root_freqs, transposed_mut_mat_row, &(*seq2_region.likelihood)[0], mutation_mat_row, seq1_region.plength_observation2node, coeff0, coeff1);
        
        coeff1 *= model.root_freqs[seq1_state];
    }
    else
    {
        // NHANLT NOTES:
        // x = seq1_state
        // l = log(1 + q_xx * t + sum(q_xy * t)
        // l' = [q_xx + sum(q_xy)]/[1 + q_xx * t + sum(q_xy * t)]
        // coeff1 = numerator = q_xx + sum(q_xy)
        assert(aln.num_states == 4);
        coeff1 += dotProduct<4>(&(*seq2_region.likelihood)[0], mutation_mat_row);
    }
    
    // NHANLT NOTES:
    // l = log(1 + q_xx * t + sum(q_xy * t)
    // l' = [q_xx + sum(q_xy)]/[1 + q_xx * t + sum(q_xy * t)]
    // coeff0 = denominator = 1 + q_xx * t + sum(q_xy * t)
    if (total_blength > 0)
        coeff0 += coeff1 * total_blength;
    
    // NHANLT NOTES:
    // l' = [q_xx + sum(q_xy)]/[1 + q_xx * t + sum(q_xy * t)] = coeff1 / coeff0;
    if (coeff1 < 0)
        coefficient += coeff1 / coeff0;
    else
        coefficient_vec.push_back(coeff0 / coeff1);
}

void Tree::estimateBlength_ACGT_RACGT(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType total_blength, const PositionType end_pos, std::vector<RealNumType> &coefficient_vec)
{
    RealNumType coeff0 = 0;
    StateType seq1_state = seq1_region.type;
    StateType seq2_state = seq2_region.type;
    if (seq2_state == TYPE_R)
        seq2_state = aln.ref_seq[end_pos];
    
    if (seq1_region.plength_observation2root >= 0)
    {
        coeff0 = model.root_freqs[seq2_state] * model.mutation_mat[model.row_index[seq2_state] + seq1_state] * seq1_region.plength_observation2node;
        RealNumType coeff1 = model.root_freqs[seq1_state] * model.mutation_mat[model.row_index[seq1_state] + seq2_state];
        
        if (total_blength > 0)
            coeff0 += coeff1 * total_blength;
        
        coeff0 /= coeff1;
    }
    // NHANLT NOTES:
    // l = log(q_xy * t)
    // l' = q_xy / (q_xy * t) = 1 / t
    else if (total_blength > 0)
        coeff0 = total_blength;
    
    coefficient_vec.push_back(coeff0);
}

RealNumType Tree::estimateBlengthFromCoeffs(RealNumType &coefficient, const std::vector<RealNumType> coefficient_vec)
{
    coefficient = -coefficient;
    PositionType num_coefficients = coefficient_vec.size();
    if (num_coefficients == 0)
        return -1;

    // Get min and max coefficients
    RealNumType min_coefficient = coefficient_vec[0];
    RealNumType max_coefficient = coefficient_vec[0];
    for (PositionType i = 1; i < num_coefficients; ++i)
    {
        RealNumType coefficient_i = coefficient_vec[i];
        if (coefficient_i < min_coefficient)
            min_coefficient = coefficient_i;
        if (coefficient_i > max_coefficient)
            max_coefficient = coefficient_i;
    }
    
    RealNumType num_coefficients_over_coefficient = num_coefficients / coefficient;
    RealNumType tDown = num_coefficients_over_coefficient - min_coefficient;
    if (tDown <= 0)
        return 0;
    RealNumType derivative_tDown = calculateDerivative(coefficient_vec, tDown);
    
    RealNumType tUp = num_coefficients_over_coefficient - max_coefficient;
    if (tUp < 0)
    {
        if (min_coefficient > 0)
            tUp = 0;
        else
            tUp = min_blength_sensitivity;
    }
    RealNumType derivative_tUp = calculateDerivative(coefficient_vec, tUp);
    
    if ((derivative_tDown > coefficient + min_blength_sensitivity) || (derivative_tUp < coefficient - min_blength_sensitivity))
        if ((derivative_tUp < coefficient - min_blength_sensitivity) && (tUp == 0))
            return 0;
    
    while (tDown - tUp > min_blength_sensitivity)
    {
        RealNumType tMiddle = (tUp + tDown) * 0.5;
        RealNumType derivative_tMiddle = calculateDerivative(coefficient_vec, tMiddle);
        
        if (derivative_tMiddle > coefficient)
            tUp = tMiddle;
        else
            tDown = tMiddle;
    }
    
    return tUp;
}

using DoubleState = uint16_t;
static constexpr DoubleState RR = (DoubleState(TYPE_R) << 8) | TYPE_R;
static constexpr DoubleState RO = (DoubleState(TYPE_R) << 8) | TYPE_O;
static constexpr DoubleState OO = (DoubleState(TYPE_O) << 8) | TYPE_O;

RealNumType Tree::estimateBranchLength(const std::unique_ptr<SeqRegions>& parent_regions, const std::unique_ptr<SeqRegions>& child_regions)
{
    // init dummy variables
    RealNumType coefficient = 0;
    vector<RealNumType> coefficient_vec;
    PositionType pos = 0;
    const SeqRegions& seq1_regions = *parent_regions;
    const SeqRegions& seq2_regions = *child_regions;
    size_t iseq1 = 0;
    size_t iseq2 = 0;
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    const RealNumType* const &cumulative_rate = model.cumulative_rate;
    
    // avoid reallocations
    coefficient_vec.reserve(parent_regions->countSharedSegments(seq2_regions, seq_length)); // avoid realloc of vector data

    while (pos < seq_length)
    {
        PositionType end_pos;
        RealNumType total_blength;
        
        // get the next shared segment in the two sequences
        SeqRegions::getNextSharedSegment(pos, seq1_regions, seq2_regions, iseq1, iseq2, end_pos);
        const auto* seq1_region = &seq1_regions[iseq1];
        const auto* seq2_region = &seq2_regions[iseq2];
        
        // 1. e1.type = N || e2.type = N
        if ((seq2_region->type == TYPE_N) + (seq1_region->type == TYPE_N))
        {
            pos = end_pos + 1;
            continue;
        }
        
        // e1.type != N && e2.type != N
        const DoubleState s1s2 = (DoubleState(seq1_region->type) << 8) | seq2_region->type;

        // total_blength will be here the total length from the root or from the upper node, down to the down node.
        if (seq1_region->plength_observation2root >= 0)
            total_blength = seq1_region->plength_observation2root;
        else if (seq1_region->plength_observation2node >= 0)
            total_blength = seq1_region->plength_observation2node;
        else
            total_blength = 0;
            
        if (seq2_region->plength_observation2node >= 0)
            total_blength = total_blength + seq2_region->plength_observation2node;
            //total_blength = (total_blength > 0 ? total_blength : 0) + seq2_region->plength_observation2node;
        
        // 2. e1.type = R
        // 2.1.e1.type = R and e2.type = R
        if (s1s2 == RR)
        {
            // NHANLT NOTES:
            // coefficient = derivative of log likelihood function wrt t
            // l = log(1 + q_xx * t) ~ q_xx * t
            // => l' = q_xx
            //if (seq2_region->type == TYPE_R)
              coefficient += cumulative_rate[end_pos + 1] - cumulative_rate[pos];
        }
        // 2.2. e1.type = R and e2.type = O
        else if (s1s2 == RO)
        {
            estimateBlength_R_O(*seq1_region, *seq2_region, total_blength, end_pos, coefficient, coefficient_vec);
        }
        // 2.3. e1.type = R and e2.type = A/C/G/T
        else if (seq1_region->type == TYPE_R)
        {
            estimateBlength_R_ACGT(*seq1_region, seq2_region->type, total_blength, end_pos, coefficient_vec);
        }
        // 3. e1.type = O
        else if (seq1_region->type == TYPE_O)
        {
            estimateBlength_O_X(*seq1_region, *seq2_region, total_blength, end_pos, coefficient, coefficient_vec);
        }
        // 4. e1.type = A/C/G/T
        // 4.1. e1.type =  e2.type
        // NHANLT NOTES:
        // coefficient = derivative of log likelihood function wrt t
        // l = log(1 + q_xx * t) ~ q_xx * t
        // => l' = q_xx
        else if (seq1_region->type == seq2_region->type)
            coefficient += model.diagonal_mut_mat[seq1_region->type];
        // e1.type = A/C/G/T and e2.type = O/A/C/G/T
        // 4.2. e1.type = A/C/G/T and e2.type = O
        else if (seq2_region->type == TYPE_O)
        {
            estimateBlength_ACGT_O(*seq1_region, *seq2_region, total_blength, coefficient, coefficient_vec);
        }
        // 4.3. e1.type = A/C/G/T and e2.type = R or A/C/G/T
        else
        {
            estimateBlength_ACGT_RACGT(*seq1_region, *seq2_region, total_blength, end_pos, coefficient_vec);
        }
        
        // update pos
        pos = end_pos + 1;
    }
    
    // now optimized branch length based on coefficients
    return estimateBlengthFromCoeffs(coefficient, coefficient_vec);
}

RealNumType Tree::calculateDerivative(const vector<RealNumType> &coefficient_vec, const RealNumType delta_t)
{
    RealNumType result = 0;
    
    for (RealNumType coefficient : coefficient_vec)
        result += 1.0 / (coefficient + delta_t);
        
    return result;
}

void Tree::handleBlengthChanged(PhyloNode& node, const Index node_index, const RealNumType best_blength)
{
    /*node->length = best_blength;
    node->neighbor->length = node->length;*/
    node.setUpperLength(best_blength);
    
    stack<Index> node_stack;
    node_stack.push(node_index);
    node_stack.push(node.getNeighborIndex(TOP));
    // node_stack.push(node->neighbor);
    updatePartialLh(node_stack);
}

void Tree::optimizeBlengthBeforeSeekingSPR(PhyloNode& node, RealNumType &best_blength, RealNumType &best_lh, bool &blength_changed, const std::unique_ptr<SeqRegions>& parent_upper_lr_lh, const std::unique_ptr<SeqRegions>& lower_lh)
{
    RealNumType original_lh = best_lh;
    
    // try different branch lengths for the current node placement (just in case branch length can be improved, in which case it counts both as tree improvment and better strategy to find a new placement).
    if (node.getUpperLength() <= 0)
    {
        best_blength = min_blength;
        best_lh = calculateSubTreePlacementCost(parent_upper_lr_lh, lower_lh, best_blength);
    }
    
    // cache best_blength
    const RealNumType cached_blength = best_blength;

    // try shorter branch lengths
    bool found_new_blength = tryShorterNewBranch<&Tree::calculateSubTreePlacementCost>(parent_upper_lr_lh, lower_lh, best_blength, best_lh, double_min_blength);
    
    // try longer branch lengths
    if (!found_new_blength)
        tryLongerNewBranch<&Tree::calculateSubTreePlacementCost>(parent_upper_lr_lh, lower_lh, best_blength, best_lh, half_max_blength);
    
    // update blength_changed
    if (cached_blength != best_blength)
        blength_changed = true;
    
    if (node.getUpperLength() <= 0 && original_lh > best_lh)
        best_lh = original_lh;
}

void Tree::checkAndApplySPR(const RealNumType best_lh_diff, const RealNumType best_blength, const RealNumType best_lh, const Index node_index, PhyloNode& node, const Index best_node_index, const Index parent_node_index, const bool is_mid_node, RealNumType& total_improvement, bool& topology_updated)
{
    const NumSeqsType parent_node_vec = parent_node_index.getVectorIndex();
    const NumSeqsType best_node_vec = best_node_index.getVectorIndex();
    PhyloNode& parent_node = nodes[parent_node_vec];
    if (best_node_vec == parent_node_vec)
        outWarning("Strange, re-placement is at same node");
    else if ((best_node_vec == parent_node.getNeighborIndex(LEFT).getVectorIndex() || best_node_vec == parent_node.getNeighborIndex(RIGHT).getVectorIndex()) && is_mid_node) // ((best_node == parent_node->next->neighbor || best_node == parent_node->next->next->neighbor) && is_mid_node)
        cout << "Re-placement is above sibling node";
    else [[likely]]
    {
        // reach the top of a multifurcation, which is the only place in a multifurcatio where placement is allowed.
        NumSeqsType top_polytomy_vec = best_node_vec;

        while (nodes[top_polytomy_vec].getUpperLength() <= 0 && top_polytomy_vec != root_vector_index)// top_polytomy->length <= 0 && top_polytomy!= root)
        {
            // top_polytomy = top_polytomy->neighbor->getTopNode();
            top_polytomy_vec = nodes[top_polytomy_vec].getNeighborIndex(TOP).getVectorIndex();
        }
        
        if (top_polytomy_vec != best_node_vec)
            outWarning("Strange, placement node not at top of polytomy");
        
        // reach the top of the multifurcation of the parent
        /*Node* parent_top_polytomy = parent_node;
        while (parent_top_polytomy->length <= 0 && parent_top_polytomy != root)
            parent_top_polytomy = parent_top_polytomy->neighbor->getTopNode();*/
        NumSeqsType parent_top_polytomy_vec = parent_node_vec;

        while (nodes[parent_top_polytomy_vec].getUpperLength() <= 0 && parent_top_polytomy_vec != root_vector_index)
            parent_top_polytomy_vec = nodes[parent_top_polytomy_vec].getNeighborIndex(TOP).getVectorIndex();
        
        if (!(parent_top_polytomy_vec == top_polytomy_vec && !is_mid_node))
        {
            total_improvement = best_lh_diff - best_lh;
            
            if (verbose_mode == VB_DEBUG)
                cout << "In improveSubTree() found SPR move with improvement " << total_improvement << endl;
            
            // apply an SPR move
            applySPR(node_index, node, best_node_index, is_mid_node, best_blength, best_lh_diff);
            
            topology_updated = true;
        }
    }
}

RealNumType Tree::improveSubTree(const Index node_index, PhyloNode& node, bool short_range_search)
{
    // dummy variables
    ASSERT(node_index.getMiniIndex() == TOP);
    NumSeqsType vec_index = node_index.getVectorIndex();
    const RealNumType threshold_prob = params->threshold_prob;
    const RealNumType thresh_placement_cost = short_range_search ? params->thresh_placement_cost_short_search : params->thresh_placement_cost;
    RealNumType total_improvement = 0;
    bool blength_changed = false; // true if a branch length has been changed
    
    // we avoid the root node since it cannot be re-placed with SPR moves
    if (vec_index != root_vector_index)
    {
        // evaluate current placement
        const std::unique_ptr<SeqRegions>& parent_upper_lr_lh = getPartialLhAtNode(node.getNeighborIndex(TOP)); // node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
        const std::unique_ptr<SeqRegions>& lower_lh = node.getPartialLh(TOP); // node->getPartialLhAtNode(aln, model, threshold_prob);
        RealNumType best_blength = node.getUpperLength(); // node->length;
        RealNumType best_lh = calculateSubTreePlacementCost(parent_upper_lr_lh, lower_lh, best_blength);
        
        // optimize branch length
        if (best_lh < thresh_placement_cost)
            optimizeBlengthBeforeSeekingSPR(node, best_blength, best_lh, blength_changed, parent_upper_lr_lh, lower_lh);
           
        // find new placement
        if (best_lh < thresh_placement_cost)
        {
            // now find the best place on the tree where to re-attach the subtree rooted at "node" but to do that we need to consider new vector probabilities after removing the node that we want to replace this is done using findBestParentTopology().
            bool topology_updated = false;
            const Index parent_index = node.getNeighborIndex(TOP);
            // PhyloNode parent_node = nodes[parent_index.getVectorIndex()]; // node->neighbor->getTopNode();
            Index best_node_index;
            RealNumType best_lh_diff = best_lh;
            bool is_mid_node = false;
            RealNumType best_up_lh_diff = MIN_NEGATIVE;
            RealNumType best_down_lh_diff = MIN_NEGATIVE;
            Index best_child_index;
            
            // seek a new placement for the subtree
            seekSubTreePlacement(best_node_index, best_lh_diff, is_mid_node, best_up_lh_diff, best_down_lh_diff, best_child_index, short_range_search, node_index, best_blength); // , true, NULL);
            
            // validate the new placement cost
            if (best_lh_diff > params->threshold_prob2)
                outError("Strange, lh cost is positive");
            else if (best_lh_diff < -1e50)
                outError("Likelihood cost is very heavy, this might mean that the reference used is not the same used to generate the input diff file");
            
            if (best_lh_diff + thresh_placement_cost > best_lh)
            {
                // check and apply SPR move
                checkAndApplySPR(best_lh_diff, best_blength, best_lh, node_index, node, best_node_index, parent_index, is_mid_node, total_improvement, topology_updated);
                
                if (!topology_updated && blength_changed)
                    handleBlengthChanged(node, node_index, best_blength);
            }
            else if (blength_changed)
                handleBlengthChanged(node, node_index, best_blength);
        }
        else if (blength_changed)
            handleBlengthChanged(node, node_index, best_blength);
    }
    
                        
    return total_improvement;
}

void calculateSubtreeCost_R_R(const SeqRegion& seq1_region, const RealNumType* const &cumulative_rate, RealNumType& total_blength, const PositionType pos, const PositionType end_pos, RealNumType& lh_cost)
{
    if (seq1_region.plength_observation2root >= 0)
      total_blength += seq1_region.plength_observation2node;

    // NHANLT NOTE:
    // approximation log(1+x)~x
    if (total_blength > 0)
      lh_cost += total_blength * (cumulative_rate[end_pos + 1] - cumulative_rate[pos]);
}

template <const StateType num_states>
void calculateSubtreeCost_R_O(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType total_blength, const StateType seq1_state, RealNumType& total_factor, const Model& model)
{
    RealNumType tot = 0;
    
    if (seq1_region.plength_observation2root >= 0)
    {
        RealNumType* transposed_mut_mat_row = model.transposed_mut_mat + model.row_index[seq1_state];
        RealNumType* mutation_mat_row = model.mutation_mat;
                            
        for (StateType i = 0; i < num_states; ++i, mutation_mat_row += num_states)
        {
            // NHANLT NOTE: UNSURE
            // tot2: likelihood that we can observe seq1_state elvoving from i at root (account for the fact that the observation might have occurred on the other side of the phylogeny with respect to the root)
            // tot2 = root_freqs[seq1_state] * (1 + mut[seq1_state,seq1_state] * plength_observation2node) + root_freqs[i] * mut[i,seq1_state] * plength_observation2node
            RealNumType tot2 = model.root_freqs[i] * transposed_mut_mat_row[i] * seq1_region.plength_observation2node + (seq1_state == i ? model.root_freqs[i] : 0);
            
            // NHANLT NOTE:
            // tot3: likelihood of i evolves to j
            // tot3 = (1 + mut[i,i] * total_blength) * lh(seq2,i) + mut[i,j] * total_blength * lh(seq2,j)
            RealNumType tot3 = total_blength > 0 ? (total_blength * dotProduct<num_states>(mutation_mat_row, &((*seq2_region.likelihood)[0]))) : 0;
            
            // NHANLT NOTE:
            // tot = tot2 * tot3
            tot += tot2 * (seq2_region.getLH(i) + tot3);
        }
        
        // NHANLT NOTE: UNCLEAR
        // why we need to divide tot by root_freqs[seq1_state]
        tot *= model.inverse_root_freqs[seq1_state];
    }
    else
    {
        // NHANLT NOTE:
        // (1 + mut[seq1_state,seq1_state] * total_blength) * lh(seq2,seq1_state) + mut[seq1_state,j] * total_blength * lh(seq2,j)
        if (total_blength > 0)
        {
            const RealNumType* mutation_mat_row = model.mutation_mat + model.row_index[seq1_state];
            tot += dotProduct<num_states>(mutation_mat_row, &((*seq2_region.likelihood)[0]));
            tot *= total_blength;
        }
        tot += seq2_region.getLH(seq1_state);
    }
    total_factor *= tot;
}

bool calculateSubtreeCost_R_ACGT(const SeqRegion& seq1_region, const RealNumType total_blength, const StateType seq1_state, const StateType seq2_state, RealNumType& total_factor, const Model& model)
{
    if (seq1_region.plength_observation2root >= 0)
    {
        if (total_blength > 0)
        {
            // NHANLT NOTE: UNSURE
            // seq1_state_evolves_seq2_state = (1) the likelihood that seq1_state evolves to seq2_state * (2) the likelihood that seq1_state unchanges from the observing position
            // (1) = model.mutation_mat[model.row_index[seq1_state] + seq2_state] * total_blength
            // (2) = (1.0 + model.diagonal_mut_mat[seq1_state] * seq1_region.plength_observation2node)
            RealNumType seq1_state_evolves_seq2_state = model.mutation_mat[model.row_index[seq1_state] + seq2_state] * total_blength * (1.0 + model.diagonal_mut_mat[seq1_state] * seq1_region.plength_observation2node);
            
            // NHANLT NOTE: UNCLEAR
            // consider the inverse process of the above
            // seq2_state_evolves_seq1_state = (1) the likelihood that seq2_state evolves to seq1_state * (2) the likelihood that seq2_state unchanges from the observing position
            // (1) = root_freqs[seq2_state] / root_freqs[seq1_state] * mutation_mat[model.row_index[seq2_state] + seq1_state] * seq1_region.plength_observation2node
            // (2) = (1.0 + model.diagonal_mut_mat[seq2_state] * total_blength)
            RealNumType seq2_state_evolves_seq1_state = model.freqi_freqj_qij[model.row_index[seq2_state] + seq1_state] * seq1_region.plength_observation2node * (1.0 + model.diagonal_mut_mat[seq2_state] * total_blength);
            
            total_factor *= seq1_state_evolves_seq2_state + seq2_state_evolves_seq1_state;
        }
        // NHANLT NOTE:
        // the same as above but total_blength = 0 then we simplify the formula to save the runtime (avoid multiplying with 0)
        else
            total_factor *= model.freqi_freqj_qij[model.row_index[seq2_state] + seq1_state] * seq1_region.plength_observation2node;
    }
    // NHANLT NOTE:
    // add the likelihood that seq1_state evoles to seq2_state = mut[seq1_state,seq2_state] * total_blength
    else if (total_blength > 0)
        total_factor *= model.mutation_mat[model.row_index[seq1_state] + seq2_state] * total_blength;
    else
        return false; // return MIN_NEGATIVE;

    // no error
    return true;
}

template <const StateType num_states>
void calculateSubtreeCost_O_O(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType total_blength, RealNumType& total_factor, const Model& model)
{
    if (total_blength > 0)
    {
      total_factor *= matrixEvolve<num_states>(&((*seq1_region.likelihood)[0]), &((*seq2_region.likelihood)[0]), model.mutation_mat, total_blength);
    }
    // NHANLT NOTE:
    // the same as above but total_blength = 0 then we simplify the formula to save the runtime (avoid multiplying with 0)
    else
    {
      total_factor *= dotProduct<num_states>(&((*seq1_region.likelihood)[0]), &((*seq2_region.likelihood)[0]));
    }
}

template <const StateType num_states>
void calculateSubtreeCost_O_RACGT(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType total_blength, const PositionType end_pos, RealNumType& total_factor, const Alignment& aln, const Model& model)
{
    StateType seq2_state = seq2_region.type;
    if (seq2_state == TYPE_R)
        seq2_state = aln.ref_seq[end_pos];
    
    if (total_blength > 0)
    {
        // NHANLT NOTE:
        // tot2: likelihood of i evolves to seq2_state
        // tot2 = (1 + mut[seq2_state,seq2_state] * total_blength) * lh(seq1,seq2_state) + lh(seq1,i) * mut[i,seq2_state] * total_blength
        RealNumType* transposed_mut_mat_row = model.transposed_mut_mat + model.row_index[seq2_state];
        RealNumType tot2 = dotProduct<num_states>(&((*seq1_region.likelihood)[0]), transposed_mut_mat_row);
        total_factor *= seq1_region.getLH(seq2_state) + total_blength * tot2;
    }
    // NHANLT NOTE:
    // the same as above but total_blength = 0 then we simplify the formula to save the runtime (avoid multiplying with 0)
    else
        total_factor *= seq1_region.getLH(seq2_state);
}

void calculateSubtreeCost_identicalACGT(const SeqRegion& seq1_region, RealNumType& total_blength, RealNumType& lh_cost, const Model& model)
{
    if (seq1_region.plength_observation2root >= 0)
        total_blength += seq1_region.plength_observation2node;
    
    // NHANLT NOTE:
    // the likelihood that seq1_state unchanges
    if (total_blength > 0)
        lh_cost += model.diagonal_mut_mat[seq1_region.type] * total_blength;
}

template <const StateType num_states>
void calculateSubtreeCost_ACGT_O(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType total_blength, RealNumType& total_factor, const Model& model)
{
    StateType seq1_state = seq1_region.type;
    if (seq1_region.plength_observation2root >= 0)
    {
        RealNumType* transposed_mut_mat_row = model.transposed_mut_mat + model.row_index[seq1_state];
        RealNumType* mutation_mat_row = model.mutation_mat;
        RealNumType tot = matrixEvolveRoot<num_states>(&((*seq2_region.likelihood)[0]), seq1_state,
          model.root_freqs, transposed_mut_mat_row, mutation_mat_row, total_blength, seq1_region.plength_observation2node);
        // NHANLT NOTE: UNCLEAR
        // why we need to divide tot by root_freqs[seq1_state]
        total_factor *= (tot * model.inverse_root_freqs[seq1_state]);
    }
    else
    {
        RealNumType* mutation_mat_row = model.mutation_mat + model.row_index[seq1_state];
        
        // NHANLT NOTE:
        // tot = the likelihood of seq1_state evolving to j
        // (1 + mut[seq1_state,seq1_state] * total_blength) * lh(seq2,seq1_state) + mut[seq1_state,j] * total_blength * lh(seq2,j)
        RealNumType tot = dotProduct<num_states>(mutation_mat_row, &((*seq2_region.likelihood)[0]));
        tot *= total_blength;
        tot += seq2_region.getLH(seq1_state);
        total_factor *= tot;
    }
}

bool calculateSubtreeCost_ACGT_RACGT(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType total_blength, const PositionType end_pos, RealNumType& total_factor, const Alignment& aln, const Model& model)
{
    StateType seq1_state = seq1_region.type;
    StateType seq2_state = seq2_region.type;
    if (seq2_state == TYPE_R)
        seq2_state = aln.ref_seq[end_pos];
    
    if (seq1_region.plength_observation2root >= 0)
    {
        if (total_blength > 0)
        {
            // NHANLT NOTE: UNSURE
            // seq1_state_evolves_seq2_state = (1) the likelihood that seq1_state evolves to seq2_state * (2) the likelihood that seq1_state unchanges from the observing position
            // (1) = model.mutation_mat[model.row_index[seq1_state] + seq2_state] * total_blength
            // (2) = (1.0 + model.diagonal_mut_mat[seq1_state] * seq1_region.plength_observation2node)
            RealNumType seq1_state_evolves_seq2_state = model.mutation_mat[model.row_index[seq1_state] + seq2_state] * total_blength * (1.0 + model.diagonal_mut_mat[seq1_state] * seq1_region.plength_observation2node);
            
            // NHANLT NOTE: UNCLEAR
            // consider the inverse process of the above
            // seq2_state_evolves_seq1_state = (1) the likelihood that seq2_state evolves to seq1_state * (2) the likelihood that seq2_state unchanges from the observing position
            // (1) = root_freqs[seq2_state] / root_freqs[seq1_state] * mutation_mat[model.row_index[seq2_state] + seq1_state] * seq1_region.plength_observation2node
            // (2) = (1.0 + model.diagonal_mut_mat[seq2_state] * total_blength)
            RealNumType seq2_state_evolves_seq1_state = model.freqi_freqj_qij[model.row_index[seq2_state] + seq1_state] * seq1_region.plength_observation2node * (1.0 + model.diagonal_mut_mat[seq2_state] * total_blength);
            
            total_factor *= seq1_state_evolves_seq2_state + seq2_state_evolves_seq1_state;
        }
        // NHANLT NOTE:
        // the same as above but total_blength = 0 then we simplify the formula to save the runtime (avoid multiplying with 0)
        else
            total_factor *= model.freqi_freqj_qij[model.row_index[seq2_state] + seq1_state] * seq1_region.plength_observation2node;
    }
    // NHANLT NOTE:
    // add the likelihood that seq1_state evoles to seq2_state = mut[seq1_state,seq2_state] * total_blength
    else if (total_blength > 0)
        total_factor *= model.mutation_mat[model.row_index[seq1_state] + seq2_state] * total_blength;
    else
        return false; // return MIN_NEGATIVE;

    // no error
    return true;
}

RealNumType Tree::calculateSubTreePlacementCost(const std::unique_ptr<SeqRegions>& parent_regions, const std::unique_ptr<SeqRegions>& child_regions, const RealNumType blength)
{
    return (this->*calculateSubTreePlacementCostPointer)(parent_regions, child_regions, blength);
}

// this implementation derives from appendProbNode
template <const StateType num_states>
RealNumType Tree::calculateSubTreePlacementCostTemplate(const std::unique_ptr<SeqRegions>& parent_regions, const std::unique_ptr<SeqRegions>& child_regions, const RealNumType blength)
{
    // 55% of runtime
    // init dummy variables
    RealNumType lh_cost = 0;
    PositionType pos = 0;
    RealNumType total_factor = 1;
    const SeqRegions& seq1_regions = *parent_regions;
    const SeqRegions& seq2_regions = *child_regions;
    size_t iseq1 = 0;
    size_t iseq2 = 0;
    const PositionType seq_length = aln.ref_seq.size();
    
    while (pos < seq_length)
    {
        PositionType end_pos;
        RealNumType total_blength;
        
        // get the next shared segment in the two sequences
        SeqRegions::getNextSharedSegment(pos, seq1_regions, seq2_regions, iseq1, iseq2, end_pos);
        const auto* const seq1_region = &seq1_regions[iseq1];
        const auto* const seq2_region = &seq2_regions[iseq2];

        // 1. e1.type = N || e2.type = N
        if ((seq2_region->type == TYPE_N) + (seq1_region->type == TYPE_N))
        {
            pos = end_pos + 1;
            continue;
        }
       
        // e1.type != N && e2.type != N
        const DoubleState s1s2 = (DoubleState(seq1_region->type) << 8) | seq2_region->type;
        
        // total_blength will be here the total length from the root or from the upper node, down to the down node.
        if (seq1_region->plength_observation2root >= 0)
            total_blength = seq1_region->plength_observation2root + (blength >= 0 ? blength : 0);
        else if (seq1_region->plength_observation2node >= 0)
            total_blength = seq1_region->plength_observation2node + (blength >= 0 ? blength : 0);
        else
            total_blength = blength;
            
        if (seq2_region->plength_observation2node >= 0)
            total_blength = (total_blength > 0 ? total_blength : 0) + seq2_region->plength_observation2node;
        
        
        //assert(total_blength >= 0); // can be -1 ..

        // 2.1. e1.type = R and e2.type = R
        if (s1s2 == RR) [[likely]]
        {
            calculateSubtreeCost_R_R(*seq1_region, model.cumulative_rate, total_blength, pos, end_pos, lh_cost);
        }
        // 2.2. e1.type = R and e2.type = O
        else if (s1s2 == RO)
        {
            calculateSubtreeCost_R_O<num_states>(*seq1_region, *seq2_region, total_blength, aln.ref_seq[end_pos], total_factor, model);
        }
        // 2.3. e1.type = R and e2.type = A/C/G/T
        else if (seq1_region->type == TYPE_R)
        {
            if (!calculateSubtreeCost_R_ACGT(*seq1_region, total_blength, aln.ref_seq[end_pos], seq2_region->type, total_factor, model)) return MIN_NEGATIVE;
        }
        // 3. e1.type = O
        // 3.1. e1.type = O and e2.type = O
        else if (s1s2 == OO)
        {
            calculateSubtreeCost_O_O<num_states>(*seq1_region, *seq2_region, total_blength, total_factor, model);
        }
        // 3.2. e1.type = O and e2.type = R or A/C/G/T
        else if (seq1_region->type == TYPE_O)
        {
            calculateSubtreeCost_O_RACGT<num_states>(*seq1_region, *seq2_region, total_blength, end_pos, total_factor, aln, model);
        }
        // 4. e1.type = A/C/G/T
        // 4.1. e1.type =  e2.type
        else if (seq1_region->type == seq2_region->type)
        {
            calculateSubtreeCost_identicalACGT(*seq1_region, total_blength, lh_cost, model);
        }
        // e1.type = A/C/G/T and e2.type = O/A/C/G/T
        // 4.2. e1.type = A/C/G/T and e2.type = O
        else if (seq2_region->type == TYPE_O)
        {
            calculateSubtreeCost_ACGT_O<num_states>(*seq1_region, *seq2_region, total_blength, total_factor, model);
        }
        // 4.3. e1.type = A/C/G/T and e2.type = R or A/C/G/T
        else
        {
            if (!calculateSubtreeCost_ACGT_RACGT(*seq1_region, *seq2_region, total_blength, end_pos, total_factor, aln, model)) return MIN_NEGATIVE;
        }
         
        // avoid underflow on total_factor
        // approximately update lh_cost and total_factor
        if (total_factor <= MIN_CARRY_OVER)
        {
            if (total_factor < MIN_POSITIVE)
                return MIN_NEGATIVE;
            
            //lh_cost += log(total_factor);
            //total_factor = 1.0;
            total_factor *= MAX_POSITIVE;
            lh_cost -= LOG_MAX_POSITIVE;
        }
        
        // update pos
        pos = end_pos + 1;
    }
    
    return lh_cost + log(total_factor);
}

void calculateSampleCost_R_R(const SeqRegion& seq1_region, const RealNumType* const &cumulative_rate, const RealNumType blength, const PositionType pos, const PositionType end_pos, RealNumType& lh_cost)
{
    if (seq1_region.plength_observation2node < 0 && seq1_region.plength_observation2root < 0)
        lh_cost += blength * (cumulative_rate[end_pos + 1] - cumulative_rate[pos]);
    else
    {
        RealNumType total_blength = blength + seq1_region.plength_observation2node;
        if (seq1_region.plength_observation2root < 0)
            lh_cost += total_blength * (cumulative_rate[end_pos + 1] - cumulative_rate[pos]);
        else
            // here contribution from root frequency gets added and subtracted so it's ignored
            lh_cost += (total_blength + seq1_region.plength_observation2root) * (cumulative_rate[end_pos + 1] - cumulative_rate[pos]);
    }
}

template <const StateType num_states>
void calculateSampleCost_R_O(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType blength, const StateType seq1_state, RealNumType& lh_cost, RealNumType& total_factor, const Model& model)
{
    if (seq1_region.plength_observation2root >= 0)
    {
        RealNumType total_blength = seq1_region.plength_observation2root + blength;
        
        if (seq2_region.getLH(seq1_state) > 0.1)
        {
            total_blength += seq1_region.plength_observation2node;
            
            // here contribution from root frequency can also be also ignored
            lh_cost += model.diagonal_mut_mat[seq1_state] * total_blength;
        }
        else
        {
            RealNumType tot = 0;
            RealNumType* freq_j_transposed_ij_row = model.freq_j_transposed_ij + model.row_index[seq1_state];
            RealNumType* mutation_mat_row = model.mutation_mat;
                                
            for (StateType i = 0; i < num_states; ++i, mutation_mat_row += num_states)
            {
                RealNumType tot2 = freq_j_transposed_ij_row[i] * seq1_region.plength_observation2node + ((seq1_state == i) ? model.root_freqs[i] : 0);
                RealNumType tot3 = ((seq2_region.getLH(i) > 0.1) ? 1 : 0) + sumMutationByLh<num_states>(&(*seq2_region.likelihood)[0], mutation_mat_row);
                
                tot += tot2 * tot3 * total_blength;
            }
            
            total_factor *= tot * model.inverse_root_freqs[seq1_state];
        }
    }
    else
    {
        if (seq2_region.getLH(seq1_state) > 0.1)
        {
            if (seq1_region.plength_observation2node >= 0)
                lh_cost += model.diagonal_mut_mat[seq1_state] * (blength + seq1_region.plength_observation2node);
            else
                lh_cost += model.diagonal_mut_mat[seq1_state] * blength;
        }
        else
        {
            RealNumType tot = 0;
            RealNumType* mutation_mat_row = model.mutation_mat + model.row_index[seq1_state];
            
            tot += sumMutationByLh<num_states>(&(*seq2_region.likelihood)[0], mutation_mat_row);
            
            if (seq1_region.plength_observation2node >= 0)
                total_factor *= tot * (blength + seq1_region.plength_observation2node);
            else
                total_factor *= tot * blength;
        }
    }
}

void calculateSampleCost_R_ACGT(const SeqRegion& seq1_region, const RealNumType blength, const StateType seq1_state, const StateType seq2_state, RealNumType& total_factor, const Model& model)
{
    if (seq1_region.plength_observation2root >= 0)
    {
        // TODO: can cache model.mutation_mat[model.row_index[seq1_state] * model.diagonal_mut_mat[seq1_state]
        // TODO: can cache  model.freqi_freqj_qij[model.row_index[seq2_state] + seq1_state] * model.diagonal_mut_mat[seq2_state]
        RealNumType seq1_state_evolves_seq2_state = model.mutation_mat[model.row_index[seq1_state] + seq2_state] * blength * (1.0 + model.diagonal_mut_mat[seq1_state] * seq1_region.plength_observation2node);
        
        RealNumType seq2_state_evolves_seq1_state = model.freqi_freqj_qij[model.row_index[seq2_state] + seq1_state] * seq1_region.plength_observation2node * (1.0 + model.diagonal_mut_mat[seq2_state] * (blength + seq1_region.plength_observation2root));
                                                                                                                                                                                    
        total_factor *= seq1_state_evolves_seq2_state + seq2_state_evolves_seq1_state;
    }
    else
    {
        total_factor *= model.mutation_mat[model.row_index[seq1_state] + seq2_state] * (blength + (seq1_region.plength_observation2node < 0 ? 0 : seq1_region.plength_observation2node));
    }
}

template <const StateType num_states>
void calculateSampleCost_O_O(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType blength, RealNumType& total_factor, const Model& model)
{
    RealNumType blength13 = blength;
    if (seq1_region.plength_observation2node >= 0)
    {
        blength13 = seq1_region.plength_observation2node;
        if (blength > 0)
            blength13 += blength;
    }
    
    RealNumType tot = 0;
    
    RealNumType* mutation_mat_row = model.mutation_mat;
                        
    for (StateType i = 0; i < num_states; ++i, mutation_mat_row += num_states)
    {
        RealNumType tot2 = blength13 * sumMutationByLh<num_states>(&(*seq2_region.likelihood)[0], mutation_mat_row);
        
        tot += (tot2 + (seq2_region.getLH(i) > 0.1 ? 1 : 0)) * seq1_region.getLH(i);
    }
        
    total_factor *= tot;
}

template <const StateType num_states>
void calculateSampleCost_O_RACGT(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType blength, const PositionType end_pos, RealNumType& total_factor, const Alignment& aln, const Model& model)
{
    RealNumType blength13 = blength;
    if (seq1_region.plength_observation2node >= 0)
    {
        blength13 = seq1_region.plength_observation2node;
        if (blength > 0)
            blength13 += blength;
    }
    
    StateType seq2_state = seq2_region.type;
    if (seq2_state == TYPE_R)
        seq2_state = aln.ref_seq[end_pos];
    
    RealNumType *transposed_mut_mat_row = model.transposed_mut_mat + model.row_index[seq2_state];
    RealNumType tot2 = dotProduct<num_states>(transposed_mut_mat_row, &((*seq1_region.likelihood)[0]));
    total_factor *= seq1_region.getLH(seq2_state) + blength13 * tot2;
}

void calculateSampleCost_identicalACGT(const SeqRegion& seq1_region, const RealNumType blength, RealNumType& lh_cost, const Model& model)
{
    RealNumType total_blength = blength;
    total_blength += (seq1_region.plength_observation2node < 0 ? 0 : seq1_region.plength_observation2node);
    total_blength += (seq1_region.plength_observation2root < 0 ? 0 : seq1_region.plength_observation2root);

    lh_cost += model.diagonal_mut_mat[seq1_region.type] * total_blength;
}

template <const StateType num_states>
void calculateSampleCost_ACGT_O(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType blength, RealNumType& lh_cost, RealNumType& total_factor, const Model& model)
{
    StateType seq1_state = seq1_region.type;
    RealNumType tot = 0.0;
    
    if (seq1_region.plength_observation2root >= 0)
    {
        RealNumType blength15 = blength + seq1_region.plength_observation2root;
        
        if (seq2_region.getLH(seq1_state) > 0.1)
            lh_cost += model.diagonal_mut_mat[seq1_state] * (blength15 + seq1_region.plength_observation2node);
        else
        {
            RealNumType* freq_j_transposed_ij_row = model.freq_j_transposed_ij + model.row_index[seq1_state];
            RealNumType* mutation_mat_row = model.mutation_mat;
                                
            for (StateType i = 0; i < num_states; ++i, mutation_mat_row += num_states)
            {
                RealNumType tot2 = freq_j_transposed_ij_row[i] * seq1_region.plength_observation2node + ((seq1_state == i) ? model.root_freqs[i] : 0);
                    
                RealNumType tot3 = sumMutationByLh<num_states>(&(*seq2_region.likelihood)[0], mutation_mat_row);
                
                tot += tot2 * blength15 * tot3 + (seq2_region.getLH(i) > 0.1 ? tot2 : 0);
            }
            
            total_factor *= (tot * model.inverse_root_freqs[seq1_state]);
        }
    }
    else
    {
        RealNumType tmp_blength = blength + (seq1_region.plength_observation2node < 0 ? 0 : seq1_region.plength_observation2node);
        if (seq2_region.getLH(seq1_state) > 0.1)
            lh_cost += model.diagonal_mut_mat[seq1_state] * tmp_blength;
        else
        {
            RealNumType* mutation_mat_row = model.mutation_mat + model.row_index[seq1_state];
            tot += sumMutationByLh<num_states>(&(*seq2_region.likelihood)[0], mutation_mat_row);
            
            total_factor *= tot * tmp_blength;
        }
    }
}

void calculateSampleCost_ACGT_RACGT(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType blength, const PositionType end_pos, RealNumType& total_factor, const Alignment& aln, const Model& model)
{
    StateType seq1_state = seq1_region.type;
    StateType seq2_state = seq2_region.type;
    if (seq2_state == TYPE_R)
        seq2_state = aln.ref_seq[end_pos];
    
    if (seq1_region.plength_observation2root >= 0)
    {
        // here we ignore contribution of non-parsimonious mutational histories
        RealNumType seq1_state_evoloves_seq2_state = model.mutation_mat[model.row_index[seq1_state] + seq2_state] * (blength + seq1_region.plength_observation2root) * (1.0 + model.diagonal_mut_mat[seq1_state] * seq1_region.plength_observation2node);
        
        RealNumType seq2_state_evolves_seq1_state = model.freqi_freqj_qij[model.row_index[seq2_state] + seq1_state] * seq1_region.plength_observation2node * (1.0 + model.diagonal_mut_mat[seq2_state] * (blength + seq1_region.plength_observation2root));
        
        total_factor *= (seq1_state_evoloves_seq2_state + seq2_state_evolves_seq1_state);
    }
    else
    {
        RealNumType tmp_blength = ((seq1_region.plength_observation2node < 0) ? blength : blength + seq1_region.plength_observation2node);
        
        total_factor *= model.mutation_mat[model.row_index[seq1_state] + seq2_state] * tmp_blength;
    }
}

RealNumType Tree::calculateSamplePlacementCost(const std::unique_ptr<SeqRegions>& parent_regions, const std::unique_ptr<SeqRegions>& child_regions, const RealNumType blength)
{
    return (this->*calculateSamplePlacementCostPointer)(parent_regions, child_regions, blength);
}

// this implementation derives from appendProb
template <const StateType num_states>
RealNumType Tree::calculateSamplePlacementCostTemplate(const std::unique_ptr<SeqRegions>& parent_regions, const std::unique_ptr<SeqRegions>& child_regions, const RealNumType input_blength)
{
    // 10% of total runtime
    // init dummy variables
    RealNumType lh_cost = 0;
    PositionType pos = 0;
    RealNumType total_factor = 1;
    const SeqRegions& seq1_regions = *parent_regions;
    const SeqRegions& seq2_regions = *child_regions;
    size_t iseq1 = 0;
    size_t iseq2 = 0;
    RealNumType blength = input_blength;
    if (blength < 0) blength = 0;
    const PositionType seq_length = aln.ref_seq.size();
    
    while (pos < seq_length)
    {
        PositionType end_pos;
        RealNumType total_blength = blength;
        
        // get the next shared segment in the two sequences
        SeqRegions::getNextSharedSegment(pos, seq1_regions, seq2_regions, iseq1, iseq2, end_pos);
        const auto* seq1_region = &seq1_regions[iseq1];
        const auto* seq2_region = &seq2_regions[iseq2]; 
        
        // 1. e1.type = N || e2.type = N
        if ((seq2_region->type == TYPE_N) + (seq1_region->type == TYPE_N))
        {
            pos = end_pos + 1;
            continue;
        }
        
        // e1.type != N && e2.type != N
        // A,C,G,T
        // R -> same as the reference
        // N -> gaps
        // O -> vector of 4 probability to observe A, C, G, T
        const DoubleState s1s2 = (DoubleState(seq1_region->type) << 8) | seq2_region->type;
        
        // 2. e1.type = R
        // 2.1. e1.type = R and e2.type = R
        if (s1s2 == RR) [[likely]]
        {
            calculateSampleCost_R_R(*seq1_region, model.cumulative_rate, blength, pos, end_pos, lh_cost);
        }
        // 2.2. e1.type = R and e2.type = O
        else if (s1s2 == RO)
        {
            calculateSampleCost_R_O<num_states>(*seq1_region, *seq2_region, blength, aln.ref_seq[end_pos], lh_cost, total_factor, model);
        }
        // 2.3. e1.type = R and e2.type = A/C/G/T
        else if (seq1_region->type == TYPE_R)
        {
            calculateSampleCost_R_ACGT(*seq1_region, blength, aln.ref_seq[end_pos], seq2_region->type, total_factor, model);
        }
        // 3. e1.type = O
        // 3.1. e1.type = O and e2.type = O
        else if (s1s2 == OO)
        {
            calculateSampleCost_O_O<4>(*seq1_region, *seq2_region, blength, total_factor, model);
        }
        // 3.2. e1.type = O and e2.type = R or A/C/G/T
        else if (seq1_region->type == TYPE_O)
        {
            calculateSampleCost_O_RACGT<num_states>(*seq1_region, *seq2_region, blength, end_pos, total_factor, aln, model);
        }
        // 4. e1.type = A/C/G/T
        // 4.1. e1.type =  e2.type
        else if (seq1_region->type == seq2_region->type)
        {
            calculateSampleCost_identicalACGT(*seq1_region, blength, lh_cost, model);
        }
        // e1.type = A/C/G/T and e2.type = O/A/C/G/T
        // 4.2. e1.type = A/C/G/T and e2.type = O
        else if (seq2_region->type == TYPE_O)
        {
            calculateSampleCost_ACGT_O<4>(*seq1_region, *seq2_region, blength, lh_cost, total_factor, model);
        }
        // 4.3. e1.type = A/C/G/T and e2.type = R or A/C/G/T
        else
        {
            calculateSampleCost_ACGT_RACGT(*seq1_region, *seq2_region, blength, end_pos, total_factor, aln, model);
        }
         
        // avoid underflow on total_factor
        // approximately update lh_cost and total_factor
        if (total_factor <= MIN_CARRY_OVER)
        {
            if (total_factor < MIN_POSITIVE)
                return MIN_NEGATIVE;
            
            //lh_cost += log(total_factor);
            //total_factor = 1.0;
            total_factor *= MAX_POSITIVE;
            lh_cost -= LOG_MAX_POSITIVE;
        }
        
        // update pos
        pos = end_pos + 1;
    }
    
    return lh_cost + log(total_factor);
}

void Tree::updateZeroBlength(const Index index, PhyloNode& node, std::stack<Index> &node_stack)
{
    // get the top node in the phylo-node
    /*Node* top_node = node->getTopNode();
    ASSERT(top_node);
    SeqRegions* upper_left_right_regions = top_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
    SeqRegions* lower_regions = top_node->getPartialLhAtNode(aln, model, threshold_prob);*/
    const NumSeqsType node_vec_index = index.getVectorIndex();
    const Index parent_index = node.getNeighborIndex(TOP);
    PhyloNode& parent_node = nodes[parent_index.getVectorIndex()];
    const std::unique_ptr<SeqRegions>& upper_left_right_regions = parent_node.getPartialLh(parent_index.getMiniIndex());
    const std::unique_ptr<SeqRegions>& lower_regions = node.getPartialLh(TOP);
    
    RealNumType best_lh = calculateSamplePlacementCost(upper_left_right_regions, lower_regions, default_blength);
    RealNumType best_length = default_blength;
    
    // try shorter lengths
    bool found_new_best_length = tryShorterNewBranch<&Tree::calculateSamplePlacementCost>(upper_left_right_regions, lower_regions, best_length, best_lh, min_blength);
    
    // try longer lengths
    if (!found_new_best_length)
        tryLongerNewBranch<&Tree::calculateSamplePlacementCost>(upper_left_right_regions, lower_regions, best_length, best_lh, max_blength);
    
    // update best_length
    /*top_node->length = best_length;
    top_node->neighbor->length = best_length;*/
    node.setUpperLength(best_length);
    
    // add current node and its parent to node_stack to for updating partials further from these nodes
    /*top_node->outdated = true;
    top_node->neighbor->getTopNode()->outdated = true;
    node_stack.push(top_node);
    node_stack.push(top_node->neighbor);*/
    node.setOutdated(true);
    parent_node.setOutdated(true);
    node_stack.push(Index(node_vec_index, TOP));
    node_stack.push(parent_index);
}

void Tree::createAnInternalNode()
{
    nodes.push_back(PhyloNode());
}

void Tree::createALeafNode(const NumSeqsType new_seq_name_index)
{
    nodes.push_back(PhyloNode(std::move(LeafNode(new_seq_name_index))));
}

std::unique_ptr<SeqRegions>& Tree::getPartialLhAtNode(const Index index)
{
    // may need assert(index.getVectorIndex() < nodes.size());
    return nodes[index.getVectorIndex()].getPartialLh(index.getMiniIndex());
}

RealNumType Tree::calculateTreeLh()
{
    // initialize the total_lh by the likelihood from root
    RealNumType total_lh = nodes[root_vector_index].getPartialLh(TOP)->computeAbsoluteLhAtRoot(aln.num_states, model);
    
    // perform a DFS to add likelihood contributions from each internal nodes
    total_lh += performDFS<&Tree::computeLhContribution>();
   
    // return total_lh
    return total_lh;
}

void Tree::updateLowerLh(RealNumType& total_lh, std::unique_ptr<SeqRegions>& new_lower_lh, PhyloNode& node, const std::unique_ptr<SeqRegions>& lower_lh_1, const std::unique_ptr<SeqRegions>& lower_lh_2, const Index neighbor_1_index, PhyloNode& neighbor_1, const Index neighbor_2_index, PhyloNode& neighbor_2, const PositionType& seq_length)
{
    lower_lh_1->mergeTwoLowers(new_lower_lh, neighbor_1.getUpperLength(), *lower_lh_2, neighbor_2.getUpperLength(), aln, model, params->threshold_prob);
    
    // NHANLT: LOGS FOR DEBUGGING
    /*if (params->debug)
        std::cout << "new_lower_lh " << new_lower_lh->size() << std::endl;*/
     
    // if new_lower_lh is NULL -> we need to update the branch lengths connecting the current node to its children
    if (!new_lower_lh)
    {
        if (neighbor_1.getUpperLength() <= 0) // next_node_1->length <= 0)
        {
            stack<Index> node_stack;
            // NHANLT: note different from original maple
            // updateBLen(nodeList,node,mutMatrix) -> the below codes update from next_node_1 instead of top_node
            updateZeroBlength(neighbor_1_index, neighbor_1, node_stack); // next_node_1->neighbor, node_stack, params->threshold_prob);
            updatePartialLh(node_stack);
        }
        else if (neighbor_2.getUpperLength() <= 0) // next_node_2->length <= 0)
        {
            stack<Index> node_stack;
            updateZeroBlength(neighbor_2_index, neighbor_2, node_stack); // updateZeroBlength(next_node_2->neighbor, node_stack, params->threshold_prob);
            updatePartialLh(node_stack);
        }
        else
            outError("Strange, branch lengths > 0 but inconsistent lower lh creation in refreshAllLowerLhs()");
    }
    // otherwise, everything is good -> update the lower lh of the current node
    else
        node.setPartialLh(TOP, std::move(new_lower_lh));
}

void Tree::computeLhContribution(RealNumType& total_lh, std::unique_ptr<SeqRegions>& new_lower_lh, PhyloNode& node, const std::unique_ptr<SeqRegions>& lower_lh_1, const std::unique_ptr<SeqRegions>& lower_lh_2, const Index neighbor_1_index, PhyloNode& neighbor_1, const Index neighbor_2_index, PhyloNode& neighbor_2, const PositionType& seq_length)
{
    total_lh += lower_lh_1->mergeTwoLowers(new_lower_lh, neighbor_1.getUpperLength(), *lower_lh_2, neighbor_2.getUpperLength(), aln, model, params->threshold_prob, true);
    
    // if new_lower_lh is NULL
    if (!new_lower_lh)
        outError("Strange, inconsistent lower genome list creation in calculateTreeLh(); old list, and children lists");
    // otherwise, everything is good -> update the lower lh of the current node
    else if (new_lower_lh->areDiffFrom(*node.getPartialLh(TOP), seq_length, aln.num_states, &params.value()))
        outError("Strange, while calculating tree likelihood encountered non-updated lower likelihood!");
}

template <void(Tree::*task)(RealNumType&, std::unique_ptr<SeqRegions>&, PhyloNode&, const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, const Index, PhyloNode&, const Index, PhyloNode&, const PositionType&)>
RealNumType Tree::performDFS()
{
    // dummy variables
    RealNumType total_lh = 0;
    const PositionType seq_length = aln.ref_seq.size();
    
    // start from root
    Index node_index = Index(root_vector_index, TOP);
    Index last_node_index;
    
    // traverse to the deepest tip, calculate the likelihoods upward from the tips
    while (node_index.getMiniIndex() != UNDEFINED) //node)
    {
        PhyloNode& node = nodes[node_index.getVectorIndex()];
        // we reach a top node by a downward traversing
        if (node_index.getMiniIndex() == TOP) //node->is_top)
        {
            // if the current node is a leaf -> we reach the deepest tip -> traversing upward to calculate the lh of its parent
            if (!node.isInternal()) // node->isLeave())
            {
                /*last_node = node;
                 node = node->neighbor;*/
                last_node_index = node_index;
                node_index = node.getNeighborIndex(TOP);
            }
            // otherwise, keep traversing downward to find the deepest tip
            else
                // node = node->next->neighbor;
                node_index = node.getNeighborIndex(RIGHT);
        }
        // we reach the current node by an upward traversing from its children
        else
        {
            // if we reach the current node by an upward traversing from its first children -> traversing downward to its second children
            if (node.getNeighborIndex(RIGHT) == last_node_index) // node->getTopNode()->next->neighbor == last_node)
            {
                // node = node->getTopNode()->next->next->neighbor;
                node_index = node.getNeighborIndex(LEFT);
            }
            // otherwise, all children of the current node are updated -> update the lower lh of the current node
            else
            {
                // calculate the new lower lh of the current node from its children
                /*Node* top_node = node->getTopNode();
                 Node* next_node_1 = top_node->next;
                 Node* next_node_2 = next_node_1->next;*/
                Index neighbor_1_index = node.getNeighborIndex(RIGHT);
                Index neighbor_2_index = node.getNeighborIndex(LEFT);
                PhyloNode& neighbor_1 = nodes[neighbor_1_index.getVectorIndex()];
                PhyloNode& neighbor_2 = nodes[neighbor_2_index.getVectorIndex()];
                
                std::unique_ptr<SeqRegions> new_lower_lh = nullptr;
                const std::unique_ptr<SeqRegions>& lower_lh_1 = neighbor_1.getPartialLh(TOP); // next_node_1->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob);
                const std::unique_ptr<SeqRegions>& lower_lh_2 = neighbor_2.getPartialLh(TOP); // next_node_2->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob);
                // lower_lh_1->mergeTwoLowers(new_lower_lh, next_node_1->length, *lower_lh_2, next_node_2->length, aln, model, params->threshold_prob);
                
                (this->*task)(total_lh, new_lower_lh, node, lower_lh_1, lower_lh_2, neighbor_1_index, neighbor_1, neighbor_2_index, neighbor_2, seq_length);
                
                last_node_index = Index(node_index.getVectorIndex(), TOP);
                node_index = node.getNeighborIndex(TOP);
            }
        }
    }
    
    return total_lh;
}
