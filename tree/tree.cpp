#include "tree.h"

#include <cassert>
#include <utils/matrix.h>

using namespace std;

Tree::Tree(Params&& n_params, Node* n_root)
{
    params = std::move(n_params);
    root = n_root;
}

Tree::~Tree()
{
    // browse tree to delete all nodes
    stack<Node*> node_stack;
    node_stack.push(root);
    
    while (!node_stack.empty())
    {
        Node* node = node_stack.top();
        node_stack.pop();
        
        // add neighbors to the stack
        stack<Node*> delete_stack;
        delete_stack.push(node);
        Node* next_node;
        FOR_NEXT(node, next_node)
        {
            // add neighbor if it's existed
            if (next_node->neighbor) node_stack.push(next_node->neighbor);
            
            // push next_node into delete_stack
            delete_stack.push(next_node);
        }
        
        // delete all nodes in delete_stack
        while (!delete_stack.empty())
        {
            Node* delete_node = delete_stack.top();
            delete_stack.pop();
            
            // delete node
            delete delete_node;
        }
    }
    
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

string Tree::exportTreeString(bool binary, Node* node)
{
    // init starting node from root
    if (!node)
        node = root;
    
    // move to its neighbor
    if (node->neighbor)
        node = node->neighbor;
    
    // do something with its neighbor
    if (node->isLeave())
        return node->exportString(binary);
        
    string output = "(";
    bool add_comma = false;
    Node* next;
    FOR_NEXT(node, next)
    {
        if (!add_comma)
            add_comma = true;
        else
            output += ",";
        output += exportTreeString(binary, next);
    }
    string length = node->length < 0 ? "0" : convertDoubleToString(node->length);
    output += "):" + length;
    
    return output;
}

void Tree::updatePartialLh(stack<Node*> &node_stack, RealNumType* cumulative_rate, RealNumType default_blength, RealNumType min_blength, RealNumType max_blength)
{
    (this->*updatePartialLhPointer)(node_stack, cumulative_rate, default_blength, min_blength, max_blength);
}

template <const StateType num_states>
void Tree::updatePartialLhTemplate(stack<Node*> &node_stack, RealNumType* cumulative_rate, RealNumType default_blength, RealNumType min_blength, RealNumType max_blength)
{
    PositionType seq_length = aln.ref_seq.size();
    
    while (!node_stack.empty())
    {
        Node* node = node_stack.top();
        node_stack.pop();
        
        // NHANLT: debug
        /*if (node->next && (node->neighbor->seq_name == "25" || (node->next->neighbor && node->next->neighbor->seq_name == "25") || (node->next->next->neighbor && node->next->next->neighbor->seq_name == "25")))
        //if (node->seq_name == "25")
            cout << "dsdas";*/
        
        bool update_blength = false;
        node->getTopNode()->outdated = true;
        
        SeqRegions* parent_upper_regions = NULL;
        bool is_non_root = root != node->getTopNode();
        if (is_non_root)
            parent_upper_regions = node->getTopNode()->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
            
        // change in likelihoods is coming from parent node
        if (node->is_top)
        {
            ASSERT(is_non_root);
            
            // if necessary, update the total probabilities at the mid node.
            if (node->length > 0)
            {
                // update vector of regions at mid-branch point
                SeqRegions* mid_branch_regions = NULL;
                SeqRegions* lower_regions = node->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
                RealNumType half_branch_length = node->length * 0.5;
                parent_upper_regions->mergeUpperLower(mid_branch_regions, half_branch_length, *lower_regions, half_branch_length, aln, model, params->threshold_prob);
                
                if (!mid_branch_regions)
                {
                    if (node->length > 1e-100)
                        outError("inside updatePartialLh(), from parent: should not have happened since node->length > 0");
                    updateZeroBlength(node, node_stack, params->threshold_prob, cumulative_rate, default_blength, min_blength, max_blength);
                    update_blength = true;
                }
                else
                {
                    // update likelihood at the mid-branch point
                    if (node->mid_branch_lh) delete node->mid_branch_lh;
                    node->mid_branch_lh = mid_branch_regions;
                    mid_branch_regions = NULL;

                    /*if (node.dist>=2*minBLenForMidNode_
                        createFurtherMidNodes(node,parent_upper_regions)*/
                }
                
                // delete mid_branch_regions
                if (mid_branch_regions) delete mid_branch_regions;
                
                // if necessary, update the total probability vector.
                if (!update_blength)
                {
                    node->computeTotalLhAtNode(aln, model, params->threshold_prob, cumulative_rate, node == root);
                    if (!node->total_lh || node->total_lh->size() == 0)
                    {
                        outError("inside updatePartialLh(), from parent 2: should not have happened since node->length > 0");
                        
                        /*updateZeroBlength(nodeList,node,mutMatrix)
                        update_blength=True
                        exit()*/
                    }
                }
            }
            
            // at valid internal node, update upLeft and upRight, and if necessary add children to node_stack.
            if (node->next && !update_blength)
            {
                Node* next_node_1 = node->next;
                Node* next_node_2 = next_node_1->next;
                
                SeqRegions* upper_left_right_regions_1 = NULL;
                SeqRegions* upper_left_right_regions_2 = NULL;
                SeqRegions* lower_regions_1 = next_node_1->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
                parent_upper_regions->mergeUpperLower(upper_left_right_regions_1, node->length, *lower_regions_1, next_node_1->length, aln, model, params->threshold_prob);
                
                if (!upper_left_right_regions_1 || upper_left_right_regions_1->size() == 0)
                {
                    if (node->length <= 0 && next_node_1->length <= 0)
                        updateZeroBlength(node, node_stack, params->threshold_prob, cumulative_rate, default_blength, min_blength, max_blength);
                    else
                        outError("Strange: None vector from non-zero distances in updatePartialLh() from parent direction.");
                    
                    update_blength = true;
                }
                
                if (!update_blength)
                {
                    SeqRegions* lower_regions_2 = next_node_2->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
                    parent_upper_regions->mergeUpperLower(upper_left_right_regions_2, node->length, *lower_regions_2, next_node_2->length, aln, model, params->threshold_prob);
                    
                    if (!upper_left_right_regions_2 || upper_left_right_regions_2->size() == 0)
                    {
                        if (node->length <= 0 && next_node_2->length <= 0)
                        {
                            updateZeroBlength(node, node_stack, params->threshold_prob, cumulative_rate, default_blength, min_blength, max_blength);
                            update_blength = true;
                        }
                        else
                            outError("Strange: None vector from non-zero distances in updatePartialLh() from parent direction, child0.");
                    }
                }
                
                if (!update_blength)
                {
                    if (next_node_1->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate)->areDiffFrom(*upper_left_right_regions_2, seq_length, num_states, &params.value()))
                    {
                        delete next_node_1->partial_lh;
                        next_node_1->partial_lh = upper_left_right_regions_2;
                        upper_left_right_regions_2 = NULL;
                        node_stack.push(next_node_1->neighbor);
                    }
                    
                    if (next_node_2->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate)->areDiffFrom(*upper_left_right_regions_1, seq_length, num_states, &params.value()))
                    {
                        delete next_node_2->partial_lh;
                        next_node_2->partial_lh = upper_left_right_regions_1;
                        upper_left_right_regions_1 = NULL;
                        
                        node_stack.push(next_node_2->neighbor);
                    }
                }
                
                // delete upper_left_right_regions_1, upper_left_right_regions_2
                if (upper_left_right_regions_1) delete upper_left_right_regions_1;
                if (upper_left_right_regions_2) delete upper_left_right_regions_2;
            }
        }
        // otherwise, change in likelihoods is coming from a child.
        else
        {
            Node* top_node = NULL;
            Node* other_next_node = NULL;
            Node* next_node = NULL;
            FOR_NEXT(node, next_node)
            {
                if (next_node->is_top)
                    top_node = next_node;
                else
                    other_next_node = next_node;
            }
            
            ASSERT(top_node && other_next_node);
            
            RealNumType this_node_distance = node->length;
            RealNumType other_next_node_distance = other_next_node->length;
            
            SeqRegions* this_node_lower_regions = node->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
            // Regions* this_node_upper_left_right_regions = node->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
            SeqRegions* next_node_upper_left_right_regions = other_next_node->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
            
            // update lower likelihoods
            SeqRegions* merged_two_lower_regions = NULL;
            SeqRegions* old_lower_regions = NULL;
            other_next_node->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate)->mergeTwoLowers(merged_two_lower_regions, other_next_node_distance, this_node_lower_regions, this_node_distance, aln, model, params->threshold_prob, cumulative_rate);
            
            if (!merged_two_lower_regions || merged_two_lower_regions->size() == 0)
            {
                if (this_node_distance <= 0 && other_next_node_distance <= 0)
                {
                    updateZeroBlength(node->neighbor, node_stack, params->threshold_prob, cumulative_rate, default_blength, min_blength, max_blength);
                    update_blength = true;
                }
                else
                {
                    outError("Strange: None vector from non-zero distances in updatePartialLh() from child direction.");
                    /*print(this_node_distance)
                    print(this_node_lower_regions)
                    print(other_next_node_distance)
                    print(otherChildVect)
                    exit()*/
                }
            }
            else
            {
                if (old_lower_regions) delete old_lower_regions;
                old_lower_regions = top_node->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
                top_node->partial_lh = merged_two_lower_regions;
                merged_two_lower_regions = NULL;
            }

            // delete merged_two_lower_regions
            if (merged_two_lower_regions) delete merged_two_lower_regions;
            
            // update total likelihood
            if (!update_blength)
            {
                if (top_node->length > 0 || top_node == root)
                {
                    SeqRegions* new_total_lh_regions = top_node->computeTotalLhAtNode(aln, model, params->threshold_prob, cumulative_rate, top_node == root, false);
                    
                    if (!new_total_lh_regions && top_node->length <= 0)
                    {
                        updateZeroBlength(top_node, node_stack, params->threshold_prob, cumulative_rate, default_blength, min_blength, max_blength);
                        update_blength = true;
                    }
                    else if (!new_total_lh_regions)
                        outError("Strange: None vector from non-zero distances in updatePartialLh() from child direction while doing overall likelihood.");
                    else
                    {
                        // init top_node->total_lh
                        if (top_node->total_lh) delete top_node->total_lh;
                        top_node->total_lh = new_total_lh_regions;
                        new_total_lh_regions = NULL;
                    }
                    
                    // delete new_total_lh_regions
                    if (new_total_lh_regions) delete new_total_lh_regions;
                }
            }
            
            // update total mid-branches likelihood
            if (!update_blength)
            {
                if (top_node->length > 0 && is_non_root)
                {
                    SeqRegions* new_mid_regions = NULL;
                    SeqRegions* tmp_lower_regions = top_node->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
                    RealNumType half_branch_length = top_node->length * 0.5;
                    parent_upper_regions->mergeUpperLower(new_mid_regions, half_branch_length, *tmp_lower_regions, half_branch_length, aln, model, params->threshold_prob);
                    
                    if (!new_mid_regions)
                    {
                        updateZeroBlength(top_node, node_stack, params->threshold_prob, cumulative_rate, default_blength, min_blength, max_blength);
                        update_blength = true;
                        // print("inside updatePartials(), from child: should not have happened since node.dist>0")
                    }
                    else
                    {
                        if (top_node->mid_branch_lh) delete top_node->mid_branch_lh;
                        top_node->mid_branch_lh = new_mid_regions;
                        new_mid_regions = NULL;
                        /*if node.dist>=2*minBLenForMidNode:
                            createFurtherMidNodes(node,parent_upper_regions)*/
                    }
                        
                    // delete new_mid_regions
                    if (new_mid_regions) delete new_mid_regions;
                }
            }
            
            if (!update_blength)
            {
                // update likelihoods at parent node
                if (top_node->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate)->areDiffFrom(*old_lower_regions, seq_length, num_states, &params.value()))
                {
                    if (root != top_node)
                        node_stack.push(top_node->neighbor);
                }

                // update likelihoods at sibling node
                SeqRegions* new_upper_regions = NULL;
                if (is_non_root)
                    parent_upper_regions->mergeUpperLower(new_upper_regions, top_node->length, *this_node_lower_regions, this_node_distance, aln, model, params->threshold_prob);
                else
                    new_upper_regions = node->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate)->computeTotalLhAtRoot(num_states, model, this_node_distance);
                
                if (!new_upper_regions || new_upper_regions->size() == 0)
                {
                    if (top_node->length <= 0 && this_node_distance <= 0)
                    {
                        updateZeroBlength(top_node, node_stack, params->threshold_prob, cumulative_rate, default_blength, min_blength, max_blength);
                        update_blength = true;
                    }
                    else
                        outError("Strange: None vector from non-zero distances in updatePartialLh() from child direction, new_upper_regions.");
                }
                else
                {
                    if (next_node_upper_left_right_regions->areDiffFrom(*new_upper_regions, seq_length, num_states, &params.value()))
                    {
                        if (other_next_node->partial_lh) delete other_next_node->partial_lh;
                        other_next_node->partial_lh = new_upper_regions;
                        new_upper_regions = NULL;
                        
                        node_stack.push(other_next_node->neighbor);
                    }
                }
                    
                    // delete new_upper_regions
                   if (new_upper_regions) delete new_upper_regions;
            }
                    
            // delete old_lower_regions
            if (old_lower_regions) delete old_lower_regions;
        }
    }
}

void Tree::seekSamplePlacement(Node* start_node, const string &seq_name, SeqRegions* sample_regions, Node* &selected_node, RealNumType &best_lh_diff , bool &is_mid_branch, RealNumType &best_up_lh_diff, RealNumType &best_down_lh_diff, Node* &best_child, RealNumType* cumulative_rate, RealNumType default_blength, RealNumType min_blength_mid)
{
    // init variables
    // output variables
    selected_node = start_node;
    // dummy variables
    RealNumType lh_diff_mid_branch = 0;
    RealNumType lh_diff_at_node = 0;
    // stack of nodes to examine positions
    stack<TraversingNode> extended_node_stack;
    extended_node_stack.push(TraversingNode(start_node, 0, MIN_NEGATIVE));
    
    // recursively examine positions for placing the new sample
    while (!extended_node_stack.empty())
    {
        TraversingNode current_extended_node = extended_node_stack.top();
        Node* current_node = current_extended_node.node;
        extended_node_stack.pop();
        
        // NHANLT: debug
        /*if (current_node->next && ((current_node->next->neighbor && current_node->next->neighbor->seq_name == "25")
                                  || (current_node->next->next->neighbor && current_node->next->next->neighbor->seq_name == "25")))
            cout << "fdsfsd";*/
    
        // if the current node is a leaf AND the new sample/sequence is strictly less informative than the current node
        // -> add the new sequence into the list of minor sequences of the current node + stop seeking the placement
        if ((!current_node->next) && (current_node->partial_lh->compareWithSample(*sample_regions, aln.ref_seq.size(), aln.num_states) == 1))
        {
            current_node->less_info_seqs.push_back(seq_name);
            selected_node = NULL;
            return;
        }
        
        // 1. try first placing as a descendant of the mid-branch point of the branch above the current node
        if (current_node != root && current_node->length > 0)
        {
            // compute the placement cost
            lh_diff_mid_branch = calculateSamplePlacementCost(cumulative_rate, current_node->mid_branch_lh, sample_regions, default_blength);
            
            // record the best_lh_diff if lh_diff_mid_branch is greater than the best_lh_diff ever
            if (lh_diff_mid_branch > best_lh_diff)
            {
                best_lh_diff = lh_diff_mid_branch;
                selected_node = current_node;
                current_extended_node.failure_count = 0;
                is_mid_branch = true;
            }
        }
        // otherwise, don't consider mid-branch point
        else
            lh_diff_mid_branch = MIN_NEGATIVE;

        // 2. try to place as descendant of the current node (this is skipped if the node has top branch length 0 and so is part of a polytomy).
        if (current_node == root || current_node->length > 0)
        {
            // compute the placement cost
            lh_diff_at_node = calculateSamplePlacementCost(cumulative_rate, current_node->total_lh, sample_regions, default_blength);
            
            // record the best_lh_diff if lh_diff_at_node is greater than the best_lh_diff ever
            if (lh_diff_at_node > best_lh_diff)
            {
                best_lh_diff = lh_diff_at_node;
                selected_node = current_node;
                current_extended_node.failure_count = 0;
                is_mid_branch = false;
                best_up_lh_diff = lh_diff_mid_branch;
            }
            else if (lh_diff_mid_branch >= (best_lh_diff - params->threshold_prob))
            {
                best_up_lh_diff = current_extended_node.likelihood_diff;
                best_down_lh_diff = lh_diff_at_node;
                best_child = current_node;
            }
            // placement at current node is considered failed if placement likelihood is not improved by a certain margin compared to best placement so far for the nodes above it.
            else if (lh_diff_at_node < (current_extended_node.likelihood_diff - params->thresh_log_lh_failure))
                ++current_extended_node.failure_count;
        }
        else
            lh_diff_at_node = current_extended_node.likelihood_diff;
        
        // keep trying to place at children nodes, unless the number of attempts has reaches the failure limit
        if ((params->strict_stop_seeking_placement_sample
             && current_extended_node.failure_count < params->failure_limit_sample
             && lh_diff_at_node > (best_lh_diff - params->thresh_log_lh_sample))
            || (!params->strict_stop_seeking_placement_sample
                && (current_extended_node.failure_count < params->failure_limit_sample
                    || lh_diff_at_node > (best_lh_diff - params->thresh_log_lh_sample))))
        {
            Node* neighbor_node;
            FOR_NEIGHBOR(current_node, neighbor_node)
            extended_node_stack.push(TraversingNode(neighbor_node, current_extended_node.failure_count, lh_diff_at_node));
        }
    }

    // exploration of the tree is finished, and we are left with the node found so far with the best appending likelihood cost. Now we explore placement just below this node for more fine-grained placement within its descendant branches.
    best_down_lh_diff = MIN_NEGATIVE;
    best_child = NULL;
    
    // if best position so far is the descendant of a node -> explore further at its children
    if (!is_mid_branch)
    {
        // current node might be part of a polytomy (represented by 0 branch lengths) so we want to explore all the children of the current node to find out if the best placement is actually in any of the branches below the current node.
        Node* neighbor_node;
        stack<Node*> node_stack;
        FOR_NEIGHBOR(selected_node, neighbor_node)
            node_stack.push(neighbor_node);
        
        while (!node_stack.empty())
        {
            Node* node = node_stack.top();
            node_stack.pop();

            if (node->length <= 0)
            {
                FOR_NEIGHBOR(node, neighbor_node)
                    node_stack.push(neighbor_node);
            }
            else
            {
                // now try to place on the current branch below the best node, at an height above the mid-branch.
                RealNumType new_blength = node->length * 0.5;
                RealNumType new_best_lh_mid_branch = MIN_NEGATIVE;
                SeqRegions* upper_lr_regions = node->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
                SeqRegions* lower_regions = node->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
                SeqRegions* mid_branch_regions = new SeqRegions(node->mid_branch_lh);

                // try to place new sample along the upper half of the current branch
                while (true)
                {
                    // compute the placement cost
                    RealNumType new_lh_mid_branch = calculateSamplePlacementCost(cumulative_rate, mid_branch_regions, sample_regions, default_blength);
                    
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
                    upper_lr_regions->mergeUpperLower(mid_branch_regions, new_blength, *lower_regions, node->length - new_blength, aln, model, params->threshold_prob);
                }
                
                // delete mid_branch_regions
                delete mid_branch_regions;
                
                //RealNumType new_best_lh_mid_branch = calculateSamplePlacementCost(cumulative_rate, node->mid_branch_lh, sample_regions, default_blength);
                
                // record new best_down_lh_diff
                if (new_best_lh_mid_branch > best_down_lh_diff)
                {
                    best_down_lh_diff = new_best_lh_mid_branch;
                    best_child = node;
                }
            }
        }
    }
}

void Tree::seekSubTreePlacement(Node* &best_node, RealNumType &best_lh_diff, bool &is_mid_branch, RealNumType &best_up_lh_diff, RealNumType &best_down_lh_diff, Node* &best_child, bool short_range_search, Node* child_node, RealNumType &removed_blength, RealNumType* cumulative_rate, RealNumType default_blength, RealNumType min_blength_mid, bool search_subtree_placement, SeqRegions* sample_regions)
{
    // init variables
    Node* node = child_node->neighbor->getTopNode();
    Node* other_child_node = child_node->neighbor->getOtherNextNode()->neighbor;
    best_node = node;
    SeqRegions* subtree_regions = NULL;
    // stack of nodes to examine positions
    stack<UpdatingNode*> node_stack;
    // dummy variables
    RealNumType threshold_prob = params->threshold_prob;
    RealNumType lh_diff_mid_branch = 0;
    RealNumType lh_diff_at_node = 0;
    SeqRegions* parent_upper_lr_regions = NULL;
    PositionType seq_length = aln.ref_seq.size();
    StateType num_states = aln.num_states;
    
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
    /*if (search_subtree_placement)
    {*/
    // get the lower regions of the child node
    subtree_regions = child_node->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
    
    // node is not the root
    if (node != root)
    {
        parent_upper_lr_regions = node->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
        SeqRegions* other_child_node_regions = other_child_node->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
       
        // add nodes (sibling and parent of the current node) into node_stack which we will need to traverse to update their regions due to the removal of the sub tree
        RealNumType branch_length = other_child_node->length;
        if (node->length > 0)
            branch_length = branch_length > 0 ? branch_length + node->length : node->length;
        
        node_stack.push(new UpdatingNode(node->neighbor, other_child_node_regions, branch_length, true, best_lh_diff, 0, false));
        node_stack.push(new UpdatingNode(other_child_node, parent_upper_lr_regions, branch_length, true, best_lh_diff, 0, false));
    }
    // node is root
    else
    {
        // there is only one sample outside of the subtree doesn't need to be considered
        if (other_child_node->next)
        {
            // add nodes (children of the sibling of the current node) into node_stack which we will need to traverse to update their regions due to the removal of the sub tree
            Node* grand_child_1 = other_child_node->next->neighbor;
            Node* grand_child_2 = other_child_node->next->next->neighbor;
            
            SeqRegions* up_lr_regions_1 = grand_child_2->computeTotalLhAtNode(aln, model, threshold_prob, cumulative_rate, true, false, grand_child_2->length);
            node_stack.push(new UpdatingNode(grand_child_1, up_lr_regions_1, grand_child_1->length, true, best_lh_diff, 0, true));
            
            SeqRegions* up_lr_regions_2 = grand_child_1->computeTotalLhAtNode(aln, model, threshold_prob, cumulative_rate, true, false, grand_child_1->length);
            node_stack.push(new UpdatingNode(grand_child_2, up_lr_regions_2, grand_child_2->length, true, best_lh_diff, 0, true));
        }
    }
    /*}
     // search a placement for a new sample
    else
    {
        // get the regions of the input sample
        subtree_regions = sample_regions;
        RealNumType down_lh = is_mid_branch ? best_down_lh_diff : best_lh_diff;
        
        // node is not the root
        if (node != root)
        {
            SeqRegions* lower_regions = new SeqRegions(node->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate), num_states);
            parent_upper_lr_regions = node->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
            
            // add the parent node of the current node into node_stack for traversing to seek the placement for the new sample
            node_stack.push(new UpdatingNode(node->neighbor, lower_regions, node->length, false, down_lh, 0));
        }
        
        // add the children nodes of the current node into node_stack for traversing to seek the placement for the new sample
        Node* neighbor_node;
        FOR_NEIGHBOR(node, neighbor_node)
        {
            SeqRegions* upper_lr_regions = new SeqRegions(neighbor_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate), num_states);
            node_stack.push(new UpdatingNode(neighbor_node, upper_lr_regions, neighbor_node->length, false, down_lh, 0));
        }
     }*/
    
    // examine each node in the node stack to seek the "best" placement
    while (!node_stack.empty())
    {
        // extract updating_node from stack
        UpdatingNode* updating_node = node_stack.top();
        node_stack.pop();
        
        // consider the case we are moving from a parent to a child
        if (updating_node->node->is_top)
        {
            if (updating_node->node->length > 0)
            {
                //  try to append mid-branch
                // avoid adding to the old position where the subtree was just been removed from
                if (updating_node->node != root && updating_node->node->neighbor->getTopNode() != node)
                {
                    SeqRegions* mid_branch_regions = NULL;
                    // get or recompute the lh regions at the mid-branch position
                    if (updating_node->need_updating)
                    {
                        SeqRegions* lower_regions = updating_node->node->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
                        RealNumType mid_branch_length = updating_node->branch_length * 0.5;
                        updating_node->incoming_regions->mergeUpperLower(mid_branch_regions, mid_branch_length, *lower_regions, mid_branch_length, aln, model, threshold_prob);
                    }
                    else
                        mid_branch_regions = updating_node->node->mid_branch_lh;
                    
                    // skip if mid_branch_regions is null (branch length == 0)
                    if (!mid_branch_regions)
                    {
                        // delete updating_node
                        delete updating_node;
                        
                        continue;
                    }
                    
                    // compute the placement cost
                    // if (search_subtree_placement)
                    lh_diff_mid_branch = calculateSubTreePlacementCost(cumulative_rate, mid_branch_regions, subtree_regions, removed_blength);
                    /*else
                        lh_diff_mid_branch = calculateSamplePlacementCost(cumulative_rate, mid_branch_regions, subtree_regions, removed_blength);*/
                    
                    // if this position is better than the best position found so far -> record it
                    if (lh_diff_mid_branch > best_lh_diff)
                    {
                        best_node = updating_node->node;
                        best_lh_diff = lh_diff_mid_branch;
                        is_mid_branch = true;
                        updating_node->failure_count = 0;
                    }
                        
                    // delete mid_branch_regions
                    if (updating_node->need_updating) delete mid_branch_regions;
                }
                // set the placement cost at the mid-branch position the most negative value if branch length is zero -> we can't place the subtree on that branch
                else
                    lh_diff_mid_branch = MIN_NEGATIVE;
                    
                // now try appending exactly at node
                SeqRegions* at_node_regions = NULL;
                bool delete_at_node_regions = false;
                if (updating_node->need_updating)
                {
                    delete_at_node_regions = true;
                    // get or recompute the lh regions at the current node position
                    SeqRegions* lower_regions = updating_node->node->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
                    updating_node->incoming_regions->mergeUpperLower(at_node_regions, updating_node->branch_length, *lower_regions, -1, aln, model, threshold_prob);
                    
                    // stop updating if the difference between the new and old regions is insignificant
                    if  (!at_node_regions->areDiffFrom(*updating_node->node->total_lh, seq_length, num_states, &params.value()))
                        updating_node->need_updating = false;
                }
                else
                    at_node_regions = updating_node->node->total_lh;

                // skip if at_node_regions is null (branch length == 0)
                if (!at_node_regions)
                {
                    // delete updating_node
                    delete updating_node;
                    
                    continue;
                }
                
                //if (search_subtree_placement)
                lh_diff_at_node = calculateSubTreePlacementCost(cumulative_rate, at_node_regions, subtree_regions, removed_blength);
                /*else
                    lh_diff_at_node = calculateSamplePlacementCost(cumulative_rate, at_node_regions, subtree_regions, removed_blength);*/
                
                // if this position is better than the best position found so far -> record it
                if (lh_diff_at_node > best_lh_diff)
                {
                    best_node = updating_node->node;
                    best_lh_diff = lh_diff_at_node;
                    is_mid_branch = false;
                    updating_node->failure_count = 0;
                    best_up_lh_diff = lh_diff_mid_branch;
                }
                else if (lh_diff_mid_branch >= (best_lh_diff - threshold_prob))
                {
                    best_up_lh_diff = updating_node->likelihood_diff;
                    best_down_lh_diff = lh_diff_at_node;
                }
                // placement at current node is considered failed if placement likelihood is not improved by a certain margin compared to best placement so far for the nodes above it.
                else if (lh_diff_at_node < (updating_node->likelihood_diff - params->thresh_log_lh_failure))
                    ++updating_node->failure_count;
                    
                // delete at_node_regions
                if (delete_at_node_regions) delete at_node_regions;
            }
            // set the placement cost at the current node position at the most negative value if branch length is zero -> we can't place the subtree on that branch
            else
                lh_diff_at_node = updating_node->likelihood_diff;
            
            // keep crawling down into children nodes unless the stop criteria for the traversal are satisfied.
            // check the stop criteria
            bool keep_traversing = false;
            /*if (search_subtree_placement)
            {*/
            if (strict_stop_seeking_placement_subtree)
            {
                if (updating_node->failure_count <= failure_limit_subtree && lh_diff_at_node > (best_lh_diff - thresh_log_lh_subtree) && updating_node->node->next)
                    keep_traversing = true;
            }
            else
            {
                if ((updating_node->failure_count <= failure_limit_subtree || lh_diff_at_node > (best_lh_diff - thresh_log_lh_subtree))
                    && updating_node->node->next)
                    keep_traversing = true;
            }
            /*}
            else
            {
                if (params->strict_stop_seeking_placement_sample)
                {
                    if (updating_node->failure_count <= params->failure_limit_sample && lh_diff_at_node > (best_lh_diff - params->thresh_log_lh_sample)
                        && updating_node->node->next)
                            keep_traversing = true;
                }
                else
                {
                    if ((updating_node->failure_count <= params->failure_limit_sample || lh_diff_at_node > (best_lh_diff - params->thresh_log_lh_sample))
                        && updating_node->node->next)
                            keep_traversing = true;
                }
            }*/
            
            // keep traversing further down to the children
            if (keep_traversing)
            {
                Node* child_1 = updating_node->node->getOtherNextNode()->neighbor;
                Node* child_2 = child_1->neighbor->getOtherNextNode()->neighbor;
                
                // add child_1 to node_stack
                SeqRegions* upper_lr_regions = NULL;
                SeqRegions* lower_regions = child_2->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
                // get or recompute the upper left/right regions of the children node
                if (updating_node->need_updating)
                    updating_node->incoming_regions->mergeUpperLower(upper_lr_regions, updating_node->branch_length, *lower_regions, child_2->length, aln, model, threshold_prob);
                else
                {
                    if (child_1->neighbor->partial_lh)
                        upper_lr_regions = child_1->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
                }
                // traverse to this child's subtree
                if (upper_lr_regions)
                    node_stack.push(new UpdatingNode(child_1, upper_lr_regions, child_1->length, updating_node->need_updating, lh_diff_at_node, updating_node->failure_count, updating_node->need_updating));
                
                // add child_2 to node_stack
                upper_lr_regions = NULL;
                lower_regions = child_1->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
                // get or recompute the upper left/right regions of the children node
                if (updating_node->need_updating)
                    updating_node->incoming_regions->mergeUpperLower(upper_lr_regions, updating_node->branch_length, *lower_regions, child_1->length, aln, model, threshold_prob);
                else
                {
                    if (child_2->neighbor->partial_lh)
                        upper_lr_regions = child_2->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
                }
                // traverse to this child's subtree
                if (upper_lr_regions)
                    node_stack.push(new UpdatingNode(child_2, upper_lr_regions, child_2->length, updating_node->need_updating, lh_diff_at_node, updating_node->failure_count, updating_node->need_updating));
            }
        }
        // case when crawling up from child to parent
        else
        {
            Node* top_node = updating_node->node->getTopNode();
            
            // append directly at the node
            if (top_node->length > 0 || top_node == root)
            {
                SeqRegions* at_node_regions = NULL;
                bool delete_at_node_regions = false;
                // get or recompute the regions when placing the subtree at the current node position
                if (updating_node->need_updating)
                {
                    delete_at_node_regions = true;
                    
                    SeqRegions* upper_lr_regions = updating_node->node->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
                    upper_lr_regions->mergeUpperLower(at_node_regions, -1, *updating_node->incoming_regions, updating_node->branch_length, aln, model, threshold_prob);
                    
                    // skip if at_node_regions is null
                    if (!at_node_regions)
                    {
                        // print("Removing a node created an inconsistency while moving up.")
                        // delete updating_node
                        delete updating_node;
                        
                        continue;
                    }
                    // stop updating if the difference between the new and old regions is insignificant
                    else if (!at_node_regions->areDiffFrom(*top_node->total_lh, seq_length, num_states, &params.value()))
                        updating_node->need_updating = false;
                }
                else
                    at_node_regions = top_node->total_lh;
                 
                // compute the placement cost
                //if (search_subtree_placement)
                lh_diff_at_node = calculateSubTreePlacementCost(cumulative_rate, at_node_regions, subtree_regions, removed_blength);
                /*else
                    lh_diff_at_node = calculateSamplePlacementCost(cumulative_rate, at_node_regions, subtree_regions, removed_blength);*/
                    
                // delete at_node_regions
                if (delete_at_node_regions) delete at_node_regions;
                
                // placement at current node is considered failed if placement likelihood is not improved by a certain margin compared to best placement so far for the nodes above it.
                if (lh_diff_at_node < (updating_node->likelihood_diff - params->thresh_log_lh_failure))
                    ++updating_node->failure_count;
                else if (lh_diff_at_node > best_lh_diff)
                {
                    best_node = top_node;
                    best_lh_diff = lh_diff_at_node;
                    is_mid_branch = false;
                    updating_node->failure_count = 0;
                }
            }
            // if placement cost at new position gets worse -> restore to the old one
            else
                lh_diff_at_node = updating_node->likelihood_diff;

            // try appending mid-branch
            Node* other_child = updating_node->node->getOtherNextNode()->neighbor;
            SeqRegions* bottom_regions = NULL;
            if (top_node->length > 0 && top_node != root)
            {
                SeqRegions* mid_branch_regions = NULL;
                // get or recompute the regions when placing the subtree at the mid-branch position
                if (updating_node->need_updating)
                {
                    SeqRegions* other_child_lower_regions = other_child->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
                    other_child_lower_regions->mergeTwoLowers(bottom_regions, other_child->length, updating_node->incoming_regions, updating_node->branch_length, aln, model, threshold_prob, cumulative_rate);
                    
                    // skip if bottom_regions is null (inconsistent)
                    if (!bottom_regions)
                    {
                        // delete updating_node
                        delete updating_node;

                        continue;
                    }
                   
                    // compute new mid-branch regions
                    SeqRegions* upper_lr_regions = top_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
                    RealNumType mid_branch_length = top_node->length * 0.5;
                    upper_lr_regions->mergeUpperLower(mid_branch_regions, mid_branch_length, *bottom_regions, mid_branch_length, aln, model, threshold_prob);
                }
                else
                    mid_branch_regions = top_node->mid_branch_lh;
                
                // skip if mid_branch_regions is null (inconsistent)
                if (!mid_branch_regions)
                {
                    // delete bottom_regions if it's existed
                    if (bottom_regions) delete bottom_regions;
                    
                    // delete updating_node
                    delete updating_node;
                    
                    continue;
                }
                
                // compute the placement cost
                //if (search_subtree_placement)
                lh_diff_mid_branch = calculateSubTreePlacementCost(cumulative_rate, mid_branch_regions, subtree_regions, removed_blength);
                /*else
                    lh_diff_mid_branch = calculateSamplePlacementCost(cumulative_rate, mid_branch_regions, subtree_regions, removed_blength);*/
                    
                // delete mid_branch_regions
                if (updating_node->need_updating) delete mid_branch_regions;
                
                if (best_node == top_node)
                    best_up_lh_diff = lh_diff_mid_branch;
                
                // if this position is better than the best position found so far -> record it
                if (lh_diff_mid_branch > best_lh_diff)
                {
                    best_node = top_node;
                    best_lh_diff = lh_diff_mid_branch;
                    is_mid_branch = true;
                    updating_node->failure_count = 0;
                    best_down_lh_diff = lh_diff_at_node;
                }
                else if (lh_diff_at_node >= (best_lh_diff - threshold_prob))
                    best_up_lh_diff = lh_diff_mid_branch;
            }
            // set the placement cost at the mid-branch position at the most negative value if branch length is zero -> we can't place the subtree on that branch
            // NHANLT: we actually don't need to do that since lh_diff_mid_branch will never be read
            // else
                // lh_diff_mid_branch = MIN_NEGATIVE;
            
            // check stop rule of the traversal process
            bool keep_traversing = false;
            /*if (search_subtree_placement)
            {*/
            if (strict_stop_seeking_placement_subtree)
            {
                if (updating_node->failure_count <= failure_limit_subtree && lh_diff_at_node > (best_lh_diff - thresh_log_lh_subtree))
                    keep_traversing = true;
            }
            else if (updating_node->failure_count <= failure_limit_subtree || lh_diff_at_node > (best_lh_diff - thresh_log_lh_subtree))
                keep_traversing = true;
            /*    }
            else
            {
                if (params->strict_stop_seeking_placement_sample)
                {
                    if (updating_node->failure_count <= params->failure_limit_sample && lh_diff_at_node > (best_lh_diff - params->thresh_log_lh_sample))
                        keep_traversing = true;
                }
                else if (updating_node->failure_count <= params->failure_limit_sample || lh_diff_at_node > (best_lh_diff - params->thresh_log_lh_sample))
                    keep_traversing = true;
            }*/
            
            // keep traversing upwards
            if (keep_traversing)
            {
                // keep crawling up into parent and sibling node
                // case the node is not the root
                if (top_node != root)
                {
                    // first pass the crawling down the other child (sibling)
                    parent_upper_lr_regions = top_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
                    
                    SeqRegions* upper_lr_regions = NULL;
                    // get or recompute the upper left/right regions of the sibling node
                    if (updating_node->need_updating)
                        parent_upper_lr_regions->mergeUpperLower(upper_lr_regions, top_node->length, *updating_node->incoming_regions, updating_node->branch_length, aln, model, threshold_prob);
                    else
                    {
                        if (updating_node->node->partial_lh)
                            upper_lr_regions = updating_node->node->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
                    }

                    // add sibling node to node_stack for traversing later; skip if upper_lr_regions is null (inconsistent)
                    if (!upper_lr_regions)
                    {
                        // delete bottom_regions if it's existed
                        if (bottom_regions) delete bottom_regions;
                        
                        // delete updating_node
                        delete updating_node;
                        
                        continue;
                    }
                    else
                        node_stack.push(new UpdatingNode(other_child, upper_lr_regions, other_child->length, updating_node->need_updating, lh_diff_at_node, updating_node->failure_count, updating_node->need_updating));
                    
                    // now pass the crawling up to the parent node
                    // get or recompute the bottom regions (comming from 2 children) of the parent node
                    if (updating_node->need_updating)
                    {
                        if (!bottom_regions)
                        {
                            SeqRegions* other_child_lower_regions = other_child->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
                            other_child_lower_regions->mergeTwoLowers(bottom_regions, other_child->length, updating_node->incoming_regions, updating_node->branch_length, aln, model, threshold_prob, cumulative_rate);
                            
                            // skip if bottom_regions is null (inconsistent)
                            if (!bottom_regions)
                            {
                                // delete updating_node
                                delete updating_node;
                                
                                continue;
                            }
                        }
                    }
                    else
                    {
                        if (bottom_regions) delete bottom_regions;
                        bottom_regions = top_node->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
                    }
                    
                    // add the parent node to node_stack for traversing later
                    node_stack.push(new UpdatingNode(top_node->neighbor, bottom_regions, top_node->length, updating_node->need_updating, lh_diff_at_node, updating_node->failure_count, updating_node->need_updating));
                }
                // now consider case of root node -> only need to care about the sibling node
                else
                {
                    SeqRegions* upper_lr_regions = NULL;
                    
                    // get or recompute the upper left/right regions of the sibling node
                    if (updating_node->need_updating)
                        upper_lr_regions = updating_node->incoming_regions->computeTotalLhAtRoot(num_states, model, updating_node->branch_length);
                    else
                        upper_lr_regions = updating_node->node->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
                    
                    // add the sibling node to node_stack for traversing later
                    node_stack.push(new UpdatingNode(other_child, upper_lr_regions, other_child->length, updating_node->need_updating, lh_diff_at_node, updating_node->failure_count, updating_node->need_updating));
                    
                    // delete bottom_regions
                    if (bottom_regions) delete bottom_regions;
                }
            }
            else
            {
                // delete bottom_regions if it's existed
                if (bottom_regions) delete bottom_regions;
            }
        }
        
        // delete updating_node
        delete updating_node;
    }
    
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
            
            best_up_lh_diff = calculateSamplePlacementCost(cumulative_rate, parent_node->total_lh, subtree_regions, removed_blength);
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
                    SeqRegions* upper_lr_regions = node->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
                    SeqRegions* lower_regions = node->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
                    SeqRegions* mid_branch_regions = new SeqRegions(node->mid_branch_lh, aln.num_states);

                    // try to place new sample along the upper half of the current branch
                    while (true)
                    {
                        // compute the placement cost
                        RealNumType new_lh_mid_branch = calculateSamplePlacementCost(cumulative_rate, mid_branch_regions, subtree_regions, removed_blength);
                        
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
                    
                    //RealNumType new_best_lh_mid_branch = calculateSamplePlacementCost(cumulative_rate, node->mid_branch_lh, sample_regions, default_blength);
                    
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
}

void Tree::applySPR(Node* subtree, Node* best_node, bool is_mid_branch, RealNumType branch_length, RealNumType best_lh_diff, RealNumType* cumulative_rate, vector< vector<PositionType> > &cumulative_base, RealNumType default_blength, RealNumType max_blength, RealNumType min_blength, RealNumType min_blength_mid)
{
    // remove subtree from the tree
    Node* parent_subtree = subtree->neighbor->getTopNode();
    Node* sibling_subtree = subtree->neighbor->getOtherNextNode()->neighbor;
    RealNumType threshold_prob = params->threshold_prob;
    StateType num_states = aln.num_states;
    
    // connect grandparent to sibling
    if (parent_subtree != root)
        parent_subtree->neighbor->neighbor = sibling_subtree;
    sibling_subtree->neighbor = parent_subtree->neighbor;
    
    // update the length of the branch connecting grandparent to sibling
    if (sibling_subtree->length > 0)
    {
        if (parent_subtree->length > 0)
            sibling_subtree->length += parent_subtree->length;
    }
    else
        sibling_subtree->length = parent_subtree->length;
    
    // update likelihood lists after subtree removal
    // case when the sibling_subtree becomes the new root
    if (!sibling_subtree->neighbor)
    {
        // update root
        root = sibling_subtree;
        
        // delete mid_branch_lh
        if (root->mid_branch_lh)
        {
            delete sibling_subtree->mid_branch_lh;
            sibling_subtree->mid_branch_lh = NULL;
        }
        
        // reset branch length (to 0) if sibling_subtree is root
        sibling_subtree->length = 0;

        // recompute the total lh regions at sibling
        sibling_subtree->computeTotalLhAtNode(aln, model, threshold_prob, cumulative_rate, true);
        
        // traverse downwards (to childrens of the sibling) to update their lh regions
        if (sibling_subtree->next)
        {
            // update upper left/right regions
            Node* next_node_1 = sibling_subtree->next;
            Node* next_node_2 = next_node_1->next;
            
            if (next_node_1->partial_lh) delete next_node_1->partial_lh;
            SeqRegions* lower_reions = next_node_2->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
            next_node_1->partial_lh = lower_reions->computeTotalLhAtRoot(num_states, model, next_node_2->length);
            
            if (next_node_2->partial_lh) delete next_node_2->partial_lh;
            lower_reions = next_node_1->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
            next_node_2->partial_lh = lower_reions->computeTotalLhAtRoot(num_states, model, next_node_1->length);
            
            // add children to node_stack for further traversing and updating likelihood regions
            stack<Node*> node_stack;
            node_stack.push(next_node_1->neighbor);
            node_stack.push(next_node_2->neighbor);
            updatePartialLh(node_stack, cumulative_rate, default_blength, min_blength, max_blength);
        }
    }
    // case when the sibling_subtree is non-root node
    else
    {
        // update branch length from the grandparent side
        sibling_subtree->neighbor->length = sibling_subtree->length;
        
        // traverse to parent and sibling node to update their likelihood regions due to subtree remova
        stack<Node*> node_stack;
        node_stack.push(sibling_subtree);
        node_stack.push(sibling_subtree->neighbor);
        updatePartialLh(node_stack, cumulative_rate, default_blength, min_blength, max_blength);
    }
    
    // replace the node and re-update the vector lists
    SeqRegions* subtree_lower_regions = subtree->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
    placeSubtree(best_node, subtree, subtree_lower_regions, is_mid_branch, branch_length, best_lh_diff, cumulative_rate, cumulative_base, default_blength, max_blength, min_blength, min_blength_mid);
}

void Tree::placeSubtree(Node* selected_node, Node* subtree, SeqRegions* subtree_regions, bool is_mid_branch, RealNumType new_branch_length, RealNumType new_lh, RealNumType* cumulative_rate, vector< vector<PositionType> > &cumulative_base, RealNumType default_blength, RealNumType max_blength, RealNumType min_blength, RealNumType min_blength_mid)
{
    // dummy variables
    RealNumType best_child_lh;
    RealNumType best_child_split = 0;
    RealNumType best_parent_lh;
    RealNumType best_parent_split = 0;
    SeqRegions* best_parent_regions = NULL;
    SeqRegions* best_child_regions = NULL;
    RealNumType best_root_blength = -1;
    StateType num_states = aln.num_states;
    RealNumType threshold_prob = params->threshold_prob;
    
    // in case of a polytomy, reach first the top of the polytomy, which is the only node at which appending is allowed.
    // NHANLT: this block seems to be unnecessary: (1) it doesn't affect the result in my tests with up to 10K sequences; (2) it is removed from new versions of MAPLE;
    /*if (selected_node)
        while (selected_node->length <= 0 && selected_node != root)
            selected_node = selected_node->neighbor->getTopNode();*/
    ASSERT(selected_node);
    
    // try to place the new sample as a descendant of a mid-branch point
    if (is_mid_branch && selected_node != root)
    {
        SeqRegions* upper_left_right_regions = selected_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
        RealNumType best_split = 0.5;
        RealNumType best_split_lh = new_lh;
        RealNumType new_split = 0.25;
        if (best_child_regions) delete best_child_regions;
        best_child_regions = new SeqRegions(selected_node->mid_branch_lh);
        SeqRegions* new_parent_regions = NULL;
        SeqRegions* lower_regions = selected_node->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
        RealNumType new_branch_length_split = selected_node->length * new_split;

        // try different positions on the existing branch
        while (new_branch_length_split > min_blength)
        {
            upper_left_right_regions->mergeUpperLower(new_parent_regions, new_branch_length_split, *lower_regions,  selected_node->length - new_branch_length_split, aln, model, threshold_prob);
            
            RealNumType placement_cost = calculateSubTreePlacementCost(cumulative_rate, new_parent_regions, subtree_regions, new_branch_length);
            
            if (placement_cost > best_split_lh)
            {
                best_split_lh = placement_cost;
                best_split = new_split;
                
                if (best_child_regions) delete best_child_regions;
                best_child_regions = new_parent_regions;
                new_parent_regions = NULL;
            }
            else
                break;
            
            new_split = best_split * 0.5;
            new_branch_length_split = selected_node->length * new_split;
        }
        
        if (best_split > 0.49)
        {
            new_split = 0.25;
            RealNumType new_branch_length_split = selected_node->length * new_split;
            while (new_branch_length_split > min_blength)
            {
                upper_left_right_regions->mergeUpperLower(new_parent_regions, selected_node->length - new_branch_length_split, *lower_regions, new_branch_length_split, aln, model, threshold_prob);
                
                RealNumType placement_cost = calculateSubTreePlacementCost(cumulative_rate, new_parent_regions, subtree_regions, new_branch_length);
                
                if (placement_cost > best_split_lh)
                {
                    best_split_lh = placement_cost;
                    best_split = new_split;
                    
                    if (best_child_regions) delete best_child_regions;
                    best_child_regions = new_parent_regions;
                    new_parent_regions = NULL;
                }
                else
                    break;
                
                new_split = best_split * 0.5;
                new_branch_length_split = selected_node->length * new_split;
            }
            if (best_split < 0.49)
                best_split = 1.0 - best_split;
        }
        
        // delete new_parent_regions
        if (new_parent_regions) delete new_parent_regions;
        
        // now try different lengths for the new branch
        RealNumType new_branch_lh = best_split_lh;
        RealNumType best_blength = new_branch_length;
        
        // change zero branch length to min branch length
        if (best_blength <= 0)
        {
            best_blength = min_blength;
            new_branch_lh = calculateSubTreePlacementCost(cumulative_rate, best_child_regions, subtree_regions, best_blength);
        }
        
        while (best_blength > min_blength)
        {
            RealNumType new_blength = best_blength * 0.5;
            RealNumType placement_cost = calculateSubTreePlacementCost(cumulative_rate, best_child_regions, subtree_regions, new_blength);
            
            if (placement_cost > new_branch_lh)
            {
                new_branch_lh = placement_cost;
                best_blength = new_blength;
            }
            else
                break;
        }
        if (new_branch_length <= 0 || best_blength > 0.7 * new_branch_length)
        {
            while (best_blength < max_blength)
            {
                RealNumType new_blength = best_blength + best_blength;
                RealNumType placement_cost = calculateSubTreePlacementCost(cumulative_rate, best_child_regions, subtree_regions, new_blength);
                if (placement_cost > new_branch_lh)
                {
                    new_branch_lh = placement_cost;
                    best_blength = new_blength;
                }
                else
                    break;
            }
        }
        
        // try zero branch length if the branch length is very short
        if (best_blength < min_blength + min_blength)
        {
            RealNumType zero_branch_lh = calculateSubTreePlacementCost(cumulative_rate, best_child_regions, subtree_regions, -1);
            if (zero_branch_lh > new_branch_lh)
                best_blength = -1;
        }
        
        // attach subtree to the branch above the selected node
        RealNumType top_distance = selected_node->length * best_split;
        RealNumType down_distance = selected_node->length - top_distance;
        
        // re-use internal nodes
        Node* next_node_1 = subtree->neighbor;
        Node* new_internal_node = next_node_1->getTopNode();
        Node* next_node_2 = next_node_1->getOtherNextNode();
        
        // NHANLT NOTES: UNNECESSARY
        // re-order next circle (not neccessary, just to make it consistent with Python code)
        new_internal_node->next = next_node_2;
        next_node_2->next = next_node_1;
        next_node_1->next = new_internal_node;
        
        // connect new_internal_node to the parent of the selected node
        new_internal_node->outdated = true;
        selected_node->neighbor->neighbor = new_internal_node;
        new_internal_node->neighbor = selected_node->neighbor;
        new_internal_node->length = top_distance;
        new_internal_node->neighbor->length = top_distance;
        
        // connect the selected_node to new_internal_node (via next_node_2)
        selected_node->neighbor = next_node_2;
        next_node_2->neighbor = selected_node;
        selected_node->length = down_distance;
        selected_node->neighbor->length = down_distance;
        
        // subtree already connected to new_internal_node (via next_node_1)
        subtree->length = best_blength;
        subtree->neighbor->length = best_blength;
                
        // update all likelihood regions
        if (next_node_1->partial_lh) delete next_node_1->partial_lh;
        next_node_1->partial_lh = best_child_regions;
        best_child_regions = NULL;
        upper_left_right_regions->mergeUpperLower(next_node_2->partial_lh, new_internal_node->length, *subtree_regions, best_blength, aln, model, threshold_prob);
        selected_node->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate)->mergeTwoLowers(new_internal_node->partial_lh, selected_node->length, subtree_regions, best_blength, aln, model, threshold_prob, cumulative_rate);
        RealNumType mid_branch_length = new_internal_node->length * 0.5;
        upper_left_right_regions->mergeUpperLower(new_internal_node->mid_branch_lh, mid_branch_length, *new_internal_node->partial_lh, mid_branch_length, aln, model, threshold_prob);
        new_internal_node->computeTotalLhAtNode(aln, model, threshold_prob, cumulative_rate, new_internal_node == root);
        
        if (!new_internal_node->total_lh || new_internal_node->total_lh->size() == 0)
            outError("Problem, None vector when re-placing sample, placing subtree at mid-branch point");
        
        /*if distTop>=2*min_blengthForMidNode:
         createFurtherMidNodes(newInternalNode,upper_left_right_regions)*/

        // iteratively traverse the tree to update partials from the current node
        stack<Node*> node_stack;
        node_stack.push(selected_node);
        node_stack.push(subtree);
        node_stack.push(new_internal_node->neighbor);
        updatePartialLh(node_stack, cumulative_rate, default_blength, min_blength, max_blength);
    }
    // otherwise, best lk so far is for appending directly to existing node
    else
    {
        // We first explore placement just below the best placement node for more fine-grained placement within its descendant branches (accounting for polytomies).
        RealNumType best_down_lh_diff = MIN_NEGATIVE;
        Node* best_child = NULL;
        
        // current node might be part of a polytomy (represented by 0 branch lengths) so we want to explore all the children of the current node to find out if the best placement is actually in any of the branches below the current node.
        stack<Node*> new_node_stack;
        RealNumType best_split = 0.5;
        RealNumType half_min_blength_mid = min_blength_mid * 0.5;
        Node* neighbor_node;
        FOR_NEIGHBOR(selected_node, neighbor_node)
            new_node_stack.push(neighbor_node);
        
        while (!new_node_stack.empty())
        {
            Node* node = new_node_stack.top();
            new_node_stack.pop();

            // add all nodes in polytomy
            if (node->length <= 0)
            {
                FOR_NEIGHBOR(node, neighbor_node)
                    new_node_stack.push(neighbor_node);
            }
            else
            {
                RealNumType new_split = 0.5;
                RealNumType new_best_split = 0.5;
                
                // now try to place on the current branch below the best node, at an height above or equal to the mid-branch.
                RealNumType new_branch_length_split;
                RealNumType tmp_best_lh_diff = MIN_NEGATIVE;
                SeqRegions* mid_branch_regions = new SeqRegions(node->mid_branch_lh);
                SeqRegions* tmp_best_child_regions = NULL;
                SeqRegions* parent_upper_lr_regions = node->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
                SeqRegions* lower_regions = node->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
                
                while (true)
                {
                    RealNumType tmp_lh_diff = calculateSubTreePlacementCost(cumulative_rate, mid_branch_regions, subtree_regions, new_branch_length);
                    
                    // if better placement found -> record it
                    if (tmp_lh_diff > tmp_best_lh_diff)
                    {
                        tmp_best_lh_diff = tmp_lh_diff;
                        new_best_split = new_split;
                        if (tmp_best_child_regions) delete tmp_best_child_regions;
                        tmp_best_child_regions = mid_branch_regions;
                        mid_branch_regions = NULL;
                    }
                    else
                        break;
                    
                    // update new_branch_length_split
                    new_split *= 0.5;
                    new_branch_length_split = node->length * new_split;
                    if (new_branch_length_split <= half_min_blength_mid)
                        break;
                    
                    // compute mid_branch_regions
                    parent_upper_lr_regions->mergeUpperLower(mid_branch_regions, new_branch_length_split, *lower_regions, node->length - new_branch_length_split, aln, model, threshold_prob);
                }
                
                // record new placement position if we found a better one
                if (tmp_best_lh_diff > best_down_lh_diff)
                {
                    best_down_lh_diff = tmp_best_lh_diff;
                    best_child = node;
                    best_split = new_best_split;
                    if (best_child_regions) delete best_child_regions;
                    best_child_regions = tmp_best_child_regions;
                    tmp_best_child_regions = NULL;
                }
                
                // delete mid_branch_regions and tmp_best_child_regions
                if (mid_branch_regions) delete mid_branch_regions;
                if (tmp_best_child_regions) delete tmp_best_child_regions;
            }
        }
        
        // place the new sample as a descendant of an existing node
        if (best_child)
        {
            RealNumType best_split_lh = best_down_lh_diff;
            SeqRegions* upper_left_right_regions = best_child->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
            SeqRegions* lower_regions = best_child->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
            RealNumType new_split = best_split * 0.5;
            RealNumType new_branch_length_split = new_split * best_child->length;
            SeqRegions* new_parent_regions = NULL;
            
            while (new_branch_length_split > min_blength)
            {
                upper_left_right_regions->mergeUpperLower(new_parent_regions, new_branch_length_split, *lower_regions, best_child->length - new_branch_length_split, aln, model, threshold_prob);
                
                RealNumType placement_cost = calculateSubTreePlacementCost(cumulative_rate, new_parent_regions, subtree_regions, new_branch_length);
                
                if (placement_cost > best_split_lh)
                {
                    best_split_lh = placement_cost;
                    best_split = new_split;
                    
                    if (best_child_regions) delete best_child_regions;
                    best_child_regions = new_parent_regions;
                    new_parent_regions = NULL;
                }
                else
                    break;
                
                new_split = best_split * 0.5;
                new_branch_length_split = new_split * best_child->length;
            }
            
            // delete new_parent_regions
            if (new_parent_regions) delete new_parent_regions;
            
            best_child_lh = best_split_lh;
            best_child_split = best_split;
        }
        else
            best_child_lh = MIN_NEGATIVE;
        
        // if node is root, try to place as sibling of the current root.
        RealNumType old_root_lh = MIN_NEGATIVE;
        if (root == selected_node)
        {
            SeqRegions* lower_regions = selected_node->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
            old_root_lh = lower_regions->computeAbsoluteLhAtRoot(aln, model, cumulative_base);
            
            // merge 2 lower vector into one
            SeqRegions* merged_root_sample_regions = NULL;
            RealNumType new_root_lh = lower_regions->mergeTwoLowers(merged_root_sample_regions, default_blength, subtree_regions, new_branch_length, aln, model, threshold_prob, cumulative_rate, true);
            
            new_root_lh += merged_root_sample_regions->computeAbsoluteLhAtRoot(aln, model, cumulative_base);
            best_parent_lh = new_root_lh - old_root_lh;
            best_root_blength = default_blength;
            
            if (best_parent_regions) delete best_parent_regions;
            best_parent_regions = merged_root_sample_regions;
            merged_root_sample_regions = NULL;
            RealNumType new_blength = 0.5 * default_blength;
            
            while (new_blength > min_blength)
            {
                // merge 2 lower vector into one
                new_root_lh = lower_regions->mergeTwoLowers(merged_root_sample_regions, new_blength, subtree_regions, new_branch_length, aln, model, threshold_prob, cumulative_rate, true);
                
                new_root_lh += merged_root_sample_regions->computeAbsoluteLhAtRoot(aln, model, cumulative_base);
                RealNumType diff_root_lh = new_root_lh - old_root_lh;
                if (diff_root_lh > best_parent_lh)
                {
                    best_parent_lh = diff_root_lh;
                    best_root_blength = new_blength;
                    
                    if (best_parent_regions) delete best_parent_regions;
                    best_parent_regions = merged_root_sample_regions;
                    merged_root_sample_regions = NULL;
                }
                else
                    break;
                
                new_blength = best_root_blength * 0.5;
            }
            
            // delete merged_root_sample_regions
            if (merged_root_sample_regions) delete merged_root_sample_regions;
        }
        // selected_node is not root
        // try to append just above node
        else
        {
            best_split = 0.5;
            SeqRegions* upper_left_right_regions = selected_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
            if (best_parent_regions) delete best_parent_regions;
            best_parent_regions = new SeqRegions(selected_node->mid_branch_lh);
            RealNumType best_split_lh = calculateSubTreePlacementCost(cumulative_rate, best_parent_regions, subtree_regions, new_branch_length);
            SeqRegions* lower_regions = selected_node->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
            SeqRegions* new_parent_regions = NULL;
            RealNumType new_split = 0.25;
            RealNumType new_branch_length_split = new_split * selected_node->length;
            
            while (new_branch_length_split > min_blength)
            {
                upper_left_right_regions->mergeUpperLower(new_parent_regions, selected_node->length - new_branch_length_split, *lower_regions, new_branch_length_split, aln, model, threshold_prob);
                
                RealNumType placement_cost = calculateSubTreePlacementCost(cumulative_rate, new_parent_regions, subtree_regions, new_branch_length);
                
                if (placement_cost > best_split_lh)
                {
                    best_split_lh = placement_cost;
                    best_split = new_split;
                    if (best_parent_regions) delete best_parent_regions;
                    best_parent_regions = new_parent_regions;
                    new_parent_regions = NULL;
                }
                else
                    break;
                
                new_split = best_split * 0.5;
                new_branch_length_split = new_split * selected_node->length;
            }
            // delete new_parent_regions
            if (new_parent_regions) delete new_parent_regions;
            
            best_parent_lh = best_split_lh;
            best_parent_split = best_split;
        }
        
        // if the best placement is below the selected_node => add an internal node below the selected_node
        /* now we have three likelihood costs,
        best_child_lh is the likelihood score of appending below node;
        best_parent_lh is the likelihood score of appending above node;
        new_lh is the likelihood cost of appending exactly at node. */
        if (best_child_lh >= best_parent_lh && best_child_lh >= new_lh)
        {
            ASSERT(best_child);
            
            SeqRegions* upper_left_right_regions = best_child->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);

            RealNumType new_branch_length_lh = best_child_lh;
            RealNumType best_length = new_branch_length;
            
            // change zero branch length to min branch length
            if (best_length <= 0)
            {
                best_length = min_blength;
                new_branch_length_lh = calculateSubTreePlacementCost(cumulative_rate, best_child_regions, subtree_regions, best_length);
            }
            
            while (best_length > min_blength)
            {
                RealNumType new_blength = best_length * 0.5;
                RealNumType placement_cost = calculateSubTreePlacementCost(cumulative_rate, best_child_regions, subtree_regions, new_blength);
                
                if (placement_cost > new_branch_length_lh)
                {
                    new_branch_length_lh = placement_cost;
                    best_length = new_blength;
                }
                else
                    break;
            }
            
            if (new_branch_length <= 0 || best_length > 0.7 * new_branch_length)
            {
                while (best_length < max_blength)
                {
                    RealNumType new_blength = best_length + best_length;
                    RealNumType placement_cost = calculateSubTreePlacementCost(cumulative_rate, best_child_regions, subtree_regions, new_blength);
                    if (placement_cost > new_branch_length_lh)
                    {
                        new_branch_length_lh = placement_cost;
                        best_length = new_blength;
                    }
                    else
                        break;
                }
            }
            
            // try zero-branch length at very short branch
            if (best_length < min_blength + min_blength)
            {
                RealNumType zero_branch_lh = calculateSubTreePlacementCost(cumulative_rate, best_child_regions, subtree_regions, -1);
                
                if (zero_branch_lh > new_branch_length_lh)
                    best_length = -1;
            }
            
            // attach subtree to the phylogenetic tree (below the selected_node ~ above the child node)
            RealNumType top_distance = best_child->length * best_child_split;
            RealNumType down_distance = best_child->length - top_distance;
            
            // re-use internal nodes
            Node* next_node_1 = subtree->neighbor;
            Node* new_internal_node = next_node_1->getTopNode();
            Node* next_node_2 = next_node_1->getOtherNextNode();
            
            // NHANLT NOTES: UNNECESSARY
            // re-order next circle (not neccessary, just to make it consistent with Python code)
            new_internal_node->next = next_node_2;
            next_node_2->next = next_node_1;
            next_node_1->next = new_internal_node;
            
            // connect new_internal_node to the parent of the selected node
            new_internal_node->outdated = true;
            best_child->neighbor->neighbor = new_internal_node;
            new_internal_node->neighbor = best_child->neighbor;
            new_internal_node->length = top_distance;
            new_internal_node->neighbor->length = top_distance;
            
            // connect the selected_node to new_internal_node (via next_node_2)
            best_child->neighbor = next_node_2;
            next_node_2->neighbor = best_child;
            best_child->length = down_distance;
            best_child->neighbor->length = down_distance;
            
            // subtree already connected to new_internal_node (via next_node_1)
            subtree->length = best_length;
            subtree->neighbor->length = best_length;
                    
            // update all likelihood regions
            if (next_node_1->partial_lh) delete next_node_1->partial_lh;
            next_node_1->partial_lh = best_child_regions;
            best_child_regions = NULL;
            upper_left_right_regions->mergeUpperLower(next_node_2->partial_lh, new_internal_node->length, *subtree_regions, best_length, aln, model, threshold_prob);
            best_child->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate)->mergeTwoLowers(new_internal_node->partial_lh, best_child->length, subtree_regions, best_length, aln, model, threshold_prob, cumulative_rate);
            RealNumType mid_branch_length = new_internal_node->length * 0.5;
            upper_left_right_regions->mergeUpperLower(new_internal_node->mid_branch_lh, mid_branch_length, *new_internal_node->partial_lh, mid_branch_length, aln, model, threshold_prob);
            new_internal_node->computeTotalLhAtNode(aln, model, threshold_prob, cumulative_rate, new_internal_node == root);
            
            if (!new_internal_node->total_lh)
                outWarning("Problem, None vector when re-placing sample, below the selected node");
                
            /*if (top_distance >= 2 * min_blengthForMidNode)
                createFurtherMidNodes(new_internal_node,this_node_upper_left_right_regions)*/

            // iteratively traverse the tree to update partials from the current node
            stack<Node*> node_stack;
            node_stack.push(best_child);
            node_stack.push(subtree);
            node_stack.push(new_internal_node->neighbor);
            updatePartialLh(node_stack, cumulative_rate, default_blength, min_blength, max_blength);
        }
        // otherwise, add new parent to the selected_node
        else
        {
            SeqRegions* lower_regions = selected_node->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
            
            // new parent is actually part of a polytomy since best placement is exactly at the node
            if (new_lh >= best_parent_lh)
            {
                best_root_blength = -1;
                best_parent_split = -1;
                best_parent_lh = new_lh;
                if (best_parent_regions) delete best_parent_regions;
                best_parent_regions = NULL;
                
                if (selected_node == root)
                    lower_regions->mergeTwoLowers(best_parent_regions, -1, subtree_regions, new_branch_length, aln, model, threshold_prob, cumulative_rate);
                else
                    best_parent_regions = new SeqRegions(selected_node->total_lh);
            }

            // add parent to the root
            if (selected_node == root)
            {
                // now try different lengths for right branch
                RealNumType best_length2 = new_branch_length;
                SeqRegions* new_root_lower_regions = NULL;
                RealNumType new_root_lh = 0;
                
                // try min branch length if branch length is zero
                if (best_length2 <= 0)
                {
                    best_length2 = min_blength;
                    new_root_lh = lower_regions->mergeTwoLowers(new_root_lower_regions, best_root_blength, subtree_regions, best_length2, aln, model, threshold_prob, cumulative_rate, true);
                    
                    new_root_lh += new_root_lower_regions->computeAbsoluteLhAtRoot(aln, model, cumulative_base);
                    
                    best_parent_lh = new_root_lh - old_root_lh;
                    if (best_parent_regions) delete best_parent_regions;
                    best_parent_regions = new_root_lower_regions;
                    new_root_lower_regions = NULL;
                }
                
                while (best_length2 > min_blength)
                {
                    RealNumType new_blength = best_length2 * 0.5;
                    
                    new_root_lh = lower_regions->mergeTwoLowers(new_root_lower_regions, best_root_blength, subtree_regions, new_blength, aln, model, threshold_prob, cumulative_rate, true);
                    
                    new_root_lh += new_root_lower_regions->computeAbsoluteLhAtRoot(aln, model, cumulative_base);
                    
                    RealNumType root_lh_diff = new_root_lh - old_root_lh;
                    if (root_lh_diff > best_parent_lh)
                    {
                        best_parent_lh = root_lh_diff;
                        best_length2 = new_blength;
                        
                        if (best_parent_regions) delete best_parent_regions;
                        best_parent_regions = new_root_lower_regions;
                        new_root_lower_regions = NULL;
                    }
                    else
                        break;
                }
                
                if (new_branch_length <= 0 || best_length2 > 0.7 * new_branch_length)
                {
                    while (best_length2 < max_blength)
                    {
                        RealNumType new_blength = best_length2 + best_length2;
                        new_root_lh = lower_regions->mergeTwoLowers(new_root_lower_regions, best_root_blength, subtree_regions, new_blength, aln, model, threshold_prob, cumulative_rate, true);
                        new_root_lh += new_root_lower_regions->computeAbsoluteLhAtRoot(aln, model, cumulative_base);
                        RealNumType root_lh_diff = new_root_lh - old_root_lh;
                        
                        if (root_lh_diff > best_parent_lh)
                        {
                            best_parent_lh = root_lh_diff;
                            best_length2 = new_blength;
                            
                            if (best_parent_regions) delete best_parent_regions;
                            best_parent_regions = new_root_lower_regions;
                            new_root_lower_regions = NULL;
                        }
                        else
                            break;
                    }
                }
                
                // try with length zero
                if (best_length2 < min_blength + min_blength)
                {
                    new_root_lh = lower_regions->mergeTwoLowers(new_root_lower_regions, best_root_blength, subtree_regions, -1, aln, model, threshold_prob, cumulative_rate, true);
                    new_root_lh += new_root_lower_regions->computeAbsoluteLhAtRoot(aln, model, cumulative_base);
                    RealNumType root_lh_diff = new_root_lh - old_root_lh;
                    if (root_lh_diff > best_parent_lh)
                    {
                        best_length2 = -1;
                        // best_parent_lh = root_lh_diff;
                        if (best_parent_regions) delete best_parent_regions;
                        best_parent_regions = new_root_lower_regions;
                        new_root_lower_regions = NULL;
                    }
                }

                // delete new_root_lower_regions
                if (new_root_lower_regions) delete new_root_lower_regions;
                
                // attach subtree to the phylogenetic tree (exactly at the seleted root node)
                // re-use internal nodes
                Node* next_node_1 = subtree->neighbor;
                Node* new_root = next_node_1->getTopNode();
                Node* next_node_2 = next_node_1->getOtherNextNode();
                
                // NHANLT NOTES: UNNECESSARY
                // re-order next circle (not neccessary, just to make it consistent with Python code)
                new_root->next = next_node_2;
                next_node_2->next = next_node_1;
                next_node_1->next = new_root;
                
                // connect new_internal_node to the parent of the selected node
                new_root->outdated = true;
                new_root->neighbor = selected_node->neighbor; // actually NULL since selected_node is root
                new_root->length = 0;
                
                // connect the selected_node to new_internal_node (via next_node_2)
                selected_node->neighbor = next_node_2;
                next_node_2->neighbor = selected_node;
                selected_node->length = best_root_blength;
                selected_node->neighbor->length = best_root_blength;
                if (best_root_blength <= 0)
                {
                    delete selected_node->total_lh;
                    selected_node->total_lh = NULL;
                    
                    if (selected_node->mid_branch_lh) delete selected_node->mid_branch_lh;
                    selected_node->mid_branch_lh = NULL;
                    //selected_node.furtherMidNodes=None
                }
                
                // subtree already connected to new_internal_node (via next_node_1)
                subtree->length = best_length2;
                subtree->neighbor->length = best_length2;
                        
                // update all likelihood regions
                if (new_root->partial_lh) delete new_root->partial_lh;
                new_root->partial_lh = best_parent_regions;
                best_parent_regions = NULL;
                if (new_root->mid_branch_lh) delete new_root->mid_branch_lh;
                new_root->mid_branch_lh = NULL;
                new_root->computeTotalLhAtNode(aln, model, threshold_prob, cumulative_rate, true);
                if (next_node_1->partial_lh) delete next_node_1->partial_lh;
                next_node_1->partial_lh = lower_regions->computeTotalLhAtRoot(num_states, model, best_root_blength);
                if (next_node_2->partial_lh) delete next_node_2->partial_lh;
                next_node_2->partial_lh = subtree_regions->computeTotalLhAtRoot(num_states, model, best_length2);
                
                if (!new_root->total_lh || new_root->total_lh->size() == 0)
                    outWarning("Problem, None vector when re-placing sample, position root");
                
                // update tree->root;
                root = new_root;
                
                // iteratively traverse the tree to update partials from the current node
                stack<Node*> node_stack;
                node_stack.push(selected_node);
                node_stack.push(subtree);
                updatePartialLh(node_stack, cumulative_rate, default_blength, min_blength, max_blength);
            }
            //add parent to non-root node (place subtree exactly at the selected non-root node)
            else
            {
                SeqRegions* upper_left_right_regions = selected_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
                
                // now try different lengths for the new branch
                RealNumType new_branch_length_lh = best_parent_lh;
                RealNumType best_length = new_branch_length;
                
                // try min branch length if the branch length is zero
                if (best_length <= 0)
                {
                    best_length = min_blength;
                    new_branch_length_lh = calculateSubTreePlacementCost(cumulative_rate, best_parent_regions, subtree_regions, best_length);
                }
                    
                while (best_length > min_blength)
                {
                    RealNumType new_blength = best_length * 0.5;
                    RealNumType placement_cost = calculateSubTreePlacementCost(cumulative_rate, best_parent_regions, subtree_regions, new_blength);
                    
                    if (placement_cost > new_branch_length_lh)
                    {
                        new_branch_length_lh = placement_cost;
                        best_length = new_blength;
                    }
                    else
                        break;
                }
                
                if (new_branch_length <= 0 || best_length > 0.7 * new_branch_length)
                {
                    RealNumType new_branch_length_10 = new_branch_length * 10;
                    while (best_length < new_branch_length_10)
                    {
                        RealNumType new_blength = best_length + best_length;
                        RealNumType placement_cost = calculateSubTreePlacementCost(cumulative_rate, best_parent_regions, subtree_regions, new_blength);
                        
                        if (placement_cost > new_branch_length_lh)
                        {
                            new_branch_length_lh = placement_cost;
                            best_length = new_blength;
                        }
                        else
                            break;
                    }
                }
                
                // try with length zero
                if (best_length < min_blength + min_blength)
                {
                    RealNumType zero_branch_length_lh = calculateSubTreePlacementCost(cumulative_rate, best_parent_regions, subtree_regions, -1);
                    
                    if (zero_branch_length_lh > new_branch_length_lh)
                        best_length = -1;
                }
                
                // attach subtree to the phylogenetic tree (exactly at the selected non-root node)
                RealNumType down_distance = selected_node->length * best_parent_split;
                RealNumType top_distance = selected_node->length - down_distance;
                if (best_parent_split <= 0)
                {
                    down_distance = -1;
                    top_distance = selected_node->length;
                    
                    if (selected_node->total_lh) delete selected_node->total_lh;
                    selected_node->total_lh = NULL;
                    
                    if (selected_node->mid_branch_lh) delete selected_node->mid_branch_lh;
                    selected_node->mid_branch_lh = NULL;
                }
                
                // re-use internal nodes
                Node* next_node_1 = subtree->neighbor;
                Node* new_internal_node = next_node_1->getTopNode();
                Node* next_node_2 = next_node_1->getOtherNextNode();
                
                // NHANLT NOTES: UNNECESSARY
                // re-order next circle (not neccessary, just to make it consistent with Python code)
                new_internal_node->next = next_node_2;
                next_node_2->next = next_node_1;
                next_node_1->next = new_internal_node;
                
                // connect new_internal_node to the parent of the selected node
                new_internal_node->outdated = true;
                selected_node->neighbor->neighbor = new_internal_node;
                new_internal_node->neighbor = selected_node->neighbor;
                new_internal_node->length = top_distance;
                new_internal_node->neighbor->length = top_distance;
                
                // connect the selected_node to new_internal_node (via next_node_2)
                selected_node->neighbor = next_node_2;
                next_node_2->neighbor = selected_node;
                selected_node->length = down_distance;
                selected_node->neighbor->length = down_distance;
                
                // subtree already connected to new_internal_node (via next_node_1)
                subtree->length = best_length;
                subtree->neighbor->length = best_length;
                        
                // update all likelihood regions
                SeqRegions* lower_regions = selected_node->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
                lower_regions->mergeTwoLowers(new_internal_node->partial_lh, selected_node->length, subtree_regions, best_length, aln, model, threshold_prob, cumulative_rate);
                if (!new_internal_node->partial_lh)
                {
                    outWarning("Problem, non lower likelihood while placing subtree -> set best branch length to min length");
                    best_length = min_blength;
                    subtree->length = best_length;
                    subtree->neighbor->length = best_length;
                    lower_regions->mergeTwoLowers(new_internal_node->partial_lh, selected_node->length, subtree_regions, best_length, aln, model, threshold_prob, cumulative_rate);
                }
                
                upper_left_right_regions->mergeUpperLower(next_node_1->partial_lh, new_internal_node->length, *lower_regions, selected_node->length, aln, model, threshold_prob);
                upper_left_right_regions->mergeUpperLower(next_node_2->partial_lh, new_internal_node->length, *subtree_regions, best_length, aln, model, threshold_prob);
                RealNumType mid_branch_length = new_internal_node->length * 0.5;
                upper_left_right_regions->mergeUpperLower(new_internal_node->mid_branch_lh, mid_branch_length, *new_internal_node->partial_lh, mid_branch_length, aln, model, threshold_prob);
                new_internal_node->computeTotalLhAtNode(aln, model, threshold_prob, cumulative_rate, new_internal_node == root);
                
                if (!new_internal_node->total_lh || new_internal_node->total_lh->size() == 0)
                    outWarning("Problem, None vector when placing sample, placing subtree at exactly non-root node");
                
                /*if (top_distance >= 2 * min_blengthForMidNode)
                    createFurtherMidNodes(new_internal_node,this_node_upper_left_right_regions)*/

                // iteratively traverse the tree to update partials from the current node
                stack<Node*> node_stack;
                node_stack.push(selected_node);
                node_stack.push(subtree);
                node_stack.push(new_internal_node->neighbor);
                updatePartialLh(node_stack, cumulative_rate, default_blength, min_blength, max_blength);
            }
        }
    }
    
    // delete best_parent_regions and best_child_regions
    if (best_parent_regions)
        delete best_parent_regions;
    if (best_child_regions)
        delete best_child_regions;
}

void Tree::placeNewSample(Node* selected_node, SeqRegions* sample, const string &seq_name, RealNumType best_lh_diff , bool is_mid_branch, RealNumType best_up_lh_diff, RealNumType best_down_lh_diff, Node* best_child, RealNumType* cumulative_rate, vector< vector<PositionType> > &cumulative_base, RealNumType default_blength, RealNumType max_blength, RealNumType min_blength)
{
    // dummy variables
    RealNumType best_child_lh;
    RealNumType best_child_split = 0;
    RealNumType best_parent_lh;
    RealNumType best_parent_split = 0;
    SeqRegions* best_parent_regions = NULL;
    SeqRegions* best_child_regions = NULL;
    RealNumType best_root_blength = -1;
    StateType num_states = aln.num_states;
    RealNumType threshold_prob = params->threshold_prob;
    
    // try to place the new sample as a descendant of a mid-branch point
    if (is_mid_branch)
    {
        SeqRegions* upper_left_right_regions = selected_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
        RealNumType best_split = 0.5;
        RealNumType best_split_lh = best_lh_diff;
        RealNumType new_split = 0.25;
        if (best_child_regions) delete best_child_regions;
        best_child_regions = new SeqRegions(selected_node->mid_branch_lh);
        //selected_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate)->mergeUpperLower(best_child_regions, selected_node->length * 0.5, selected_node->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate), selected_node->length * 0.5, aln, model, threshold_prob);
        SeqRegions* new_parent_regions = NULL;
        SeqRegions* lower_regions = selected_node->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
        RealNumType new_branch_length_split = new_split * selected_node->length;

        // try different positions on the existing branch
        while (new_branch_length_split > min_blength)
        {
            upper_left_right_regions->mergeUpperLower(new_parent_regions, new_branch_length_split, *lower_regions,  selected_node->length - new_branch_length_split, aln, model, threshold_prob);
            RealNumType placement_cost = calculateSamplePlacementCost(cumulative_rate, new_parent_regions, sample, default_blength);
            
            if (placement_cost > best_split_lh)
            {
                best_split_lh = placement_cost;
                best_split = new_split;
                
                if (best_child_regions) delete best_child_regions;
                best_child_regions = new_parent_regions;
                new_parent_regions = NULL;
            }
            else
                break;
            
            new_split = best_split * 0.5;
            new_branch_length_split = new_split * selected_node->length;
        }
        
        if (best_split > 0.49)
        {
            new_split = 0.25;
            new_branch_length_split = new_split * selected_node->length;
            while (new_branch_length_split > min_blength)
            {
                upper_left_right_regions->mergeUpperLower(new_parent_regions, selected_node->length - new_branch_length_split, *lower_regions, new_branch_length_split, aln, model, threshold_prob);
                
                RealNumType placement_cost = calculateSamplePlacementCost(cumulative_rate, new_parent_regions, sample, default_blength);
                if (placement_cost > best_split_lh)
                {
                    best_split_lh = placement_cost;
                    best_split = new_split;
                    
                    if (best_child_regions) delete best_child_regions;
                    best_child_regions = new_parent_regions;
                    new_parent_regions = NULL;
                }
                else
                    break;
                
                new_split = best_split * 0.5;
                new_branch_length_split = new_split * selected_node->length;
            }
            if (best_split < 0.49)
                best_split = 1.0 - best_split;
        }
        
        // delete new_parent_regions
        if (new_parent_regions) delete new_parent_regions;
        
        // now try different lengths for the new branch
        RealNumType new_branch_lh = best_split_lh;
        RealNumType best_blength = default_blength;
        while (best_blength > min_blength)
        {
            RealNumType new_blength = best_blength * 0.5;
            RealNumType placement_cost = calculateSamplePlacementCost(cumulative_rate, best_child_regions, sample, new_blength);
            
            if (placement_cost > new_branch_lh)
            {
                new_branch_lh = placement_cost;
                best_blength = new_blength;
            }
            else
                break;
        }
        if (best_blength > 0.7 * default_blength)
        {
            while (best_blength < max_blength)
            {
                RealNumType new_blength = best_blength * 2;
                RealNumType placement_cost = calculateSamplePlacementCost(cumulative_rate, best_child_regions, sample, new_blength);
                if (placement_cost > new_branch_lh)
                {
                    new_branch_lh = placement_cost;
                    best_blength = new_blength;
                }
                else
                    break;
            }
        }
        if (best_blength < min_blength)
        {
            RealNumType zero_branch_lh = calculateSamplePlacementCost(cumulative_rate, best_child_regions, sample, -1);
            if (zero_branch_lh > new_branch_lh)
                best_blength = -1;
        }
        
        // create new internal node and append child to it
        Node* new_internal_node = new Node(true);
        Node* next_node_1 = new Node();
        Node* next_node_2 = new Node();
        Node* new_sample_node = new Node(seq_name);
        
        new_internal_node->next = next_node_2;
        next_node_2->next = next_node_1;
        next_node_1->next = new_internal_node;
        
        new_internal_node->neighbor = selected_node->neighbor;
        selected_node->neighbor->neighbor = new_internal_node;
        RealNumType top_distance = selected_node->length * best_split;
        new_internal_node->length = top_distance;
        new_internal_node->neighbor->length = top_distance;
        
        selected_node->neighbor = next_node_2;
        next_node_2->neighbor = selected_node;
        RealNumType down_distance = selected_node->length - top_distance;
        selected_node->length = down_distance;
        selected_node->neighbor->length = down_distance;
        
        new_sample_node->neighbor = next_node_1;
        next_node_1->neighbor = new_sample_node;
        new_sample_node->length = best_blength;
        new_sample_node->neighbor->length = best_blength;
        
        new_sample_node->partial_lh = sample;
        next_node_1->partial_lh = best_child_regions;
        best_child_regions = NULL;
        upper_left_right_regions->mergeUpperLower(next_node_2->partial_lh, new_internal_node->length, *sample, best_blength, aln, model, threshold_prob);
        selected_node->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate)->mergeTwoLowers(new_internal_node->partial_lh, selected_node->length, sample, best_blength, aln, model, threshold_prob, cumulative_rate);
        RealNumType half_branch_length = new_internal_node->length * 0.5;
        upper_left_right_regions->mergeUpperLower(new_internal_node->mid_branch_lh, half_branch_length, *new_internal_node->partial_lh, half_branch_length, aln, model, threshold_prob);
        new_internal_node->computeTotalLhAtNode(aln, model, threshold_prob, cumulative_rate, new_internal_node == root);
        
        if (!new_internal_node->total_lh || new_internal_node->total_lh->empty())
            outError("Problem, None vector when placing sample, below node");
        
        if (best_blength > 0)
        {
            new_sample_node->computeTotalLhAtNode(aln, model, threshold_prob, cumulative_rate, new_sample_node == root);
            RealNumType half_branch_length = new_sample_node->length * 0.5;
            next_node_1->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate)->mergeUpperLower(new_sample_node->mid_branch_lh, half_branch_length, *sample, half_branch_length, aln, model, threshold_prob);
        }
        
        // update pseudo_count
        model.updatePesudoCount(aln, *next_node_1->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate), *sample);

        // iteratively traverse the tree to update partials from the current node
        stack<Node*> node_stack;
        node_stack.push(selected_node);
        node_stack.push(new_internal_node->neighbor);
        updatePartialLh(node_stack, cumulative_rate, default_blength, min_blength, max_blength);
    }
    // otherwise, best lk so far is for appending directly to existing node
    else
    {
        // place the new sample as a descendant of an existing node
        if (best_child)
        {
            RealNumType best_split = 0.5;
            RealNumType best_split_lh = best_down_lh_diff;
            SeqRegions* upper_left_right_regions = best_child->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
            SeqRegions* lower_regions = best_child->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
            if (best_child_regions) delete best_child_regions;
            best_child_regions = new SeqRegions(best_child->mid_branch_lh);
            //upper_left_right_regions->mergeUpperLower(best_child_regions, best_child->length * 0.5, lower_regions, best_child->length * 0.5, aln, model, threshold_prob);
            RealNumType new_split = 0.25;
            RealNumType new_branch_length_split = new_split * best_child->length;
            
            SeqRegions* new_parent_regions = NULL;
            while (new_branch_length_split > min_blength)
            {
                upper_left_right_regions->mergeUpperLower(new_parent_regions, new_branch_length_split, *lower_regions, best_child->length - new_branch_length_split, aln, model, threshold_prob);
                
                RealNumType placement_cost = calculateSamplePlacementCost(cumulative_rate, new_parent_regions, sample, default_blength);
                
                if (placement_cost > best_split_lh)
                {
                    best_split_lh = placement_cost;
                    best_split = new_split;
                    
                    if (best_child_regions) delete best_child_regions;
                    best_child_regions = new_parent_regions;
                    new_parent_regions = NULL;
                }
                else
                    break;
                
                new_split = best_split * 0.5;
                new_branch_length_split = new_split * best_child->length;
            }
            // delete new_parent_regions
            if (new_parent_regions) delete new_parent_regions;
            
            best_child_lh = best_split_lh;
            best_child_split = best_split;
        }
        else
            best_child_lh = MIN_NEGATIVE;
        
        // if node is root, try to place as sibling of the current root.
        RealNumType old_root_lh = MIN_NEGATIVE;
        if (root == selected_node)
        {
            old_root_lh = selected_node->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate)->computeAbsoluteLhAtRoot(aln, model, cumulative_base);
            RealNumType new_root_lh;
            SeqRegions* merged_root_sample_regions = NULL;
            SeqRegions* lower_regions = selected_node->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
            
            // merge 2 lower vector into one
            new_root_lh = lower_regions->mergeTwoLowers(merged_root_sample_regions, default_blength, sample, default_blength, aln, model, threshold_prob, cumulative_rate, true);
            
            new_root_lh += merged_root_sample_regions->computeAbsoluteLhAtRoot(aln, model, cumulative_base);
            best_parent_lh = new_root_lh - old_root_lh;
            best_root_blength = default_blength;
            
            if (best_parent_regions) delete best_parent_regions;
            best_parent_regions = merged_root_sample_regions;
            merged_root_sample_regions = NULL;
            
            RealNumType new_blength = 0.5 * default_blength;
            
            while (new_blength > min_blength)
            {
                // merge 2 lower vector into one
                new_root_lh = lower_regions->mergeTwoLowers(merged_root_sample_regions, new_blength, sample, default_blength, aln, model, threshold_prob, cumulative_rate, true);
                
                new_root_lh += merged_root_sample_regions->computeAbsoluteLhAtRoot(aln, model, cumulative_base);
                RealNumType diff_root_lh = new_root_lh-old_root_lh;
                if (diff_root_lh > best_parent_lh)
                {
                    best_parent_lh = diff_root_lh;
                    best_root_blength = new_blength;
                    
                    if (best_parent_regions) delete best_parent_regions;
                    best_parent_regions = merged_root_sample_regions;
                    merged_root_sample_regions = NULL;
                }
                else
                    break;
                
                new_blength = best_root_blength * 0.5;
            }
            
            // delete merged_root_sample_regions
            if (merged_root_sample_regions) delete merged_root_sample_regions;
        }
        // selected_node is not root
        else
        {
            RealNumType best_split = 0.5;
            RealNumType best_split_lh = best_up_lh_diff;
            SeqRegions* upper_left_right_regions = selected_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
            SeqRegions* lower_regions = selected_node->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
            if (best_parent_regions) delete best_parent_regions;
            best_parent_regions = new SeqRegions(selected_node->mid_branch_lh);
            // upper_left_right_regions->mergeUpperLower(best_parent_regions, selected_node->length * 0.5, lower_regions, selected_node->length * 0.5, aln, model, threshold_prob);
            RealNumType new_split = 0.25;
            RealNumType new_branch_length_split = new_split * selected_node->length;
            
            SeqRegions* new_parent_regions = NULL;
            while (new_branch_length_split > min_blength)
            {
                upper_left_right_regions->mergeUpperLower(new_parent_regions, selected_node->length - new_branch_length_split, *lower_regions, new_branch_length_split, aln, model, threshold_prob);
                
                RealNumType placement_cost = calculateSamplePlacementCost(cumulative_rate, new_parent_regions, sample, default_blength);
                
                if (placement_cost > best_split_lh)
                {
                    best_split_lh = placement_cost;
                    best_split = new_split;
                    if (best_parent_regions) delete best_parent_regions;
                    best_parent_regions = new_parent_regions;
                    new_parent_regions = NULL;
                }
                else
                    break;
                
                new_split = best_split * 0.5;
                new_branch_length_split = new_split * selected_node->length;
            }
            // delete new_parent_regions
            if (new_parent_regions) delete new_parent_regions;
            
            best_parent_lh = best_split_lh;
            best_parent_split = best_split;
        }
        
        // if the best placement is below the selected_node => add an internal node below the selected_node
        if (best_child_lh >= best_parent_lh && best_child_lh >= best_lh_diff)
        {
            ASSERT(best_child);
            
            SeqRegions* upper_left_right_regions = best_child->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);

            RealNumType new_branch_length_lh = best_child_lh;
            RealNumType best_length = default_blength;
            
            while (best_length > min_blength)
            {
                RealNumType new_blength = best_length * 0.5;
                RealNumType placement_cost = calculateSamplePlacementCost(cumulative_rate, best_child_regions, sample, new_blength);
                
                if (placement_cost > new_branch_length_lh)
                {
                    new_branch_length_lh = placement_cost;
                    best_length = new_blength;
                }
                else
                    break;
            }
            
            if (best_length > 0.7 * default_blength)
            {
                while (best_length < max_blength)
                {
                    RealNumType new_blength = best_length * 2;
                    RealNumType placement_cost = calculateSamplePlacementCost(cumulative_rate, best_child_regions, sample, new_blength);
                    if (placement_cost > new_branch_length_lh)
                    {
                        new_branch_length_lh = placement_cost;
                        best_length = new_blength;
                    }
                    else
                        break;
                }
            }
            
            if (best_length < min_blength)
            {
                RealNumType tmp_lh = calculateSamplePlacementCost(cumulative_rate, best_child_regions, sample, -1);
                
                if (tmp_lh > new_branch_length_lh)
                    best_length = -1;
            }
            
            // create new internal node and append child to it
            Node* new_internal_node = new Node(true);
            Node* next_node_1 = new Node();
            Node* next_node_2 = new Node();
            Node* new_sample_node = new Node(seq_name);
            
            new_internal_node->next = next_node_2;
            next_node_2->next = next_node_1;
            next_node_1->next = new_internal_node;
            
            new_internal_node->neighbor = best_child->neighbor;
            best_child->neighbor->neighbor = new_internal_node;
            RealNumType top_distance = best_child->length * best_child_split;
            new_internal_node->length = top_distance;
            new_internal_node->neighbor->length = top_distance;
                
            best_child->neighbor = next_node_2;
            next_node_2->neighbor = best_child;
            RealNumType down_distance = best_child->length * (1 - best_child_split);
            best_child->length = down_distance;
            best_child->neighbor->length = down_distance;
            
            new_sample_node->neighbor = next_node_1;
            next_node_1->neighbor = new_sample_node;
            new_sample_node->length = best_length;
            new_sample_node->neighbor->length = best_length;
            
            new_sample_node->partial_lh = sample;
            next_node_1->partial_lh = best_child_regions;
            best_child_regions = NULL;
            upper_left_right_regions->mergeUpperLower(next_node_2->partial_lh, new_internal_node->length, *sample, best_length, aln, model, threshold_prob);
            best_child->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate)->mergeTwoLowers(new_internal_node->partial_lh, best_child->length, sample, best_length, aln, model, threshold_prob, cumulative_rate);
            RealNumType half_branch_length = new_internal_node->length * 0.5;
            upper_left_right_regions->mergeUpperLower(new_internal_node->mid_branch_lh, half_branch_length, *new_internal_node->partial_lh, half_branch_length, aln, model, threshold_prob);
            new_internal_node->computeTotalLhAtNode(aln, model, threshold_prob, cumulative_rate, new_internal_node == root);
            
            /*if (top_distance >= 2 * min_blengthForMidNode)
                createFurtherMidNodes(new_internal_node,this_node_upper_left_right_regions)*/
            if (best_length > 0)
            {
                new_sample_node->computeTotalLhAtNode(aln, model, threshold_prob, cumulative_rate, new_sample_node == root);
                RealNumType half_branch_length = new_sample_node->length * 0.5;
                next_node_1->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate)->mergeUpperLower(new_sample_node->mid_branch_lh, half_branch_length, *sample, half_branch_length, aln, model, threshold_prob);
                /*if best_length>=2*min_blengthForMidNode:
                    createFurtherMidNodes(new_sample_node,best_child_regions)*/
            }
            
            // update pseudo_count
            model.updatePesudoCount(aln, *next_node_1->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate), *sample);

            // iteratively traverse the tree to update partials from the current node
            stack<Node*> node_stack;
            node_stack.push(best_child);
            node_stack.push(new_internal_node->neighbor);
            updatePartialLh(node_stack, cumulative_rate, default_blength, min_blength, max_blength);
        }
        // otherwise, add new parent to the selected_node
        else
        {
            // new parent is actually part of a polytomy since best placement is exactly at the node
            if (best_lh_diff >= best_parent_lh)
            {
                best_root_blength = -1;
                best_parent_split = -1;
                best_parent_lh = best_lh_diff;
                if (best_parent_regions) delete best_parent_regions;
                best_parent_regions = NULL;
                
                if (selected_node == root)
                    selected_node->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate)->mergeTwoLowers(best_parent_regions, -1, sample, default_blength, aln, model, threshold_prob, cumulative_rate);
                else
                    best_parent_regions = new SeqRegions(selected_node->total_lh);
            }

            // add parent to the root
            if (selected_node == root)
            {
                // now try different lengths for right branch
                RealNumType best_length2 = default_blength;
                SeqRegions* new_root_lower_regions = NULL;
                RealNumType new_root_lh = 0;
                
                while (best_length2 > min_blength)
                {
                    RealNumType new_blength = best_length2 * 0.5;
                    
                    new_root_lh = selected_node->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate)->mergeTwoLowers(new_root_lower_regions, best_root_blength, sample, new_blength, aln, model, threshold_prob, cumulative_rate, true);
                    
                    new_root_lh += new_root_lower_regions->computeAbsoluteLhAtRoot(aln, model, cumulative_base);
                    
                    RealNumType root_lh_diff = new_root_lh - old_root_lh;
                    if (root_lh_diff > best_parent_lh)
                    {
                        best_parent_lh = root_lh_diff;
                        best_length2 = new_blength;
                        
                        if (best_parent_regions) delete best_parent_regions;
                        best_parent_regions = new_root_lower_regions;
                        new_root_lower_regions = NULL;
                    }
                    else
                        break;
                }
                
                if (best_length2 > 0.7 * default_blength)
                {
                    while (best_length2 < max_blength)
                    {
                        RealNumType new_blength = best_length2 * 2;
                        new_root_lh = selected_node->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate)->mergeTwoLowers(new_root_lower_regions, best_root_blength, sample, new_blength, aln, model, threshold_prob, cumulative_rate, true);
                        new_root_lh += new_root_lower_regions->computeAbsoluteLhAtRoot(aln, model, cumulative_base);
                        RealNumType root_lh_diff = new_root_lh - old_root_lh;
                        
                        if (root_lh_diff > best_parent_lh)
                        {
                            best_parent_lh = root_lh_diff;
                            best_length2 = new_blength;
                            
                            if (best_parent_regions) delete best_parent_regions;
                            best_parent_regions = new_root_lower_regions;
                            new_root_lower_regions = NULL;
                        }
                        else
                            break;
                    }
                }
                
                // try with length zero
                if (best_length2 < min_blength)
                {
                    new_root_lh = selected_node->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate)->mergeTwoLowers(new_root_lower_regions, best_root_blength, sample, -1, aln, model, threshold_prob, cumulative_rate, true);
                    new_root_lh += new_root_lower_regions->computeAbsoluteLhAtRoot(aln, model, cumulative_base);
                    RealNumType root_lh_diff = new_root_lh - old_root_lh;
                    if (root_lh_diff > best_parent_lh)
                    {
                        best_length2 = -1;
                        // best_parent_lh = root_lh_diff; // best_parent_lh will never be read
                        if (best_parent_regions) delete best_parent_regions;
                        best_parent_regions = new_root_lower_regions;
                        new_root_lower_regions = NULL;
                    }
                }

                // delete new_root_lower_regions
                if (new_root_lower_regions) delete new_root_lower_regions;
                
                // add new root node into tree
                Node* new_root = new Node(true);
                Node* next_node_1 = new Node();
                Node* next_node_2 = new Node();
                Node* new_sample_node = new Node(seq_name);
                
                new_root->next = next_node_2;
                next_node_2->next = next_node_1;
                next_node_1->next = new_root;
                
                // attach the left child
                selected_node->neighbor = next_node_2;
                next_node_2->neighbor = selected_node;
                selected_node->length = best_root_blength;
                selected_node->neighbor->length = best_root_blength;
                
                if (best_root_blength <= 0)
                {
                    delete selected_node->total_lh;
                    selected_node->total_lh = NULL;
                    
                    if (selected_node->mid_branch_lh) delete selected_node->mid_branch_lh;
                    selected_node->mid_branch_lh = NULL;
                    //selected_node.furtherMidNodes=None
                }
                
                // attach the right child
                new_sample_node->neighbor = next_node_1;
                next_node_1->neighbor = new_sample_node;
                new_sample_node->length = best_length2;
                new_sample_node->neighbor->length = best_length2;
                
                new_root->partial_lh = best_parent_regions;
                best_parent_regions = NULL;
                new_root->total_lh = new_root->computeTotalLhAtNode(aln, model, threshold_prob, cumulative_rate, true);

                next_node_1->partial_lh = selected_node->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate)->computeTotalLhAtRoot(num_states, model, best_root_blength);
                next_node_2->partial_lh = sample->computeTotalLhAtRoot(num_states, model, best_length2);
                
                new_sample_node->partial_lh = sample;
                
                if (!new_root->total_lh || new_root->total_lh->size() == 0)
                {
                    outError("Problem, None vector when placing sample, new root");
                    /*print(merged_root_sample_regions)
                    print(node.probVect)
                    print(sample)
                    print(best_length2)
                    print(best_root_blength)*/
                }
                
                if (best_root_blength < 0)
                {
                    if (selected_node->total_lh)
                        delete selected_node->total_lh;
                    selected_node->total_lh = NULL;
                    
                    if (selected_node->mid_branch_lh)
                        delete selected_node->mid_branch_lh;
                    selected_node->mid_branch_lh = NULL;
                    /*node.furtherMidNodes=None*/
                }
                
                if (best_length2 > 0)
                {
                    new_sample_node->computeTotalLhAtNode(aln, model, threshold_prob, cumulative_rate, new_sample_node == root);
                    RealNumType half_branch_length = new_sample_node->length * 0.5;
                    next_node_1->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate)->mergeUpperLower(new_sample_node->mid_branch_lh, half_branch_length, *sample, half_branch_length, aln, model, threshold_prob);
                    
                    /*if best_length2>=2*min_blengthForMidNode:
                        createFurtherMidNodes(new_root.children[1],new_root.probVectUpLeft)*/
                }
                
                // update tree->root;
                root = new_root;
                
                // iteratively traverse the tree to update partials from the current node
                stack<Node*> node_stack;
                node_stack.push(selected_node);
                updatePartialLh(node_stack, cumulative_rate, default_blength, min_blength, max_blength);
            }
            //add parent to non-root node
            else
            {
                SeqRegions* upper_left_right_regions = selected_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
                
                // now try different lengths for the new branch
                RealNumType new_branch_length_lh = best_parent_lh;
                RealNumType best_length = default_blength;
                while (best_length > min_blength)
                {
                    RealNumType new_blength = best_length * 0.5;
                    RealNumType placement_cost = calculateSamplePlacementCost(cumulative_rate, best_parent_regions, sample, new_blength);
                    
                    if (placement_cost > new_branch_length_lh)
                    {
                        new_branch_length_lh = placement_cost;
                        best_length = new_blength;
                    }
                    else
                        break;
                }
                
                if (best_length > 0.7 * default_blength)
                {
                    while (best_length < max_blength)
                    {
                        RealNumType new_blength = best_length * 2;
                        RealNumType placement_cost = calculateSamplePlacementCost(cumulative_rate, best_parent_regions, sample, new_blength);
                        
                        if (placement_cost > new_branch_length_lh)
                        {
                            new_branch_length_lh = placement_cost;
                            best_length = new_blength;
                        }
                        else
                            break;
                    }
                }
                
                // try with length zero
                if (best_length < min_blength)
                {
                    RealNumType placement_cost = calculateSamplePlacementCost(cumulative_rate, best_parent_regions, sample, -1);
                    
                    if (placement_cost > new_branch_length_lh)
                        best_length = -1;
                }
                
                // now create new internal node and append child to it
                Node* new_internal_node = new Node(true);
                Node* next_node_1 = new Node();
                Node* next_node_2 = new Node();
                Node* new_sample_node = new Node(seq_name);
                
                new_internal_node->next = next_node_2;
                next_node_2->next = next_node_1;
                next_node_1->next = new_internal_node;
                
                RealNumType down_distance = selected_node->length * best_parent_split;
                RealNumType top_distance = selected_node->length - down_distance;
                if (best_parent_split < 0)
                {
                    down_distance = -1;
                    top_distance = selected_node->length;
                    
                    if (selected_node->total_lh) delete selected_node->total_lh;
                    selected_node->total_lh = NULL;
                    
                    if (selected_node->mid_branch_lh) delete selected_node->mid_branch_lh;
                    selected_node->mid_branch_lh = NULL;
                    
                    /*
                    node.furtherMidNodes=None*/
                }
                
                // attach to the parent node
                new_internal_node->neighbor = selected_node->neighbor;
                selected_node->neighbor->neighbor = new_internal_node;
                new_internal_node->length = top_distance;
                new_internal_node->neighbor->length = top_distance;
                
                // attach to the right child
                selected_node->neighbor = next_node_2;
                next_node_2->neighbor = selected_node;
                selected_node->length = down_distance;
                selected_node->neighbor->length = down_distance;
                
                // attach to the left child
                new_sample_node->neighbor = next_node_1;
                next_node_1->neighbor = new_sample_node;
                new_sample_node->length = best_length;
                new_sample_node->neighbor->length = best_length;
                
                new_sample_node->partial_lh = sample;
                next_node_1->partial_lh = best_parent_regions;
                best_parent_regions = NULL;
                upper_left_right_regions->mergeUpperLower(next_node_2->partial_lh, new_internal_node->length, *sample, best_length, aln, model, threshold_prob);
                selected_node->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate)->mergeTwoLowers(new_internal_node->partial_lh, selected_node->length, sample, best_length, aln, model, threshold_prob, cumulative_rate);
                RealNumType half_branch_length = new_internal_node->length * 0.5;
                upper_left_right_regions->mergeUpperLower(new_internal_node->mid_branch_lh, half_branch_length, *new_internal_node->partial_lh, half_branch_length, aln, model, threshold_prob);
                new_internal_node->computeTotalLhAtNode(aln, model, threshold_prob, cumulative_rate, new_internal_node == root);
                //mergeLhUpDown(new_internal_node->total_lh, best_parent_regions, -1, sample, best_length);
                
                if (!new_internal_node->total_lh || new_internal_node->total_lh->size() == 0)
                {
                    outError("Problem, None vector when placing sample, new parent");
                    /*print(best_parent_regions)
                    print(upper_left_right_regions)
                    print(sample)
                    print(best_length)
                    print(top_distance)
                    print(down_distance)*/
                }
                
                /*if (top_distance >= 2 * min_blengthForMidNode)
                    createFurtherMidNodes(new_internal_node,this_node_upper_left_right_regions)*/
                if (best_length > 0)
                {
                    new_sample_node->computeTotalLhAtNode(aln, model, threshold_prob, cumulative_rate, new_sample_node == root);
                    RealNumType half_branch_length = new_sample_node->length * 0.5;
                    next_node_1->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate)->mergeUpperLower(new_sample_node->mid_branch_lh, half_branch_length, *sample, half_branch_length, aln, model, threshold_prob);
                    /*if best_length>=2*min_blengthForMidNode:
                        createFurtherMidNodes(new_sample_node,best_parent_regions)*/
                }
                
                // update pseudo_count
                model.updatePesudoCount(aln, *next_node_1->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate), *sample);

                // iteratively traverse the tree to update partials from the current node
                stack<Node*> node_stack;
                node_stack.push(selected_node);
                node_stack.push(new_internal_node->neighbor);
                updatePartialLh(node_stack, cumulative_rate, default_blength, min_blength, max_blength);
            }
        }
    }
    
    // delete best_parent_regions and best_child_regions
    if (best_parent_regions)
        delete best_parent_regions;
    if (best_child_regions)
        delete best_child_regions;
}

void Tree::refreshAllLhs(RealNumType *cumulative_rate, RealNumType default_blength, RealNumType max_blength, RealNumType min_blength)
{
    // 1. update all the lower lhs along the tree
    refreshAllLowerLhs(cumulative_rate, default_blength, max_blength, min_blength);
    
    // 2. update all the non-lower lhs along the tree
    refreshAllNonLowerLhs(cumulative_rate, default_blength, max_blength, min_blength);
}

void Tree::refreshAllLowerLhs(RealNumType *cumulative_rate, RealNumType default_blength, RealNumType max_blength, RealNumType min_blength)
{
    // start from root
    Node* node = root;
    Node* last_node = NULL;
    
    // traverse to the deepest tip, update the lower lhs upward from the tips
    while (node)
    {
        // we reach a top node by a downward traversing
        if (node->is_top)
        {
            // if the current node is a leaf -> we reach the deepest tip -> traversing upward to update the lower lh of it parent
            if (node->isLeave())
            {
                last_node = node;
                node = node->neighbor;
            }
            // otherwise, keep traversing downward to find the deepest tip
            else
                node = node->next->neighbor;
        }
        // we reach the current node by an upward traversing from its children
        else
        {
            // if we reach the current node by an upward traversing from its first children -> traversing downward to its second children
            if (node->getTopNode()->next->neighbor == last_node)
                node = node->getTopNode()->next->next->neighbor;
            // otherwise, all children of the current node are updated -> update the lower lh of the current node
            else
            {
                // calculate the new lower lh of the current node from its children
                Node* top_node = node->getTopNode();
                Node* next_node_1 = top_node->next;
                Node* next_node_2 = next_node_1->next;
                
                SeqRegions* new_lower_lh = NULL;
                SeqRegions* lower_lh_1 = next_node_1->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
                SeqRegions* lower_lh_2 = next_node_2->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
                lower_lh_1->mergeTwoLowers(new_lower_lh, next_node_1->length, lower_lh_2, next_node_2->length, aln, model, params->threshold_prob, cumulative_rate);
                 
                // if new_lower_lh is NULL -> we need to update the branch lengths connecting the current node to its children
                if (!new_lower_lh)
                {
                    if (next_node_1->length <= 0)
                    {
                        stack<Node*> node_stack;
                        // NHANLT: note different from original maple
                        // updateBLen(nodeList,node,mutMatrix) -> the below codes update from next_node_1 instead of top_node
                        updateZeroBlength(next_node_1->neighbor, node_stack, params->threshold_prob, cumulative_rate, default_blength, min_blength, max_blength);
                        updatePartialLh(node_stack, cumulative_rate, default_blength, min_blength, max_blength);
                    }
                    else if (next_node_2->length <= 0)
                    {
                        stack<Node*> node_stack;
                        updateZeroBlength(next_node_2->neighbor, node_stack, params->threshold_prob, cumulative_rate, default_blength, min_blength, max_blength);
                        updatePartialLh(node_stack, cumulative_rate, default_blength, min_blength, max_blength);
                    }
                    else
                        outError("Strange, branch lengths > 0 but inconsistent lower lh creation in refreshAllLowerLhs()");
                }
                // otherwise, everything is good -> update the lower lh of the current node
                else
                {
                    delete top_node->partial_lh;
                    top_node->partial_lh = new_lower_lh;
                    new_lower_lh = NULL;
                }

                // delete new_lower_lh
                if (new_lower_lh)
                    delete new_lower_lh;
                
                // traverse upward to the parent of the current node
                last_node = top_node;
                node = top_node->neighbor;
            }
        }
    }
}

void Tree::refreshAllNonLowerLhs(RealNumType *cumulative_rate, RealNumType default_blength, RealNumType max_blength, RealNumType min_blength)
{
    // dummy variables
    StateType num_states = aln.num_states;
    RealNumType threshold_prob = params->threshold_prob;
    
    // start from the root
    Node* node = root;
    
    // update the total lh at root
    delete node->total_lh;
    node->total_lh = NULL;
    node->computeTotalLhAtNode(aln, model, params->threshold_prob, cumulative_rate, true);
    
    // if the root has children -> update its upper left/right lh then traverse downward to update non-lower lhs of other nodes
    if (!node->isLeave())
    {
        // update upper left/right lh of the root
        Node* next_node_1 = node->next;
        Node* next_node_2 = next_node_1->next;
        delete next_node_1->partial_lh;
        next_node_1->partial_lh = next_node_2->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate)->computeTotalLhAtRoot(num_states, model, next_node_2->length);
        delete next_node_2->partial_lh;
        next_node_2->partial_lh = next_node_1->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate)->computeTotalLhAtRoot(num_states, model, next_node_1->length);
        
        // traverse the tree downward and update the non-lower genome lists for all other nodes of the tree.
        Node* last_node = NULL;
        node = next_node_1->neighbor;
        while (node)
        {
            // we reach a top node by a downward traversing
            if (node->is_top)
            {
                SeqRegions* parent_upper_lr_lh = node->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
                
                // update the total lh, total lh at the mid-branch point of the current node
                if (node->length > 0)
                {
                    // update the total lh
                    if (node->total_lh)
                    {
                        delete node->total_lh;
                        node->total_lh = NULL;
                    }
                    node->computeTotalLhAtNode(aln, model, threshold_prob, cumulative_rate, node == root);
                    
                    if (!node->total_lh)
                        outError("Strange, inconsistent total lh creation in refreshAllNonLowerLhs()");
                    else
                    {
                        // update total lh at the mid-branch point
                        if (node->mid_branch_lh)
                        {
                            delete node->mid_branch_lh;
                            node->mid_branch_lh = NULL;
                        }
                        SeqRegions* lower_lh = node->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
                        RealNumType half_branch_length = node->length * 0.5;
                        parent_upper_lr_lh->mergeUpperLower(node->mid_branch_lh, half_branch_length, *lower_lh, half_branch_length, aln, model, threshold_prob);
                        
                        /*if node.dist>=2*minBLenForMidNode:
                            createFurtherMidNodes(node,upper_lr_lh)*/
                    }
                }
                
                // if the current node is an internal node (~having children) -> update its upper left/right lh then traverse downward to update non-lower lhs of other nodes
                if (!node->isLeave())
                {
                    Node* next_node_1 = node->next;
                    Node* next_node_2 = next_node_1->next;
                    
                    // recalculate the FIRST upper left/right lh of the current node
                    SeqRegions* new_upper_lr_lh = NULL;
                    SeqRegions* lower_lh = next_node_2->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
                    parent_upper_lr_lh->mergeUpperLower(new_upper_lr_lh, node->length, *lower_lh, next_node_2->length, aln, model, threshold_prob);
                    
                    // if the upper left/right lh is null -> try to increase the branch length
                    if (!new_upper_lr_lh)
                    {
                        if (next_node_2->length <= 0)
                        {
                            stack<Node*> node_stack;
                            updateZeroBlength(next_node_2->neighbor, node_stack, params->threshold_prob, cumulative_rate, default_blength, min_blength, max_blength);
                            updatePartialLh(node_stack, cumulative_rate, default_blength, min_blength, max_blength);
                        }
                        else if (node->length <= 0)
                        {
                            stack<Node*> node_stack;
                            updateZeroBlength(node, node_stack, params->threshold_prob, cumulative_rate, default_blength, min_blength, max_blength);
                            updatePartialLh(node_stack, cumulative_rate, default_blength, min_blength, max_blength);
                        }
                        else
                            outError("Strange, inconsistent upper left/right lh creation in refreshAllNonLowerLhs()");
                    }
                    // otherwise, everything is good -> update upper left/right lh of the current node
                    else
                    {
                        if (next_node_1->partial_lh) delete next_node_1->partial_lh;
                        next_node_1->partial_lh = new_upper_lr_lh;
                        new_upper_lr_lh = NULL;
                    }
                    
                    // recalculate the SECOND upper left/right lh of the current node
                    lower_lh = next_node_1->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
                    parent_upper_lr_lh->mergeUpperLower(new_upper_lr_lh, node->length, *lower_lh, next_node_1->length, aln, model, threshold_prob);
                    
                    // if the upper left/right lh is null -> try to increase the branch length
                    if (!new_upper_lr_lh)
                    {
                        if (next_node_1->length <= 0)
                        {
                            stack<Node*> node_stack;
                            updateZeroBlength(next_node_1->neighbor, node_stack, params->threshold_prob, cumulative_rate, default_blength, min_blength, max_blength);
                            updatePartialLh(node_stack, cumulative_rate, default_blength, min_blength, max_blength);
                        }
                        else if (node->length <= 0)
                        {
                            stack<Node*> node_stack;
                            updateZeroBlength(node, node_stack, params->threshold_prob, cumulative_rate, default_blength, min_blength, max_blength);
                            updatePartialLh(node_stack, cumulative_rate, default_blength, min_blength, max_blength);
                        }
                        else
                            outError("Strange, inconsistent upper left/right lh creation in refreshAllNonLowerLhs()");
                    }
                    // otherwise, everything is good -> update upper left/right lh of the current node
                    else
                    {
                        if (next_node_2->partial_lh) delete next_node_2->partial_lh;
                        next_node_2->partial_lh = new_upper_lr_lh;
                        new_upper_lr_lh = NULL;
                    }
                    
                    // delete new_upper_lr_lh
                    if (new_upper_lr_lh) delete new_upper_lr_lh;
                    
                    // keep traversing downward to its firt child
                    node = next_node_1->neighbor;
                }
                // if the current node is a leaf -> traverse upward to its parent
                else
                {
                    last_node = node;
                    node = node->neighbor;
                }
            }
            // we reach the current node by an upward traversing from its children
            else
            {
                Node* top_node = node->getTopNode();
                Node* next_node_1 = top_node->next;
                Node* next_node_2 = next_node_1->next;
                
                // if we reach the current node by an upward traversing from its first children -> traversing downward to its second children
                if (last_node == next_node_1->neighbor)
                    node = next_node_2->neighbor;
                // otherwise, all children of the current node are updated -> update the lower lh of the current node
                else
                {
                    last_node = top_node;
                    node = top_node->neighbor;
                }
            }
        }
    }
}

void Tree::setAllNodeOutdated()
{
    // start from the root
    stack<Node*> node_stack;
    node_stack.push(root);
    
    // traverse downward to set all descentdant outdated
    while (!node_stack.empty())
    {
        // pick the top node from the stack
        Node* node = node_stack.top();
        node_stack.pop();
        
        // set the current node outdated
        node->outdated = true;
        
        // traverse downward
        Node* neighbor_node;
        FOR_NEIGHBOR(node, neighbor_node)
            node_stack.push(neighbor_node);
    }
}

RealNumType Tree::improveEntireTree(bool short_range_search, RealNumType *cumulative_rate, vector< vector<PositionType> > &cumulative_base, RealNumType default_blength, RealNumType max_blength, RealNumType min_blength, RealNumType min_blength_mid)
{
    // start from the root
    stack<Node*> node_stack;
    node_stack.push(root);
    
    // dummy variables
    RealNumType total_improvement = 0;
    PositionType num_nodes = 0;
    
    // traverse downward the tree
    while (!node_stack.empty())
    {
        // pick the top node from the stack
        Node* node = node_stack.top();
        node_stack.pop();
        
        // add all children of the current nodes to the stack for further traversing later
        Node* neighbor_node = NULL;
        FOR_NEIGHBOR(node, neighbor_node)
            node_stack.push(neighbor_node);
        
        // only process outdated node to avoid traversing the same part of the tree multiple times
        if (node->outdated)
        {
            node->outdated = false;
            
            /*if checkEachSPR:
                root=node
                while root.up!=None:
                    root=root.up
                #print("Pre-SPR tree: "+createBinaryNewick(root))
                oldTreeLK=calculateTreeLikelihood(root,mutMatrix,checkCorrectness=True)
                #print("Pre-SPR tree likelihood: "+str(oldTreeLK))
                reCalculateAllGenomeLists(root,mutMatrix, checkExistingAreCorrect=True)*/
            
            // do SPR moves to improve the tree
            RealNumType improvement = improveSubTree(node, short_range_search, cumulative_rate, cumulative_base, default_blength, max_blength, min_blength, min_blength_mid);
            
            /*if checkEachSPR:
                #print(" apparent improvement "+str(improvement))
                root=node
                while root.up!=None:
                    root=root.up
                #print("Post-SPR tree: "+createBinaryNewick(root))
                newTreeLK=calculateTreeLikelihood(root,mutMatrix)
                reCalculateAllGenomeLists(root,mutMatrix, checkExistingAreCorrect=True)
                #print("Post-SPR tree likelihood: "+str(newTreeLK))
                if newTreeLK-oldTreeLK < improvement-1.0:
                    print("In startTopologyUpdates, LK score of improvement "+str(newTreeLK)+" - "+str(oldTreeLK)+" = "+str(newTreeLK-oldTreeLK)+" is less than what is supposd to be "+str(improvement))
                    exit()
             */
            
            // update total_improvement
            total_improvement += improvement;
            
            // Show log every 1000 nodes
            num_nodes += 1;
            if (num_nodes % 1000 == 0)
                cout << "Processed topology for " << convertIntToString(num_nodes) << " nodes." << endl;
        }
    }
    
    return total_improvement;
}

PositionType Tree::optimizeBranchLengths(RealNumType *cumulative_rate, RealNumType default_blength, RealNumType max_blength, RealNumType min_blength, RealNumType min_blength_sensitivity)
{
    // start from the root's children
    stack<Node*> node_stack;
    if (!root || !root->next)
        return 0;
    Node* neighbor_node = NULL;
    FOR_NEIGHBOR(root, neighbor_node)
        node_stack.push(neighbor_node);
    
    // dummy variables
    PositionType num_improvement = 0;
    RealNumType threshold_prob = params->threshold_prob;
    
    // traverse downward the tree
    while (!node_stack.empty())
    {
        // pick the top node from the stack
        Node* node = node_stack.top();
        node_stack.pop();
        
        SeqRegions* upper_lr_regions = node->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
        SeqRegions* lower_regions = node->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
        
        // add all children of the current nodes to the stack for further traversing later
        neighbor_node = NULL;
        FOR_NEIGHBOR(node, neighbor_node)
            node_stack.push(neighbor_node);
        
        // only process outdated node to avoid traversing the same part of the tree multiple times
        if (node->outdated)
        {
            // estimate the branch length
            RealNumType best_length = estimateBranchLength(upper_lr_regions, lower_regions, cumulative_rate, min_blength_sensitivity);
            
            if (best_length > 0 || node->length > 0)
            {
                if (best_length <= 0 || node->length <= 0 || (node->length > 1.01 * best_length) || (node->length < 0.99 * best_length))
                {
                    node->length = best_length;
                    node->neighbor->length = node->length;
                    ++num_improvement;
                    
                    // update partial likelihood regions
                    stack<Node*> new_node_stack;
                    new_node_stack.push(node);
                    new_node_stack.push(node->neighbor);
                    updatePartialLh(new_node_stack, cumulative_rate, default_blength, min_blength, max_blength);
                }
            }
        }
    }
    
    return num_improvement;
}

RealNumType Tree::estimateBranchLength(const SeqRegions* const parent_regions, const SeqRegions* const child_regions, const RealNumType* const cumulative_rate, RealNumType min_blength_sensitivity)
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
        if (seq2_region->type == TYPE_N || seq1_region->type == TYPE_N)
        {
            pos = end_pos + 1;
            continue;
        }
        
        // e1.type != N && e2.type != N
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
        if (seq1_region->type == TYPE_R)
        {
            // 2.1. e1.type = R and e2.type = R
            // NHANLT NOTES:
            // coefficient = derivative of log likelihood function wrt t
            // l = log(1 + q_xx * t) ~ q_xx * t
            // => l' = q_xx
            if (seq2_region->type == TYPE_R)
                coefficient += cumulative_rate[end_pos + 1] - cumulative_rate[pos];
            // 2.2. e1.type = R and e2.type = O
            else if (seq2_region->type == TYPE_O)
            {
                StateType seq1_state = aln.ref_seq[end_pos];
                RealNumType* mutation_mat_row = model.mutation_mat + model.row_index[seq1_state];
                
                RealNumType coeff0 = seq2_region->getLH(seq1_state);
                RealNumType coeff1 = 0;
                
                if (seq1_region->plength_observation2root >= 0)
                {
                    coeff0 *= model.root_freqs[seq1_state];

                    RealNumType* transposed_mut_mat_row = model.transposed_mut_mat + model.row_index[seq1_state];
                                        
                    for (StateType i = 0; i < num_states; ++i)
                    {
                        coeff0 += model.root_freqs[i] * transposed_mut_mat_row[i] * seq1_region->plength_observation2node * seq2_region->getLH(i);
                        coeff1 += mutation_mat_row[i] * seq2_region->getLH(i);
                    }
                    
                    coeff1 *=  model.root_freqs[seq1_state];
                }
                else
                {
                    // NHANLT NOTES:
                    // x = seq1_state
                    // l = log(1 + q_xx * t + sum(q_xy * t)
                    // l' = [q_xx + sum(q_xy)]/[1 + q_xx * t + sum(q_xy * t)]
                    // coeff1 = numerator = q_xx + sum(q_xy)
                    for (StateType j = 0; j < num_states; ++j)
                        coeff1 += mutation_mat_row[j] * seq2_region->getLH(j);
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
            // 2.3. e1.type = R and e2.type = A/C/G/T
            else
            {
                if (seq1_region->plength_observation2root >= 0)
                {
                    StateType seq1_state = aln.ref_seq[end_pos];
                    StateType seq2_state = seq2_region->type;
                    
                    RealNumType coeff1 = model.root_freqs[seq1_state] * model.mutation_mat[model.row_index[seq1_state] + seq2_state];
                    RealNumType coeff0 = model.root_freqs[seq2_state] * model.mutation_mat[model.row_index[seq2_state] + seq1_state] * seq1_region->plength_observation2node;
                    
                    if (total_blength > 0)
                        coeff0 += coeff1 * total_blength;
                    
                    coefficient_vec.push_back(coeff0 / coeff1);
                }
                // NHANLT NOTES:
                // l = log(q_xy * t)
                // l' = q_xy / (q_xy * t) = 1 / t
                else if (total_blength > 0)
                    coefficient_vec.push_back(total_blength);
                else
                    coefficient_vec.push_back(0);
            }
        }
        // 3. e1.type = O
        else if (seq1_region->type == TYPE_O)
        {
            RealNumType coeff0 = 0;
            RealNumType coeff1 = 0;
            
            // 3.1. e1.type = O and e2.type = O
            if (seq2_region->type == TYPE_O)
            {
                RealNumType* mutation_mat_row = model.mutation_mat;
                
                // NHANLT NOTES:
                // l = log(sum_x(1 + q_xx * t + sum_y(q_xy * t)))
                // l' = [sum_x(q_xx + sum_y(q_xy))]/[sum_x(1 + q_xx * t + sum_y(q_xy * t))]
                // coeff1 = numerator = sum_x(q_xx + sum_y(q_xy))
                // coeff0 = denominator = sum_x(1 + q_xx * t + sum_y(q_xy * t))
                for (StateType i = 0; i < num_states; ++i, mutation_mat_row += num_states)
                {
                    RealNumType seq1_lh_i = seq1_region->getLH(i);
                    coeff0 += seq1_lh_i * seq2_region->getLH(i);
                    
                    for (StateType j = 0; j < num_states; ++j)
                        coeff1 += seq1_lh_i * seq2_region->getLH(j) * mutation_mat_row[j];
                }
            }
            // 3.2. e1.type = O and e2.type = R or A/C/G/T
            else
            {
                StateType seq2_state = seq2_region->type;
                if (seq2_state == TYPE_R)
                    seq2_state = aln.ref_seq[end_pos];
                
                coeff0 = seq1_region->getLH(seq2_state);

                // NHANLT NOTES:
                // y = seq2_state
                // l = log(1 + q_yy * t + sum_x(q_xy * t)
                // l' = [q_yy + sum_x(q_xy))]/[1 + q_xx * t + sum_y(q_xy * t)]
                // coeff1 = numerator = q_yy + sum_x(q_xy))
                // coeff0 = denominator = 1 + q_xx * t + sum_y(q_xy * t)
                RealNumType* transposed_mut_mat_row = model.transposed_mut_mat + model.row_index[seq2_state];
                for (StateType i = 0; i < num_states; ++i)
                    coeff1 += seq1_region->getLH(i) * transposed_mut_mat_row[i];
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
        // 4. e1.type = A/C/G/T
        else
        {
            int seq1_state = seq1_region->type;
            
            // 4.1. e1.type =  e2.type
            // NHANLT NOTES:
            // coefficient = derivative of log likelihood function wrt t
            // l = log(1 + q_xx * t) ~ q_xx * t
            // => l' = q_xx
            if (seq1_region->type == seq2_region->type)
                coefficient += model.diagonal_mut_mat[seq1_state];
            // e1.type = A/C/G/T and e2.type = O/A/C/G/T
            else
            {
                // 4.2. e1.type = A/C/G/T and e2.type = O
                if (seq2_region->type == TYPE_O)
                {
                    RealNumType coeff0 = seq2_region->getLH(seq1_state);
                    RealNumType coeff1 = 0;
                    
                    RealNumType* mutation_mat_row = model.mutation_mat + model.row_index[seq1_state];
                    
                    if (seq1_region->plength_observation2root >= 0)
                    {
                        coeff0 *= model.root_freqs[seq1_state];

                        RealNumType* transposed_mut_mat_row = model.transposed_mut_mat + model.row_index[seq1_state];
                                            
                        for (StateType i = 0; i < num_states; ++i)
                        {
                            coeff0 += model.root_freqs[i] * transposed_mut_mat_row[i] * seq1_region->plength_observation2node * seq2_region->getLH(i);
                            coeff1 += mutation_mat_row[i] * seq2_region->getLH(i);
                        }
                        
                        coeff1 *= model.root_freqs[seq1_state];
                    }
                    else
                    {
                        // NHANLT NOTES:
                        // x = seq1_state
                        // l = log(1 + q_xx * t + sum(q_xy * t)
                        // l' = [q_xx + sum(q_xy)]/[1 + q_xx * t + sum(q_xy * t)]
                        // coeff1 = numerator = q_xx + sum(q_xy)
                        for (StateType j = 0; j < num_states; ++j)
                            coeff1 += mutation_mat_row[j] * seq2_region->getLH(j);
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
                // 4.3. e1.type = A/C/G/T and e2.type = R or A/C/G/T
                else
                {
                    RealNumType coeff0 = 0;
                    
                    StateType seq2_state = seq2_region->type;
                    if (seq2_state == TYPE_R)
                        seq2_state = aln.ref_seq[end_pos];
                    
                    if (seq1_region->plength_observation2root >= 0)
                    {
                        coeff0 = model.root_freqs[seq2_state] * model.mutation_mat[model.row_index[seq2_state] + seq1_state] * seq1_region->plength_observation2node;
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
            }
        }
        
        // update pos
        pos = end_pos + 1;
    }
    
    // now optimized branch length based on coefficients
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

RealNumType Tree::calculateDerivative(vector<RealNumType> &coefficient_vec, RealNumType delta_t)
{
    RealNumType result = 0;
    
    for (RealNumType coefficient : coefficient_vec)
        result += 1.0 / (coefficient + delta_t);
        
    return result;
}

RealNumType Tree::improveSubTree(Node* node, bool short_range_search, RealNumType *cumulative_rate, vector< vector<PositionType> > &cumulative_base, RealNumType default_blength, RealNumType max_blength, RealNumType min_blength, RealNumType min_blength_mid)
{
    // dummy variables
    RealNumType threshold_prob = params->threshold_prob;
    RealNumType threshold_prob2 = params->threshold_prob2;
    RealNumType thresh_placement_cost = short_range_search ? params->thresh_placement_cost_short_search : params->thresh_placement_cost;
    RealNumType total_improvement = 0;
    bool blength_changed = false; // true if a branch length has been changed
    
    // we avoid the root node since it cannot be re-placed with SPR moves
    if (node != root)
    {
        // evaluate current placement
        Node* parent_node = node->neighbor->getTopNode();
        SeqRegions* parent_upper_lr_lh = node->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
        
        // score of current tree
        RealNumType best_blength = node->length;
        SeqRegions* lower_lh = node->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
        RealNumType best_lh = calculateSubTreePlacementCost(cumulative_rate, parent_upper_lr_lh, lower_lh, best_blength);
        RealNumType original_lh = best_lh;
        
        // optimize branch length
        if (best_lh < thresh_placement_cost)
        {
            // try different branch lengths for the current node placement (just in case branch length can be improved, in which case it counts both as tree improvment and better strategy to find a new placement).
            if (node->length <= 0)
            {
                best_blength = min_blength;
                best_lh = calculateSubTreePlacementCost(cumulative_rate, parent_upper_lr_lh, lower_lh, best_blength);
            }
            
            RealNumType best_split = 1;
            RealNumType new_split = 0.5;
            RealNumType new_bl_split = new_split * best_blength;
            while (new_bl_split > min_blength)
            {
                RealNumType new_lh = calculateSubTreePlacementCost(cumulative_rate, parent_upper_lr_lh, lower_lh, new_bl_split);
                
                if (new_lh > best_lh)
                {
                    best_lh = new_lh;
                    best_split = new_split;
                    blength_changed = true;
                }
                else
                    break;
                
                new_split = best_split * 0.5;
                new_bl_split = new_split * best_blength;
            }
            
            if (best_split > 0.7)
            {
                new_split = 2;
                RealNumType new_bl_split = new_split * best_blength;
                while (new_bl_split < max_blength)
                {
                    RealNumType new_lh = calculateSubTreePlacementCost(cumulative_rate, parent_upper_lr_lh, lower_lh, new_bl_split);
                    
                    if (new_lh > best_lh)
                    {
                        best_lh = new_lh;
                        best_split = new_split;
                        blength_changed = true;
                    }
                    else
                        break;
                    
                    new_split = best_split * 2;
                    new_bl_split = new_split * best_blength;
                }
            }
            
            best_blength = best_blength * best_split;
            
            if (node->length <= 0)
                if (original_lh > best_lh)
                    best_lh = original_lh;
        }
           
        // find new placement
        if (best_lh < thresh_placement_cost)
        {
            // now find the best place on the tree where to re-attach the subtree rooted at "node" but to do that we need to consider new vector probabilities after removing the node that we want to replace this is done using findBestParentTopology().
            bool topology_updated = false;
            
            Node* best_node = NULL;
            RealNumType best_lh_diff = best_lh;
            bool is_mid_node = false;
            RealNumType best_up_lh_diff = MIN_NEGATIVE;
            RealNumType best_down_lh_diff = MIN_NEGATIVE;
            Node* best_child = NULL;
            
            // seek a new placement for the subtree
            seekSubTreePlacement(best_node, best_lh_diff, is_mid_node, best_up_lh_diff, best_down_lh_diff, best_child, short_range_search, node, best_blength, cumulative_rate, default_blength, min_blength_mid, true, NULL);
            
            // validate the new placement cost
            if (best_lh_diff > threshold_prob2)
                outError("Strange, lh cost is positive");
            else if (best_lh_diff < -1e50)
                outError("Likelihood cost is very heavy, this might mean that the reference used is not the same used to generate the input diff file");
            
            if (best_lh_diff + thresh_placement_cost > best_lh)
            {
                if (best_node == parent_node)
                    outWarning("Strange, re-placement is at same node");
                else if ((best_node == parent_node->next->neighbor || best_node == parent_node->next->next->neighbor) && is_mid_node)
                    cout << "Re-placement is above sibling node";
                else
                {
                    // reach the top of a multifurcation, which is the only place in a multifurcatio where placement is allowed.
                    Node* top_polytomy = best_node;
                    while (top_polytomy->length <= 0 && top_polytomy!= root)
                        top_polytomy = top_polytomy->neighbor->getTopNode();
                    
                    if (top_polytomy != best_node)
                        outWarning("Strange, placement node not at top of polytomy");
                    
                    // reach the top of the multifurcation of the parent
                    Node* parent_top_polytomy = parent_node;
                    PositionType parent_top_count = 0;
                    while (parent_top_polytomy->length <= 0 && parent_top_polytomy != root)
                    {
                        parent_top_polytomy = parent_top_polytomy->neighbor->getTopNode();
                        parent_top_count += 1;
                    }
                    
                    if (!(parent_top_polytomy == top_polytomy && !is_mid_node))
                    {
                        total_improvement = best_lh_diff - best_lh;
                        
                        if (verbose_mode == VB_DEBUG)
                            cout << "In improveSubTree() found SPR move with improvement " << total_improvement << endl;
                        
                        // apply an SPR move
                        applySPR(node, best_node, is_mid_node, best_blength, best_lh_diff, cumulative_rate, cumulative_base, default_blength, max_blength, min_blength, min_blength_mid);
                        
                        topology_updated = true;
                    }
                }
                if (!topology_updated && blength_changed)
                {
                    node->length = best_blength;
                    node->neighbor->length = node->length;
                    
                    stack<Node*> node_stack;
                    node_stack.push(node);
                    node_stack.push(node->neighbor);
                    updatePartialLh(node_stack, cumulative_rate, default_blength, min_blength, max_blength);
                }
            }
            else if (blength_changed)
            {
                node->length = best_blength;
                node->neighbor->length = node->length;
                
                stack<Node*> node_stack;
                node_stack.push(node);
                node_stack.push(node->neighbor);
                updatePartialLh(node_stack, cumulative_rate, default_blength, min_blength, max_blength);
            }
        }
        else if (blength_changed)
        {
            node->length = best_blength;
            node->neighbor->length = node->length;
            
            stack<Node*> node_stack;
            node_stack.push(node);
            node_stack.push(node->neighbor);
            updatePartialLh(node_stack, cumulative_rate, default_blength, min_blength, max_blength);
        }
    }
                        
    return total_improvement;
}

RealNumType Tree::calculateSubTreePlacementCost(RealNumType* cumulative_rate, const SeqRegions* const parent_regions, const SeqRegions* const child_regions, RealNumType blength)
{
    return (this->*calculateSubTreePlacementCostPointer)(cumulative_rate, parent_regions, child_regions, blength);
}

// this implementation derives from appendProbNode
template <const StateType num_states>
RealNumType Tree::calculateSubTreePlacementCostTemplate(
  RealNumType* cumulative_rate, 
  const SeqRegions* const parent_regions, 
  const SeqRegions* const child_regions, 
  RealNumType blength)
{  // 55% of runtime
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

        // 2. e1.type = R
        if (seq1_region->type == TYPE_R)
        {
            // 2.1. e1.type = R and e2.type = R
            if (seq2_region->type == TYPE_R)
            {
                if (seq1_region->plength_observation2root >= 0)
                    total_blength += seq1_region->plength_observation2node;
                
                // NHANLT NOTE:
                // approximation log(1+x)~x
                if (total_blength > 0)
                    lh_cost += total_blength * (cumulative_rate[end_pos + 1] - cumulative_rate[pos]);
            }
            // 2.2. e1.type = R and e2.type = O
            else if (seq2_region->type == TYPE_O)
            {
                RealNumType tot = 0;
                StateType seq1_state = aln.ref_seq[end_pos];
                
                if (seq1_region->plength_observation2root >= 0)
                {
                    RealNumType* transposed_mut_mat_row = model.transposed_mut_mat + model.row_index[seq1_state];
                    RealNumType* mutation_mat_row = model.mutation_mat;
                                        
                    for (StateType i = 0; i < num_states; ++i, mutation_mat_row += num_states)
                    {
                        // NHANLT NOTE: UNSURE
                        // tot2: likelihood that we can observe seq1_state elvoving from i at root (account for the fact that the observation might have occurred on the other side of the phylogeny with respect to the root)
                        // tot2 = root_freqs[seq1_state] * (1 + mut[seq1_state,seq1_state] * plength_observation2node) + root_freqs[i] * mut[i,seq1_state] * plength_observation2node
                        RealNumType tot2;
                        
                        if (seq1_state == i)
                            tot2 = model.root_freqs[i] * (1.0 + transposed_mut_mat_row[i] * seq1_region->plength_observation2node);
                        else
                            tot2 = model.root_freqs[i] * (transposed_mut_mat_row[i] * seq1_region->plength_observation2node);
                        
                        // NHANLT NOTE:
                        // tot3: likelihood of i evolves to j
                        // tot3 = (1 + mut[i,i] * total_blength) * lh(seq2,i) + mut[i,j] * total_blength * lh(seq2,j)
                        RealNumType tot3 = 0;
                        if (total_blength > 0)
                        {
                            for (StateType j = 0; j < num_states; ++j)
                                tot3 += mutation_mat_row[j] * seq2_region->getLH(j);
                            
                            tot3 *= total_blength;
                        }
                        
                        // NHANLT NOTE:
                        // tot = tot2 * tot3
                        tot += tot2 * (seq2_region->getLH(i) + tot3);
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
                        tot += dotProduct<num_states>(mutation_mat_row, &((*seq2_region->likelihood)[0]));
                        tot *= total_blength;
                    }
                    tot += seq2_region->getLH(seq1_state);
                }
                total_factor *= tot;
            }
            // 2.3. e1.type = R and e2.type = A/C/G/T
            else
            {
                const StateType seq1_state = aln.ref_seq[end_pos];
                const StateType seq2_state = seq2_region->type;
                
                if (seq1_region->plength_observation2root >= 0)
                {
                    if (total_blength > 0)
                    {
                        // NHANLT NOTE: UNSURE
                        // seq1_state_evolves_seq2_state = (1) the likelihood that seq1_state evolves to seq2_state * (2) the likelihood that seq1_state unchanges from the observing position
                        // (1) = model.mutation_mat[model.row_index[seq1_state] + seq2_state] * total_blength
                        // (2) = (1.0 + model.diagonal_mut_mat[seq1_state] * seq1_region->plength_observation2node)
                        RealNumType seq1_state_evolves_seq2_state = model.mutation_mat[model.row_index[seq1_state] + seq2_state] * total_blength * (1.0 + model.diagonal_mut_mat[seq1_state] * seq1_region->plength_observation2node);
                        
                        // NHANLT NOTE: UNCLEAR
                        // consider the inverse process of the above
                        // seq2_state_evolves_seq1_state = (1) the likelihood that seq2_state evolves to seq1_state * (2) the likelihood that seq2_state unchanges from the observing position
                        // (1) = root_freqs[seq2_state] / root_freqs[seq1_state] * mutation_mat[model.row_index[seq2_state] + seq1_state] * seq1_region->plength_observation2node
                        // (2) = (1.0 + model.diagonal_mut_mat[seq2_state] * total_blength)
                        RealNumType seq2_state_evolves_seq1_state = model.freqi_freqj_qij[model.row_index[seq2_state] + seq1_state] * seq1_region->plength_observation2node * (1.0 + model.diagonal_mut_mat[seq2_state] * total_blength);
                        
                        total_factor *= seq1_state_evolves_seq2_state + seq2_state_evolves_seq1_state;
                    }
                    // NHANLT NOTE:
                    // the same as above but total_blength = 0 then we simplify the formula to save the runtime (avoid multiplying with 0)
                    else
                        total_factor *= model.freqi_freqj_qij[model.row_index[seq2_state] + seq1_state] * seq1_region->plength_observation2node;
                }
                // NHANLT NOTE:
                // add the likelihood that seq1_state evoles to seq2_state = mut[seq1_state,seq2_state] * total_blength
                else if (total_blength > 0)
                    total_factor *= model.mutation_mat[model.row_index[seq1_state] + seq2_state] * total_blength;
                else
                    return MIN_NEGATIVE;
            }
        }
        // 3. e1.type = O
        else if (seq1_region->type == TYPE_O)
        {
            // 3.1. e1.type = O and e2.type = O
            if (seq2_region->type == TYPE_O)
            {
                if (total_blength > 0)
                {
                  total_factor *= matrixEvolve<num_states>(&((*seq1_region->likelihood)[0]), &((*seq2_region->likelihood)[0]), model.mutation_mat, total_blength);
                }
                // NHANLT NOTE:
                // the same as above but total_blength = 0 then we simplify the formula to save the runtime (avoid multiplying with 0)
                else
                {
                  total_factor *= dotProduct<num_states>(&((*seq1_region->likelihood)[0]), &((*seq2_region->likelihood)[0]));
                }
            }
            // 3.2. e1.type = O and e2.type = R or A/C/G/T
            else
            {
                StateType seq2_state = seq2_region->type;
                if (seq2_state == TYPE_R)
                    seq2_state = aln.ref_seq[end_pos];
                
                if (total_blength > 0)
                {
                    // NHANLT NOTE:
                    // tot2: likelihood of i evolves to seq2_state
                    // tot2 = (1 + mut[seq2_state,seq2_state] * total_blength) * lh(seq1,seq2_state) + lh(seq1,i) * mut[i,seq2_state] * total_blength
                    RealNumType* transposed_mut_mat_row = model.transposed_mut_mat + model.row_index[seq2_state];
                    RealNumType tot2 = dotProduct<num_states>(&((*seq1_region->likelihood)[0]), transposed_mut_mat_row);
                    total_factor *= seq1_region->getLH(seq2_state) + total_blength * tot2;
                }
                // NHANLT NOTE:
                // the same as above but total_blength = 0 then we simplify the formula to save the runtime (avoid multiplying with 0)
                else
                    total_factor *= seq1_region->getLH(seq2_state);
            }
        }
        // 4. e1.type = A/C/G/T
        else
        {
            int seq1_state = seq1_region->type;
            
            // 4.1. e1.type =  e2.type
            if (seq1_region->type == seq2_region->type)
            {
                if (seq1_region->plength_observation2root >= 0)
                    total_blength += seq1_region->plength_observation2node;
                
                // NHANLT NOTE:
                // the likelihood that seq1_state unchanges
                if (total_blength > 0)
                    lh_cost += model.diagonal_mut_mat[seq1_state] * total_blength;
            }
            // e1.type = A/C/G/T and e2.type = O/A/C/G/T
            else
            {
                // 4.2. e1.type = A/C/G/T and e2.type = O
                if (seq2_region->type == TYPE_O)
                {
                    if (seq1_region->plength_observation2root >= 0)
                    {
                        RealNumType* transposed_mut_mat_row = model.transposed_mut_mat + model.row_index[seq1_state];
                        RealNumType* mutation_mat_row = model.mutation_mat;
                        RealNumType tot = matrixEvolveRoot<num_states>(&((*seq2_region->likelihood)[0]), seq1_state,
                          model.root_freqs, transposed_mut_mat_row, mutation_mat_row, total_blength, seq1_region->plength_observation2node);
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
                        RealNumType tot = dotProduct<num_states>(mutation_mat_row, &((*seq2_region->likelihood)[0]));
                        tot *= total_blength;
                        tot += seq2_region->getLH(seq1_state);
                        total_factor *= tot;
                    }
                }
                // 4.3. e1.type = A/C/G/T and e2.type = R or A/C/G/T
                else
                {
                    StateType seq2_state = seq2_region->type;
                    if (seq2_state == TYPE_R)
                        seq2_state = aln.ref_seq[end_pos];
                    
                    if (seq1_region->plength_observation2root >= 0)
                    {
                        if (total_blength > 0)
                        {
                            // NHANLT NOTE: UNSURE
                            // seq1_state_evolves_seq2_state = (1) the likelihood that seq1_state evolves to seq2_state * (2) the likelihood that seq1_state unchanges from the observing position
                            // (1) = model.mutation_mat[model.row_index[seq1_state] + seq2_state] * total_blength
                            // (2) = (1.0 + model.diagonal_mut_mat[seq1_state] * seq1_region->plength_observation2node)
                            RealNumType seq1_state_evolves_seq2_state = model.mutation_mat[model.row_index[seq1_state] + seq2_state] * total_blength * (1.0 + model.diagonal_mut_mat[seq1_state] * seq1_region->plength_observation2node);
                            
                            // NHANLT NOTE: UNCLEAR
                            // consider the inverse process of the above
                            // seq2_state_evolves_seq1_state = (1) the likelihood that seq2_state evolves to seq1_state * (2) the likelihood that seq2_state unchanges from the observing position
                            // (1) = root_freqs[seq2_state] / root_freqs[seq1_state] * mutation_mat[model.row_index[seq2_state] + seq1_state] * seq1_region->plength_observation2node
                            // (2) = (1.0 + model.diagonal_mut_mat[seq2_state] * total_blength)
                            RealNumType seq2_state_evolves_seq1_state = model.freqi_freqj_qij[model.row_index[seq2_state] + seq1_state] * seq1_region->plength_observation2node * (1.0 + model.diagonal_mut_mat[seq2_state] * total_blength);
                            
                            total_factor *= seq1_state_evolves_seq2_state + seq2_state_evolves_seq1_state;
                        }
                        // NHANLT NOTE:
                        // the same as above but total_blength = 0 then we simplify the formula to save the runtime (avoid multiplying with 0)
                        else
                            total_factor *= model.freqi_freqj_qij[model.row_index[seq2_state] + seq1_state] * seq1_region->plength_observation2node;
                    }
                    // NHANLT NOTE:
                    // add the likelihood that seq1_state evoles to seq2_state = mut[seq1_state,seq2_state] * total_blength
                    else if (total_blength > 0)
                        total_factor *= model.mutation_mat[model.row_index[seq1_state] + seq2_state] * total_blength;
                    else
                        return MIN_NEGATIVE;
                }
            }
        }
         
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

RealNumType Tree::calculateSamplePlacementCost(RealNumType* cumulative_rate, const SeqRegions* const parent_regions, const SeqRegions* const child_regions, RealNumType blength)
{
    return (this->*calculateSamplePlacementCostPointer)(cumulative_rate, parent_regions, child_regions, blength);
}

// this implementation derives from appendProb
template <const StateType num_states>
RealNumType Tree::calculateSamplePlacementCostTemplate(RealNumType* cumulative_rate, const SeqRegions* const parent_regions, const SeqRegions* const child_regions, RealNumType blength)
{     // 10% of total runtime
    // init dummy variables
    RealNumType lh_cost = 0;
    PositionType pos = 0;
    RealNumType total_factor = 1;
    const SeqRegions& seq1_regions = *parent_regions;
    const SeqRegions& seq2_regions = *child_regions;
    size_t iseq1 = 0;
    size_t iseq2 = 0;
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
        if (seq2_region->type == TYPE_N || seq1_region->type == TYPE_N)
        {
            pos = end_pos + 1;
            continue;
        }
        // e1.type != N && e2.type != N
        else {
            // A,C,G,T
            // R -> same as the reference
            // N -> gaps
            // O -> vector of 4 probability to observe A, C, G, T
            // 2. e1.type = R
            if (seq1_region->type == TYPE_R)
            {
                // 2.1. e1.type = R and e2.type = R
                if (seq2_region->type == TYPE_R)
                {
                    if (seq1_region->plength_observation2node < 0 && seq1_region->plength_observation2root < 0)
                        lh_cost += blength * (cumulative_rate[end_pos + 1] - cumulative_rate[pos]);
                    else
                    {
                        total_blength = blength + seq1_region->plength_observation2node;
                        if (seq1_region->plength_observation2root < 0)
                            lh_cost += total_blength * (cumulative_rate[end_pos + 1] - cumulative_rate[pos]);
                        else
                            // here contribution from root frequency gets added and subtracted so it's ignored
                            lh_cost += (total_blength + seq1_region->plength_observation2root) * (cumulative_rate[end_pos + 1] - cumulative_rate[pos]);
                    }
                }
                // 2.2. e1.type = R and e2.type = O
                else if (seq2_region->type == TYPE_O)
                {
                    StateType seq1_state = aln.ref_seq[end_pos];
                    if (seq1_region->plength_observation2root >= 0)
                    {
                        total_blength = seq1_region->plength_observation2root + blength;
                        
                        if (seq2_region->getLH(seq1_state) > 0.1)
                        {
                            total_blength += seq1_region->plength_observation2node;
                            
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
                                RealNumType tot2 = freq_j_transposed_ij_row[i] * seq1_region->plength_observation2node + ((seq1_state == i) ? model.root_freqs[i] : 0);
                                RealNumType tot3 = (seq2_region->getLH(i) > 0.1) ? 1 : 0;
                                
                                for (StateType j = 0; j < num_states; ++j)
                                    if (seq2_region->getLH(j) > 0.1)
                                        tot3 += mutation_mat_row[j];
                                tot3 *= total_blength;
                                
                                tot += tot2 * tot3;
                            }
                            
                            total_factor *= tot * model.inverse_root_freqs[seq1_state];
                        }
                    }
                    else
                    {
                        if (seq2_region->getLH(seq1_state) > 0.1)
                        {
                            if (seq1_region->plength_observation2node >= 0)
                                lh_cost += model.diagonal_mut_mat[seq1_state] * (blength + seq1_region->plength_observation2node);
                            else
                                lh_cost += model.diagonal_mut_mat[seq1_state] * blength;
                        }
                        else
                        {
                            RealNumType tot = 0;
                            RealNumType* mutation_mat_row = model.mutation_mat + model.row_index[seq1_state];
                            for (StateType i = 0; i < num_states; ++i)
                                if (seq2_region->getLH(i) > 0.1)
                                    tot += mutation_mat_row[i];
                            
                            if (seq1_region->plength_observation2node >= 0)
                                total_factor *= tot * (blength + seq1_region->plength_observation2node);
                            else
                                total_factor *= tot * blength;
                        }
                    }
                }
                // 2.3. e1.type = R and e2.type = A/C/G/T
                else
                {
                    StateType seq1_state = aln.ref_seq[end_pos];
                    StateType seq2_state = seq2_region->type;
                    
                    if (seq1_region->plength_observation2root >= 0)
                    {
                        // TODO: can cache model.mutation_mat[model.row_index[seq1_state] * model.diagonal_mut_mat[seq1_state]
                        // TODO: can cache  model.freqi_freqj_qij[model.row_index[seq2_state] + seq1_state] * model.diagonal_mut_mat[seq2_state]
                        RealNumType seq1_state_evolves_seq2_state = model.mutation_mat[model.row_index[seq1_state] + seq2_state] * blength * (1.0 + model.diagonal_mut_mat[seq1_state] * seq1_region->plength_observation2node);
                        
                        RealNumType seq2_state_evolves_seq1_state = model.freqi_freqj_qij[model.row_index[seq2_state] + seq1_state] * seq1_region->plength_observation2node * (1.0 + model.diagonal_mut_mat[seq2_state] * (blength + seq1_region->plength_observation2root));
                                                                                                                                                                                                    
                        total_factor *= seq1_state_evolves_seq2_state + seq2_state_evolves_seq1_state;
                    }
                    else
                    {
                        total_factor *= model.mutation_mat[model.row_index[seq1_state] + seq2_state] * (blength + (seq1_region->plength_observation2node < 0 ? 0 : seq1_region->plength_observation2node));
                    }
                }
            }
            // 3. e1.type = O
            else if (seq1_region->type == TYPE_O)
            {
                RealNumType blength13 = blength;
                if (seq1_region->plength_observation2node >= 0)
                {
                    blength13 = seq1_region->plength_observation2node;
                    if (blength > 0)
                        blength13 += blength;
                }
                    
                // 3.1. e1.type = O and e2.type = O
                if (seq2_region->type == TYPE_O)
                {
                    RealNumType tot = 0;
                    
                    RealNumType* mutation_mat_row = model.mutation_mat;
                                        
                    for (StateType i = 0; i < num_states; ++i, mutation_mat_row += num_states)
                    {
                        RealNumType tot2 = 0;
                        
                        for (StateType j = 0; j < num_states; ++j)
                            tot2 += (seq2_region->getLH(j) > 0.1) ? mutation_mat_row[j] : 0;
                        
                        tot2 *= blength13;
                        
                        if (seq2_region->getLH(i) > 0.1)
                            tot2 += 1;
                        
                        tot += tot2 * seq1_region->getLH(i);
                    }
                        
                    total_factor *= tot;
                }
                // 3.2. e1.type = O and e2.type = R or A/C/G/T
                else
                {
                    StateType seq2_state = seq2_region->type;
                    if (seq2_state == TYPE_R)
                        seq2_state = aln.ref_seq[end_pos];
                    
                    RealNumType tot2 = 0;
                    RealNumType *transposed_mut_mat_row = model.transposed_mut_mat + model.row_index[seq2_state];
                    for (StateType j = 0; j < num_states; ++j)
                        tot2 += transposed_mut_mat_row[j] * seq1_region->getLH(j);
                    
                    total_factor *= seq1_region->getLH(seq2_state) + blength13 * tot2;
                }
            }
            // 4. e1.type = A/C/G/T
            else
            {
                int seq1_state = seq1_region->type;
                
                // 4.1. e1.type =  e2.type
                if (seq1_region->type == seq2_region->type)
                {
                    total_blength = blength;
                    total_blength += (seq1_region->plength_observation2node < 0 ? 0 : seq1_region->plength_observation2node);
                    total_blength += (seq1_region->plength_observation2root < 0 ? 0 : seq1_region->plength_observation2root);

                    lh_cost += model.diagonal_mut_mat[seq1_state] * total_blength;
                }
                // e1.type = A/C/G/T and e2.type = O/A/C/G/T
                else
                {
                    // 4.2. e1.type = A/C/G/T and e2.type = O
                    if (seq2_region->type == TYPE_O)
                    {
                        RealNumType tot = 0.0;
                        
                        if (seq1_region->plength_observation2root >= 0)
                        {
                            RealNumType blength15 = blength + seq1_region->plength_observation2root;
                            
                            if (seq2_region->getLH(seq1_state) > 0.1)
                                lh_cost += model.diagonal_mut_mat[seq1_state] * (blength15 + seq1_region->plength_observation2node);
                            else
                            {
                                RealNumType* freq_j_transposed_ij_row = model.freq_j_transposed_ij + model.row_index[seq1_state];
                                RealNumType* mutation_mat_row = model.mutation_mat;
                                                    
                                for (StateType i = 0; i < num_states; ++i, mutation_mat_row += num_states)
                                {
                                    RealNumType tot2 = freq_j_transposed_ij_row[i] * seq1_region->plength_observation2node + ((seq1_state == i) ? model.root_freqs[i] : 0);
                                        
                                    RealNumType tot3 = 0;
                                    for (StateType j = 0; j < num_states; ++j)
                                        if (seq2_region->getLH(j) > 0.1)
                                            tot3 += mutation_mat_row[j];
                                    
                                    if (seq2_region->getLH(i) > 0.1)
                                        tot += tot2 * (1.0 + blength15 * tot3);
                                    else
                                        tot += tot2 * blength15 * tot3;
                                }
                                
                                total_factor *= (tot * model.inverse_root_freqs[seq1_state]);
                            }
                        }
                        else
                        {
                            RealNumType tmp_blength = blength + (seq1_region->plength_observation2node < 0 ? 0 : seq1_region->plength_observation2node);
                            if (seq2_region->getLH(seq1_state) > 0.1)
                                lh_cost += model.diagonal_mut_mat[seq1_state] * tmp_blength;
                            else
                            {
                                RealNumType* mutation_mat_row = model.mutation_mat + model.row_index[seq1_state];
                                for (StateType j = 0; j < num_states; ++j)
                                    if (seq2_region->getLH(j) > 0.1)
                                        tot += mutation_mat_row[j];
                                
                                total_factor *= tot * tmp_blength;
                            }
                        }
                    }
                    // 4.3. e1.type = A/C/G/T and e2.type = R or A/C/G/T
                    else
                    {
                        StateType seq2_state = seq2_region->type;
                        if (seq2_state == TYPE_R)
                            seq2_state = aln.ref_seq[end_pos];
                        
                        if (seq1_region->plength_observation2root >= 0)
                        {
                            // here we ignore contribution of non-parsimonious mutational histories
                            RealNumType seq1_state_evoloves_seq2_state = model.mutation_mat[model.row_index[seq1_state] + seq2_state] * (blength + seq1_region->plength_observation2root) * (1.0 + model.diagonal_mut_mat[seq1_state] * seq1_region->plength_observation2node);
                            
                            RealNumType seq2_state_evolves_seq1_state = model.freqi_freqj_qij[model.row_index[seq2_state] + seq1_state] * seq1_region->plength_observation2node * (1.0 + model.diagonal_mut_mat[seq2_state] * (blength + seq1_region->plength_observation2root));
                            
                            total_factor *= (seq1_state_evoloves_seq2_state + seq2_state_evolves_seq1_state);
                        }
                        else
                        {
                            RealNumType tmp_blength = ((seq1_region->plength_observation2node < 0) ? blength : blength + seq1_region->plength_observation2node);
                            
                            total_factor *= model.mutation_mat[model.row_index[seq1_state] + seq2_state] * tmp_blength;
                        }
                    }
                }
            }
        }
         
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

void Tree::updateZeroBlength(Node* node, stack<Node*> &node_stack, RealNumType threshold_prob, RealNumType* cumulative_rate, RealNumType default_blength, RealNumType min_blength, RealNumType max_blength)
{
    // get the top node in the phylo-node
    Node* top_node = node->getTopNode();
    ASSERT(top_node);
    SeqRegions* upper_left_right_regions = top_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
    SeqRegions* lower_regions = top_node->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
    
    RealNumType best_lh = calculateSamplePlacementCost(cumulative_rate, upper_left_right_regions, lower_regions, default_blength);
    RealNumType best_length = default_blength;
    
    while (best_length > min_blength)
    {
        RealNumType new_blength = best_length * 0.5;
        RealNumType new_lh = calculateSamplePlacementCost(cumulative_rate, upper_left_right_regions, lower_regions, new_blength);
        
        if (new_lh > best_lh)
        {
            best_lh = new_lh;
            best_length = new_blength;
        }
        else
            break;
    }
    
    if (best_length > 0.7 * default_blength)
    {
        while (best_length < max_blength)
        {
            RealNumType new_blength = best_length + best_length;
            RealNumType new_lh = calculateSamplePlacementCost(cumulative_rate, upper_left_right_regions, lower_regions, new_blength);
            if (new_lh > best_lh)
            {
                best_lh = new_lh;
                best_length = new_blength;
            }
            else
                break;
        }
    }
    
    // update best_length
    top_node->length = best_length;
    top_node->neighbor->length = best_length;
    
    // add current node and its parent to node_stack to for updating partials further from these nodes
    top_node->outdated = true;
    top_node->neighbor->getTopNode()->outdated = true;
    node_stack.push(top_node);
    node_stack.push(top_node->neighbor);
}
