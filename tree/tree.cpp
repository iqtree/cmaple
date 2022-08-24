#include "tree.h"

Tree::Tree()
{
    params = NULL;
    aln = new Alignment();
    model = new Model();
    root = NULL;
}

Tree::Tree(Params *n_params, Node* n_root)
{
    params = n_params;
    aln = new Alignment();
    model = new Model();
    root = n_root;
}

Tree::~Tree()
{
    if (aln)
    {
        delete aln;
        aln = NULL;
    }
    
    if (model)
    {
        delete model;
        model = NULL;
    }
    
    // browse tree to delete all nodes
    // NHANLT: TODO
}

void Tree::setupFunctionPointers()
{
    switch (aln->num_states) {
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

string Tree::exportTreeString(Node* node)
{
    // init starting node from root
    if (!node)
        node = root;
    
    // move to its neighbor
    if (node->neighbor)
        node = node->neighbor;
    
    // do something with its neighbor
    if (node->isLeave())
        return node->exportString();
        
    string output = "(";
    bool add_comma = false;
    Node* next;
    FOR_NEXT(node, next)
    {
        if (!add_comma)
            add_comma = true;
        else
            output += ",";
        output += exportTreeString(next);
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
    PositionType seq_length = aln->ref_seq.size();
    
    while (!node_stack.empty())
    {
        Node* node = node_stack.top();
        node_stack.pop();
        
        // NHANLT: debug
       /* if (node->next && node->next->neighbor && node->next->neighbor->seq_name == "711")
            cout << "dsdas";*/
        
        bool update_blength = false;
        node->outdated = true;
        
        SeqRegions* parent_upper_regions = NULL;
        if (root != node->getTopNode())
            parent_upper_regions = node->getTopNode()->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
            
        // change in likelihoods is coming from parent node
        if (node->is_top)
        {
            // if necessary, update the total probabilities at the mid node.
            if (node->length > 0)
            {
                // update vector of regions at mid-branch point
                SeqRegions* mid_branch_regions = NULL;
                SeqRegions* lower_regions = node->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
                parent_upper_regions->mergeUpperLower(mid_branch_regions, node->length * 0.5, lower_regions, node->length * 0.5, aln, model, params->threshold_prob);
                
                if (!mid_branch_regions)
                {
                    if (node->length > 1e-100)
                        outError("inside updatePartialLh(), from parent: should not have happened since node->length > 0");
                    updateZeroBlength(node, node_stack, aln, model, params->threshold_prob, cumulative_rate, default_blength, min_blength, max_blength);
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
                parent_upper_regions->mergeUpperLower(upper_left_right_regions_1, node->length, lower_regions_1, next_node_1->length, aln, model, params->threshold_prob);
                
                if (!upper_left_right_regions_1 || upper_left_right_regions_1->size() == 0)
                {
                    if (node->length <= 0 && next_node_1->length <= 0)
                        updateZeroBlength(node, node_stack, aln, model, params->threshold_prob, cumulative_rate, default_blength, min_blength, max_blength);
                    else
                        outError("Strange: None vector from non-zero distances in updatePartialLh() from parent direction.");
                    
                    update_blength = true;
                }
                
                if (!update_blength)
                {
                    SeqRegions* lower_regions_2 = next_node_2->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
                    parent_upper_regions->mergeUpperLower(upper_left_right_regions_2, node->length, lower_regions_2, next_node_2->length, aln, model, params->threshold_prob);
                    
                    if (!upper_left_right_regions_2 || upper_left_right_regions_2->size() == 0)
                    {
                        if (node->length <= 0 && next_node_2->length <= 0)
                        {
                            updateZeroBlength(node, node_stack, aln, model, params->threshold_prob, cumulative_rate, default_blength, min_blength, max_blength);
                            update_blength = true;
                        }
                        else
                            outError("Strange: None vector from non-zero distances in updatePartialLh() from parent direction, child0.");
                    }
                }
                
                if (!update_blength)
                {
                    if (next_node_1->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate)->areDiffFrom(upper_left_right_regions_2, seq_length, num_states, params))
                    {
                        delete next_node_1->partial_lh;
                        next_node_1->partial_lh = upper_left_right_regions_2;
                        upper_left_right_regions_2 = NULL;
                        node_stack.push(next_node_1->neighbor);
                    }
                    
                    if (next_node_2->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate)->areDiffFrom(upper_left_right_regions_1, seq_length, num_states, params))
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
                    updateZeroBlength(node->neighbor, node_stack, aln, model, params->threshold_prob, cumulative_rate, default_blength, min_blength, max_blength);
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
                        updateZeroBlength(top_node, node_stack, aln, model, params->threshold_prob, cumulative_rate, default_blength, min_blength, max_blength);
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
                if (top_node->length > 0 && top_node != root)
                {
                    SeqRegions* new_mid_regions = NULL;
                    SeqRegions* tmp_lower_regions = top_node->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
                    parent_upper_regions->mergeUpperLower(new_mid_regions, top_node->length * 0.5, tmp_lower_regions, top_node->length * 0.5, aln, model, params->threshold_prob);
                    
                    if (!new_mid_regions)
                    {
                        updateZeroBlength(top_node, node_stack, aln, model, params->threshold_prob, cumulative_rate, default_blength, min_blength, max_blength);
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
                if (top_node->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate)->areDiffFrom(old_lower_regions, seq_length, num_states, params))
                {
                    if (root != top_node)
                        node_stack.push(top_node->neighbor);
                }

                // update likelihoods at sibling node
                SeqRegions* new_upper_regions = NULL;
                if (root != top_node)
                    parent_upper_regions->mergeUpperLower(new_upper_regions, top_node->length, this_node_lower_regions, this_node_distance, aln, model, params->threshold_prob);
                else
                    new_upper_regions = node->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate)->computeTotalLhAtRoot(num_states, model, this_node_distance);
                
                if (!new_upper_regions || new_upper_regions->size() == 0)
                {
                    if (top_node->length <= 0 && this_node_distance <= 0)
                    {
                        updateZeroBlength(top_node, node_stack, aln, model, params->threshold_prob, cumulative_rate, default_blength, min_blength, max_blength);
                        update_blength = true;
                    }
                    else
                        outError("Strange: None vector from non-zero distances in updatePartialLh() from child direction, new_upper_regions.");
                }
                else
                {
                    if (next_node_upper_left_right_regions->areDiffFrom(new_upper_regions, seq_length, num_states, params))
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

void Tree::seekSamplePlacement(Node* start_node, string seq_name, SeqRegions* sample_regions, Node* &selected_node, RealNumType &best_lh_diff , bool &is_mid_branch, RealNumType &best_up_lh_diff, RealNumType &best_down_lh_diff, Node* &best_child, RealNumType* cumulative_rate, RealNumType default_blength, RealNumType min_blength_mid)
{
    // init variables
    // output variables
    selected_node = start_node;
    best_lh_diff = -DBL_MAX;
    is_mid_branch = false;
    best_up_lh_diff = -DBL_MAX;
    best_down_lh_diff = -DBL_MAX;
    best_child = NULL;
    // dummy variables
    RealNumType lh_diff_mid_branch = 0;
    RealNumType lh_diff_at_node = 0;
    // stack of nodes to examine positions
    stack<ExtendedNode> extended_node_stack;
    extended_node_stack.push(ExtendedNode(start_node, 0, -DBL_MAX));
    
    // recursively examine positions for placing the new sample
    while (!extended_node_stack.empty())
    {
        ExtendedNode current_extended_node = extended_node_stack.top();
        Node* current_node = current_extended_node.node;
        extended_node_stack.pop();
        
        // NHANLT: debug
       /* if (current_node->next && ((current_node->next->neighbor && current_node->next->neighbor->seq_name == "25")
                                   || (current_node->next->next->neighbor && current_node->next->next->neighbor->seq_name == "25")))
            cout << "fdsfsd";*/
        /*if (current_node->seq_name == "639")
            cout << "fdsfsd";
        if (current_node->seq_name == "635")
            cout << "fdsfsd";*/
    
        // if the current node is a leaf AND the new sample/sequence is strictly less informative than the current node
        // -> add the new sequence into the list of minor sequences of the current node + stop seeking the placement
        if ((!current_node->next) && (current_node->partial_lh->compareWithSample(sample_regions, aln->ref_seq.size(), aln->num_states) == 1))
        {
            current_node->less_info_seqs.push_back(seq_name);
            selected_node = NULL;
            return;
        }
        
        // 1. try first placing as a descendant of the mid-branch point of the branch above the current node
        if (current_node != root && current_node->length > 0)
        {
            // compute the placement cost
            lh_diff_mid_branch = calculateSamplePlacementCost(aln, model, cumulative_rate, current_node->mid_branch_lh, sample_regions, default_blength);
            
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
            lh_diff_mid_branch = -DBL_MAX;

        // 2. try to place as descendant of the current node (this is skipped if the node has top branch length 0 and so is part of a polytomy).
        if (current_node == root || current_node->length > 0)
        {
            // compute the placement cost
            lh_diff_at_node = calculateSamplePlacementCost(aln, model, cumulative_rate, current_node->total_lh, sample_regions, default_blength);
            
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
        if ((params->strict_stop_seeking_placement
             && current_extended_node.failure_count < params->failure_limit
             && lh_diff_at_node > (best_lh_diff - params->thresh_log_lh_subtree_explore))
            || (!params->strict_stop_seeking_placement
                && (current_extended_node.failure_count < params->failure_limit
                    || lh_diff_at_node > (best_lh_diff - params->thresh_log_lh_subtree_explore))))
        {
            Node* neighbor_node;
            FOR_NEIGHBOR(current_node, neighbor_node)
            extended_node_stack.push(ExtendedNode(neighbor_node, current_extended_node.failure_count, lh_diff_at_node));
        }
    }

    // exploration of the tree is finished, and we are left with the node found so far with the best appending likelihood cost. Now we explore placement just below this node for more fine-grained placement within its descendant branches.
    best_down_lh_diff = -DBL_MAX;
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
                // NHANLT: only try at the mid-branch position to save time
                /*// now try to place on the current branch below the best node, at an height above the mid-branch.
                RealNumType new_blength = node->length * 0.5;
                RealNumType new_best_lh_mid_branch = -DBL_MAX;

                // try to place new sample along the upper half of the current branch
                while (true)
                {
                    // compute the placement cost
                    RealNumType new_lh_mid_branch = node->mid_branch_lh->calculatePlacementCost(aln, model, cumulative_rate, sample_regions, default_blength);
                    
                    // record new_best_lh_mid_branch
                    if (new_lh_mid_branch > new_best_lh_mid_branch)
                        new_best_lh_mid_branch = new_lh_mid_branch;
                    // otherwise, stop trying along the current branch
                    else
                        break;
                    
                    // try at different position along the current branch
                    new_blength *= 0.5;
                    
                    // stop trying if reaching the minimum branch length
                    if (new_blength <= min_blength_mid * 0.5)
                        break;
                }*/
                
                RealNumType new_best_lh_mid_branch = calculateSamplePlacementCost(aln, model, cumulative_rate, node->mid_branch_lh, sample_regions, default_blength);
                
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

void Tree::seekPlacement(Node* &best_node, RealNumType &best_lh_diff, bool &is_mid_branch, RealNumType &best_up_lh_diff, RealNumType &best_down_lh_diff, Node* &best_child, Node* child_node, RealNumType &removed_blength, RealNumType* cumulative_rate, RealNumType default_blength, RealNumType min_blength_mid, bool search_subtree_placement, SeqRegions* sample_regions)
{
    // init variables
    Node* node = child_node->neighbor->getTopNode();
    best_node = node;
    best_lh_diff = -DBL_MAX;
    is_mid_branch = false;
    best_up_lh_diff = -DBL_MAX;
    best_down_lh_diff = -DBL_MAX;
    best_child = NULL;
    SeqRegions* subtree_regions = NULL;
    // stack of nodes to examine positions
    stack<ExtendedNode> node_stack;
    // dummy variables
    RealNumType threshold_prob = params->threshold_prob;
    RealNumType lh_diff_mid_branch = 0;
    RealNumType lh_diff_at_node = 0;
    SeqRegions* parent_upper_lr_regions = NULL;

    /*if (search_subtree_placement)
    {
        subtree_regions = child_node->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
        
        if (node != root)
        {
            parent_upper_lr_regions = node->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
            #the list node_stack keeps trak of the nodes of the tree we wil need to traverse, (first element of each entry),
            # the direction from where we are visitng them (0=from parent, and 1,2=from one of the children),
            # an updated genome list from the direction where we come from (taking into account the removal of the given subtree),
            # a branch length value separating the node from this updated genome list (useful for the fact that the removal of the subtree changes the branch length at the removal node),
            # a flag that says if the updated genome list passed needs still updating, or if it has become identical to the pre-existing genome list in the tree (which usually happens after a while),
            # the likelihood cost of appending at the last node encountered in this direction,
            # a number of consecutively failed traversal steps since the last improvement found (if this number goes beyond a threshold, traversal in the considered direction might be halted).
            node_stack.append((node.up,childUp,node.children[1-child].probVect,node.children[1-child].dist+node.dist,True,bestLKdiff,0))
            node_stack.append((node.children[1-child],0,parent_upper_lr_regions,node.children[1-child].dist+node.dist,True,bestLKdiff,0))
        }
        else
        {
            // node is root
            if node.children[1-child].children: # case there is only one sample outside of the subtree doesn't need to be considered
            {
                child1=node.children[1-child].children[0]
                child2=node.children[1-child].children[1]
                vectUp1=rootVector(child2.probVect,child2.dist,mutMatrix)
                node_stack.append((child1,0,vectUp1,child1.dist,True,bestLKdiff,0))
                vectUp2=rootVector(child1.probVect,child1.dist,mutMatrix)
                node_stack.append((child2,0,vectUp2,child2.dist,True,bestLKdiff,0))
            }
        }
    }
    else:
        subtree_regions=newSamplePartials
        if is_mid_branch:
            downLK=best_down_lh_diff
        else:
            downLK=bestLKdiff
        if node.up!=None:
            if node.up.children[0]==node:
                childUp=1
                parent_upper_lr_regions=node.up.probVectUpRight
            else:
                childUp=2
                parent_upper_lr_regions=node.up.probVectUpLeft
            node_stack.append((node.up,childUp,node.probVect,node.dist,False,downLK,0))
        if node.children:
            node_stack.append((node.children[0],0,node.probVectUpRight,node.children[0].dist,False,downLK,0))
            node_stack.append((node.children[1],0,node.probVectUpLeft,node.children[1].dist,False,downLK,0))

    while node_stack:
        t1,direction,passedPartials,distance,needsUpdating,lastLK,failedPasses=node_stack.pop()
        if direction==0:
            #consider the case we are moving from a parent to a child
            if t1.dist:
                if not (t1.up==node or t1.up==None): #try to append mid-branch
                    if needsUpdating:
                        midTot=mergeVectorsUpDown(passedPartials,distance * 0.5,t1.probVect,distance * 0.5, mutMatrix)
                    else:
                        midTot=t1.probVectTotUp
                    if midTot==None:
                        continue
                    if search_subtree_placement:
                        midProb=appendProbNode(midTot,subtree_regions,removedBLen,mutMatrix)
                    else:
                        midProb=appendProb(midTot,subtree_regions,removedBLen,mutMatrix)
                    if midProb>bestLKdiff:
                        best_node=t1
                        bestLKdiff=midProb
                        is_mid_branch=True
                        failedPasses=0
                else:
                    midProb=float("-inf")
                #now try appending exactly at node
                if needsUpdating:
                    nodeTot=mergeVectorsUpDown(passedPartials,distance,t1.probVect,False,mutMatrix)
                    if not areVectorsDifferent(nodeTot,t1.probVectTot):
                        needsUpdating=False
                else:
                    nodeTot=t1.probVectTot

                if nodeTot==None:
                    continue
                if search_subtree_placement:
                    nodeProb=appendProbNode(nodeTot,subtree_regions,removedBLen,mutMatrix)
                else:
                    nodeProb=appendProb(nodeTot,subtree_regions,removedBLen,mutMatrix)
                if nodeProb>bestLKdiff:
                    best_node=t1
                    bestLKdiff=nodeProb
                    is_mid_branch=False
                    failedPasses=0
                    best_up_lh_diff=midProb
                elif midProb>=(bestLKdiff-thresholdProb):
                    best_up_lh_diff=lastLK
                    best_down_lh_diff=nodeProb
                elif nodeProb<(lastLK-thresholdLogLKconsecutivePlacement): #placement at current node is considered failed if placement likelihood is not improved by a certain margin compared to best placement so far for the nodes above it.
                    failedPasses+=1
            else:
                nodeProb=lastLK
            
            #keep crawling down into children nodes unless the stop criteria for the traversal are satisfied.
            traverseChildren=False
            if search_subtree_placement:
                if strictTopologyStopRules:
                    if failedPasses<=allowedFailsTopology and nodeProb>(bestLKdiff-thresholdLogLKtopology) and t1.children:
                        traverseChildren=True
                else:
                    if failedPasses<=allowedFailsTopology or nodeProb>(bestLKdiff-thresholdLogLKtopology):
                        if t1.children:
                            traverseChildren=True
            else:
                if strictInitialStopRules:
                    if failedPasses<=allowedFails and nodeProb>(bestLKdiff-thresholdLogLK):
                        if t1.children:
                            traverseChildren=True
                else:
                    if failedPasses<=allowedFails or nodeProb>(bestLKdiff-thresholdLogLK):
                        if t1.children:
                            traverseChildren=True
            if traverseChildren:
                #if (failedPasses<=allowedFailsTopology or nodeProb>(bestLKdiff-thresholdLogLKtopology) ) and len(t1.children)==2:
                child=t1.children[0]
                otherChild=t1.children[1]
                if needsUpdating:
                    vectUpRight=mergeVectorsUpDown(passedPartials,distance,otherChild.probVect,otherChild.dist,mutMatrix)
                else:
                    vectUpRight=t1.probVectUpRight
                if vectUpRight!=None:
                    node_stack.append((child,0,vectUpRight,child.dist,needsUpdating,nodeProb,failedPasses))
                child=t1.children[1]
                otherChild=t1.children[0]
                if needsUpdating:
                    vectUpLeft=mergeVectorsUpDown(passedPartials,distance,otherChild.probVect,otherChild.dist,mutMatrix)
                else:
                    vectUpLeft=t1.probVectUpLeft
                if vectUpLeft!=None:
                    node_stack.append((child,0,vectUpLeft,child.dist,needsUpdating,nodeProb,failedPasses))

        else: #case when crawling up from child to parent
            if t1.dist or t1.up==None: #append directly at the node
                if needsUpdating:
                    if direction==1:
                        nodeTot=mergeVectorsUpDown(t1.probVectUpRight,False,passedPartials,distance,mutMatrix)
                    else:
                        nodeTot=mergeVectorsUpDown(t1.probVectUpLeft,False,passedPartials,distance,mutMatrix)
                    if nodeTot==None:
                        #print("Removing a node created an inconsistency while moving up.")
                        continue
                    elif not areVectorsDifferent(nodeTot,t1.probVectTot):
                        needsUpdating=False
                else:
                    nodeTot=t1.probVectTot
                if search_subtree_placement:
                    nodeProb=appendProbNode(nodeTot,subtree_regions,removedBLen,mutMatrix)
                else:
                    nodeProb=appendProb(nodeTot,subtree_regions,removedBLen,mutMatrix)
                if nodeProb<(lastLK-thresholdLogLKconsecutivePlacement): #placement at current node is considered failed if placement likelihood is not improved by a certain margin compared to best placement so far for the nodes above it.
                    failedPasses+=1
                elif nodeProb>bestLKdiff:
                    best_node=t1
                    bestLKdiff=nodeProb
                    is_mid_branch=False
                    failedPasses=0
            else:
                nodeProb=lastLK

            otherChild=t1.children[2-direction]
            midBottom=None
            if t1.dist and t1.up!=None: #try appending mid-branch
                if needsUpdating:
                    midBottom=mergeVectors(otherChild.probVect,otherChild.dist,passedPartials,distance,mutMatrix)
                    if midBottom==None:
                        continue
                    if t1==t1.up.children[0]:
                        vectUp=t1.up.probVectUpRight
                    else:
                        vectUp=t1.up.probVectUpLeft
                    midTot=mergeVectorsUpDown(vectUp,t1.dist * 0.5,midBottom,t1.dist * 0.5,mutMatrix)
                else:
                    midTot=t1.probVectTotUp
                if midTot==None:
                    continue
                if search_subtree_placement:
                    midProb=appendProbNode(midTot,subtree_regions,removedBLen,mutMatrix)
                else:
                    midProb=appendProb(midTot,subtree_regions,removedBLen,mutMatrix)
                if best_node==t1:
                    best_up_lh_diff=midProb
                if midProb>bestLKdiff:
                    best_node=t1
                    bestLKdiff=midProb
                    is_mid_branch=True
                    failedPasses=0
                    best_down_lh_diff=nodeProb
                elif nodeProb>=(bestLKdiff-thresholdProb):
                    best_up_lh_diff=midProb
            else:
                midProb=float("-inf")
            
            #testing for the stop rule of the traversal process
            keepTraversing=False
            if search_subtree_placement:
                if strictTopologyStopRules:
                    if failedPasses<=allowedFailsTopology and nodeProb>(bestLKdiff-thresholdLogLKtopology):
                        keepTraversing=True
                else:
                    if failedPasses<=allowedFailsTopology or nodeProb>(bestLKdiff-thresholdLogLKtopology):
                        keepTraversing=True
            else:
                if strictInitialStopRules:
                    if failedPasses<=allowedFails and nodeProb>(bestLKdiff-thresholdLogLK):
                        keepTraversing=True
                else:
                    if failedPasses<=allowedFails or nodeProb>(bestLKdiff-thresholdLogLK):
                        keepTraversing=True
            if keepTraversing:
                #if failedPasses<=allowedFailsTopology or nodeProb>(bestLKdiff-thresholdLogLKtopology) :
                # keep crawling up into parent and sibling node
                if t1.up!=None: #case the node is not the root
                    #first pass the crawling down the other child
                    if t1==t1.up.children[0]:
                        upChild=0
                        if needsUpdating:
                            parent_upper_lr_regions=t1.up.probVectUpRight
                    else:
                        upChild=1
                        if needsUpdating:
                            parent_upper_lr_regions=t1.up.probVectUpLeft
                    if needsUpdating:
                        vectUp=mergeVectorsUpDown(parent_upper_lr_regions,t1.dist,passedPartials,distance,mutMatrix)
                    else:
                        if direction==1:
                            vectUp=t1.probVectUpLeft
                        else:
                            vectUp=t1.probVectUpRight

                    if vectUp==None:
                        continue
                    else:
                        node_stack.append((otherChild,0,vectUp,otherChild.dist,needsUpdating,nodeProb,failedPasses))
                    #now pass the crawling up to the parent node
                    if needsUpdating:
                        if midBottom==None:
                            midBottom=mergeVectors(otherChild.probVect,otherChild.dist,passedPartials,distance,mutMatrix)
                            if midBottom==None:
                                continue
                    else:
                        midBottom=t1.probVect
                    node_stack.append((t1.up,upChild+1,midBottom,t1.dist,needsUpdating,nodeProb,failedPasses))
                #now consider case of root node
                else:
                    if needsUpdating:
                        vectUp=rootVector(passedPartials,distance,mutMatrix)
                    else:
                        if direction==1:
                            vectUp=t1.probVectUpLeft
                        else:
                            vectUp=t1.probVectUpRight
                    node_stack.append((otherChild,0,vectUp,otherChild.dist,needsUpdating,nodeProb,failedPasses))
    if search_subtree_placement:
        return best_node , bestLKdiff , is_mid_branch
    else:
        #exploration of the tree is finished, and we are left with the node found so far with the best appending likelihood cost.
        #Now we explore placement just below this node for more fine-grained placement within its descendant branches.
        bestOfChildLK=float("-inf")
        bestChild=None
        if is_mid_branch:
            nodeUp=best_node.up
            while (not nodeUp.dist) and (nodeUp.up!=None):
                nodeUp=nodeUp.up
            best_up_lh_diff=appendProb(nodeUp.probVectTot,subtree_regions,removedBLen,mutMatrix)
            return best_node , bestLKdiff , is_mid_branch, best_up_lh_diff, best_down_lh_diff, best_node
        else:
            #current node might be part of a polytomy (represented by 0 branch lengths) so we want to explore all the children of the current node
            #to find out if the best placement is actually in any of the branches below the current node.
            node_stack=[]
            for c in best_node.children:
                node_stack.append(c)
            while node_stack:
                t1=node_stack.pop()
                if not t1.dist:
                    for c in t1.children:
                        node_stack.append(c)
                else:
                    #now try to place on the current branch below the best node, at an height above the mid-branch.
                    newBLen2=t1.dist * 0.5
                    bestLKdiff2=float("-inf")
                    furtherNode=-1
                    newProbVect2=t1.probVectTotUp
                    while True:
                        newLKdiff2=appendProb(newProbVect2,subtree_regions,removedBLen,mutMatrix)
                        if newLKdiff2>bestLKdiff2:
                            bestLKdiff2=newLKdiff2
                        else:
                            break
                        newBLen2=newBLen2 * 0.5
                        if newBLen2<=(minBLenForMidNode * 0.5):
                            break
                        furtherNode+=1
                        newProbVect2=t1.furtherMidNodes[furtherNode]
                    
                    if bestLKdiff2>bestOfChildLK:
                        bestOfChildLK=bestLKdiff2
                        bestChild=t1
            #pass on the best child found to the next function, which will place the new sample somewhere between the best node and the best child.
            return best_node , bestLKdiff , is_mid_branch, best_up_lh_diff, bestOfChildLK, bestChild*/
}


void Tree::placeNewSample(Node* selected_node, SeqRegions* sample, string seq_name, RealNumType best_lh_diff , bool is_mid_branch, RealNumType best_up_lh_diff, RealNumType best_down_lh_diff, Node* best_child, RealNumType* cumulative_rate, vector< vector<PositionType> > &cumulative_base, RealNumType default_blength, RealNumType max_blength, RealNumType min_blength)
{
    // dummy variables
    RealNumType best_child_lh;
    RealNumType best_child_split = 0;
    RealNumType best_parent_lh;
    RealNumType best_parent_split = 0;
    SeqRegions* best_parent_regions = NULL;
    SeqRegions* best_child_regions = NULL;
    RealNumType best_root_blength = -1;
    StateType num_states = aln->num_states;
    
    // try to place the new sample as a descendant of a mid-branch point
    if (is_mid_branch)
    {
        SeqRegions* upper_left_right_regions = selected_node->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
        RealNumType best_split = 0.5;
        RealNumType best_split_lh = best_lh_diff;
        RealNumType new_split = 0.25;
        selected_node->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate)->mergeUpperLower(best_child_regions, selected_node->length * 0.5, selected_node->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate), selected_node->length * 0.5, aln, model, params->threshold_prob);
        SeqRegions* new_parent_regions = NULL;

        // try different positions on the existing branch
        while (new_split * selected_node->length > min_blength)
        {
            upper_left_right_regions->mergeUpperLower(new_parent_regions, selected_node->length * new_split, selected_node->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate),  selected_node->length * (1 - new_split), aln, model, params->threshold_prob);
            RealNumType placement_cost = calculateSamplePlacementCost(aln, model, cumulative_rate, new_parent_regions, sample, default_blength);
            
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
        }
        
        if (best_split > 0.49)
        {
            new_split = 0.25;
            while (new_split * selected_node->length > min_blength)
            {
                upper_left_right_regions->mergeUpperLower(new_parent_regions, selected_node->length * (1.0 - new_split), selected_node->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate),selected_node->length * new_split, aln, model, params->threshold_prob);
                
                RealNumType placement_cost = calculateSamplePlacementCost(aln, model, cumulative_rate, new_parent_regions, sample, default_blength);
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
            RealNumType placement_cost = calculateSamplePlacementCost(aln, model, cumulative_rate, best_child_regions, sample, new_blength);
            
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
                RealNumType placement_cost = calculateSamplePlacementCost(aln, model, cumulative_rate, best_child_regions, sample, new_blength);
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
            RealNumType zero_branch_lh = calculateSamplePlacementCost(aln, model, cumulative_rate, best_child_regions, sample, -1);
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
        RealNumType down_distance = selected_node->length * (1 - best_split);
        selected_node->length = down_distance;
        selected_node->neighbor->length = down_distance;
        
        new_sample_node->neighbor = next_node_1;
        next_node_1->neighbor = new_sample_node;
        new_sample_node->length = best_blength;
        new_sample_node->neighbor->length = best_blength;
        
        new_sample_node->partial_lh = sample;
        next_node_1->partial_lh = best_child_regions;
        best_child_regions = NULL;
        upper_left_right_regions->mergeUpperLower(next_node_2->partial_lh, new_internal_node->length, sample, best_blength, aln, model, params->threshold_prob);
        selected_node->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate)->mergeTwoLowers(new_internal_node->partial_lh, selected_node->length, sample, best_blength, aln, model, params->threshold_prob, cumulative_rate);
        upper_left_right_regions->mergeUpperLower(new_internal_node->mid_branch_lh, new_internal_node->length * 0.5, new_internal_node->partial_lh, new_internal_node->length * 0.5, aln, model, params->threshold_prob);
        new_internal_node->computeTotalLhAtNode(aln, model, params->threshold_prob, cumulative_rate, new_internal_node == root);
        
        if (!new_internal_node->total_lh || new_internal_node->total_lh->size() == 0)
            outError("Problem, None vector when placing sample, below node");
        
        /*if distTop>=2*min_blengthForMidNode:
         createFurtherMidNodes(newInternalNode,upper_left_right_regions)*/
        if (best_blength > 0)
        {
            new_sample_node->computeTotalLhAtNode(aln, model, params->threshold_prob, cumulative_rate, new_sample_node == root);
            next_node_1->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate)->mergeUpperLower(new_sample_node->mid_branch_lh, new_sample_node->length * 0.5, sample, new_sample_node->length * 0.5, aln, model, params->threshold_prob);
            /*if best_blength>=2*min_blengthForMidNode:
                createFurtherMidNodes(newInternalNode.children[1],best_child_regions)*/
        }
        
        // update pseudo_count
        model->updatePesudoCount(aln, next_node_1->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate), sample);

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
            SeqRegions* upper_left_right_regions = best_child->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
            SeqRegions* lower_regions = best_child->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
            best_child_regions = NULL;
            upper_left_right_regions->mergeUpperLower(best_child_regions, best_child->length * 0.5, lower_regions, best_child->length * 0.5, aln, model, params->threshold_prob);
            RealNumType new_split = 0.25;
            
            while (new_split * best_child->length > min_blength)
            {
                SeqRegions* new_parent_regions = NULL;
                upper_left_right_regions->mergeUpperLower(new_parent_regions, best_child->length * new_split, lower_regions, best_child->length * (1 - new_split), aln, model, params->threshold_prob);
                
                RealNumType placement_cost = calculateSamplePlacementCost(aln, model, cumulative_rate, new_parent_regions, sample, default_blength);
                
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
                
                // delete new_parent_regions
                if (new_parent_regions) delete new_parent_regions;
            }
            
            best_child_lh = best_split_lh;
            best_child_split = best_split;
        }
        else
            best_child_lh = -DBL_MAX;
        
        // if node is root, try to place as sibling of the current root.
        RealNumType old_root_lh = -DBL_MAX;
        if (root == selected_node)
        {
            old_root_lh = selected_node->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate)->computeAbsoluteLhAtRoot(aln, model, cumulative_base);
            RealNumType new_root_lh;
            SeqRegions* merged_root_sample_regions = NULL;
            SeqRegions* lower_regions = selected_node->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
            
            // merge 2 lower vector into one
            new_root_lh = lower_regions->mergeTwoLowers(merged_root_sample_regions, default_blength, sample, default_blength, aln, model, params->threshold_prob, cumulative_rate, true);
            
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
                new_root_lh = lower_regions->mergeTwoLowers(merged_root_sample_regions, new_blength, sample, default_blength, aln, model, params->threshold_prob, cumulative_rate, true);
                
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
            SeqRegions* upper_left_right_regions = selected_node->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
            SeqRegions* lower_regions = selected_node->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
            upper_left_right_regions->mergeUpperLower(best_parent_regions, selected_node->length * 0.5, lower_regions, selected_node->length * 0.5, aln, model, params->threshold_prob);
            RealNumType new_split = 0.25;
            
            SeqRegions* new_parent_regions = NULL;
            while (new_split * selected_node->length > min_blength)
            {
                upper_left_right_regions->mergeUpperLower(new_parent_regions, selected_node->length * (1 - new_split), lower_regions, selected_node->length * new_split, aln, model, params->threshold_prob);
                
                RealNumType placement_cost = calculateSamplePlacementCost(aln, model, cumulative_rate, new_parent_regions, sample, default_blength);
                
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
            }
            // delete new_parent_regions
            if (new_parent_regions) delete new_parent_regions;
            
            best_parent_lh = best_split_lh;
            best_parent_split = best_split;
        }
        
        // if the best placement is below the selected_node => add an internal node below the selected_node
        if (best_child_lh >= best_parent_lh && best_child_lh >= best_lh_diff)
        {
            SeqRegions* upper_left_right_regions = best_child->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);

            RealNumType new_branch_length_lh = best_child_lh;
            RealNumType best_length = default_blength;
            
            while (best_length > min_blength)
            {
                RealNumType new_blength = best_length * 0.5;
                RealNumType placement_cost = calculateSamplePlacementCost(aln, model, cumulative_rate, best_child_regions, sample, new_blength);
                
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
                    RealNumType placement_cost = calculateSamplePlacementCost(aln, model, cumulative_rate, best_child_regions, sample, new_blength);
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
                RealNumType tmp_lh = calculateSamplePlacementCost(aln, model, cumulative_rate, best_child_regions, sample, -1);
                
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
            upper_left_right_regions->mergeUpperLower(next_node_2->partial_lh, new_internal_node->length, sample, best_length, aln, model, params->threshold_prob);
            best_child->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate)->mergeTwoLowers(new_internal_node->partial_lh, best_child->length, sample, best_length, aln, model, params->threshold_prob, cumulative_rate);
            upper_left_right_regions->mergeUpperLower(new_internal_node->mid_branch_lh, new_internal_node->length * 0.5, new_internal_node->partial_lh, new_internal_node->length * 0.5, aln, model, params->threshold_prob);
            new_internal_node->computeTotalLhAtNode(aln, model, params->threshold_prob, cumulative_rate, new_internal_node == root);
            
            /*if (top_distance >= 2 * min_blengthForMidNode)
                createFurtherMidNodes(new_internal_node,this_node_upper_left_right_regions)*/
            if (best_length > 0)
            {
                new_sample_node->computeTotalLhAtNode(aln, model, params->threshold_prob, cumulative_rate, new_sample_node == root);
                next_node_1->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate)->mergeUpperLower(new_sample_node->mid_branch_lh, new_sample_node->length * 0.5, sample, new_sample_node->length * 0.5, aln, model, params->threshold_prob);
                /*if best_length>=2*min_blengthForMidNode:
                    createFurtherMidNodes(new_sample_node,best_child_regions)*/
            }
            
            // update pseudo_count
            model->updatePesudoCount(aln, next_node_1->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate), sample);

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
                if (!best_parent_regions) best_parent_regions = new SeqRegions();
                best_parent_regions->copyRegions(selected_node->total_lh, num_states);
                if (selected_node == root)
                {
                    delete best_parent_regions;
                    best_parent_regions = NULL;
                    
                    selected_node->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate)->mergeTwoLowers(best_parent_regions, -1, sample, default_blength, aln, model, params->threshold_prob, cumulative_rate);
                }
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
                    
                    new_root_lh = selected_node->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate)->mergeTwoLowers(new_root_lower_regions, best_root_blength, sample, new_blength, aln, model, params->threshold_prob, cumulative_rate, true);
                    
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
                        new_root_lh = selected_node->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate)->mergeTwoLowers(new_root_lower_regions, best_root_blength, sample, new_blength, aln, model, params->threshold_prob, cumulative_rate, true);
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
                    new_root_lh = selected_node->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate)->mergeTwoLowers(new_root_lower_regions, best_root_blength, sample, -1, aln, model, params->threshold_prob, cumulative_rate, true);
                    new_root_lh += new_root_lower_regions->computeAbsoluteLhAtRoot(aln, model, cumulative_base);
                    RealNumType root_lh_diff = new_root_lh - old_root_lh;
                    if (root_lh_diff > best_parent_lh)
                    {
                        best_length2 = -1;
                        best_parent_lh = root_lh_diff;
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
                new_root->total_lh = new_root->computeTotalLhAtNode(aln, model, params->threshold_prob, cumulative_rate, true);

                next_node_1->partial_lh = selected_node->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate)->computeTotalLhAtRoot(num_states, model, best_root_blength);
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
                    new_sample_node->computeTotalLhAtNode(aln, model, params->threshold_prob, cumulative_rate, new_sample_node == root);
                    
                    new_root->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate)->mergeUpperLower(new_sample_node->mid_branch_lh, new_sample_node->length * 0.5, sample, new_sample_node->length * 0.5, aln, model, params->threshold_prob);
                    
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
                SeqRegions* upper_left_right_regions = selected_node->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
                
                // now try different lengths for the new branch
                RealNumType new_branch_length_lh = best_parent_lh;
                RealNumType best_length = default_blength;
                while (best_length > min_blength)
                {
                    RealNumType new_blength = best_length * 0.5;
                    RealNumType placement_cost = calculateSamplePlacementCost(aln, model, cumulative_rate, best_parent_regions, sample, new_blength);
                    
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
                        RealNumType placement_cost = calculateSamplePlacementCost(aln, model, cumulative_rate, best_parent_regions, sample, new_blength);
                        
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
                    RealNumType placement_cost = calculateSamplePlacementCost(aln, model, cumulative_rate, best_parent_regions, sample, -1);
                    
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
                RealNumType top_distance = selected_node->length * (1.0 - best_parent_split);
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
                upper_left_right_regions->mergeUpperLower(next_node_2->partial_lh, new_internal_node->length, sample, best_length, aln, model, params->threshold_prob);
                selected_node->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate)->mergeTwoLowers(new_internal_node->partial_lh, selected_node->length, sample, best_length, aln, model, params->threshold_prob, cumulative_rate);
                upper_left_right_regions->mergeUpperLower(new_internal_node->mid_branch_lh, new_internal_node->length * 0.5, new_internal_node->partial_lh, new_internal_node->length * 0.5, aln, model, params->threshold_prob);
                new_internal_node->computeTotalLhAtNode(aln, model, params->threshold_prob, cumulative_rate, new_internal_node == root);
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
                    new_sample_node->computeTotalLhAtNode(aln, model, params->threshold_prob, cumulative_rate, new_sample_node == root);
                    
                    next_node_1->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate)->mergeUpperLower(new_sample_node->mid_branch_lh, new_sample_node->length * 0.5, sample, new_sample_node->length * 0.5, aln, model, params->threshold_prob);
                    /*if best_length>=2*min_blengthForMidNode:
                        createFurtherMidNodes(new_sample_node,best_parent_regions)*/
                }
                
                // update pseudo_count
                model->updatePesudoCount(aln, next_node_1->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate), sample);

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
                        updateZeroBlength(next_node_1->neighbor, node_stack, aln, model, params->threshold_prob, cumulative_rate, default_blength, min_blength, max_blength);
                        updatePartialLh(node_stack, cumulative_rate, default_blength, min_blength, max_blength);
                    }
                    else if (next_node_2->length <= 0)
                    {
                        stack<Node*> node_stack;
                        updateZeroBlength(next_node_2->neighbor, node_stack, aln, model, params->threshold_prob, cumulative_rate, default_blength, min_blength, max_blength);
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
    StateType num_states = aln->num_states;
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
                        parent_upper_lr_lh->mergeUpperLower(node->mid_branch_lh, node->length * 0.5, lower_lh, node->length * 0.5, aln, model, threshold_prob);
                        
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
                    parent_upper_lr_lh->mergeUpperLower(new_upper_lr_lh, node->length, lower_lh, next_node_2->length, aln, model, threshold_prob);
                    
                    // if the upper left/right lh is null -> try to increase the branch length
                    if (!new_upper_lr_lh)
                    {
                        if (next_node_2->length <= 0)
                        {
                            stack<Node*> node_stack;
                            updateZeroBlength(next_node_2->neighbor, node_stack, aln, model, params->threshold_prob, cumulative_rate, default_blength, min_blength, max_blength);
                            updatePartialLh(node_stack, cumulative_rate, default_blength, min_blength, max_blength);
                        }
                        else if (node->length <= 0)
                        {
                            stack<Node*> node_stack;
                            updateZeroBlength(node, node_stack, aln, model, params->threshold_prob, cumulative_rate, default_blength, min_blength, max_blength);
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
                    parent_upper_lr_lh->mergeUpperLower(new_upper_lr_lh, node->length, lower_lh, next_node_1->length, aln, model, threshold_prob);
                    
                    // if the upper left/right lh is null -> try to increase the branch length
                    if (!new_upper_lr_lh)
                    {
                        if (next_node_1->length <= 0)
                        {
                            stack<Node*> node_stack;
                            updateZeroBlength(next_node_1->neighbor, node_stack, aln, model, params->threshold_prob, cumulative_rate, default_blength, min_blength, max_blength);
                            updatePartialLh(node_stack, cumulative_rate, default_blength, min_blength, max_blength);
                        }
                        else if (node->length <= 0)
                        {
                            stack<Node*> node_stack;
                            updateZeroBlength(node, node_stack, aln, model, params->threshold_prob, cumulative_rate, default_blength, min_blength, max_blength);
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

RealNumType Tree::improveEntireTree(RealNumType *cumulative_rate, RealNumType default_blength, RealNumType max_blength, RealNumType min_blength)
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
            RealNumType improvement = improveSubTree(node, cumulative_rate, default_blength, max_blength, min_blength);
            
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
                cout << "Processed topology for " + convertIntToString(num_nodes) + " nodes." << endl;
        }
    }
    
    return total_improvement;
}

RealNumType Tree::improveSubTree(Node* node, RealNumType *cumulative_rate, RealNumType default_blength, RealNumType max_blength, RealNumType min_blength)
{
    // dummy variables
    RealNumType threshold_prob = params->threshold_prob;
    RealNumType threshold_prob2 = threshold_prob * threshold_prob;
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
        RealNumType best_lh = calculateSubTreePlacementCost(aln, model, cumulative_rate, parent_upper_lr_lh, lower_lh, best_blength);
        RealNumType original_lh = best_lh;
        
        // optimize branch length
        if (best_lh < params->thresh_placement_cost)
        {
            // try different branch lengths for the current node placement (just in case branch length can be improved, in which case it counts both as tree improvment and better strategy to find a new placement).
            if (node->length <= 0)
            {
                best_blength = min_blength;
                best_lh = calculateSubTreePlacementCost(aln, model, cumulative_rate, parent_upper_lr_lh, lower_lh, best_blength);
            }
            
            RealNumType best_split = 1;
            RealNumType new_split = 0.5;
            while (new_split * best_blength > min_blength)
            {
                RealNumType new_lh = calculateSubTreePlacementCost(aln, model, cumulative_rate, parent_upper_lr_lh, lower_lh, new_split * best_blength);
                
                if (new_lh > best_lh)
                {
                    best_lh = new_lh;
                    best_split = new_split;
                    blength_changed = true;
                }
                else
                    break;
                
                new_split = best_split * 0.5;
            }
            
            if (best_split > 0.7)
            {
                new_split = 2;
                while (new_split * best_blength < max_blength)
                {
                    RealNumType new_lh = calculateSubTreePlacementCost(aln, model, cumulative_rate, parent_upper_lr_lh, lower_lh, new_split * best_blength);
                    
                    if (new_lh > best_lh)
                    {
                        best_lh = new_lh;
                        best_split = new_split;
                        blength_changed = true;
                    }
                    else
                        break;
                    
                    new_split = best_split * 2;
                }
            }
            
            best_blength = best_blength * best_split;
            
            if (node->length <= 0)
                if (original_lh > best_lh)
                    best_lh = original_lh;
        }
           
        // find new placement
        if (best_lh < params->thresh_placement_cost)
        {
            // now find the best place on the tree where to re-attach the subtree rooted at "node" but to do that we need to consider new vector probabilities after removing the node that we want to replace this is done using findBestParentTopology().
            bool topology_updated = false;
            
            Node* best_node = NULL;
            RealNumType best_lh_diff = -DBL_MAX;
            bool is_mid_node = false;
            
            // NHANLT: ongoing work
            //findBestParent(parent_node,child,best_lh,best_blength,mutMatrix,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology);
            
            // validate the new placement cost
            if (best_lh_diff > threshold_prob2)
                outError("Strange, lh cost is positive");
            else if (best_lh_diff < -1e50)
                outError("Likelihood cost is very heavy, this might mean that the reference used is not the same used to generate the input diff file");
            
            if (best_lh_diff + params->thresh_placement_cost > best_lh)
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
                        // NHANLT: ongoing work
                        // newRoot = cutAndPasteNode(node,best_node,is_mid_node,best_blength,best_lh_diff,mutMatrix)
                        
                        topology_updated = true;
                    }
                }
                if (!topology_updated && blength_changed)
                {
                    node->length = best_blength;
                    
                    stack<Node*> node_stack;
                    node_stack.push(node);
                    node_stack.push(node->neighbor);
                    updatePartialLh(node_stack, cumulative_rate, default_blength, min_blength, max_blength);
                }
            }
            else if (blength_changed)
            {
                node->length = best_blength;
                
                stack<Node*> node_stack;
                node_stack.push(node);
                node_stack.push(node->neighbor);
                updatePartialLh(node_stack, cumulative_rate, default_blength, min_blength, max_blength);
            }
        }
        else if (blength_changed)
        {
            node->length = best_blength;
            
            stack<Node*> node_stack;
            node_stack.push(node);
            node_stack.push(node->neighbor);
            updatePartialLh(node_stack, cumulative_rate, default_blength, min_blength, max_blength);
        }
    }
                        
    return total_improvement;
}

RealNumType Tree::calculateSubTreePlacementCost(Alignment* aln, Model* model, RealNumType* cumulative_rate, SeqRegions* parent_regions, SeqRegions* child_regions, RealNumType blength)
{
    return (this->*calculateSubTreePlacementCostPointer)(aln, model, cumulative_rate, parent_regions, child_regions, blength);
}

// this implementation derives from appendProbNode
template <const StateType num_states>
RealNumType Tree::calculateSubTreePlacementCostTemplate(Alignment* aln, Model* model, RealNumType* cumulative_rate, SeqRegions* parent_regions, SeqRegions* child_regions, RealNumType blength)
{
    // init dummy variables
    RealNumType lh_cost = 0;
    PositionType seq1_index = -1;
    PositionType seq2_index = -1;
    PositionType pos = 0;
    RealNumType total_factor = 1;
    SeqRegion *seq1_region, *seq2_region;
    PositionType end_pos;
    RealNumType minimum_carry_over = DBL_MIN * 1e50;
    RealNumType total_blength = blength;
    PositionType seq_length = aln->ref_seq.size();
    
    while (pos < seq_length)
    {
        // get the next shared segment in the two sequences
        SeqRegions::getNextSharedSegment(pos, parent_regions, child_regions, seq1_index, seq2_index, seq1_region, seq2_region, end_pos);
        
        // 1. e1.type = N || e2.type = N
        if (seq2_region->type == TYPE_N || seq1_region->type == TYPE_N)
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
        
        // 2. e1.type = R
        if (seq1_region->type == TYPE_R)
        {
            // 2.1. e1.type = R and e2.type = R
            if (seq2_region->type == TYPE_R)
            {
                if (seq1_region->plength_observation2root >= 0)
                    total_blength += seq1_region->plength_observation2node;
                
                if (total_blength > 0)
                    lh_cost += total_blength * (cumulative_rate[end_pos + 1] - cumulative_rate[pos]);
            }
            // 2.2. e1.type = R and e2.type = O
            else if (seq2_region->type == TYPE_O)
            {
                RealNumType tot = 0;
                StateType seq1_state = aln->ref_seq[end_pos];
                
                if (seq1_region->plength_observation2root >= 0)
                {
                    StateType mutation_index = model->row_index[seq1_state];
                    for (StateType i = 0; i < num_states; ++i)
                    {
                        RealNumType tot2;
                        
                        if (seq1_state == i)
                            tot2 = model->root_freqs[i] * (1.0 + model->transposed_mut_mat[mutation_index + i] * seq1_region->plength_observation2node);
                        else
                            tot2 = model->root_freqs[i] * (model->transposed_mut_mat[mutation_index + i] * seq1_region->plength_observation2node);
                            
                        RealNumType tot3 = 0;
                        if (total_blength > 0)
                        {
                            for (StateType j = 0; j < num_states; ++j)
                                tot3 += model->mutation_mat[model->row_index[i] + j] * seq2_region->likelihood[j];
                        }
                        
                        tot += tot2 * (seq2_region->likelihood[i] + total_blength * tot3);
                    }
                    
                    tot *= model->inverse_root_freqs[seq1_state];
                }
                else
                {
                    if (total_blength > 0)
                    {
                        StateType mutation_index = model->row_index[seq1_state];
                        for (StateType j = 0; j < num_states; ++j)
                            tot += model->mutation_mat[mutation_index + j] * seq2_region->likelihood[j];
                        
                        tot *= total_blength;
                    }
                    
                    tot += seq2_region->likelihood[seq1_state];
                }
                    
                total_factor *= tot;
            }
            // 2.3. e1.type = R and e2.type = A/C/G/T
            else
            {
                StateType seq1_state = aln->ref_seq[end_pos];
                StateType seq2_state = seq2_region->type;
                
                if (seq1_region->plength_observation2root >= 0)
                {
                    if (total_blength > 0)
                    {
                        RealNumType seq1_state_evolves_seq2_state = model->mutation_mat[model->row_index[seq1_state] + seq2_state] * total_blength * (1.0 + model->diagonal_mut_mat[seq1_state] * seq1_region->plength_observation2node);
                        
                        RealNumType seq2_state_evolves_seq1_state = model->freqi_freqj_qij[model->row_index[seq2_state] + seq1_state] * seq1_region->plength_observation2node * (1.0 + model->diagonal_mut_mat[seq2_state] * total_blength);
                        
                        total_factor *= seq1_state_evolves_seq2_state + seq2_state_evolves_seq1_state;
                    }
                    else
                        total_factor *= model->freqi_freqj_qij[model->row_index[seq2_state] + seq1_state] * seq1_region->plength_observation2node;
                }
                else if (total_blength > 0)
                    total_factor *= model->mutation_mat[model->row_index[seq1_state] + seq2_state] * total_blength;
                else
                    return -DBL_MAX;
            }
        }
        // 3. e1.type = O
        else if (seq1_region->type == TYPE_O)
        {
            // 3.1. e1.type = O and e2.type = O
            if (seq2_region->type == TYPE_O)
            {
                RealNumType tot = 0;
                
                if (total_blength > 0)
                {
                    for (StateType i = 0; i < num_states; ++i)
                    {
                        RealNumType tot2 = 0;
                        
                        for (StateType j = 0; j < num_states; ++j)
                            tot2 += model->mutation_mat[model->row_index[i] + j] * seq2_region->likelihood[j];
                        
                        tot += seq1_region->likelihood[i] * (seq2_region->likelihood[i] + total_blength * tot2);
                    }
                }
                else
                    for (StateType i = 0; i < num_states; ++i)
                        tot += seq1_region->likelihood[i] * seq2_region->likelihood[i];
                
                total_factor *= tot;
            }
            // 3.2. e1.type = O and e2.type = R or A/C/G/T
            else
            {
                StateType seq2_state = seq2_region->type;
                if (seq2_state == TYPE_R)
                    seq2_state = aln->ref_seq[end_pos];
                
                if (total_blength > 0)
                {
                    RealNumType tot2 = 0;
                    StateType mutation_index = model->row_index[seq2_state];
                    for (StateType j = 0; j < num_states; ++j)
                        tot2 += seq1_region->likelihood[j] * model->transposed_mut_mat[mutation_index + j];
                    
                    total_factor *= seq1_region->likelihood[seq2_state] + total_blength * tot2;
                }
                else
                    total_factor *= seq1_region->likelihood[seq2_state];
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
                
                if (total_blength > 0)
                    lh_cost += model->diagonal_mut_mat[seq1_state] * total_blength;
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
                        StateType mutation_index = model->row_index[seq1_state];
                        for (StateType i = 0; i < num_states; ++i)
                        {
                            RealNumType tot2;
                            
                            if (seq1_state == i)
                                tot2 = model->root_freqs[i] * (1.0 + model->transposed_mut_mat[mutation_index + i] * seq1_region->plength_observation2node);
                            else
                                tot2 = model->root_freqs[i] * (model->transposed_mut_mat[mutation_index + i] * seq1_region->plength_observation2node);
                                
                            RealNumType tot3 = 0;
                            for (StateType j = 0; j < num_states; ++j)
                                tot3 += model->mutation_mat[model->row_index[i] + j] * seq2_region->likelihood[j];
                            tot += tot2 * (seq2_region->likelihood[i] + total_blength * tot3);
                        }
                        
                        total_factor *= (tot * model->inverse_root_freqs[seq1_state]);
                    }
                    else
                    {
                        StateType mutation_index = model->row_index[seq1_state];
                        for (StateType j = 0; j < num_states; ++j)
                            tot += model->mutation_mat[mutation_index + j] * seq2_region->likelihood[j];
                        
                        tot *= total_blength;
                        tot += seq2_region->likelihood[seq1_state];
                        total_factor *= tot;
                    }
                }
                // 4.3. e1.type = A/C/G/T and e2.type = R or A/C/G/T
                else
                {
                    StateType seq2_state = seq2_region->type;
                    if (seq2_state == TYPE_R)
                        seq2_state = aln->ref_seq[end_pos];
                    
                    if (seq1_region->plength_observation2root >= 0)
                    {
                        if (total_blength > 0)
                        {
                            RealNumType seq1_state_evolves_seq2_state = model->mutation_mat[model->row_index[seq1_state] + seq2_state] * total_blength * (1.0 + model->diagonal_mut_mat[seq1_state] * seq1_region->plength_observation2node);
                            
                            RealNumType seq2_state_evolves_seq1_state = model->freqi_freqj_qij[model->row_index[seq2_state] + seq1_state] * seq1_region->plength_observation2node * (1.0 + model->diagonal_mut_mat[seq2_state] * total_blength);
                            
                            total_factor *= seq1_state_evolves_seq2_state + seq2_state_evolves_seq1_state;
                        }
                        else
                            total_factor *= model->freqi_freqj_qij[model->row_index[seq2_state] + seq1_state] * seq1_region->plength_observation2node;
                    }
                    else if (total_blength > 0)
                        total_factor *= model->mutation_mat[model->row_index[seq1_state] + seq2_state] * total_blength;
                    else
                        return -DBL_MAX;
                }
            }
        }
         
        // approximately update lh_cost and total_factor
        if (total_factor <= minimum_carry_over)
        {
            if (total_factor < DBL_MIN)
                return -DBL_MAX;
            lh_cost += log(total_factor);
            total_factor = 1.0;
        }
        
        // update pos
        pos = end_pos + 1;
    }
    
    return lh_cost + log(total_factor);
}

RealNumType Tree::calculateSamplePlacementCost(Alignment* aln, Model* model, RealNumType* cumulative_rate, SeqRegions* parent_regions, SeqRegions* child_regions, RealNumType blength)
{
    return (this->*calculateSamplePlacementCostPointer)(aln, model, cumulative_rate, parent_regions, child_regions, blength);
}

// this implementation derives from appendProb
#define MIN_CARRY_OVER 1e-250
template <const StateType num_states>
RealNumType Tree::calculateSamplePlacementCostTemplate(Alignment* aln, Model* model, RealNumType* cumulative_rate, SeqRegions* parent_regions, SeqRegions* child_regions, RealNumType blength)
{
    // init dummy variables
    RealNumType lh_cost = 0;
    PositionType seq1_index = -1;
    PositionType seq2_index = -1;
    PositionType pos = 0;
    RealNumType total_factor = 1;
    SeqRegion *seq1_region, *seq2_region;
    PositionType end_pos;
    //RealNumType minimum_carry_over = DBL_MIN * 1e50;
    if (blength < 0) blength = 0;
    RealNumType total_blength = blength;
    PositionType seq_length = aln->ref_seq.size();
    
    while (pos < seq_length)
    {
        // get the next shared segment in the two sequences
        SeqRegions::getNextSharedSegment(pos, parent_regions, child_regions, seq1_index, seq2_index, seq1_region, seq2_region, end_pos);
        
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
                    StateType seq1_state = aln->ref_seq[end_pos];
                    if (seq1_region->plength_observation2root >= 0)
                    {
                        total_blength = seq1_region->plength_observation2root + blength;
                        
                        if (seq2_region->likelihood[seq1_state] > 0.1)
                        {
                            total_blength += seq1_region->plength_observation2node;
                            
                            // here contribution from root frequency can also be also ignored
                            lh_cost += model->diagonal_mut_mat[seq1_state] * total_blength;
                        }
                        else
                        {
                            RealNumType tot = 0;
                            StateType mutation_index = model->row_index[seq1_state];
                            for (StateType i = 0; i < num_states; ++i)
                            {
                                RealNumType tot2;
                                
                                if (seq1_state == i)
                                    tot2 = model->root_freqs[i] * (1.0 + model->transposed_mut_mat[mutation_index + i] * seq1_region->plength_observation2node);
                                else
                                    tot2 = model->root_freqs[i] * (model->transposed_mut_mat[mutation_index + i] * seq1_region->plength_observation2node);
                                    
                                RealNumType tot3 = 0;
                                if (seq2_region->likelihood[i] > 0.1)
                                    tot3 = 1;
                                
                                for (StateType j = 0; j < num_states; ++j)
                                    if (seq2_region->likelihood[j] > 0.1)
                                        tot3 += model->mutation_mat[model->row_index[i] + j];
                                tot3 *= total_blength;
                                
                                tot += tot2 * tot3;
                            }
                            
                            total_factor *= tot * model->inverse_root_freqs[seq1_state];
                        }
                    }
                    else
                    {
                        if (seq2_region->likelihood[seq1_state] > 0.1)
                        {
                            if (seq1_region->plength_observation2node >= 0)
                                lh_cost += model->diagonal_mut_mat[seq1_state] * (blength + seq1_region->plength_observation2node);
                            else
                                lh_cost += model->diagonal_mut_mat[seq1_state] * blength;
                        }
                        else
                        {
                            RealNumType tot = 0;
                            StateType mutation_index = model->row_index[seq1_state];
                            for (StateType i = 0; i < num_states; ++i)
                                if (seq2_region->likelihood[i] > 0.1)
                                    tot += model->mutation_mat[mutation_index + i];
                            
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
                    StateType seq1_state = aln->ref_seq[end_pos];
                    StateType seq2_state = seq2_region->type;
                    
                    if (seq1_region->plength_observation2root >= 0)
                    {
                        RealNumType seq1_state_evolves_seq2_state = model->mutation_mat[model->row_index[seq1_state] + seq2_state] * blength * (1.0 + model->diagonal_mut_mat[seq1_state] * seq1_region->plength_observation2node);
                        
                        RealNumType seq2_state_evolves_seq1_state = model->freqi_freqj_qij[model->row_index[seq2_state] + seq1_state] * seq1_region->plength_observation2node * (1.0 + model->diagonal_mut_mat[seq2_state] * (blength + seq1_region->plength_observation2root));
                                                                                                                                                                                                    
                        total_factor *= seq1_state_evolves_seq2_state + seq2_state_evolves_seq1_state;
                    }
                    else
                    {
                        total_factor *= model->mutation_mat[model->row_index[seq1_state] + seq2_state] * (blength + (seq1_region->plength_observation2node < 0 ? 0 : seq1_region->plength_observation2node));
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
                    
                    for (StateType i = 0; i < num_states; ++i)
                    {
                        RealNumType tot2 = 0;
                        
                        for (StateType j = 0; j < num_states; ++j)
                            if (seq2_region->likelihood[j] > 0.1)
                                tot2 += model->mutation_mat[model->row_index[i] + j];
                        
                        tot2 *= blength13;
                        
                        if (seq2_region->likelihood[i] > 0.1)
                            tot2 += 1;
                        
                        tot += tot2 * seq1_region->likelihood[i];
                    }
                        
                    total_factor *= tot;
                }
                // 3.2. e1.type = O and e2.type = R or A/C/G/T
                else
                {
                    StateType seq2_state = seq2_region->type;
                    if (seq2_state == TYPE_R)
                        seq2_state = aln->ref_seq[end_pos];
                    
                    RealNumType tot2 = 0;
                    StateType mutation_index = model->row_index[seq2_state];
                    for (StateType j = 0; j < num_states; ++j)
                        tot2 += model->transposed_mut_mat[mutation_index + j] * seq1_region->likelihood[j];
                    
                    total_factor *= seq1_region->likelihood[seq2_state] + blength13 * tot2;
                }
            }
            // 4. e1.type = A/C/G/T
            else
            {
                int seq1_state = seq1_region->type;
                
                // 4.1. e1.type =  e2.type
                if (seq1_region->type == seq2_region->type)
                {
                    RealNumType total_blength = blength;
                    total_blength += (seq1_region->plength_observation2node < 0 ? 0 : seq1_region->plength_observation2node);
                    total_blength += (seq1_region->plength_observation2root < 0 ? 0 : seq1_region->plength_observation2root);

                    lh_cost += model->diagonal_mut_mat[seq1_state] * total_blength;
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
                            
                            if (seq2_region->likelihood[seq1_state] > 0.1)
                                lh_cost += model->diagonal_mut_mat[seq1_state] * (blength15 + seq1_region->plength_observation2node);
                            else
                            {
                                StateType mutation_index = model->row_index[seq1_state];

                                for (StateType i = 0; i < num_states; ++i)
                                {
                                    RealNumType tot2;
                                    if (seq1_state == i)
                                        tot2 = model->root_freqs[i] * (1.0 + model->transposed_mut_mat[mutation_index + i] * seq1_region->plength_observation2node);
                                    else
                                        tot2 = model->root_freqs[i] * (model->transposed_mut_mat[mutation_index + i] * seq1_region->plength_observation2node);
                                        
                                    RealNumType tot3 = 0;
                                    for (StateType j = 0; j < num_states; ++j)
                                        if (seq2_region->likelihood[j] > 0.1)
                                            tot3 += model->mutation_mat[model->row_index[i] + j];
                                    
                                    if (seq2_region->likelihood[i] > 0.1)
                                        tot += tot2 * (1.0 + blength15 * tot3);
                                    else
                                        tot += tot2 * blength15 * tot3;
                                }
                                
                                total_factor *= (tot * model->inverse_root_freqs[seq1_state]);
                            }
                        }
                        else
                        {
                            RealNumType tmp_blength = blength + (seq1_region->plength_observation2node < 0 ? 0 : seq1_region->plength_observation2node);
                            if (seq2_region->likelihood[seq1_state] > 0.1)
                                lh_cost += model->diagonal_mut_mat[seq1_state] * tmp_blength;
                            else
                            {
                                StateType start_index = model->row_index[seq1_state];
                                for (StateType j = 0; j < num_states; ++j)
                                    if (seq2_region->likelihood[j] > 0.1)
                                        tot += model->mutation_mat[start_index + j];
                                
                                total_factor *= tot * tmp_blength;
                            }
                        }
                    }
                    // 4.3. e1.type = A/C/G/T and e2.type = R or A/C/G/T
                    else
                    {
                        StateType seq2_state = seq2_region->type;
                        if (seq2_state == TYPE_R)
                            seq2_state = aln->ref_seq[end_pos];
                        
                        if (seq1_region->plength_observation2root >= 0)
                        {
                            // here we ignore contribution of non-parsimonious mutational histories
                            RealNumType seq1_state_evoloves_seq2_state = model->mutation_mat[model->row_index[seq1_state] + seq2_state] * (blength + seq1_region->plength_observation2root) * (1.0 + model->diagonal_mut_mat[seq1_state] * seq1_region->plength_observation2node);
                            
                            RealNumType seq2_state_evolves_seq2_state = model->freqi_freqj_qij[model->row_index[seq2_state] + seq1_state] * seq1_region->plength_observation2node * (1.0 + model->diagonal_mut_mat[seq2_state] * (blength + seq1_region->plength_observation2root));
                            
                            total_factor *= seq1_state_evoloves_seq2_state + seq2_state_evolves_seq2_state;
                        }
                        else
                        {
                            RealNumType tmp_blength = blength + (seq1_region->plength_observation2node < 0 ? 0 : seq1_region->plength_observation2node);
                            
                            total_factor *= model->mutation_mat[model->row_index[seq1_state] + seq2_state] * tmp_blength;
                        }
                    }
                }
            }
        }
         
        // approximately update lh_cost and total_factor
        if (total_factor <= MIN_CARRY_OVER)
        {
            if (total_factor < DBL_MIN)
                return -DBL_MAX;
            lh_cost += log(total_factor);
            total_factor = 1.0;
        }
        
        // update pos
        pos = end_pos + 1;
    }
    
    return lh_cost + log(total_factor);
}

void Tree::updateZeroBlength(Node* node, stack<Node*> &node_stack, Alignment* aln, Model* model, RealNumType threshold_prob, RealNumType* cumulative_rate, RealNumType default_blength, RealNumType min_blength, RealNumType max_blength)
{
    // get the top node in the phylo-node
    Node* top_node = node->getTopNode();
    ASSERT(top_node);
    SeqRegions* upper_left_right_regions = top_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
    SeqRegions* lower_regions = top_node->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
    
    RealNumType best_lh = calculateSamplePlacementCost(aln, model, cumulative_rate, upper_left_right_regions, lower_regions, default_blength);
    RealNumType best_length = default_blength;
    
    while (best_length > min_blength)
    {
        RealNumType new_blength = best_length * 0.5;
        RealNumType new_lh = calculateSamplePlacementCost(aln, model, cumulative_rate, upper_left_right_regions, lower_regions, new_blength);
        
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
            RealNumType new_blength = best_length * 2;
            RealNumType new_lh = calculateSamplePlacementCost(aln, model, cumulative_rate, upper_left_right_regions, lower_regions, new_blength);
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
    top_node->neighbor->outdated = true;
    node_stack.push(top_node);
    node_stack.push(top_node->neighbor);
}
