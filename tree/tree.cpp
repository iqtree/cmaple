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

void Tree::addNode(Node* new_node, Node* current_node, int &new_id, bool insert_to_right)
{
    Node* new_node_1 = new Node(new_id++);
    Node* new_node_2 = new Node(new_id++);
    Node* new_node_3 = new Node(new_id++);
    new_node_1->next = new_node_2;
    new_node_2->next = new_node_3;
    new_node_3->next = new_node_1;
    
    new_node_1->neighbor = current_node->neighbor;
    if (current_node == root)
        root = new_node_1;
    else
        new_node_1->neighbor->neighbor = new_node_1;
    
    Node *left_tip, *right_tip;
    if (insert_to_right)
    {
        left_tip = current_node;
        right_tip = new_node;
    }
    else
    {
        left_tip = new_node;
        right_tip = current_node;
    }
    
    new_node_2->neighbor = right_tip;
    right_tip->neighbor = new_node_2;
    
    new_node_3->neighbor = left_tip;
    left_tip->neighbor = new_node_3;
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
    StateType num_states = aln->num_states;
    PositionType seq_length = aln->ref_seq.size();
    
    while (!node_stack.empty())
    {
        Node* node = node_stack.top();
        node_stack.pop();
        
        // NHANLT: debug
       /* if (node->next && node->next->neighbor && node->next->neighbor->seq_name == "711")
            cout << "dsdas";*/
        
        bool update_blength = false;
        node->dirty = true;
        
        Regions* parent_upper_regions = NULL;
        if (root != node->getTopNode())
            parent_upper_regions = node->getTopNode()->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
            
        // change in likelihoods is coming from parent node
        if (node->real_number_attributes.find(IS_TOP_NODE) != node->real_number_attributes.end())
        {
            // if necessary, update the total probabilities at the mid node.
            if (node->length > 0)
            {
                // update vector of regions at mid-branch point
                Regions* mid_branch_regions = NULL;
                Regions* lower_regions = node->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
                parent_upper_regions->mergeUpperLower(mid_branch_regions, node->length / 2, lower_regions, node->length / 2, aln, model, params->threshold_prob);
                
                if (!mid_branch_regions)
                {
                    if (node->length > 1e-100)
                        outError("inside updatePartialLh(), from parent: should not have happened since node->length > 0");
                    node->updateZeroBlength(node_stack, aln, model, params->threshold_prob, cumulative_rate, default_blength, min_blength, max_blength);
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
                
                Regions* upper_left_right_regions_1 = NULL;
                Regions* upper_left_right_regions_2 = NULL;
                Regions* lower_regions_1 = next_node_1->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
                parent_upper_regions->mergeUpperLower(upper_left_right_regions_1, node->length, lower_regions_1, next_node_1->length, aln, model, params->threshold_prob);
                
                if (!upper_left_right_regions_1 || upper_left_right_regions_1->size() == 0)
                {
                    if (node->length <= 0 && next_node_1->length <= 0)
                        node->updateZeroBlength(node_stack, aln, model, params->threshold_prob, cumulative_rate, default_blength, min_blength, max_blength);
                    else
                        outError("Strange: None vector from non-zero distances in updatePartialLh() from parent direction.");
                    
                    update_blength = true;
                }
                
                if (!update_blength)
                {
                    Regions* lower_regions_2 = next_node_2->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
                    parent_upper_regions->mergeUpperLower(upper_left_right_regions_2, node->length, lower_regions_2, next_node_2->length, aln, model, params->threshold_prob);
                    
                    if (!upper_left_right_regions_2 || upper_left_right_regions_2->size() == 0)
                    {
                        if (node->length <= 0 && next_node_2->length <= 0)
                        {
                            node->updateZeroBlength(node_stack, aln, model, params->threshold_prob, cumulative_rate, default_blength, min_blength, max_blength);
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
            Node* top_node;
            Node* other_next_node;
            Node* next_node;
            FOR_NEXT(node, next_node)
            {
                if (next_node->real_number_attributes.find(IS_TOP_NODE) != next_node->real_number_attributes.end())
                    top_node = next_node;
                else
                    other_next_node = next_node;
            }
            
            ASSERT(top_node && other_next_node);
            
            RealNumType this_node_distance = node->length;
            RealNumType other_next_node_distance = other_next_node->length;
            
            Regions* this_node_lower_regions = node->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
            Regions* this_node_upper_left_right_regions = node->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
            Regions* next_node_upper_left_right_regions = other_next_node->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
            
            // update lower likelihoods
            Regions* merged_two_lower_regions = NULL;
            Regions* old_lower_regions = NULL;
            other_next_node->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate)->mergeTwoLowers(merged_two_lower_regions, other_next_node_distance, this_node_lower_regions, this_node_distance, aln, model, params->threshold_prob, cumulative_rate);
            
            if (!merged_two_lower_regions || merged_two_lower_regions->size() == 0)
            {
                if (this_node_distance <= 0 && other_next_node_distance <= 0)
                {
                    node->neighbor->updateZeroBlength(node_stack, aln, model, params->threshold_prob, cumulative_rate, default_blength, min_blength, max_blength);
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
                    Regions* new_total_lh_regions = top_node->computeTotalLhAtNode(aln, model, params->threshold_prob, cumulative_rate, top_node == root, false);
                    
                    if (!new_total_lh_regions && top_node->length <= 0)
                    {
                        top_node->updateZeroBlength(node_stack, aln, model, params->threshold_prob, cumulative_rate, default_blength, min_blength, max_blength);
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
                    Regions* new_mid_regions = NULL;
                    Regions* tmp_lower_regions = top_node->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
                    parent_upper_regions->mergeUpperLower(new_mid_regions, top_node->length / 2, tmp_lower_regions, top_node->length / 2, aln, model, params->threshold_prob);
                    
                    if (!new_mid_regions)
                    {
                        top_node->updateZeroBlength(node_stack, aln, model, params->threshold_prob, cumulative_rate, default_blength, min_blength, max_blength);
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
                Regions* new_upper_regions = NULL;
                if (root != top_node)
                    parent_upper_regions->mergeUpperLower(new_upper_regions, top_node->length, this_node_lower_regions, this_node_distance, aln, model, params->threshold_prob);
                else
                    new_upper_regions = node->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate)->computeTotalLhAtRoot(num_states, model, this_node_distance);
                
                if (!new_upper_regions || new_upper_regions->size() == 0)
                {
                    if (top_node->length <= 0 && this_node_distance <= 0)
                    {
                        top_node->updateZeroBlength(node_stack, aln, model, params->threshold_prob, cumulative_rate, default_blength, min_blength, max_blength);
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

void Tree::seekPlacement(Node* start_node, string seq_name, Regions* sample_regions, Node* &selected_node, RealNumType &best_lh_diff , bool &is_mid_branch, RealNumType &best_up_lh_diff, RealNumType &best_down_lh_diff, Node* &best_child, RealNumType* cumulative_rate, RealNumType default_blength, RealNumType min_blength_mid)
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
    stack<Node*> node_stack;
    start_node->real_number_attributes[LH_DIFF] = -DBL_MAX;
    start_node->real_number_attributes[FAILURE_COUNT] = 0;
    node_stack.push(start_node);
    
    // recursively examine positions for placing the new sample
    while (!node_stack.empty())
    {
        Node* current_node = node_stack.top();
        node_stack.pop();
        
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
            lh_diff_mid_branch = current_node->mid_branch_lh->calculatePlacementCost(aln, model, cumulative_rate, sample_regions, default_blength);
            
            // record the best_lh_diff if lh_diff_mid_branch is greater than the best_lh_diff ever
            if (lh_diff_mid_branch > best_lh_diff)
            {
                best_lh_diff = lh_diff_mid_branch;
                selected_node = current_node;
                current_node->real_number_attributes[FAILURE_COUNT] = 0;
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
            lh_diff_at_node = current_node->total_lh->calculatePlacementCost(aln, model, cumulative_rate, sample_regions, default_blength);
            
            // record the best_lh_diff if lh_diff_at_node is greater than the best_lh_diff ever
            if (lh_diff_at_node > best_lh_diff)
            {
                best_lh_diff = lh_diff_at_node;
                selected_node = current_node;
                current_node->real_number_attributes[FAILURE_COUNT] = 0;
                is_mid_branch = false;
                best_up_lh_diff = lh_diff_mid_branch;
            }
            else if (lh_diff_mid_branch >= (best_lh_diff - params->threshold_prob))
            {
                best_up_lh_diff = current_node->real_number_attributes[LH_DIFF];
                best_down_lh_diff = lh_diff_at_node;
                best_child = current_node;
            }
            // placement at current node is considered failed if placement likelihood is not improved by a certain margin compared to best placement so far for the nodes above it.
            else if (lh_diff_at_node < (current_node->real_number_attributes[LH_DIFF] - params->thresh_log_lh_failure))
                current_node->real_number_attributes[FAILURE_COUNT] = current_node->real_number_attributes[FAILURE_COUNT] + 1;
        }
        else
            lh_diff_at_node = current_node->real_number_attributes[LH_DIFF];
        
        // keep trying to place at children nodes, unless the number of attempts has reaches the failure limit
        if ((params->strict_stop_seeking_placement
             && current_node->real_number_attributes[FAILURE_COUNT]  < params->failure_limit
             && lh_diff_at_node > (best_lh_diff - params->thresh_log_lh_subtree_explore))
            || (!params->strict_stop_seeking_placement
                && (current_node->real_number_attributes[FAILURE_COUNT]  < params->failure_limit
                    || lh_diff_at_node > (best_lh_diff - params->thresh_log_lh_subtree_explore))))
        {
            Node* neighbor_node;
            FOR_NEIGHBOR(current_node, neighbor_node)
            {
                neighbor_node->real_number_attributes[LH_DIFF] = lh_diff_at_node;
                neighbor_node->real_number_attributes[FAILURE_COUNT] = current_node->real_number_attributes[FAILURE_COUNT];
                node_stack.push(neighbor_node);
            }
        }
        
        // remove real_number_attributes of the current node
        current_node->real_number_attributes.erase(FAILURE_COUNT);
        current_node->real_number_attributes.erase(LH_DIFF);
    }

    // exploration of the tree is finished, and we are left with the node found so far with the best appending likelihood cost. Now we explore placement just below this node for more fine-grained placement within its descendant branches.
    best_down_lh_diff = -DBL_MAX;
    best_child = NULL;
    
    // if best position so far is the descendant of a node -> explore further at its children
    if (!is_mid_branch)
    {
        // current node might be part of a polytomy (represented by 0 branch lengths) so we want to explore all the children of the current node to find out if the best placement is actually in any of the branches below the current node.
        Node* neighbor_node;
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
                RealNumType new_blength = node->length / 2;
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
                    new_blength /= 2;
                    
                    // stop trying if reaching the minimum branch length
                    if (new_blength <= min_blength_mid / 2)
                        break;
                }*/
                
                RealNumType new_best_lh_mid_branch = node->mid_branch_lh->calculatePlacementCost(aln, model, cumulative_rate, sample_regions, default_blength);
                
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

void Tree::placeNewSample(Node* selected_node, Regions* sample, string seq_name, RealNumType best_lh_diff , bool is_mid_branch, RealNumType best_up_lh_diff, RealNumType best_down_lh_diff, Node* best_child, RealNumType* cumulative_rate, vector<vector<PositionType>> &cumulative_base, RealNumType default_blength, RealNumType max_blength, RealNumType min_blength)
{
    // dummy variables
    RealNumType best_child_lh;
    RealNumType best_child_split;
    RealNumType best_parent_lh;
    RealNumType best_parent_split;
    Regions* best_parent_regions = NULL;
    Regions* best_child_regions = NULL;
    RealNumType best_root_blength;
    StateType num_states = aln->num_states;
    
    // try to place the new sample as a descendant of a mid-branch point
    if (is_mid_branch)
    {
        Regions* upper_left_right_regions = selected_node->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
        RealNumType best_split = 0.5;
        RealNumType best_split_lh = best_lh_diff;
        RealNumType new_split = 0.25;
        selected_node->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate)->mergeUpperLower(best_child_regions, selected_node->length / 2, selected_node->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate), selected_node->length / 2, aln, model, params->threshold_prob);
        Regions* new_parent_regions = NULL;

        // try different positions on the existing branch
        while (new_split * selected_node->length > min_blength)
        {
            upper_left_right_regions->mergeUpperLower(new_parent_regions, selected_node->length * new_split, selected_node->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate),  selected_node->length * (1 - new_split), aln, model, params->threshold_prob);
            RealNumType placement_cost = new_parent_regions->calculatePlacementCost(aln, model, cumulative_rate, sample, default_blength);
            
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
            new_split = best_split / 2;
        }
        
        if (best_split > 0.49)
        {
            new_split = 0.25;
            while (new_split * selected_node->length > min_blength)
            {
                upper_left_right_regions->mergeUpperLower(new_parent_regions, selected_node->length * (1.0 - new_split), selected_node->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate),selected_node->length * new_split, aln, model, params->threshold_prob);
                
                RealNumType placement_cost = new_parent_regions->calculatePlacementCost(aln, model, cumulative_rate, sample, default_blength);
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
                
                new_split = best_split / 2;
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
            RealNumType new_blength = best_blength / 2;
            RealNumType placement_cost = best_child_regions->calculatePlacementCost(aln, model, cumulative_rate, sample, new_blength);
            
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
                RealNumType placement_cost = best_child_regions->calculatePlacementCost(aln, model, cumulative_rate, sample, new_blength);
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
            RealNumType zero_branch_lh = best_child_regions->calculatePlacementCost(aln, model, cumulative_rate, sample, -1);
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
        upper_left_right_regions->mergeUpperLower(new_internal_node->mid_branch_lh, new_internal_node->length / 2, new_internal_node->partial_lh, new_internal_node->length / 2, aln, model, params->threshold_prob);
        new_internal_node->computeTotalLhAtNode(aln, model, params->threshold_prob, cumulative_rate, new_internal_node == root);
        
        if (!new_internal_node->total_lh || new_internal_node->total_lh->size() == 0)
            outError("Problem, None vector when placing sample, below node");
        
        /*if distTop>=2*min_blengthForMidNode:
         createFurtherMidNodes(newInternalNode,upper_left_right_regions)*/
        if (best_blength > 0)
        {
            new_sample_node->computeTotalLhAtNode(aln, model, params->threshold_prob, cumulative_rate, new_sample_node == root);
            next_node_1->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate)->mergeUpperLower(new_sample_node->mid_branch_lh, new_sample_node->length / 2, sample, new_sample_node->length / 2, aln, model, params->threshold_prob);
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
            Regions* upper_left_right_regions = best_child->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
            Regions* lower_regions = best_child->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
            best_child_regions = NULL;
            upper_left_right_regions->mergeUpperLower(best_child_regions, best_child->length / 2, lower_regions, best_child->length / 2, aln, model, params->threshold_prob);
            RealNumType new_split = 0.25;
            
            while (new_split * best_child->length > min_blength)
            {
                Regions* new_parent_regions = NULL;
                upper_left_right_regions->mergeUpperLower(new_parent_regions, best_child->length * new_split, lower_regions, best_child->length * (1 - new_split), aln, model, params->threshold_prob);
                
                RealNumType placement_cost = new_parent_regions->calculatePlacementCost(aln, model, cumulative_rate, sample, default_blength);
                
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
                
                new_split = best_split / 2;
                
                // delete new_parent_regions
                if (new_parent_regions) delete new_parent_regions;
            }
            
            best_child_lh = best_split_lh;
            best_child_split = best_split;
        }
        else
            best_child_lh = -DBL_MAX;
        
        // if node is root, try to place as sibling of the current root.
        RealNumType old_root_lh;
        if (root == selected_node)
        {
            old_root_lh = selected_node->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate)->computeAbsoluteLhAtRoot(aln, model, cumulative_base);
            RealNumType new_root_lh;
            Regions* merged_root_sample_regions = NULL;
            Regions* lower_regions = selected_node->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
            
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
                
                new_blength = best_root_blength / 2;
            }
            
            // delete merged_root_sample_regions
            if (merged_root_sample_regions) delete merged_root_sample_regions;
        }
        // selected_node is not root
        else
        {
            RealNumType best_split = 0.5;
            RealNumType best_split_lh = best_up_lh_diff;
            Regions* upper_left_right_regions = selected_node->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
            Regions* lower_regions = selected_node->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
            upper_left_right_regions->mergeUpperLower(best_parent_regions, selected_node->length / 2, lower_regions, selected_node->length / 2, aln, model, params->threshold_prob);
            RealNumType new_split = 0.25;
            
            Regions* new_parent_regions = NULL;
            while (new_split * selected_node->length > min_blength)
            {
                upper_left_right_regions->mergeUpperLower(new_parent_regions, selected_node->length * (1 - new_split), lower_regions, selected_node->length * new_split, aln, model, params->threshold_prob);
                
                RealNumType placement_cost = new_parent_regions->calculatePlacementCost(aln, model, cumulative_rate, sample, default_blength);
                
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
                
                new_split = best_split / 2;
            }
            // delete new_parent_regions
            if (new_parent_regions) delete new_parent_regions;
            
            best_parent_lh = best_split_lh;
            best_parent_split = best_split;
        }
        
        // if the best placement is below the selected_node => add an internal node below the selected_node
        if (best_child_lh >= best_parent_lh && best_child_lh >= best_lh_diff)
        {
            Regions* upper_left_right_regions = best_child->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);

            RealNumType new_branch_length_lh = best_child_lh;
            RealNumType best_length = default_blength;
            
            while (best_length > min_blength)
            {
                RealNumType new_blength = best_length / 2;
                RealNumType placement_cost = best_child_regions->calculatePlacementCost(aln, model, cumulative_rate, sample, new_blength);
                
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
                    RealNumType placement_cost = best_child_regions->calculatePlacementCost(aln, model, cumulative_rate, sample, new_blength);
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
                RealNumType tmp_lh = best_child_regions->calculatePlacementCost(aln, model, cumulative_rate, sample, -1);
                
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
            upper_left_right_regions->mergeUpperLower(new_internal_node->mid_branch_lh, new_internal_node->length / 2, new_internal_node->partial_lh, new_internal_node->length / 2, aln, model, params->threshold_prob);
            new_internal_node->computeTotalLhAtNode(aln, model, params->threshold_prob, cumulative_rate, new_internal_node == root);
            
            /*if (top_distance >= 2 * min_blengthForMidNode)
                createFurtherMidNodes(new_internal_node,this_node_upper_left_right_regions)*/
            if (best_length > 0)
            {
                new_sample_node->computeTotalLhAtNode(aln, model, params->threshold_prob, cumulative_rate, new_sample_node == root);
                next_node_1->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate)->mergeUpperLower(new_sample_node->mid_branch_lh, new_sample_node->length / 2, sample, new_sample_node->length / 2, aln, model, params->threshold_prob);
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
                if (!best_parent_regions) best_parent_regions = new Regions();
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
                Regions* new_root_lower_regions = NULL;
                RealNumType new_root_lh = 0;
                
                while (best_length2 > min_blength)
                {
                    RealNumType new_blength = best_length2 / 2;
                    
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
                new_root->total_lh = new_root->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate)->computeTotalLhAtRoot(num_states, model);

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
                    
                    new_root->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate)->mergeUpperLower(new_sample_node->mid_branch_lh, new_sample_node->length / 2, sample, new_sample_node->length / 2, aln, model, params->threshold_prob);
                    
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
                Regions* upper_left_right_regions = selected_node->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate);
                
                // now try different lengths for the new branch
                RealNumType new_branch_length_lh = best_parent_lh;
                RealNumType best_length = default_blength;
                while (best_length > min_blength)
                {
                    RealNumType new_blength = best_length / 2;
                    RealNumType placement_cost = best_parent_regions->calculatePlacementCost(aln, model, cumulative_rate, sample, new_blength);
                    
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
                        RealNumType placement_cost = best_parent_regions->calculatePlacementCost(aln, model, cumulative_rate, sample, new_blength);
                        
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
                    RealNumType placement_cost = best_parent_regions->calculatePlacementCost(aln, model, cumulative_rate, sample, -1);
                    
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
                upper_left_right_regions->mergeUpperLower(new_internal_node->mid_branch_lh, new_internal_node->length / 2, new_internal_node->partial_lh, new_internal_node->length / 2, aln, model, params->threshold_prob);
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
                    
                    next_node_1->getPartialLhAtNode(aln, model, params->threshold_prob, cumulative_rate)->mergeUpperLower(new_sample_node->mid_branch_lh, new_sample_node->length / 2, sample, new_sample_node->length / 2, aln, model, params->threshold_prob);
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
