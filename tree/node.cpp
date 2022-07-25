#include "node.h"

Node::Node(bool is_top_node)
{
    id = -1;
    seq_name = "";
    length = 0;
    if (is_top_node)
        double_attributes[IS_TOP_NODE] = 1;
    next = NULL;
    neighbor = NULL;
    partial_lh = NULL;
    total_lh = NULL;
    dirty = false;
}

Node::Node(int n_id, string n_seq_name)
{
    id = n_id;
    seq_name = n_seq_name;
    length = 0;
    next = NULL;
    neighbor = NULL;
    partial_lh = NULL;
    total_lh = NULL;
    dirty = false;
}

Node::Node(string n_seq_name)
{
    id = -1;
    seq_name = n_seq_name;
    length = 0;
    double_attributes[IS_TOP_NODE] = 1;
    next = NULL;
    neighbor = NULL;
    partial_lh = NULL;
    total_lh = NULL;
    dirty = false;
}

Node::~Node()
{
    // delete other neighbors
    // NHANLT: TODO
}

bool Node::isLeave()
{
    return next == NULL;
}

Node* Node::getTopNode()
{
    Node* next_node;
    Node* node = this;
    
    if (node->double_attributes.find(IS_TOP_NODE) != node->double_attributes.end())
        return node;
    
    FOR_NEXT(node, next_node)
    {
        if (next_node->double_attributes.find(IS_TOP_NODE) != next_node->double_attributes.end())
            return next_node;
    }
    
    return NULL;
}

string Node::exportString()
{
    if (isLeave())
    {
        string length_str = length < 0 ? "0" : convertDoubleToString(length);
        // without minor sequences -> simply return node's name and its branch length
        if (less_info_seqs.size() == 0)
            return seq_name + ":" + length_str;
        // with minor sequences -> return minor sequences' names with zero branch lengths
        else
        {
            string output = "(" + seq_name + ":0";
            for (string minor_seq_name : less_info_seqs)
                output += "," + minor_seq_name + ":0";
            output += "):" + length_str;
            
            return output;
        }
    }
    
    return "";
}

Regions* Node::getPartialLhAtNode(Alignment* aln, Model* model, double threshold_prob, double* cumulative_rate)
{
    // if partial_lh has not yet computed (~NULL) -> compute it from next nodes
    if (!partial_lh)
    {
        // init partial_lh
        partial_lh = new Regions();
        
        // the phylonode is an internal node
        if (next)
        {
            // if node is a top node -> partial_lh is the lower lh regions
            if (double_attributes.find(IS_TOP_NODE) != double_attributes.end())
            {
                // extract the two lower vectors of regions
                Node* next_node_1 = next;
                Regions* regions1 = next_node_1->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
                Node* next_node_2 = next_node_1->next;
                Regions* regions2 = next_node_2->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
                
                // compute partial_lh
                regions1->mergeTwoLowers(partial_lh, next_node_1->length, regions2, next_node_2->length, aln, model, threshold_prob, cumulative_rate);
            }
            // otherwise -> partial_lh is the upper left/right regions
            else
            {
                // extract the upper and the lower vectors of regions
                Regions* upper_regions, *lower_regions;
                double upper_blength, lower_blength;
                
                Node* next_node_1 = next;
                Node* next_node_2 = next_node_1->next;
                
                if (next_node_1->double_attributes.find(IS_TOP_NODE) != next_node_1->double_attributes.end())
                {
                    upper_regions = next_node_1->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
                    upper_blength = next_node_1->length;
                    lower_regions = next_node_2->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
                    lower_blength = next_node_2->length;
                }
                else
                {
                    upper_regions = next_node_2->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
                    upper_blength = next_node_2->length;
                    lower_regions = next_node_1->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
                    lower_blength = next_node_1->length;
                }
                
                // compute partial_lh
                upper_regions->mergeUpperLower(partial_lh, upper_blength, lower_regions, lower_blength, aln, model, threshold_prob);
            }
        }
        // the phylonode is a tip, partial_lh must be already computed
        else
            outError("Something went wrong! Lower likelihood regions has not been computed at tip!");
    }
    
    // return
    return partial_lh;
}

void Node::updateZeroBlength(stack<Node*> &node_stack, Alignment* aln, Model* model, double threshold_prob, double* cumulative_rate, double default_blength, double min_blength, double max_blength)
{
    // get the top node in the phylo-node
    Node* top_node = getTopNode();
    ASSERT(top_node);
    Regions* upper_left_right_regions = top_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
    Regions* lower_regions = top_node->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
    
    double best_lh = upper_left_right_regions->calculatePlacementCost(aln, model, cumulative_rate, lower_regions, default_blength);
    double best_length = default_blength;
    
    while (best_length > min_blength)
    {
        double new_blength = best_length/2;
        double new_lh = upper_left_right_regions->calculatePlacementCost(aln, model, cumulative_rate, lower_regions, new_blength);
        
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
            double new_blength = best_length * 2;
            double new_lh = upper_left_right_regions->calculatePlacementCost(aln, model, cumulative_rate, lower_regions, new_blength);
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
    top_node->dirty = true;
    top_node->neighbor->dirty = true;
    node_stack.push(top_node);
    node_stack.push(top_node->neighbor);
}

Regions* Node::computeTotalLhAtNode(Alignment* aln, Model* model, double threshold_prob, double* cumulative_rate, bool is_root, bool update)
{
    Regions* new_regions;
    
    // if node is root
    if (is_root)
        new_regions = getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate)->computeTotalLhAtRoot(aln->num_states, model);
    // if not is normal nodes
    else
    {
        new_regions = new Regions();
        Regions* lower_regions = getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
        neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate)->mergeUpperLower(new_regions, length, lower_regions, -1, aln, model, threshold_prob);
    }
    
    // update if necessary
    if (update)
    {
        if (total_lh) delete total_lh;
        total_lh = new_regions;
    }
    
    return new_regions;
}
