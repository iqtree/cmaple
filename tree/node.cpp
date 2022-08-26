#include "node.h"

Node::Node(bool is_top_node)
{
    id = -1;
    seq_name = "";
    length = 0;
    is_top = is_top_node;
    next = NULL;
    neighbor = NULL;
    partial_lh = NULL;
    mid_branch_lh = NULL;
    total_lh = NULL;
    outdated = false;
}

Node::Node(int n_id, string n_seq_name)
{
    id = n_id;
    seq_name = n_seq_name;
    length = 0;
    is_top = false;
    next = NULL;
    neighbor = NULL;
    partial_lh = NULL;
    mid_branch_lh = NULL;
    total_lh = NULL;
    outdated = false;
}

Node::Node(string n_seq_name)
{
    id = -1;
    seq_name = n_seq_name;
    length = 0;
    is_top = true;
    next = NULL;
    neighbor = NULL;
    partial_lh = NULL;
    mid_branch_lh = NULL;
    total_lh = NULL;
    outdated = false;
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
    
    if (node->is_top)
        return node;
    
    FOR_NEXT(node, next_node)
    {
        if (next_node->is_top)
            return next_node;
    }
    
    return NULL;
}

Node* Node::getOtherNextNode()
{
    Node* next_node;
    Node* node = this;
    
    FOR_NEXT(node, next_node)
    {
        if (!next_node->is_top)
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

SeqRegions* Node::getPartialLhAtNode(Alignment* aln, Model* model, RealNumType threshold_prob, RealNumType* cumulative_rate)
{
    // if partial_lh has not yet computed (~NULL) -> compute it from next nodes
    if (!partial_lh)
    {
        // init partial_lh
        partial_lh = NULL;
        
        // the phylonode is an internal node
        if (next)
        {
            // if node is a top node -> partial_lh is the lower lh regions
            if (is_top)
            {
                // extract the two lower vectors of regions
                Node* next_node_1 = next;
                SeqRegions* regions1 = next_node_1->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
                Node* next_node_2 = next_node_1->next;
                SeqRegions* regions2 = next_node_2->neighbor->getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
                
                // compute partial_lh
                regions1->mergeTwoLowers(partial_lh, next_node_1->length, regions2, next_node_2->length, aln, model, threshold_prob, cumulative_rate);
            }
            // otherwise -> partial_lh is the upper left/right regions
            else
            {
                // extract the upper and the lower vectors of regions
                SeqRegions* upper_regions, *lower_regions;
                RealNumType upper_blength, lower_blength;
                
                Node* next_node_1 = next;
                Node* next_node_2 = next_node_1->next;
                
                if (next_node_1->is_top)
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

SeqRegions* Node::computeTotalLhAtNode(Alignment* aln, Model* model, RealNumType threshold_prob, RealNumType* cumulative_rate, bool is_root, bool update)
{
    SeqRegions* new_regions = NULL;
    
    // if node is root
    if (is_root)
        new_regions = getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate)->computeTotalLhAtRoot(aln->num_states, model);
    // if not is normal nodes
    else
    {
        SeqRegions* lower_regions = getPartialLhAtNode(aln, model, threshold_prob, cumulative_rate);
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

TraversingNode::TraversingNode()
{
    node = NULL;
    failure_count = 0;
    likelihood_diff = 0;
}

TraversingNode::TraversingNode(Node* n_node, short int n_failure_count, RealNumType n_lh_diff)
{
    node = n_node;
    failure_count = n_failure_count;
    likelihood_diff = n_lh_diff;
}

TraversingNode::~TraversingNode()
{
    // do nothing
}

UpdatingNode::UpdatingNode():TraversingNode()
{
    incoming_regions = NULL;
    branch_length = 0;
    need_updating = false;
}

UpdatingNode::UpdatingNode(Node* n_node, SeqRegions* n_incoming_regions, RealNumType n_branch_length, bool n_need_updating, short int n_failure_count, RealNumType n_lh_diff):TraversingNode(n_node, n_failure_count, n_lh_diff)
{
    incoming_regions = n_incoming_regions;
    branch_length = n_branch_length;
    need_updating = n_need_updating;
}

UpdatingNode::~UpdatingNode()
{
    // do nothing
}
