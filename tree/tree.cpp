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

vector<Region*> Tree::getOverallLhSeqAtRoot(vector<Region*> input_regions, double blength, double threshold_prob4)
{
    // init dummy variables
    vector<Region*> regions;
    StateType num_states = aln->num_states;
    
    // process regions one by one
    for (Region* input_region:input_regions)
    {
        Region* region = new Region(input_region, aln->num_states, false);
        if (region->type == TYPE_O)
        {
            // delete and re-initialize region->likelihood
            reinitDoubleArr(region->likelihood, num_states);
            
            // compute new likelihoods
            if (region->plength_upward > threshold_prob4)
            {
                for (StateType i = 0; i < num_states; i++)
                {
                    double total_lh = 0;
                    for (StateType j = 0; j < num_states; j++)
                    {
                        if (i == j)
                            total_lh += model->root_freqs[i] * (1.0 + model->mutation_mat[i * num_states + j] * (region->plength_upward + blength)) * input_region->likelihood[j];
                        else
                            total_lh += model->root_freqs[i] * model->mutation_mat[i * num_states + j] * (region->plength_upward + blength) * input_region->likelihood[j];
                    }
                    region->likelihood[i] = total_lh;
                }
            }
            else
            {
                for (StateType i = 0; i < num_states; i++)
                    region->likelihood[i] = model->root_freqs[i] * input_region->likelihood[i];
            }
                
            // set the distance from root of the upward branch at 0
            region->plength_upward = 0;
        }
        else
        {
            region->likelihood = NULL;
            if (region->type != TYPE_N)
            {
                // update the distance from root of the upward branch
                region->plength_upward += blength;
                
                // set the distance from root on the downward branch at 0
                region->plength_downward = 0;
            }
        }
        
        // add the new region into regions
        regions.push_back(region);
    }
    
    // return the overall likelihood sequence
    return regions;
}

vector<Region*> Tree::getLowerLhSeqAtTip(vector<Mutation*> mutations, double blength)
{
    PositionType sequence_length = aln->ref_seq.size();
    PositionType pos = 0;
    vector<Region*> regions;
    
    for (Mutation* mutation: mutations)
    {
        // insert Region of type R (if necessary)
        if (mutation->position > pos)
            regions.push_back(new Region(TYPE_R, pos, aln->seq_type, aln->num_states));
        
        // convert the current mutation
        pos = mutation->position + mutation->getLength();
        regions.push_back(new Region(mutation, aln->seq_type, aln->num_states));
    }
    
    // insert the last Region of type R (if necessary)
    if (pos < sequence_length)
        regions.push_back(new Region(TYPE_R, pos, aln->seq_type, aln->num_states));
    
    
    return regions;
}

void Tree::addNode(Node* new_node, Node* current_node, int &new_id, bool insert_to_right)
{
    Node* new_node_1 = new Node(new_id++);
    Node* new_node_2 = new Node(new_id++);
    Node* new_node_3 = new Node(new_id++);
    new_node_1->relative = new_node_2;
    new_node_2->relative = new_node_3;
    new_node_3->relative = new_node_1;
    
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
