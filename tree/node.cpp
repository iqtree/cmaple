#include "node.h"

Node::Node()
{
    id = -1;
    seq_name = "";
    length = 0;
    next = NULL;
    neighbor = NULL;
    partial_lh = NULL;
    total_lh = NULL;
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

void Node::computeTotalLhForRoot(Model* model, StateType num_states, double blength)
{
    // init total_lh
    if (total_lh)
        delete total_lh;
    total_lh = new Regions();
    
    for (Region* region : (*partial_lh))
    {
        // type N
        if (region->type == TYPE_N)
        {
            Region* new_region = new Region(region, num_states, false);
            total_lh->push_back(new_region);
        }
        else
            // type O
            if (region->type == TYPE_O)
            {
                double total_plength_observation = region->plength_observation + blength;
                double* new_partial_likelihood = new double[num_states];
                double sum_partial_likelihood = 0;
                
                for (StateType i = 0; i < num_states; i++)
                {
                    double tot = 0.0;
                    
                    if (total_plength_observation > 0)
                    {
                        for (StateType j = 0; j < num_states; j++)
                            tot += model->mutation_mat[i * num_states + j] * region->likelihood[j];
                        tot *= total_plength_observation;
                    }
                    
                    tot += region->likelihood[i];
                    new_partial_likelihood[i] = tot * model->root_freqs[i];
                    sum_partial_likelihood += new_partial_likelihood[i];
                }
                
                // normalize partial likelihood
                for (StateType i = 0; i < num_states; i++)
                    new_partial_likelihood[i] /= sum_partial_likelihood;
                
                // add new region to the total_lh_regions
                Region* new_region = new Region(region, num_states, false);
                new_region->likelihood = new_partial_likelihood;
                total_lh->push_back(new_region);
            }
            // other types: R or A/C/G/T
            else
            {
                // add new region to the total_lh_regions
                Region* new_region = new Region(region, num_states, false);
                new_region->plength_observation += blength;
                total_lh->push_back(new_region);
            }
    }
}
