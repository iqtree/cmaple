#include "model.h"

Model::Model()
{
    model_name = "";
    mutation_mat = NULL;
    root_freqs = NULL;
    root_log_freqs = NULL;
}

Model::~Model()
{
    if (mutation_mat)
    {
        delete [] mutation_mat;
        mutation_mat = NULL;
    }
    
    if (root_freqs)
    {
        delete [] root_freqs;
        root_freqs = NULL;
    }
    
    if (root_log_freqs)
    {
        delete [] root_log_freqs;
        root_log_freqs = NULL;
    }
}

void Model::extractRefInfo(vector<StateType> ref_seq, StateType num_states)
{
    ASSERT(ref_seq.size() > 0);
    
    // init variables
    PositionType seq_length = ref_seq.size();
    root_freqs = new double[num_states];
    root_log_freqs = new double[num_states];
    
    // init root_freqs
    for (StateType i = 0; i < num_states; i++)
        root_freqs[i] = 0;
    
    // browse all sites in the ref one by one to count bases
    for (PositionType i = 0; i < seq_length; i++)
        root_freqs[ref_seq[i]]++;
    
    // update root_freqs and root_log_freqs
    double inverse_seq_length = 1.0/seq_length;
    for (StateType i = 0; i < num_states; i++)
    {
        root_freqs[i] *= inverse_seq_length;
        root_log_freqs[i] = log(root_freqs[i]);
    }
}

void Model::updateMutationMat(double* pseu_mutation_count, StateType num_states)
{
    // init the "zero" mutation matrix
    string model_rates = "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
    convert_doubles(mutation_mat, model_rates);
    
    // init UNREST model
    if (model_name.compare("UNREST") == 0 || model_name.compare("unrest") == 0)
    {
        for (int i = 0; i <  num_states; i++)
        {
            int start_index = i * num_states;
            double inverse_roota_freq = 1.0 / root_freqs[i];
            double sum_rate = 0;
            
            for (int j = 0; j <  num_states; j++)
                if (i != j)
                {
                    mutation_mat[start_index + j] = pseu_mutation_count[start_index + j] * inverse_roota_freq;
                    sum_rate += mutation_mat[start_index + j];
                }
            
            // update the diagonal entry
            mutation_mat[start_index + i] = -sum_rate;
        }
    }
    // init GTR model
    else if (model_name.compare("GTR") == 0 || model_name.compare("gtr") == 0)
    {
        for (int i = 0; i <  num_states; i++)
        {
            int start_index = i * num_states;
            double sum_rate = 0;
            
            for (int j = 0; j <  num_states; j++)
                if (i != j)
                {
                    mutation_mat[start_index + j] = (pseu_mutation_count[start_index + j] + pseu_mutation_count[j * num_states + i]) / root_freqs[i];
                    sum_rate += mutation_mat[start_index + j];
                }
            
            // update the diagonal entry
            mutation_mat[start_index + i] = -sum_rate;
        }
    }
    // handle other model names
    else
        outError("Unsupported model! Please check and try again!");
    
    // compute the total rate regarding the root freqs
    double total_rate = 0;
    for (int i = 0; i < num_states; i++)
        total_rate -= root_freqs[i] * mutation_mat[i * (num_states + 1)];
    
    // normalize the mutation_mat
    for (int i = 0; i <  num_states; i++)
    {
        int start_index = i * num_states;
        
        for (int j = 0; j <  num_states; j++)
            mutation_mat[start_index + j] /= total_rate;
    }
}

void Model::initMutationMat(string n_model_name, StateType num_states, double* &pseu_mutation_count)
{
    model_name = n_model_name;
    if (model_name.compare("JC") == 0 || model_name.compare("jc") == 0)
    {
        // init mutation_mat
        string model_rates = "0.0 0.333 0.333 0.333 0.333 0.0 0.333 0.333 0.333 0.333 0.0 0.333 0.333 0.333 0.333 0.0";
        convert_doubles(mutation_mat, model_rates);
        
        // update root_freqs
        double log_freq = log(0.25);
        for (int i = 0; i < num_states; i++)
        {
            root_freqs[i] = 0.25;
            root_log_freqs[i] = log_freq;
        }
    }
    else
    {
        // init pseu_mutation_counts
        string model_rates = "0.0 1.0 5.0 2.0 2.0 0.0 1.0 40.0 5.0 2.0 0.0 20.0 2.0 3.0 1.0 0.0";
        convert_doubles(pseu_mutation_count, model_rates);
        
        updateMutationMat(pseu_mutation_count, num_states);
    }
}
