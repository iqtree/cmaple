#include "model.h"

Model::Model()
{
    model_name = "";
    mutation_mat = NULL;
    root_freqs = NULL;
    root_log_freqs = NULL;
    pseu_mutation_count = NULL;
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
    
    if (pseu_mutation_count)
    {
        delete[] pseu_mutation_count;
        pseu_mutation_count = NULL;
    }
}

void Model::extractRefInfo(vector<StateType> ref_seq, StateType num_states)
{
    ASSERT(ref_seq.size() > 0);
    
    // init variables
    PositionType seq_length = ref_seq.size();
    root_freqs = new RealNumType[num_states];
    root_log_freqs = new RealNumType[num_states];
    
    // init root_freqs
    for (StateType i = 0; i < num_states; i++)
        root_freqs[i] = 0;
    
    // browse all sites in the ref one by one to count bases
    for (PositionType i = 0; i < seq_length; i++)
        root_freqs[ref_seq[i]]++;
    
    // update root_freqs and root_log_freqs
    for (StateType i = 0; i < num_states; i++)
    {
        root_freqs[i] /= seq_length;
        root_log_freqs[i] = log(root_freqs[i]);
    }
}

void Model::updateMutationMat(StateType num_states)
{
    // init the "zero" mutation matrix
    string model_rates = "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
    convert_real_numbers(mutation_mat, model_rates);
    
    // init UNREST model
    if (model_name.compare("UNREST") == 0 || model_name.compare("unrest") == 0)
    {
        for (int i = 0; i <  num_states; i++)
        {
            int start_index = i * num_states;
            RealNumType sum_rate = 0;
            RealNumType inverse_root_freq = 1 / root_freqs[i];
            
            for (int j = 0; j <  num_states; j++)
                if (i != j)
                {
                    mutation_mat[start_index + j] = pseu_mutation_count[start_index + j] * inverse_root_freq;
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
            RealNumType sum_rate = 0;
            RealNumType inverse_root_freq = 1 / root_freqs[i];
            
            for (int j = 0; j <  num_states; j++)
                if (i != j)
                {
                    mutation_mat[start_index + j] = (pseu_mutation_count[start_index + j] + pseu_mutation_count[j * num_states + i]) * inverse_root_freq;
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
    RealNumType total_rate = 0;
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

void Model::initMutationMat(string n_model_name, StateType num_states)
{
    model_name = n_model_name;
    if (model_name.compare("JC") == 0 || model_name.compare("jc") == 0)
    {
        // init mutation_mat
        string model_rates = "0.0 0.333 0.333 0.333 0.333 0.0 0.333 0.333 0.333 0.333 0.0 0.333 0.333 0.333 0.333 0.0";
        convert_real_numbers(mutation_mat, model_rates);
        
        // update root_freqs
        RealNumType log_freq = log(0.25);
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
        convert_real_numbers(pseu_mutation_count, model_rates);
        
        updateMutationMat(num_states);
    }
}

void Model::computeCumulativeRate(RealNumType *&cumulative_rate, vector< vector<PositionType> > &cumulative_base, Alignment* aln)
{
    PositionType sequence_length = aln->ref_seq.size();
    ASSERT(sequence_length > 0);
    
    // init cumulative_rate
    if (!cumulative_rate)
        cumulative_rate = new RealNumType[sequence_length + 1];
    
    // init cumulative_base and cumulative_rate
    cumulative_base.resize(sequence_length + 1);
    cumulative_rate[0] = 0;
    cumulative_base[0].resize(aln->num_states, 0);
    
    // compute cumulative_base and cumulative_rate
    for (PositionType i = 0; i < sequence_length; i++)
    {
        StateType state = aln->ref_seq[i];
        cumulative_rate[i + 1] = cumulative_rate[i] + mutation_mat[state * (aln->num_states + 1)];
        
        cumulative_base[i + 1] =  cumulative_base[i];
        cumulative_base[i + 1][state] = cumulative_base[i][state] + 1;
    }
        
}

void Model::updateMutationMatEmpirical(RealNumType *&cumulative_rate, vector< vector<PositionType> > &cumulative_base, Alignment* aln)
{
    StateType num_states = aln->num_states;
    
    // backup the current mutation matrix
    RealNumType* tmp_mutation_mat = new RealNumType[num_states * num_states];
    memcpy(tmp_mutation_mat, mutation_mat, num_states * num_states * sizeof(RealNumType));
    
    // update the mutation matrix regarding the pseu_mutation_count
    updateMutationMat(num_states);
    
    // update cumulative_rate if the mutation matrix changes more than a threshold
    RealNumType change_thresh = 1e-3;
    bool update = false;
    for (StateType j = 0; j < num_states; j++)
    {
        StateType index = j * (num_states + 1);
        if (fabs(tmp_mutation_mat[index] - mutation_mat[index]) > change_thresh)
        {
            update = true;
            break;
        }
    }
    
    // update the cumulative_rate
    if (update)
        computeCumulativeRate(cumulative_rate, cumulative_base, aln);
    
    // delete tmp_mutation_mat
    delete[] tmp_mutation_mat;
}

void Model::updatePesudoCount(Alignment* aln, Regions* regions1, Regions* regions2)
{
    if (model_name != "JC")
    {
        // init variables
        PositionType seq1_index = -1;
        PositionType seq2_index = -1;
        PositionType pos = 0;
        StateType num_states = aln->num_states;
        Region *seq1_region, *seq2_region;
        PositionType seq1_end = -1, seq2_end = -1;
        PositionType length;
        PositionType seq_length = aln->ref_seq.size();
                    
        while (pos < seq_length)
        {
            // get the next shared segment in the two sequences
            Regions::getNextSharedSegment(pos, seq_length, regions1, regions2, seq1_index, seq2_index, seq1_region, seq2_region, seq1_end, seq2_end, length);
        
            if (seq1_region->type != seq2_region->type && (seq1_region->type < num_states || seq1_region->type == TYPE_R) && (seq2_region->type < num_states || seq2_region->type == TYPE_R))
            {
                if (seq1_region->type == TYPE_R)
                    pseu_mutation_count[aln->ref_seq[pos] * num_states + seq2_region->type] += 1;
                else if (seq2_region->type == TYPE_R)
                    pseu_mutation_count[seq1_region->type * num_states + aln->ref_seq[pos]] += 1;
                else
                    pseu_mutation_count[seq1_region->type * num_states + seq2_region->type] += 1;
            }

            // update pos
            pos += length;
        }
    }
}
