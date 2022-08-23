#include "model.h"

Model::Model()
{
    model_name = "";
    mutation_mat = NULL;
    diagonal_mut_mat = NULL;
    transposed_mut_mat = NULL;
    freqi_freqj_qij = NULL;
    root_freqs = NULL;
    root_log_freqs = NULL;
    inverse_root_freqs = NULL;
    pseu_mutation_count = NULL;
    row_index = NULL;
}

Model::~Model()
{
    if (mutation_mat)
    {
        delete [] mutation_mat;
        mutation_mat = NULL;
    }
    
    if (diagonal_mut_mat)
    {
        delete [] diagonal_mut_mat;
        diagonal_mut_mat = NULL;
    }
    
    if (transposed_mut_mat)
    {
        delete [] transposed_mut_mat;
        transposed_mut_mat = NULL;
    }
    
    if (freqi_freqj_qij)
    {
        delete [] freqi_freqj_qij;
        freqi_freqj_qij = NULL;
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
    
    if (inverse_root_freqs)
    {
        delete [] inverse_root_freqs;
        inverse_root_freqs = NULL;
    }
    
    if (pseu_mutation_count)
    {
        delete[] pseu_mutation_count;
        pseu_mutation_count = NULL;
    }
    
    if (row_index)
    {
        delete[] row_index;
        row_index = NULL;
    }
}

void Model::extractRefInfo(vector<StateType> ref_seq, StateType num_states)
{
    ASSERT(ref_seq.size() > 0);
    
    // init variables
    PositionType seq_length = ref_seq.size();
    root_freqs = new RealNumType[num_states];
    root_log_freqs = new RealNumType[num_states];
    inverse_root_freqs = new RealNumType[num_states];
    
    // init root_freqs
    for (StateType i = 0; i < num_states; ++i)
        root_freqs[i] = 0;
    
    // browse all sites in the ref one by one to count bases
    for (PositionType i = 0; i < seq_length; ++i)
        ++root_freqs[ref_seq[i]];
    
    // update root_freqs and root_log_freqs
    RealNumType inverse_seq_length = 1.0 / seq_length;
    for (StateType i = 0; i < num_states; ++i)
    {
        root_freqs[i] *= inverse_seq_length;
        inverse_root_freqs[i] = 1.0 / root_freqs[i];
        root_log_freqs[i] = log(root_freqs[i]);
    }
}

void Model::updateMutationMat(StateType num_states)
{
    /*// init the "zero" mutation matrix
    string model_rates = "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
    convert_real_numbers(mutation_mat, model_rates);*/
    
    // init mutation_mat, diagonal_mut_mat, and transposed_mut_mat
    StateType matrix_size = row_index[num_states];
    if (!mutation_mat) mutation_mat = new RealNumType[matrix_size];
    if (!transposed_mut_mat) transposed_mut_mat = new RealNumType[matrix_size];
    if (!freqi_freqj_qij) freqi_freqj_qij = new RealNumType[matrix_size];
    if (!diagonal_mut_mat) diagonal_mut_mat = new RealNumType[num_states];
    
    // init UNREST model
    if (model_name.compare("UNREST") == 0 || model_name.compare("unrest") == 0)
    {
        for (StateType i = 0; i <  num_states; ++i)
        {
            RealNumType sum_rate = 0;
            
            for (StateType j = 0; j <  num_states; ++j)
            {
                if (i != j)
                {
                    StateType index = row_index[i] + j;
                    RealNumType new_rate = pseu_mutation_count[index] * inverse_root_freqs[i];
                    mutation_mat[index] = new_rate;
                    sum_rate += new_rate;
                }
            }
            
            // update the diagonal entry
            mutation_mat[row_index[i] + i] = -sum_rate;
            diagonal_mut_mat[i] = -sum_rate;
        }
    }
    // init GTR model
    else if (model_name.compare("GTR") == 0 || model_name.compare("gtr") == 0)
    {
        for (StateType i = 0; i <  num_states; ++i)
        {
            RealNumType sum_rate = 0;
            
            for (StateType j = 0; j <  num_states; ++j)
                if (i != j)
                {
                    StateType index = row_index[i] + j;
                    mutation_mat[index] = (pseu_mutation_count[index] + pseu_mutation_count[row_index[j] + i]) * inverse_root_freqs[i];
                    sum_rate += mutation_mat[index];
                }
            
            // update the diagonal entry
            mutation_mat[row_index[i] + i] = -sum_rate;
            diagonal_mut_mat[i] = -sum_rate;
        }
    }
    // handle other model names
    else
        outError("Unsupported model! Please check and try again!");
    
    // compute the total rate regarding the root freqs
    RealNumType total_rate = 0;
    for (StateType i = 0; i < num_states; ++i)
        total_rate -= root_freqs[i] * diagonal_mut_mat[i];
    // inverse total_rate
    total_rate = 1.0 / total_rate;
    
    // normalize the mutation_mat
    for (StateType i = 0; i <  num_states; ++i)
    {
        for (StateType j = 0; j <  num_states; ++j)
        {
            StateType index = row_index[i] + j;
            
            mutation_mat[index] *= total_rate;
            
            // update freqi_freqj_qij
            if (i != j)
                freqi_freqj_qij[index] = root_freqs[i] * inverse_root_freqs[j] * mutation_mat[index];
            else
                freqi_freqj_qij[index] = mutation_mat[index];
            
            // update the transposed mutation matrix
            transposed_mut_mat[row_index[j] + i] = mutation_mat[index];
        }
        
        // update diagonal
        diagonal_mut_mat[i] = mutation_mat[row_index[i] + i];
    }
}

void Model::initMutationMat(string n_model_name, StateType num_states)
{
    model_name = n_model_name;
    
    // init row_index
    row_index = new StateType[num_states + 1];
    StateType start_index = 0;
    for (StateType i = 0; i < num_states + 1; i++)
    {
        row_index[i] = start_index;
        start_index += num_states;
    }
        
    if (model_name.compare("JC") == 0 || model_name.compare("jc") == 0)
    {
        // update root_freqs
        RealNumType log_freq = log(0.25);
        for (StateType i = 0; i < num_states; ++i)
        {
            root_freqs[i] = 0.25;
            inverse_root_freqs[i] = 4;
            root_log_freqs[i] = log_freq;
        }
        
        // init mutation_mat, transposed_mut_mat, and diagonal_mut_mat
        StateType mat_size = row_index[num_states];
        mutation_mat = new RealNumType[mat_size];
        transposed_mut_mat = new RealNumType[mat_size];
        diagonal_mut_mat = new RealNumType[num_states];
        RealNumType jc_rate = 1.0 / 3.0;
        
        for (StateType i = 0; i < num_states; ++i)
        {
            for (StateType j = 0; j < num_states; ++j)
            {
                StateType index = row_index[i] + j;
                
                if (i == j)
                {
                    mutation_mat[index] = -1;
                    transposed_mut_mat[index] = -1;
                    // update freqi_freqj_qij
                    freqi_freqj_qij[index] = -1;
                }
                else
                {
                    mutation_mat[index] = jc_rate;
                    transposed_mut_mat[row_index[j] + i] = jc_rate;
                    // update freqi_freqj_qij
                    freqi_freqj_qij[index] = root_freqs[i] * inverse_root_freqs[j] * jc_rate;
                }
            }
            
            // update diagonal entries
            diagonal_mut_mat[i] = -1;
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
    for (PositionType i = 0; i < sequence_length; ++i)
    {
        StateType state = aln->ref_seq[i];
        cumulative_rate[i + 1] = cumulative_rate[i] + diagonal_mut_mat[state];
        
        cumulative_base[i + 1] =  cumulative_base[i];
        cumulative_base[i + 1][state] = cumulative_base[i][state] + 1;
    }
        
}

void Model::updateMutationMatEmpirical(RealNumType *&cumulative_rate, vector< vector<PositionType> > &cumulative_base, Alignment* aln)
{
    // don't update JC model parameters
    if (model_name == "JC" || model_name == "jc") return;
    
    StateType num_states = aln->num_states;
    
    // clone the current mutation matrix
    RealNumType* tmp_diagonal_mut_mat = new RealNumType[num_states];
    memcpy(tmp_diagonal_mut_mat, diagonal_mut_mat, num_states * sizeof(RealNumType));
    
    // update the mutation matrix regarding the pseu_mutation_count
    updateMutationMat(num_states);
    
    // update cumulative_rate if the mutation matrix changes more than a threshold
    RealNumType change_thresh = 1e-3;
    bool update = false;
    for (StateType j = 0; j < num_states; ++j)
    {
        if (fabs(tmp_diagonal_mut_mat[j] - diagonal_mut_mat[j]) > change_thresh)
        {
            update = true;
            break;
        }
    }
    
    // update the cumulative_rate
    if (update)
        computeCumulativeRate(cumulative_rate, cumulative_base, aln);
    
    // delete tmp_diagonal_mutation_mat
    delete[] tmp_diagonal_mut_mat;
}

void Model::updatePesudoCount(Alignment* aln, SeqRegions* regions1, SeqRegions* regions2)
{
    if (model_name != "JC" && model_name != "jc")
    {
        // init variables
        PositionType seq1_index = -1;
        PositionType seq2_index = -1;
        PositionType pos = 0;
        StateType num_states = aln->num_states;
        SeqRegion *seq1_region, *seq2_region;
        PositionType seq1_end = -1, seq2_end = -1;
        PositionType length;
        PositionType seq_length = aln->ref_seq.size();
                    
        while (pos < seq_length)
        {
            // get the next shared segment in the two sequences
            SeqRegions::getNextSharedSegment(pos, seq_length, regions1, regions2, seq1_index, seq2_index, seq1_region, seq2_region, seq1_end, seq2_end, length);
        
            if (seq1_region->type != seq2_region->type && (seq1_region->type < num_states || seq1_region->type == TYPE_R) && (seq2_region->type < num_states || seq2_region->type == TYPE_R))
            {
                if (seq1_region->type == TYPE_R)
                    pseu_mutation_count[row_index[aln->ref_seq[pos]] + seq2_region->type] += 1;
                else if (seq2_region->type == TYPE_R)
                    pseu_mutation_count[row_index[seq1_region->type] + aln->ref_seq[pos]] += 1;
                else
                    pseu_mutation_count[row_index[seq1_region->type] + seq2_region->type] += 1;
            }

            // update pos
            pos += length;
        }
    }
}
