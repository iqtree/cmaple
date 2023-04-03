#include "model_dna.h"
using namespace std;

ModelDNA::~ModelDNA()
{
    if (pseu_mutation_count)
    {
        delete[] pseu_mutation_count;
        pseu_mutation_count = NULL;
    }
}

void ModelDNA::updateMutMatbyMutCount()
{
    RealNumType* pseu_mutation_count_row = pseu_mutation_count;
    RealNumType* mutation_mat_row = mutation_mat;
    
    // init UNREST model
    if (model_name.compare("UNREST") == 0 || model_name.compare("unrest") == 0)
    {
        for (StateType i = 0; i <  num_states_; ++i, pseu_mutation_count_row += num_states_, mutation_mat_row += num_states_)
        {
            RealNumType sum_rate = 0;
            
            for (StateType j = 0; j <  num_states_; ++j)
            {
                if (i != j)
                {
                    RealNumType new_rate = pseu_mutation_count_row[j] * inverse_root_freqs[i];
                    mutation_mat_row[j] = new_rate;
                    sum_rate += new_rate;
                }
            }
            
            // update the diagonal entry
            mutation_mat_row[i] = -sum_rate;
            diagonal_mut_mat[i] = -sum_rate;
        }
    }
    // init GTR model
    else if (model_name.compare("GTR") == 0 || model_name.compare("gtr") == 0)
    {
        for (StateType i = 0; i <  num_states_; ++i, pseu_mutation_count_row += num_states_, mutation_mat_row += num_states_)
        {
            RealNumType sum_rate = 0;
            
            for (StateType j = 0; j <  num_states_; ++j)
                if (i != j)
                {
                    mutation_mat_row[j] = (pseu_mutation_count_row[j] + pseu_mutation_count[row_index[j] + i]) * inverse_root_freqs[i];
                    sum_rate += mutation_mat_row[j];
                }
            
            // update the diagonal entry
            mutation_mat_row[i] = -sum_rate;
            diagonal_mut_mat[i] = -sum_rate;
        }
    }
    // handle other model names
    else
        outError("Unsupported model! Please check and try again!");
}

void ModelDNA::updateMutationMat()
{
    // update Mutation matrix regarding the pseudo muation count
    updateMutMatbyMutCount();
    
    // compute the total rate regarding the root freqs
    RealNumType total_rate = 0;
    assert(num_states_ == 4);
    total_rate -= dotProduct<4>(root_freqs, diagonal_mut_mat);
    
    // inverse total_rate
    total_rate = 1.0 / total_rate;
    
    // normalize the mutation_mat
    RealNumType* mutation_mat_row = mutation_mat;
    RealNumType* freqi_freqj_qij_row = freqi_freqj_qij;
    for (StateType i = 0; i <  num_states_; ++i, mutation_mat_row += num_states_, freqi_freqj_qij_row += num_states_)
    {
        for (StateType j = 0; j <  num_states_; ++j)
        {
            mutation_mat_row[j] *= total_rate;
            
            // update freqi_freqj_qij
            if (i != j)
                freqi_freqj_qij_row[j] = root_freqs[i] * inverse_root_freqs[j] * mutation_mat_row[j];
            else
                freqi_freqj_qij_row[j] = mutation_mat_row[j];
            
            // update the transposed mutation matrix
            transposed_mut_mat[row_index[j] + i] = mutation_mat_row[j];
        }
        
        // update diagonal
        diagonal_mut_mat[i] = mutation_mat_row[i];
    }
    
    // pre-compute matrix to speedup
    RealNumType* transposed_mut_mat_row = transposed_mut_mat;
    RealNumType* freq_j_transposed_ij_row = freq_j_transposed_ij;
    
    assert(num_states_ == 4);
    for (StateType i = 0; i < num_states_; ++i, transposed_mut_mat_row += num_states_, freq_j_transposed_ij_row += num_states_)
        setVecByProduct<4>(freq_j_transposed_ij_row, root_freqs, transposed_mut_mat_row);
}

void ModelDNA::initMutationMatJC()
{
    // update root_freqs
    const RealNumType log_freq = log(0.25);
    
    for (StateType i = 0; i < num_states_; ++i)
    {
        root_freqs[i] = 0.25;
        inverse_root_freqs[i] = 4;
        root_log_freqs[i] = log_freq;
    }
    
    // compute some fixed value for JC model
    const RealNumType jc_rate = 1.0 / 3.0;
    const RealNumType freq_j_jc_rate = 0.25 * jc_rate;
    
    RealNumType* mutation_mat_row = mutation_mat;
    RealNumType* freqi_freqj_qij_row = freqi_freqj_qij;
    RealNumType* freq_j_transposed_ij_row = freq_j_transposed_ij;
    for (StateType i = 0; i < num_states_; ++i, mutation_mat_row += num_states_, freqi_freqj_qij_row += num_states_, freq_j_transposed_ij_row += num_states_)
    {
        for (StateType j = 0; j < num_states_; ++j)
        {
            if (i == j)
            {
                mutation_mat_row[j] = -1;
                transposed_mut_mat[row_index[j] + i] = -1;
                // update freqi_freqj_qij
                freqi_freqj_qij_row[j] = -1;
                // update freq_j_transposed_ij_row
                freq_j_transposed_ij_row[j] = -0.25; // 0.25 * -1 (freq(j) * transposed_ij)
            }
            else
            {
                mutation_mat_row[j] = jc_rate;
                transposed_mut_mat[row_index[j] + i] = jc_rate;
                // update freqi_freqj_qij
                freqi_freqj_qij_row[j] = root_freqs[i] * inverse_root_freqs[j] * jc_rate;
                // update freq_j_transposed_ij_row
                freq_j_transposed_ij_row[j] = freq_j_jc_rate; // 0.25 * 0.333 (freq(j) * transposed_ij)
            }
        }
        
        // update diagonal entries
        diagonal_mut_mat[i] = -1;
    }
}

void ModelDNA::initMutationMat()
{
    // init variable pointers
    initPointers();
    
    // for JC model
    if (model_name.compare("JC") == 0 || model_name.compare("jc") == 0)
        initMutationMatJC();
    // for other models
    else
    {
        // convert model name to uppercase
        transform(model_name.begin(), model_name.end(), model_name.begin(), ::toupper);
        if (model_name.compare("UNREST") != 0 && model_name.compare("GTR") != 0)
            outError("Invalid or unsupported DNA model: " + model_name + ". Please check and try again!");
            
        // init pseu_mutation_counts
        string model_rates = "0.0 1.0 5.0 2.0 2.0 0.0 1.0 40.0 5.0 2.0 0.0 20.0 2.0 3.0 1.0 0.0";
        convert_real_numbers(pseu_mutation_count, model_rates);
        
        updateMutationMat();
    }
}

void ModelDNA::updateMutationMatEmpirical(const Alignment& aln)
{
    // don't update JC model parameters
    if (model_name == "JC" || model_name == "jc") return;
    
    // clone the current mutation matrix
    RealNumType* tmp_diagonal_mut_mat = new RealNumType[num_states_];
    memcpy(tmp_diagonal_mut_mat, diagonal_mut_mat, num_states_ * sizeof(RealNumType));
    
    // update the mutation matrix regarding the pseu_mutation_count
    updateMutationMat();
    
    // update cumulative_rate if the mutation matrix changes more than a threshold
    RealNumType change_thresh = 1e-3;
    bool update = false;
    for (StateType j = 0; j < num_states_; ++j)
    {
        if (fabs(tmp_diagonal_mut_mat[j] - diagonal_mut_mat[j]) > change_thresh)
        {
            update = true;
            break;
        }
    }
    
    // update the cumulative_rate
    if (update)
        computeCumulativeRate(aln);
    
    // delete tmp_diagonal_mutation_mat
    delete[] tmp_diagonal_mut_mat;
}

void ModelDNA::updatePesudoCount(const Alignment& aln, const SeqRegions& regions1, const SeqRegions& regions2)
{
    if (model_name != "JC" && model_name != "jc")
    {
        // init variables
        PositionType pos = 0;
        const SeqRegions& seq1_regions = regions1;
        const SeqRegions& seq2_regions = regions2;
        size_t iseq1 = 0;
        size_t iseq2 = 0;
        const PositionType seq_length = aln.ref_seq.size();
                    
        while (pos < seq_length)
        {
            PositionType end_pos;
            
            // get the next shared segment in the two sequences
            SeqRegions::getNextSharedSegment(pos, seq1_regions, seq2_regions, iseq1, iseq2, end_pos);
            const auto* const seq1_region = &seq1_regions[iseq1];
            const auto* const seq2_region = &seq2_regions[iseq2];

            if (seq1_region->type != seq2_region->type && (seq1_region->type < num_states_ || seq1_region->type == TYPE_R) && (seq2_region->type < num_states_ || seq2_region->type == TYPE_R))
            {
                if (seq1_region->type == TYPE_R)
                    pseu_mutation_count[row_index[aln.ref_seq[end_pos]] + seq2_region->type] += 1;
                else if (seq2_region->type == TYPE_R)
                    pseu_mutation_count[row_index[seq1_region->type] + aln.ref_seq[end_pos]] += 1;
                else
                    pseu_mutation_count[row_index[seq1_region->type] + seq2_region->type] += 1;
            }

            // update pos
            pos = end_pos + 1;
        }
    }
}
