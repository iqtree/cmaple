#include "model.h"

using namespace std;
using namespace cmaple;

// explicit instantiation of templates
template void Model::updateMutationMat<4>();
template void Model::updateMutationMat<20>();
template void Model::updateMutationMatEmpiricalTemplate<4>(const Alignment&);
template void Model::updateMutationMatEmpiricalTemplate<20>(const Alignment&);

Model::Model(const std::string n_model_name)
{
    model_name = std::move(n_model_name);
    mutation_mat = NULL;
    diagonal_mut_mat = NULL;
    transposed_mut_mat = NULL;
    freqi_freqj_qij = NULL;
    freq_j_transposed_ij = NULL;
    root_freqs = NULL;
    root_log_freqs = NULL;
    inverse_root_freqs = NULL;
    row_index = NULL;
    model_block = NULL;
    pseu_mutation_count = NULL;
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
    
    if (freq_j_transposed_ij)
    {
        delete [] freq_j_transposed_ij;
        freq_j_transposed_ij = NULL;
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
    
    if (row_index)
    {
        delete[] row_index;
        row_index = NULL;
    }
    
    if (cumulative_rate)
    {
        delete[] cumulative_rate;
    }
    
    if (pseu_mutation_count)
    {
        delete[] pseu_mutation_count;
        pseu_mutation_count = NULL;
    }
}

void Model::extractRefInfo(const Alignment& aln)
{
    // init variables
    root_freqs = new RealNumType[num_states_];
    root_log_freqs = new RealNumType[num_states_];
    inverse_root_freqs = new RealNumType[num_states_];
    
    // init root_freqs
    switch (aln.getSeqType()) {
        case SEQ_PROTEIN:
            resetVec<20>(root_freqs);
            break;
            
        default: // dna
            resetVec<4>(root_freqs);
            break;
    }
    
    // extract root freqs from the ref_seq
    extractRootFreqs(aln);
}

void Model::extractRootFreqs(const Alignment& aln)
{
    // init variables
    const vector<StateType>& ref_seq = aln.ref_seq;
    ASSERT(ref_seq.size() > 0);
    const PositionType seq_length = ref_seq.size();
    
    // browse all sites in the ref one by one to count bases
    for (PositionType i = 0; i < seq_length; ++i)
        ++root_freqs[ref_seq[i]];
    
    // update root_freqs and root_log_freqs
    RealNumType inverse_seq_length = 1.0 / seq_length;
    for (StateType i = 0; i < num_states_; ++i)
    {
        root_freqs[i] *= inverse_seq_length;
        inverse_root_freqs[i] = 1.0 / root_freqs[i];
        root_log_freqs[i] = log(root_freqs[i]);
    }
}

void Model::computeCumulativeRate(const Alignment& aln)
{
    const PositionType sequence_length = aln.ref_seq.size();
    ASSERT(sequence_length > 0);
    
    // init cumulative_rate
    if (!cumulative_rate)
        cumulative_rate = new RealNumType[sequence_length + 1];
    
    // init cumulative_base and cumulative_rate
    cumulative_base.resize(sequence_length + 1);
    cumulative_rate[0] = 0;
    cumulative_base[0].resize(num_states_, 0);
    
    // compute cumulative_base and cumulative_rate
    for (PositionType i = 0; i < sequence_length; ++i)
    {
        StateType state = aln.ref_seq[i];
        cumulative_rate[i + 1] = cumulative_rate[i] + diagonal_mut_mat[state];
        
        cumulative_base[i + 1] =  cumulative_base[i];
        cumulative_base[i + 1][state] = cumulative_base[i][state] + 1;
    }
}

std::string Model::exportRootFrequenciesStr(Alignment& aln)
{
    string output{};
    string header{};
    
    for (StateType i = 0; i < num_states_; ++i)
    {
        header += aln.convertState2Char(i);
        header += "\t\t\t";
        output += convertDoubleToString(root_freqs[i]) + "\t";
    }
    
    return header + "\n" + output + "\n";
}

std::string Model::exportQMatrixStr(Alignment& aln)
{
    string output{};
    
    // generate header
    output += "\t";
    for (StateType i = 0; i < num_states_; ++i)
    {
        output += aln.convertState2Char(i);
        output += "\t\t\t";
    }
    output += "\n";
    
    RealNumType* mut_mat_row = mutation_mat;
    for (StateType i = 0; i < num_states_; ++i, mut_mat_row += num_states_)
    {
        output += aln.convertState2Char(i);
        output += "\t";
        
        for (StateType j = 0; j < num_states_; ++j)
            output += convertDoubleToString(mut_mat_row[j]) + "\t";
        
        output += "\n";
    }
    
    return output;
}

std::string Model::exportString(Alignment& aln)
{
    string output{};
    
    // root frequencies
    output += "\nROOT FREQUENCIES\n";
    output += exportRootFrequenciesStr(aln);
    
    // Q matrix
    output += "\nMUTATION MATRIX\n";
    output += exportQMatrixStr(aln);
    
    return output;
}

ModelsBlock* Model::readModelsDefinition(const char* builtin_models) {

    ModelsBlock *models_block = new ModelsBlock;

    /*try
    {
        // loading internal model definitions
        stringstream in(builtin_mixmodels_definition);
        ASSERT(in && "stringstream is OK");
        NxsReader nexus;
        nexus.Add(models_block);
        MyToken token(in);
        nexus.Execute(token);
    } catch (...) {
        ASSERT(0 && "predefined mixture models not initialized");
    }*/

    try
    {
        // loading internal protei model definitions
        stringstream in(builtin_models);
        ASSERT(in && "stringstream is OK");
        NxsReader nexus;
        nexus.Add(models_block);
        MyToken token(in);
        nexus.Execute(token);
    } catch (...) {
        ASSERT(0 && "predefined protein models not initialized");
    }

    /*if (params.model_def_file) {
        cout << "Reading model definition file " << params.model_def_file << " ... ";
        MyReader nexus(params.model_def_file);
        nexus.Add(models_block);
        MyToken token(nexus.inf);
        nexus.Execute(token);
        int num_model = 0, num_freq = 0;
        for (ModelsBlock::iterator it = models_block->begin(); it != models_block->end(); it++)
            if (it->second.flag & NM_FREQ) num_freq++; else num_model++;
        cout << num_model << " models and " << num_freq << " frequency vectors loaded" << endl;
    }*/
    return models_block;
}

bool Model::readParametersString(string& model_str) {

    // if detect if reading full matrix or half matrix by the first entry
    PositionType end_pos;
    RealNumType d = convert_real_number(model_str.c_str(), end_pos);
    const bool is_reversible = (d >= 0);
    try {
        stringstream in(model_str);
        readRates(in, is_reversible);
        readStateFreq(in);
    }
    catch (const char *str) {
        outError(str);
    }
    return is_reversible;
}

void Model::readStateFreq(istream &in)
{
    StateType i;
    for (i = 0; i < num_states_; i++) {
        string tmp_value;
        in >> tmp_value;
        if (tmp_value.length() == 0)
            throw "State frequencies could not be read";
        root_freqs[i] = convert_real_number(tmp_value.c_str());
        if (root_freqs[i] < 0.0)
            throw "Negative state frequencies found";
    }
    
    RealNumType sum = 0.0;
    for (i = 0; i < num_states_; i++) sum += root_freqs[i];
    if (fabs(sum-1.0) >= 1e-7)
    {
        outWarning("Normalizing state frequencies so that sum of them equals to 1");
        sum = 1.0/sum;
        for (i = 0; i < num_states_; i++)
            root_freqs[i] *= sum;
    }
}

void Model::normalizeQMatrix()
{
    ASSERT(root_freqs && mutation_mat);
    
    RealNumType sum = 0.0;
    RealNumType* mutation_mat_row = mutation_mat;
    for (StateType i = 0; i < num_states_; ++i, mutation_mat_row += num_states_)
        sum -= mutation_mat_row[i] * root_freqs[i];
    
    if (sum == 0.0) throw "Empty Q matrix";
    
    double delta = 1.0 / sum;
    
    mutation_mat_row = mutation_mat;
    for (StateType i = 0; i < num_states_; ++i, mutation_mat_row += num_states_)
        for (StateType j = 0; j < num_states_; ++j)
            mutation_mat_row[j] *= delta;
}

void Model::initPointers()
{
    // init row_index
    row_index = new StateType[num_states_ + 1];
    uint16_t start_index = 0;
    for (StateType i = 0; i < num_states_ + 1; i++, start_index += num_states_)
        row_index[i] = start_index;
    
    // init root_freqs, root_log_freqs, inverse_root_freqs if they have not yet been initialized
    if (!root_freqs)
    {
        root_freqs = new RealNumType[num_states_];
        root_log_freqs = new RealNumType[num_states_];
        inverse_root_freqs = new RealNumType[num_states_];
    }
    
    // init mutation_mat, transposed_mut_mat, and diagonal_mut_mat
    uint16_t mat_size = row_index[num_states_];
    mutation_mat = new RealNumType[mat_size];
    transposed_mut_mat = new RealNumType[mat_size];
    diagonal_mut_mat = new RealNumType[num_states_];
    freqi_freqj_qij = new RealNumType[mat_size];
    freq_j_transposed_ij = new RealNumType[mat_size];
}

void Model::updateMutMatbyMutCount()
{
    RealNumType* pseu_mutation_count_row = pseu_mutation_count;
    RealNumType* mutation_mat_row = mutation_mat;
    
    // init UNREST model
    if (model_name.compare("UNREST") == 0 || model_name.compare("unrest") == 0 || model_name.compare("NONREV") == 0 || model_name.compare("nonrev") == 0)
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

template <StateType num_states>
void Model::updateMutationMat()
{
    // update Mutation matrix regarding the pseudo muation count
    updateMutMatbyMutCount();
    
    // compute the total rate regarding the root freqs
    RealNumType total_rate = 0;
    total_rate -= dotProduct<num_states>(root_freqs, diagonal_mut_mat);
    
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
    
    for (StateType i = 0; i < num_states_; ++i, transposed_mut_mat_row += num_states_, freq_j_transposed_ij_row += num_states_)
        setVecByProduct<num_states>(freq_j_transposed_ij_row, root_freqs, transposed_mut_mat_row);
}

template <StateType num_states>
void Model::updateMutationMatEmpiricalTemplate(const Alignment& aln)
{
    // clone the current mutation matrix
    RealNumType* tmp_diagonal_mut_mat = new RealNumType[num_states_];
    memcpy(tmp_diagonal_mut_mat, diagonal_mut_mat, num_states_ * sizeof(RealNumType));
    
    // update the mutation matrix regarding the pseu_mutation_count
    updateMutationMat<num_states>();
    
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

void Model::updatePesudoCount(const Alignment& aln, const SeqRegions& regions1, const SeqRegions& regions2)
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
