#include "modelbase.h"

using namespace std;
using namespace cmaple;

// explicit instantiation of templates
template void cmaple::ModelBase::updateMutationMat<4>();
template void cmaple::ModelBase::updateMutationMat<20>();
template void cmaple::ModelBase::updateMutationMatEmpiricalTemplate<4>(const Alignment*);
template void cmaple::ModelBase::updateMutationMatEmpiricalTemplate<20>(const Alignment*);

cmaple::ModelBase::ModelBase(const SubModel n_sub_model):sub_model(n_sub_model), mutation_mat(nullptr), diagonal_mut_mat(nullptr), transposed_mut_mat(nullptr), freqi_freqj_qij(nullptr),freq_j_transposed_ij(nullptr), root_freqs(nullptr), root_log_freqs(nullptr), inverse_root_freqs(nullptr), row_index(nullptr), model_block(nullptr), pseu_mutation_count(nullptr){}

cmaple::ModelBase::~ModelBase()
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

void cmaple::ModelBase::initEqualStateFreqs()
{
    // init variables
    if (!root_freqs) root_freqs = new RealNumType[num_states_];
    if (!root_log_freqs) root_log_freqs = new RealNumType[num_states_];
    if (!inverse_root_freqs) inverse_root_freqs = new RealNumType[num_states_];
    
    RealNumType state_freq = 1.0 / num_states_;
    RealNumType log_state_freq = log(state_freq);
    RealNumType inverse_state_freq = num_states_;
    for (auto i = 0; i < num_states_; ++i)
    {
        root_freqs[i] = state_freq;
        root_log_freqs[i] = log_state_freq;
        inverse_root_freqs[i] = inverse_state_freq;
    }
}

void cmaple::ModelBase::init()
{
    // init equal state_freqs
    initEqualStateFreqs();
    
    // init the mutation matrix
    initMutationMat();
}

void cmaple::ModelBase::extractRefInfo(const Alignment* aln)
{
    // init variables
    if (!root_freqs) root_freqs = new RealNumType[num_states_];
    if (!root_log_freqs) root_log_freqs = new RealNumType[num_states_];
    if (!inverse_root_freqs) inverse_root_freqs = new RealNumType[num_states_];
    
    // extract root freqs from the ref_seq
    extractRootFreqs(aln);
}

void cmaple::ModelBase::extractRootFreqs(const Alignment* aln)
{
    // init variables
    const vector<StateType>& ref_seq = aln->ref_seq;
    if(!ref_seq.size())
        throw std::logic_error("The reference genome is empty!");
    const PositionType seq_length = ref_seq.size();
    
    // init root_freqs
    switch (aln->getSeqType()) {
        case SEQ_PROTEIN:
            resetVec<20>(root_freqs);
            break;
            
        default: // dna
            resetVec<4>(root_freqs);
            break;
    }
    
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

void cmaple::ModelBase::computeCumulativeRate(const Alignment* aln)
{
    const PositionType sequence_length = aln->ref_seq.size();
    
    if (sequence_length <= 0)
        throw std::logic_error("Reference genome is empty");
    
    // init cumulative_rate
    if (!cumulative_rate)
        cumulative_rate = new RealNumType[sequence_length + 1];
    
    // init cumulative_base and cumulative_rate
    cumulative_base.resize(sequence_length + 1);
    cumulative_rate[0] = 0;
    cumulative_base[0].resize(num_states_, 0);
    
    // compute cumulative_base and cumulative_rate
    const std::vector<cmaple::StateType>& ref_seq = aln->ref_seq;
    for (PositionType i = 0; i < sequence_length; ++i)
    {
        StateType state = ref_seq[i];
        cumulative_rate[i + 1] = cumulative_rate[i] + diagonal_mut_mat[state];
        
        cumulative_base[i + 1] =  cumulative_base[i];
        cumulative_base[i + 1][state] = cumulative_base[i][state] + 1;
    }
}

std::string cmaple::ModelBase::exportRootFrequenciesStr()
{
    // Handle cases when root_freqs is null
    if (!root_freqs) return "\n \n";
    
    string output{};
    string header{};
    
    // Init a dummy alignment to get states (in readable character)
    Alignment aln;
    aln.setSeqType(getSeqType());
    
    for (StateType i = 0; i < num_states_; ++i)
    {
        header += aln.convertState2Char(i);
        header += "\t\t\t";
        output += convertDoubleToString(root_freqs[i]) + "\t";
    }
    
    return header + "\n" + output + "\n";
}

std::string cmaple::ModelBase::exportQMatrixStr()
{
    // Handle cases when mutation_mat is null
    if (!mutation_mat) return "\n";
    
    string output{};
    
    // Init a dummy alignment to get states (in readable character)
    Alignment aln;
    aln.setSeqType(getSeqType());
    
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

cmaple::ModelParams cmaple::ModelBase::exportModelParams()
{
    cmaple::ModelParams model_params;
    
    // model_name
    model_params.model_name = getModelName();
    
    // root frequencies
    model_params.state_freqs = exportRootFrequenciesStr();
    
    // Q matrix
    model_params.mut_rates = exportQMatrixStr();
    
    return model_params;
}

ModelsBlock* cmaple::ModelBase::readModelsDefinition(const char* builtin_models) {

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

bool cmaple::ModelBase::readParametersString(string& model_str) {

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
        throw std::logic_error(str);
    }
    return is_reversible;
}

void cmaple::ModelBase::readStateFreq(istream &in)
{
    StateType i;
    for (i = 0; i < num_states_; i++) {
        string tmp_value;
        in >> tmp_value;
        if (!tmp_value.length())
            throw "State frequencies could not be read";
        root_freqs[i] = convert_real_number(tmp_value.c_str());
        if (root_freqs[i] < 0.0)
            throw "Negative state frequencies found";
    }
    
    RealNumType sum = 0.0;
    for (i = 0; i < num_states_; i++) sum += root_freqs[i];
    if (fabs(sum-1.0) >= 1e-7)
    {
        if (cmaple::verbose_mode >= cmaple::VB_MED)
            outWarning("Normalizing state frequencies so that sum of them equals to 1");
        sum = 1.0/sum;
        for (i = 0; i < num_states_; i++)
            root_freqs[i] *= sum;
    }
}

void cmaple::ModelBase::normalizeQMatrix()
{
    ASSERT(root_freqs && mutation_mat);
    
    RealNumType sum = 0.0;
    RealNumType* mutation_mat_row = mutation_mat;
    for (StateType i = 0; i < num_states_; ++i, mutation_mat_row += num_states_)
        sum -= mutation_mat_row[i] * root_freqs[i];
    
    if (sum == 0.0) throw std::logic_error("Empty Q matrix");
    
    double delta = 1.0 / sum;
    
    mutation_mat_row = mutation_mat;
    for (StateType i = 0; i < num_states_; ++i, mutation_mat_row += num_states_)
        for (StateType j = 0; j < num_states_; ++j)
            mutation_mat_row[j] *= delta;
}

void cmaple::ModelBase::initPointers()
{
    // init row_index
    row_index = new StateType[num_states_ + 1];
    uint16_t start_index = 0;
    for (StateType i = 0; i < num_states_ + 1; i++, start_index += num_states_)
        row_index[i] = start_index;
    
    // init root_freqs, root_log_freqs, inverse_root_freqs if they have not yet been initialized
    if (!root_freqs) root_freqs = new RealNumType[num_states_];
    if (!root_log_freqs) root_log_freqs = new RealNumType[num_states_];
    if (!inverse_root_freqs) inverse_root_freqs = new RealNumType[num_states_];
    
    // init mutation_mat, transposed_mut_mat, and diagonal_mut_mat
    uint16_t mat_size = row_index[num_states_];
    mutation_mat = new RealNumType[mat_size];
    transposed_mut_mat = new RealNumType[mat_size];
    diagonal_mut_mat = new RealNumType[num_states_];
    freqi_freqj_qij = new RealNumType[mat_size];
    freq_j_transposed_ij = new RealNumType[mat_size];
}

void cmaple::ModelBase::updateMutMatbyMutCount()
{
    RealNumType* pseu_mutation_count_row = pseu_mutation_count;
    RealNumType* mutation_mat_row = mutation_mat;
    
    // init UNREST model
    if (sub_model == UNREST || sub_model == NONREV)
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
    else if (sub_model == GTR || sub_model == GTR20)
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
        throw std::logic_error("Unsupported model! Please check and try again!");
}

template <StateType num_states>
void cmaple::ModelBase::updateMutationMat()
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
void cmaple::ModelBase::updateMutationMatEmpiricalTemplate(const Alignment* aln)
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

void cmaple::ModelBase::updatePesudoCount(const Alignment* aln, const SeqRegions& regions1, const SeqRegions& regions2)
{
    // init variables
    PositionType pos = 0;
    const SeqRegions& seq1_regions = regions1;
    const SeqRegions& seq2_regions = regions2;
    size_t iseq1 = 0;
    size_t iseq2 = 0;
    const std::vector<cmaple::StateType>& ref_seq = aln->ref_seq;
    const PositionType seq_length = ref_seq.size();
                
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
                pseu_mutation_count[row_index[ref_seq[end_pos]] + seq2_region->type] += 1;
            else if (seq2_region->type == TYPE_R)
                pseu_mutation_count[row_index[seq1_region->type] + ref_seq[end_pos]] += 1;
            else
                pseu_mutation_count[row_index[seq1_region->type] + seq2_region->type] += 1;
        }

        // update pos
        pos = end_pos + 1;
    }
}

SeqType cmaple::ModelBase::getSeqType()
{
    switch (num_states_) {
        case 20:
            return SEQ_PROTEIN;
            break;
        case 4:
            return SEQ_DNA;
            break;
            
        default: // unkown
            return SEQ_UNKNOWN;
            break;
    }
}

SeqType cmaple::ModelBase::detectSeqType(const SubModel sub_model)
{
    // search in the list of dna models
    for(auto &it : dna_models_mapping)
        if(it.second == sub_model)
            return SEQ_DNA;
    
    // search in the list of protein models
    for(auto &it : aa_models_mapping)
        if(it.second == sub_model)
            return SEQ_PROTEIN;
    
    // not found
    return SEQ_UNKNOWN;
}
