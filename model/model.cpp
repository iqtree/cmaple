#include "model.h"

using namespace std;
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
