
#pragma once

#include <string>
#include "utils/tools.h"
#include "alignment/alignment.h"
#include "utils/matrix.h"
#include "nclextra/modelsblock.h"
#include "ncl/ncl.h"
#include "nclextra/myreader.h"
#include <cassert>

class SeqRegions;

/** Class of evolutionary models */
class Model
{
private:
    
    /**
        Export state frequencies at root
     */
    std::string exportRootFrequenciesStr(Alignment& aln);
    
    /**
        Export Q matrix
     */
    std::string exportQMatrixStr(Alignment& aln);
    
protected:
    // NHANLT: we can change to use unique_ptr(s) instead of normal pointers
    /** Model definitions*/
    ModelsBlock* model_block;
    
    /**
        Read model definitions from string/file
     */
    ModelsBlock* readModelsDefinition(const char* builtin_models);
    
    /**
        Read model parameters
        return TRUE if the model is reversible
     */
    bool readParametersString(string& model_str);
    
    /**
        Read model's rates from string/file
    */
    virtual void readRates(istream &in, const bool is_reversible) {};
    
    /**
        Read root state frequencies from string/file
     */
    void readStateFreqWithNumStates(istream &in, const StateType num_states);
    
    /**
        Read root state frequencies from string/file
     */
    virtual void readStateFreq(istream &in) {};
    
    /**
        Normalize the Q matrix so that the expected number of subtitution is 1
     */
    void normalizeQMatrix(const StateType num_states);
    
public:
    // NHANLT: we can change to use unique_ptr(s) instead of normal pointers in the following
    /** Name of the model */
    std::string model_name;
    
    /** State frequencies*/
    RealNumType *root_freqs = nullptr;
    
    /** Mutation matrix */
    RealNumType *mutation_mat;
    
    /** cumulative rates/bases*/
    RealNumType *cumulative_rate = nullptr;
    std::vector< std::vector<PositionType> > cumulative_base;
    
    /**
        Caches to reduce runtime
     */
    RealNumType *root_log_freqs = nullptr; // log of state frequencies
    RealNumType *diagonal_mut_mat; // diagonal of the mutation matrix
    RealNumType *transposed_mut_mat; // the transposed matrix of the mutation matrix
    RealNumType *inverse_root_freqs = nullptr; // the inversed values of state frequencies
    RealNumType *freqi_freqj_qij; // freq(i) / freq(j) * Qij
    RealNumType *freq_j_transposed_ij; // freq[j] * transposed[i][j]
    StateType *row_index; // the starting index of row i: i * num_states
    
    /**
        Constructor
    */
    Model() = default;
    
    /**
		Constructor
	*/
    Model(const std::string n_model_name);
    
    /**
        Destructor
    */
    ~Model();
    
    /**
        Extract reference-related info (freqs, log_freqs)
     */
    void extractRefInfo(const Alignment& aln);
    
    /**
        Init the mutation rate matrix from a model
        @param n_model_name: name of the model; num_states: the number of states
     */
    virtual void initMutationMat() {};
    
    /**
        Compute cumulative rate of the ref genome
     */
    void computeCumulativeRate(const Alignment& aln);
    
    /**
        Update the mutation matrix periodically from the empirical count of mutations
     */
    virtual void updateMutationMatEmpirical(const Alignment& aln) {};
    
    /**
        Update pseudocounts from new sample to improve the estimate of the substitution rates
        @param node_regions the genome list at the node where the appending happens;
        @param sample_regions the genome list for the new sample.
     */
    virtual void updatePesudoCount(const Alignment& aln, const SeqRegions& node_regions, const SeqRegions& sample_regions) {};
    
    /**
        Export model parameters to string
     */
    std::string exportString(Alignment& aln);
};
