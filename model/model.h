#include <string>
#include "utils/tools.h"
#include "alignment/alignment.h"
#include "utils/matrix.h"
#include <cassert>

#ifndef MODEL_H
#define MODEL_H

/** Class of evolutionary models */
class Model
{
public:
    /** Name of the model */
    std::string model_name;
    
    /** State frequencies*/
    RealNumType *root_freqs;
    
    /** Mutation matrix */
    RealNumType *mutation_mat;
    
    /** Pseudo mutation count */
    RealNumType* pseu_mutation_count;
    
    /**
        Caches to reduce runtime
     */
    RealNumType *root_log_freqs; // log of state frequencies
    RealNumType *diagonal_mut_mat; // diagonal of the mutation matrix
    RealNumType *transposed_mut_mat; // the transposed matrix of the mutation matrix
    RealNumType *inverse_root_freqs; // the inversed values of state frequencies
    RealNumType *freqi_freqj_qij; // freq(i) / freq(j) * Qij
    RealNumType *freq_j_transposed_ij; // freq[j] * transposed[i][j]
    StateType *row_index; // the starting index of row i: i * num_states
    
	/**
		Constructor
	*/
    Model();
    
    /**
        Destructor
    */
    ~Model();
    
    /**
        Extract reference-related info (freqs, log_freqs)
     */
    void extractRefInfo(std::vector<StateType> ref_seq, StateType num_states);
    
    /**
        Init the mutation rate matrix from a model
        @param n_model_name: name of the model; num_states: the number of states
     */
    void initMutationMat(std::string n_model_name, StateType num_states);
    
    /**
        Update the mutation rate matrix regarding the pseu_mutation_count
        @param num_states: the number of states
     */
    void updateMutationMat(StateType num_states);
    
    /**
        Compute cumulative rate of the ref genome
     */
    void computeCumulativeRate(RealNumType *&cumulative_rate, std::vector< std::vector<PositionType> > &cumulative_base, const Alignment& aln);
    
    /**
        Update the mutation matrix periodically from the empirical count of mutations
     */
    void updateMutationMatEmpirical(RealNumType *&cumulative_rate, std::vector< std::vector<PositionType> > &cumulative_base, const Alignment& aln);
    
    /**
        Update pseudocounts from new sample to improve the estimate of the substitution rates
        @param node_regions the genome list at the node where the appending happens;
        @param sample_regions the genome list for the new sample.
     */
    void updatePesudoCount(const Alignment& aln, const SeqRegions& node_regions, const SeqRegions& sample_regions);
};
#endif
