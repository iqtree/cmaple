#include <string>
#include "utils/tools.h"
#include "alignment/alignment.h"

using namespace std;

#ifndef MODEL_H
#define MODEL_H

class Model
{

public:
    string model_name;
    RealNumType *root_freqs;
    RealNumType *mutation_mat;
    RealNumType* pseu_mutation_count;
    
    /**
     caches to reduce runtime
     */
    RealNumType *root_log_freqs;
    RealNumType *diagonal_mut_mat;
    RealNumType *transposed_mut_mat;
    RealNumType *inverse_root_freqs;
    RealNumType *freqi_freqj_qij; // freq(i) / freq(j) * Qij
    StateType *row_index; // the starting index of row i: i * num_states
    
	/**
		constructor
	*/
    Model();
    
    /**
        deconstructor
    */
    ~Model();
    
    /**
        Extract reference-related info (freqs, log_freqs)
     */
    void extractRefInfo(vector<StateType> ref_seq, StateType num_states);
    
    /**
        Init the mutation rate matrix from a model
        @param n_model_name: name of the model; num_states: the number of states
     */
    void initMutationMat(string n_model_name, StateType num_states);
    
    /**
        Update the mutation rate matrix regarding the pseu_mutation_count
        @param num_states: the number of states
     */
    void updateMutationMat(StateType num_states);
    
    /**
        compute cumulative rate of the ref genome
     */
    void computeCumulativeRate(RealNumType *&cumulative_rate, vector< vector<PositionType> > &cumulative_base, Alignment* aln);
    
    /**
        Update the mutation matrix periodically from the empirical count of mutations
     */
    void updateMutationMatEmpirical(RealNumType *&cumulative_rate, vector< vector<PositionType> > &cumulative_base, Alignment* aln);
    
    /**
        update pseudocounts from new sample to improve the estimate of the substitution rates
        @param node_regions the genome list at the node where the appending happens;
        @param sample_regions the genome list for the new sample.
     */
    void updatePesudoCount(Alignment* aln, SeqRegions* node_regions, SeqRegions* sample_regions);
};
#endif
