
#pragma once

#include "model.h"

/** Class of DNA evolutionary models */
class ModelDNA: public Model
{
private:
    
    /**
        Init the mutation rate matrix from JC model
        @param n_model_name: name of the model; num_states: the number of states
     */
    void initMutationMatJC(const StateType num_states);
    
    /**
        Update the mutation rate matrix regarding the pseu_mutation_count
        @param num_states: the number of states
     */
    void updateMutMatbyMutCount(const StateType num_states);
    
    /**
        Update the mutation rate matrix
        @param num_states: the number of states
     */
    void updateMutationMat(const StateType num_states);
    
public:
    /** Pseudo mutation count */
    RealNumType* pseu_mutation_count;
    
	/**
		Constructor
	*/
    ModelDNA();
    
    /**
        Destructor
    */
    ~ModelDNA();
    
    /**
        Init the mutation rate matrix from a model
        @param n_model_name: name of the model; num_states: the number of states
     */
    virtual void initMutationMat(const std::string n_model_name, const StateType num_states);
    
    /**
        Update the mutation matrix periodically from the empirical count of mutations
     */
    virtual void updateMutationMatEmpirical(const Alignment& aln);
    
    /**
        Update pseudocounts from new sample to improve the estimate of the substitution rates
        @param node_regions the genome list at the node where the appending happens;
        @param sample_regions the genome list for the new sample.
     */
    virtual void updatePesudoCount(const Alignment& aln, const SeqRegions& node_regions, const SeqRegions& sample_regions);
};
