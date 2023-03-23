
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
    
protected:
    
    /**
        Read root state frequencies from string/file
     */
    inline virtual void readStateFreq(istream& ist) {
        readStateFreqWithNumStates(ist, 4);
    };
    
public:
    // NHANLT: we can change to use unique_ptr(s) instead of normal pointers
    /** Pseudo mutation count */
    RealNumType* pseu_mutation_count;
    
	/**
		Constructor
	*/
    ModelDNA(const string n_model_name):Model(n_model_name), pseu_mutation_count(nullptr) {};
    
    /**
        Destructor
    */
    ~ModelDNA();
    
    /**
        Init the mutation rate matrix from a model
        @param n_model_name: name of the model; num_states: the number of states
     */
    virtual void initMutationMat();
    
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
