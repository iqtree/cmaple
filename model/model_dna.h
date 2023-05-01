
#pragma once

#include "model.h"

/** Class of DNA evolutionary models */
class ModelDNA: public Model
{
private:
    
    /**
        Init the mutation rate matrix from JC model
     */
    void initMutationMatJC();
    
public:    
	/**
		Constructor
	*/
    ModelDNA(const string n_model_name):Model(n_model_name) { num_states_ = 4; };
    
    /**
        Init the mutation rate matrix from a model
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
