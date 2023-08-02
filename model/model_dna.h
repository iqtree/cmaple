
#pragma once

#include "modelbase.h"

namespace cmaple
{
    /** Class of DNA evolutionary models */
    class ModelDNA: public ModelBase
    {
    private:
        
        /**
         Init the mutation rate matrix from JC model
         */
        void initMutationMatJC();
        
        /**
         extract root freqs from the reference sequence
         */
        virtual void extractRootFreqs(const AlignmentBase* aln);
        
        /**
         Get the model name
         */
        virtual std::string getModelName() const;
        
    public:    
        /**
         Constructor
         @throw std::invalid\_argument if sub\_model is unknown/unsupported
         */
        ModelDNA(const SubModel sub_model);
        
        /**
         Init the mutation rate matrix from a model
         @throw std::logic\_error if the substitution model is unknown/unsupported
         */
        virtual void initMutationMat();
        
        /**
         Update the mutation matrix periodically from the empirical count of mutations
         @throw  std::logic\_error if any of the following situations occur.
         - the substitution model is unknown/unsupported
         - the reference genome is empty
         */
        virtual void updateMutationMatEmpirical(const AlignmentBase* aln);
        
        /**
         Update pseudocounts from new sample to improve the estimate of the substitution rates
         @param node_regions the genome list at the node where the appending happens;
         @param sample_regions the genome list for the new sample.
         */
        virtual void updatePesudoCount(const AlignmentBase* aln, const SeqRegions& node_regions, const SeqRegions& sample_regions);
    };
}
