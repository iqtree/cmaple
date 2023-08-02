
#pragma once

#include "modelbase.h"

namespace cmaple
{
    extern const char* builtin_prot_models;

    /** Class of AA evolutionary models */
    class ModelAA: public ModelBase
    {
    private:
        
        /**
         Read model's rates from string/file
         */
        virtual void readRates(istream &in, const bool is_reversible);
        
        /**
         Rescale the lower diagonal rates
         */
        void rescaleLowerDiagonalRates();
        
        /**
         Rescale all rates
         */
        void rescaleAllRates();
        
        /**
         extract root freqs from the reference sequence
         for AA models, we directly get root_freqs from predefined state_freqs
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
        ModelAA(const SubModel sub_model);
        
        /**
         Destructor
         */
        ~ModelAA();
        
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
