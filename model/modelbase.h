
#pragma once

#include <string>
#include "../utils/tools.h"
#include "../alignment/alignmentbase.h"
#include "../utils/matrix.h"
#include "../libraries/nclextra/modelsblock.h"
#include "../libraries/ncl/ncl.h"
#include "../libraries/nclextra/myreader.h"
#include <cassert>

namespace cmaple
{
    class SeqRegions;
    
    /**
     A structure to store model parameters
     */
    struct ModelParams
    {
        /**
         Name of the model in string
         */
        std::string model_name;
        
        /**
         State frequencies in string
         */
        std::string state_freqs;
        
        /**
         Mutation rates in string
         */
        std::string mut_rates;
    };

    /** Base class of evolutionary models */
    class ModelBase
    {
    private:        
        /**
         Initialize equal state frequencies
         */
        void initEqualStateFreqs();
        
        /**
         Get SeqType from num_states
         */
        SeqType getSeqType();
        
        /**
         Export state frequencies at root
         */
        std::string exportRootFrequenciesStr();
        
        /**
         Export Q matrix
         */
        std::string exportQMatrixStr();
        
        /**
         Read root state frequencies from string/file
         */
        void readStateFreq(istream &in);
        
        /**
         Update the mutation rate matrix regarding the pseu_mutation_count
         
         @throw  std::logic\_error if the substitution model is unknown/unsupported
         */
        void updateMutMatbyMutCount();
        
    protected:
        // NHANLT: we can change to use unique_ptr(s) instead of normal pointers
        /** Model definitions*/
        ModelsBlock* model_block;
        
        /**
         Initialize the model
         @throw std::logic\_error if the substitution model is unknown/unsupported
         */
        void init();
        
        /**
         Read model definitions from string/file
         */
        ModelsBlock* readModelsDefinition(const char* builtin_models);
        
        /**
         Read model parameters
         @return TRUE if the model is reversible
         @throw std::logic\_error if failing to read the parameters
         */
        bool readParametersString(string& model_str);
        
        /**
         Read model's rates from string/file
         */
        virtual void readRates(istream &in, const bool is_reversible) {};
        
        /**
         Normalize the Q matrix so that the expected number of subtitution is 1
         @throw std::logic\_error if the Q matrix is empty
         */
        void normalizeQMatrix();
        
        /**
         Extract root freqs from the reference sequence
         */
        virtual void extractRootFreqs(const AlignmentBase* aln);
        
        /**
         Init pointers
         */
        void initPointers();
        
        /**
         Update the mutation rate matrix
         @throw  std::logic\_error if the substitution model is unknown/unsupported
         */
        template <cmaple::StateType num_states>
        void updateMutationMat();
        
        /**
         Update the mutation matrix periodically from the empirical count of mutations (template)
         @throw  std::logic\_error if the substitution model is unknown/unsupported
         */
        template <cmaple::StateType num_states>
        void updateMutationMatEmpiricalTemplate(const AlignmentBase* aln);
        
        /**
         Get the model name
         */
        virtual std::string getModelName() const {return "";};
        
    public:
        // NHANLT: we can change to use unique_ptr(s) instead of normal pointers in the following
        /** Substitution model */
        SubModel sub_model;
        
        /** Number of states */
        cmaple::StateType num_states_;
        
        /** Pseudo mutation count */
        cmaple::RealNumType* pseu_mutation_count = nullptr;
        
        /** State frequencies*/
        cmaple::RealNumType *root_freqs = nullptr;
        
        /** Mutation matrix */
        cmaple::RealNumType *mutation_mat = nullptr;
        
        /** cumulative rates/bases*/
        cmaple::RealNumType *cumulative_rate = nullptr;
        std::vector< std::vector<cmaple::PositionType> > cumulative_base;
        
        /**
         Caches to reduce runtime
         */
        cmaple::RealNumType *root_log_freqs = nullptr; // log of state frequencies
        cmaple::RealNumType *diagonal_mut_mat; // diagonal of the mutation matrix
        cmaple::RealNumType *transposed_mut_mat; // the transposed matrix of the mutation matrix
        cmaple::RealNumType *inverse_root_freqs = nullptr; // the inversed values of state frequencies
        cmaple::RealNumType *freqi_freqj_qij; // freq(i) / freq(j) * Qij
        cmaple::RealNumType *freq_j_transposed_ij; // freq[j] * transposed[i][j]
        cmaple::StateType *row_index; // the starting index of row i: i * num_states
        
        /**
         Constructor
         */
        ModelBase() = default;
        
        /**
         Constructor
         */
        ModelBase(const SubModel sub_model);
        
        /**
         Destructor
         */
        ~ModelBase();
        
        /**
         Extract reference-related info (freqs, log_freqs)
         */
        void extractRefInfo(const AlignmentBase* aln);
        
        /**
         Init the mutation rate matrix from a model
         @throw std::logic\_error if the substitution model is unknown/unsupported
         */
        virtual void initMutationMat() {};
        
        /**
         Compute cumulative rate of the ref genome
         */
        void computeCumulativeRate(const AlignmentBase* aln);
        
        /**
         Update the mutation matrix periodically from the empirical count of mutations
         @throw std::logic\_error if the substitution model is unknown/unsupported
         */
        virtual void updateMutationMatEmpirical(const AlignmentBase* aln) {};
        
        /**
         Update pseudocounts from new sample to improve the estimate of the substitution rates
         @param node_regions the genome list at the node where the appending happens;
         @param sample_regions the genome list for the new sample.
         */
        virtual void updatePesudoCount(const AlignmentBase* aln, const SeqRegions& node_regions, const SeqRegions& sample_regions);
        
        /**
         Export model parameters to a dictionary
         */
        cmaple::ModelParams exportModelParams();
        
        /**
         Detect SeqType from a SubModel enum
         @param[in] sub_model SubModel enum
         */
        static SeqType detectSeqType(const SubModel sub_model)
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
    };
}
