
#pragma once

#include <string>
#include "../utils/tools.h"
#include "../alignment/alignment.h"
#include "../utils/matrix.h"
#include "../libraries/nclextra/modelsblock.h"
#include "../libraries/ncl/ncl.h"
#include "../libraries/nclextra/myreader.h"
#include <cassert>

namespace cmaple
{
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
        
        /**
         Read root state frequencies from string/file
         */
        void readStateFreq(istream &in);
        
        /**
         Update the mutation rate matrix regarding the pseu_mutation_count
         */
        void updateMutMatbyMutCount();
        
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
         Normalize the Q matrix so that the expected number of subtitution is 1
         */
        void normalizeQMatrix();
        
        /**
         Extract root freqs from the reference sequence
         */
        virtual void extractRootFreqs(const Alignment& aln);
        
        /**
         Init pointers
         */
        void initPointers();
        
        /**
         Update the mutation rate matrix
         */
        template <cmaple::StateType num_states>
        void updateMutationMat();
        
        /**
         Update the mutation matrix periodically from the empirical count of mutations (template)
         */
        template <cmaple::StateType num_states>
        void updateMutationMatEmpiricalTemplate(const Alignment& aln);
        
    public:
        // NHANLT: we can change to use unique_ptr(s) instead of normal pointers in the following
        /** Name of the model */
        std::string model_name;
        
        /** Number of states */
        cmaple::StateType num_states_;
        
        /** Pseudo mutation count */
        cmaple::RealNumType* pseu_mutation_count = nullptr;
        
        /** State frequencies*/
        cmaple::RealNumType *root_freqs = nullptr;
        
        /** Mutation matrix */
        cmaple::RealNumType *mutation_mat;
        
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
        virtual void updatePesudoCount(const Alignment& aln, const SeqRegions& node_regions, const SeqRegions& sample_regions);
        
        /**
         Export model parameters to string
         */
        std::string exportString(Alignment& aln);
    };
}
