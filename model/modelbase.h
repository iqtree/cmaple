
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
        cmaple::SeqRegion::SeqType getSeqType();
        
        /**
         Read root state frequencies from string/file
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        void readStateFreq(istream &in);
        
        /**
         Update the mutation rate matrix regarding the pseu_mutation_count
         
         @throw  std::logic\_error if the substitution model is unknown/unsupported
         */
        void updateMutMatbyMutCount();
        
    protected:
        // NHANLT: we can change to use unique_ptr(s) instead of normal pointers
        /**
         @private
         Model definitions
         */
        ModelsBlock* model_block;
        
        /**
         @private
         Initialize the model
         @throw std::logic\_error if the substitution model is unknown/unsupported
         */
        void init();
        
        /**
         @private
         Read model definitions from string/file
         */
        ModelsBlock* readModelsDefinition(const char* builtin_models);
        
        /**
         @private
         Read model parameters
         @return TRUE if the model is reversible
         @throw std::logic\_error if failing to read the parameters
         */
        bool readParametersString(string& model_str);
        
        /**
         @private
         Read model's rates from string/file
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        virtual void readRates(istream &in, const bool is_reversible) {};
        
        /**
         @private
         Normalize the Q matrix so that the expected number of subtitution is 1
         @throw std::logic\_error if the Q matrix is empty
         */
        void normalizeQMatrix();
        
        /**
         @private
         Extract root freqs from the reference sequence
         @throw std::logic\_error if the reference genome is empty
         */
        virtual void extractRootFreqs(const Alignment* aln);
        
        /**
         @private
         Init pointers
         */
        void initPointers();
        
        /**
         @private
         Update the mutation rate matrix
         @throw  std::logic\_error if the substitution model is unknown/unsupported
         */
        template <cmaple::StateType num_states>
        void updateMutationMat();
        
        /**
         @private
         Update the mutation matrix periodically from the empirical count of mutations (template)
         @throw  std::logic\_error if the substitution model is unknown/unsupported
         */
        template <cmaple::StateType num_states>
        void updateMutationMatEmpiricalTemplate(const Alignment* aln);
        
    public:
        /*!
         * List of substitution models. See [Substitution models](http://www.iqtree.org/doc/Substitution-Models) for references of those models.
         */
        enum SubModel {
            ///@{
            /// See [DNA models](http://www.iqtree.org/doc/Substitution-Models#dna-model)
            JC,
            GTR,
            UNREST,
            ///@}
            
            ///@{
            /// See [Protein models](http://www.iqtree.org/doc/Substitution-Models#protein-models)
            GTR20, NONREV, LG, WAG, JTT, Q_PFAM, Q_BIRD, Q_MAMMAL, Q_INSECT, Q_PLANT, Q_YEAST, JTTDCMUT, DCMUT, VT, PMB, BLOSUM62, DAYHOFF, MTREV, MTART, MTZOA, MTMET, MTVER, MTINV, MTMAM, FLAVI, HIVB, HIVW, FLU, RTREV, CPREV, NQ_PFAM, NQ_BIRD, NQ_MAMMAL, NQ_INSECT, NQ_PLANT, NQ_YEAST,
            ///@}
            
            UNKNOWN_MODEL, ///<  Unknown model
        };
        
        /**
         * @private
         * mapping between DNA model names and their enums
         */
        const static std::map<std::string, SubModel> dna_models_mapping;
        
        /**
         * @private
         * mapping between Protein model names and their enums
         */
        const static std::map<std::string, SubModel> aa_models_mapping;
        
        // NHANLT: we can change to use unique_ptr(s) instead of normal pointers in the following
        /**
         @private
         Substitution model
         */
        SubModel sub_model;
        
        /**
         @private
         Number of states
         */
        cmaple::StateType num_states_;
        
        /**
         @private
         Pseudo mutation count
         */
        cmaple::RealNumType* pseu_mutation_count = nullptr;
        
        /**
         @private
         State frequencies
         */
        cmaple::RealNumType *root_freqs = nullptr;
        
        /**
         @private
         Mutation matrix
         */
        cmaple::RealNumType *mutation_mat = nullptr;
        
        /**
         @private
         cumulative rates
         */
        cmaple::RealNumType *cumulative_rate = nullptr;
        
        /**
         @private
         cumulative bases
         */
        std::vector< std::vector<cmaple::PositionType> > cumulative_base;
        
        /**
         @private
         log of state frequencies
         */
        cmaple::RealNumType *root_log_freqs = nullptr;
        
        /**
         @private
         diagonal of the mutation matrix
         */
        cmaple::RealNumType *diagonal_mut_mat;
        
        /**
         @private
         the transposed matrix of the mutation matrix
         */
        cmaple::RealNumType *transposed_mut_mat;
        
        /**
         @private
         the inversed values of state frequencies
         */
        cmaple::RealNumType *inverse_root_freqs = nullptr;
        
        /**
         @private
         freq(i) / freq(j) * Qij
         */
        cmaple::RealNumType *freqi_freqj_qij;
        
        /**
         @private
         freq[j] * transposed[i][j]
         */
        cmaple::RealNumType *freq_j_transposed_ij;
        
        /**
         @private
         the starting index of row i: i * num_states
         */
        cmaple::StateType *row_index;
        
        /**
         @private
         Constructor
         */
        ModelBase() = default;
        
        /**
         @private
         Constructor
         */
        ModelBase(const SubModel sub_model);
        
        /**
         @private
         Destructor
         */
        ~ModelBase();
        
        /**
         @private
         Get the model name
         */
        virtual std::string getModelName() const {return "";};
        
        /**
         @private
         Export state frequencies at root
         */
        std::string exportRootFrequenciesStr();
        
        /**
         @private
         Export Q matrix
         */
        std::string exportQMatrixStr();
        
        /**
         @private
         Extract reference-related info (freqs, log_freqs)
         
         @throw std::logic\_error if the reference genome is empty
         */
        void extractRefInfo(const Alignment* aln);
        
        /**
         @private
         Init the mutation rate matrix from a model
         @throw std::logic\_error if the substitution model is unknown/unsupported
         */
        virtual void initMutationMat() {};
        
        /**
         @private
         Compute cumulative rate of the ref genome
         @throw std::logic\_error if the reference genome is empty
         */
        void computeCumulativeRate(const Alignment* aln);
        
        /**
         @private
         Update the mutation matrix periodically from the empirical count of mutations
         @throw std::logic\_error if any of the following situations occur.
         - the substitution model is unknown/unsupported
         - the reference genome is empty
         */
        virtual void updateMutationMatEmpirical(const Alignment* aln) {};
        
        /**
         @private
         Update pseudocounts from new sample to improve the estimate of the substitution rates
         @param node_regions the genome list at the node where the appending happens;
         @param sample_regions the genome list for the new sample.
         */
        virtual void updatePesudoCount(const Alignment* aln, const SeqRegions& node_regions, const SeqRegions& sample_regions);
        
        /**
         @private
         Detect SeqType from a SubModel enum
         @param[in] sub_model SubModel enum
         */
        static cmaple::SeqRegion::SeqType detectSeqType(const SubModel sub_model);
        
        /**
        * @private
        * Parse model from its name in a string
        * @param n_seqtype_str a sequence type in string
        */
        static cmaple::ModelBase::SubModel parseModel(const std::string& model_name);
    };
}
