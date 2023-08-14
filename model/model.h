
#pragma once

#include "modelbase.h"

namespace cmaple
{
    /** Class represents a substitution model */
    class Model
    {
    public:
        /*!
         A structure to store model parameters
         */
        struct ModelParams
        {
            /*!
             Name of the model in string
             */
            std::string model_name;
            
            /*!
             State frequencies in string
             */
            std::string state_freqs;
            
            /*!
             Mutation rates in string
             */
            std::string mut_rates;
        };
        
        /*! \brief Constructor from a model name
         * @param[in] sub_model a substitution model. Default: MODEL_AUTO - auto select GTR for DNA, and LG for Protein data. List of supported models:
         * <br>**DNA models**: JC, GTR, UNREST;
         * <br>**Protein models**: GTR20, NONREV, LG, WAG, JTT, Q_PFAM, Q_BIRD, Q_MAMMAL, Q_INSECT, Q_PLANT, Q_YEAST, JTTDCMUT, DCMUT, VT, PMB, BLOSUM62, DAYHOFF, MTREV, MTART, MTZOA, MTMET, MTVER, MTINV, MTMAM, FLAVI, HIVB, HIVW, FLU, RTREV, CPREV, NQ_PFAM, NQ_BIRD, NQ_MAMMAL, NQ_INSECT, NQ_PLANT, NQ_YEAST;
         * <br> <em> See [**Substitution models**](http://www.iqtree.org/doc/Substitution-Models) for references of those models.</em>
         * @param[in] seqtype Data type of sequences (optional): SEQ_DNA (nucleotide data), SEQ_PROTEIN (amino acid data), or SEQ_AUTO (auto detection)
         * @throw std::invalid\_argument if any of the following situations occur.
         * - sub\_model is unknown/unsupported
         * - both sub_model and seqtype are specified as AUTO
         */
        Model(const cmaple::ModelBase::SubModel sub_model = cmaple::ModelBase::MODEL_AUTO, const cmaple::SeqRegion::SeqType seqtype = cmaple::SeqRegion::SEQ_AUTO);
        
        /*! \brief Destructor
         */
        ~Model();
        
        /*! \brief Keep the model parameters unchanged. This API is useful when using one model for multiple trees - after the model parameters are estimated according to a tree, users can keep them unchanged, then use the model for other trees.
         * @param[in] fixed_model_params TRUE to keep the model parameters unchanged
         */
        void fixParameters(const bool& fixed_model_params);
        
        /*! \brief Export the substitution model and its parameters to ModelParams structure
         * @return A ModelParams structure.
         */
        ModelParams getParams();
        
        // TODO: allow users to specify model parameters
        
        // declare Tree as a friend class
        friend class Tree;
    
    private:
        /**
         A base instance of Model
         */
        ModelBase* model_base;
    };
}
