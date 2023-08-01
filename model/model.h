
#pragma once

#include "modelbase.h"

namespace cmaple
{
    /** Class represents a substitution model */
    class Model
    {
    public:
        /*! \brief Constructor from a model name
         * @param[in] sub_model a substitution model. List of supported models:
         * <br>**DNA models**: JC, GTR, UNREST;
         * <br>**Protein models**: GTR20, NONREV, LG, WAG, JTT, Q_PFAM, Q_BIRD, Q_MAMMAL, Q_INSECT, Q_PLANT, Q_YEAST, JTTDCMUT, DCMUT, VT, PMB, BLOSUM62, DAYHOFF, MTREV, MTART, MTZOA, MTMET, MTVER, MTINV, MTMAM, FLAVI, HIVB, HIVW, FLU, RTREV, CPREV, NQ_PFAM, NQ_BIRD, NQ_MAMMAL, NQ_INSECT, NQ_PLANT, NQ_YEAST;
         * <br> <em> See [**Substitution models**](http://www.iqtree.org/doc/Substitution-Models) for references of those models.</em>
         * @param[in] seqtype Data type of sequences (optional): SEQ_DNA (nucleotide data), SEQ_PROTEIN (amino acid data), or SEQ_UNKNOWN (auto detection)
         * @throw std::invalid\_argument if sub\_model is unknown/unsupported
         */
        Model(const SubModel sub_model = GTR, const SeqType seqtype = SEQ_UNKNOWN);
        
        /*! \brief Destructor
         */
        ~Model();
        
        /*! \brief Export the substitution model and its parameters to ModelParams structure
         * @return A ModelParams structure contains the following members.
         * <br><em>model_name</em>: the name of the model in string
         * <br><em>state_freqs</em>: the state frequencies in string
         * <br><em>mut_rates</em>: the mutation matrix in string
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
