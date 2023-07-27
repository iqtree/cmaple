
#pragma once

#include "modelbase.h"

namespace cmaple
{
    /** Class represents a substitution model */
    class Model
    {
    public:
        /*! \brief Constructor from a model name
         * @param[in] model_name Name of a substitution model. List of supported models:
         * <br>**DNA models**: "JC", "GTR", "UNREST";
         * <br>**Protein models**: "GTR20", "NONREV", "LG", "WAG", "JTT", "Q.PFAM", "Q.BIRD", "Q.MAMMAL", "Q.INSECT", "Q.PLANT", "Q.YEAST", "JTTDCMUT", "DCMUT", "VT", "PMB", "BLOSUM62", "DAYHOFF", "MTREV", "MTART", "MTZOA", "MTMET" , "MTVER" , "MTINV", "MTMAM", "FLAVI", "HIVB", "HIVW", "FLU", "RTREV", "CPREV", "NQ.PFAM", "NQ.BIRD", "NQ.MAMMAL", "NQ.INSECT", "NQ.PLANT", "NQ.YEAST";
         * <br> <em> See [**Substitution models**](http://www.iqtree.org/doc/Substitution-Models) for references of those models.</em>
         * @param[in] seqtype Data type of sequences (optional): "" (auto detect), "DNA" (nucleotide data), "AA" (amino acid data)
         * @param[in] seqtype Data type of sequences (optional): SEQ_DNA (nucleotide data), SEQ_PROTEIN (amino acid data), or SEQ_UNKNOWN (auto detection)
         */
        Model(const std::string& model_name = "GTR", const SeqType seqtype = SEQ_UNKNOWN);
        
        /*! \brief Destructor
         */
        ~Model();
        
        /*! \brief Export the substitution model and its parameters to a dictionary
         * @return A dictionary (std::map<std::string, std::string>) with the keys and values as in the following.
         * <br><em>"model_name"</em>: the name of the model in string
         * <br><em>"model_freqs"</em>: the state frequencies in string
         * <br><em>"model_rates"</em>: the mutation matrix in string
         */
        std::map<std::string, std::string> getParams();
        
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
