
#pragma once

#include "modelbase.h"

namespace cmaple
{
    /** Class represents a substitutionmodel */
    class Model
    {
    public:
        /*! \brief Constructor from a model name
         * @param[in] model_name Name of a substitution model
         * @param[in] seqtype Data type of sequences (optional): "", "DNA", "AA"
         */
        Model(const std::string& model_name = "GTR", const std::string& seqtype = "");
        
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
