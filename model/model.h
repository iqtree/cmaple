
#pragma once

#include "modelbase.h"

namespace cmaple
{
    /** Class of evolutionary models */
    class Model
    {
    public:
        /*! \brief Model constructor
         *
         * Model constructor with a model name
         * @param[in] model_name Name of a substitution model
         * @param[in] seqtype Data type of sequences (optional): "", "DNA", "AA"
         */
        Model(const std::string& model_name = "GTR", const std::string& seqtype = "");
        
        /*! \brief Model destructor
         *
         * Model destructor
         */
        ~Model();
        
        /*! \brief Export the substitution model and its parameters in a dictionary
         *
         * Export the substitution model and its parameters in a dictionary
         * @return a dictionary (std::map<std::string, std::string>) with the keys and values as in the following.
         * key: "model_name", value: the name of the model in string
         * key: "model_freqs", value: the state frequencies in string
         * key: "model_rates", value: the mutation matrix in string
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
