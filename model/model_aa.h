
#pragma once

#include "model.h"

extern const char* builtin_prot_models;

/** Class of AA evolutionary models */
class ModelAA: public Model
{
private:
    
    /**
        Read model's rates from string/file
    */
    virtual void readRates(istream &in, const bool is_reversible);
    
    /**
        Read root state frequencies from string/file
     */
    inline virtual void readStateFreq(istream& ist) {
        readStateFreqWithNumStates(ist, 20);
    };
    
    /**
        Rescale the lower diagonal rates
     */
    void rescaleLowerDiagonalRates();
    
    /**
        Rescale all rates
     */
    void rescaleAllRates();
    
public:
    
    /**
        Constructor
    */
    ModelAA(const std::string n_model_name);
    
    /**
        Destructor
    */
    ~ModelAA();
    
    /// Move CTor
    ModelAA(ModelAA&& model) noexcept
    {
        model_name = std::move(model.model_name);
        
        model_block = model.model_block;
        model.model_block = nullptr;
    };
    
    /**
        Init the mutation rate matrix from a model
        @param n_model_name: name of the model; num_states: the number of states
     */
    virtual void initMutationMat();
};
