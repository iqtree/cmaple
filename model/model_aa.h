
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
    virtual void extractRootFreqs(const Alignment& aln) {};
    
public:
    
    /**
        Constructor
    */
    ModelAA(const std::string n_model_name):Model(n_model_name) { num_states_ = 20; };
    
    /**
        Destructor
    */
    ~ModelAA();
    
    /**
        Init the mutation rate matrix from a model
     */
    virtual void initMutationMat();
};
