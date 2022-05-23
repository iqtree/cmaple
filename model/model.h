#include <string>
#include "utils/tools.h"

using namespace std;

#ifndef MODEL_H
#define MODEL_H

class Model
{

public:
    string model_name;
    double *root_freqs, *root_log_freqs;
    double *mutation_mat;
    
	/**
		constructor
	*/
    Model();
    
    /**
        deconstructor
    */
    ~Model();
    
    /**
        Extract reference-related info (freqs, log_freqs)
     */
    void extractRefInfo(vector<StateType> ref_seq, StateType num_states);
    
    /**
        Init the mutation rate matrix from a model
        @param n_model_name: name of the model; num_states: the number of states
     */
    void initMutationMat(string n_model_name, StateType num_states, double* &pseu_mutation_count);
    
    /**
        Update the mutation rate matrix regarding the pseu_mutation_count
        @param pseu_mutation_count: pseudo mutation count; num_states: the number of states
     */
    void updateMutationMat(double* pseu_mutation_count, StateType num_states);
};
#endif
