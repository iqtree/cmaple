//
//  region.h
//  alignment
//
//  Created by Nhan Ly-Trong on 24/01/2022.
//

#include "mutation.h"

#ifndef REGION_H
#define REGION_H

class Region: public Mutation {
private:
    /**
    *  compute the relative likelihood for an ambiguious state
    */
    void computeLhAmbiguity(IntVector entries);
    
    /**
    *  Convert Ambiguious state into typical states: nucleotides/amino-acids...; N; O; R
    */
    void convertAmbiguiousState(SeqType seq_type, int max_num_states);
    
    /**
    *  Convert Ambiguious state (of DNA data) into typical states: A/C/G/T; N; O; R
    */
    void convertAmbiguiousStateDNA(int max_num_states);
    
public:
    // distance from root of the upward branch (~length of the path between the current phylo node and the node where the likelihood is calculated)
    double plength_upward;
    
    // distance from root of the downward branch
    double plength_downward;
    
    // the relative partial likelihood
    double* likelihood;
    
    /**
    *  Region constructor
    */
    Region();
    
    /**
    *  Region constructor
    */
    Region(StateType n_type, PositionType n_position, SeqType seq_type, int max_num_states);
    
    /**
    *  Region constructor
    */
    Region(Mutation* n_mutation, SeqType seq_type, int max_num_states);
    
    /**
    *  Region constructor
    */
    Region(Region* region, StateType num_states, bool copy_likelihood = true);
};
#endif
