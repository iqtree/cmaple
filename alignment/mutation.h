//
//  mutation.h
//  alignment
//
//  Created by Nhan Ly-Trong on 24/01/2022.
//

#include "utils/tools.h"

#ifndef MUTATION_H
#define MUTATION_H

class Mutation {
public:
    // type: A, C, G, T, N, O, -, R
    StateType type;
    // (starting) position
    PositionType position;
    
    /**
    *  Mutation constructor
    */
    Mutation();
    
    /**
    *  Mutation constructor
    */
    Mutation(StateType n_type, PositionType n_position);
    
    /**
    *  return Length of the current mutation, default: 1
    */
    virtual PositionType getLength() { return 1;};
};
#endif
