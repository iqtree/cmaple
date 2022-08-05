//
//  mutation.h
//  alignment
//
//  Created by Nhan Ly-Trong on 24/01/2022.
//

#include "utils/tools.h"

#ifndef MUTATION_H
#define MUTATION_H

/** A mutation compared with a reference sequence */
class Mutation {
public:
    /** alternative allele, for DNA, it is A, C, G, T, N, O, - */
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
    *  Mutation deconstructor
    */
    virtual ~Mutation();
    
    /**
    *  return Length of the current mutation, default: 1
    */
    virtual PositionType getLength() { return 1;};
};
#endif
