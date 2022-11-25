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
    /**
        Alternative allele, for DNA, it is A, C, G, T, N, O, -
     */
    StateType type = TYPE_N;
    
    /**
        (starting) position
     */
    PositionType position = 0;

    PositionType length = 1;

    /**
    *  Mutation constructor
    */
    Mutation() = default;
    
    /**
    *  Mutation constructor
    */
    Mutation(StateType n_type, PositionType n_position)
      : type(n_type),
      position(n_position)
    {}

    Mutation(StateType n_type, PositionType n_position, PositionType n_length)
     : type(n_type),
       position(n_position),
       length (n_length)
    {
      // validate the data
      if (type != TYPE_N && type != TYPE_DEL && type != TYPE_R)
        outError("Invalid mutation. Only mutation type N, -, or R can have length greater than 1.");
    }

    PositionType getLength() const { return length;}
};
#endif
