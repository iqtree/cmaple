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
    LengthType length_ = 1;
public:
    /**
        Alternative allele, for DNA, it is A, C, G, T, N, O, -
     */
    StateType type = TYPE_N;

    /**
        (starting) position
     */
    PositionType position = 0;

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

    Mutation(StateType n_type, PositionType n_position, LengthTypeLarge n_length)
     : type(n_type),
       position(n_position),
       length_(n_length)
    {
      // validate the data
      if (type != TYPE_N && type != TYPE_DEL && type != TYPE_R)
        outError("Invalid mutation. Only mutation type N, -, or R can have length greater than 1.");
      if (n_length > (std::numeric_limits<LengthType>::max)())
        outError("Invalid mutation. Length is larger than 2^15. Recompile with larger 'LengthType' at the cost of higher memory consumption.");
    }

   LengthType getLength() const { return length_;}
};
#endif
