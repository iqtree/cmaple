//
//  mutation.h
//  alignment
//
//  Created by Nhan Ly-Trong on 24/01/2022.
//

#include "utils/tools.h"

#pragma once

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
    Mutation(StateType n_type, PositionType n_position);

    /**
    *  Mutation constructor
    */
    Mutation(StateType n_type, PositionType n_position, LengthTypeLarge n_length);

    /**
    *  Return the length of the mutation
    */
    LengthType getLength() const;
};
