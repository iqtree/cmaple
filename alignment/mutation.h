#include "../utils/tools.h"

#pragma once

namespace cmaple
{
    /** A mutation compared with a reference sequence */
    class Mutation {
        cmaple::LengthType length_ = 1;
    public:
        /**
         Alternative allele, for DNA, it is A, C, G, T, N, O, -
         */
        cmaple::StateType type = cmaple::TYPE_N;
        
        /**
         (starting) position
         */
        cmaple::PositionType position = 0;
        
        /**
         *  Mutation constructor
         */
        Mutation() = default;
        
        /**
         *  Mutation constructor
         */
        Mutation(cmaple::StateType n_type, cmaple::PositionType n_position);
        
        /**
         *  Mutation constructor
         *  @throw throw std::invalid\_argument any of the following situations occur.
         *  - n_type is invalid
         *  - n_length is greater than the maximum value of the LengthType
         */
        Mutation(cmaple::StateType n_type, cmaple::PositionType n_position, cmaple::LengthTypeLarge n_length);
        
        /**
         *  Return the length of the mutation
         */
        cmaple::LengthType getLength() const;
    };
}
