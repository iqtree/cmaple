//
//  mutation.cpp
//  alignment
//
//  Created by NhanLT on 31/3/2022.
//

#include "mutation.h"
using namespace cmaple;

cmaple::Mutation::Mutation(StateType n_type, PositionType n_position)
    : type(n_type),
    position(n_position)
{
    // do nothing else
}

cmaple::Mutation::Mutation(StateType n_type, PositionType n_position, LengthTypeLarge n_length)
 : type(n_type),
   position(n_position),
   length_(n_length)
{
  // validate the data
  if (n_length > 1 && type != TYPE_N && type != TYPE_DEL && type != TYPE_R)
    outError("Invalid mutation. Only mutation type N, -, or R can have length greater than 1.");
  if (n_length > (std::numeric_limits<LengthType>::max)())
    outError("Invalid mutation. Length is larger than 2^15. Recompile with larger 'LengthType' at the cost of higher memory consumption.");
}

LengthType cmaple::Mutation::getLength() const
{
    return length_;
}


