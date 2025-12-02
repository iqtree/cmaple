//
//  mutation.cpp
//  alignment
//
//  Created by NhanLT on 31/3/2022.
//

#include "mutation.h"
using namespace cmaple;

cmaple::Mutation::Mutation(StateType n_type, PositionType n_position, cmaple::StateType n_prev_state)
    : type(n_type),
    position(n_position),
    prev_state(n_prev_state)
{
    // do nothing else
}

cmaple::Mutation::Mutation(StateType n_type, PositionType n_position, LengthTypeLarge n_length, cmaple::StateType n_prev_state)
 : type(n_type),
   position(n_position),
   length_(n_length),
   prev_state(n_prev_state)
{
  assert(n_length > 0);
  assert(n_position >= 0);
  
  // validate the data
  if (n_length > 1 && type != TYPE_N && type != TYPE_DEL && type != TYPE_R) {
    throw std::invalid_argument("Invalid mutation. Only mutation type N, -, or "
                                "R can have length greater than 1.");
  }
  if (n_length > (std::numeric_limits<LengthType>::max)()) {
    throw std::invalid_argument(
        "The sequence is longer than the maximum supported by the current build."
        " Recompile CMAPLE with a larger `LengthType` to handle longer sequences"
        " (will increase memory use).");
  }
}

auto cmaple::Mutation::getLength() const -> LengthType { return length_; }
