//
//  mutation.cpp
//  alignment
//
//  Created by NhanLT on 31/3/2022.
//

#include "mutation.h"

Mutation::Mutation()
{
    type = TYPE_N;
    position = 0;
}

Mutation::Mutation(StateType n_type, PositionType n_position)
{
    type = n_type;
    position = n_position;
}
