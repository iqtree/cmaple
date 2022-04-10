//
//  mutationindel.cpp
//  alignment
//
//  Created by NhanLT on 31/3/2022.
//

#include "mutationindel.h"

MutationIndel::MutationIndel():Mutation()
{
    length = 0;
}


MutationIndel::MutationIndel(StateType n_type, PositionType n_position, PositionType n_length):Mutation(n_type, n_position)
{
    length = n_length;
    
    // validate the data
    if (type != TYPE_N && type != TYPE_DEL && type != TYPE_R)
        outError("Invalid mutation. Only mutation type N, -, or R can have length greater than 1.");
}
