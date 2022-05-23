//
//  sequence.cpp
//  alignment
//
//  Created by NhanLT on 31/3/2022.
//

#include "sequence.h"

Sequence::Sequence(string n_seq_name)
{
    seq_name = n_seq_name;
    mutations.resize(0);
}

Sequence::Sequence(string n_seq_name, vector<Mutation*> n_mutations)
{
    seq_name = n_seq_name;
    mutations = n_mutations;
}

Sequence::~Sequence()
{
    for (Mutation *mutation:mutations)
        delete mutation;
}
