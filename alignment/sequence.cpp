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

vector<Region*> Sequence::getRegionsFromMutations(PositionType sequence_length, SeqType seq_type, int max_num_states)
{
    PositionType pos = 0;
    vector<Region*> regions;
    
    for (Mutation* mutation: mutations)
    {
        // insert Region of type R (if necessary)
        if (mutation->position > pos)
            regions.push_back(new Region(TYPE_R, pos, seq_type, max_num_states));
        
        // convert the current mutation
        pos = mutation->position + mutation->getLength();
        regions.push_back(new Region(mutation, seq_type, max_num_states));
    }
    
    // insert the last Region of type R (if necessary)
    if (pos < sequence_length)
        regions.push_back(new Region(TYPE_R, pos, seq_type, max_num_states));
    
    
    return regions;
}
