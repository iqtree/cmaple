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

SeqRegions* Sequence::getLowerLhVector(PositionType sequence_length, StateType num_states, SeqType seq_type)
{
    SeqRegions* regions = new SeqRegions();
    PositionType pos = 0;
    
    for (Mutation* mutation: mutations)
    {
        // insert Region of type R (if necessary)
        if (mutation->position > pos)
            regions->push_back(new SeqRegion(TYPE_R, pos));
        
        // convert the current mutation
        pos = mutation->position + mutation->getLength();
        regions->push_back(new SeqRegion(mutation, seq_type, num_states));
    }
    
    // insert the last Region of type R (if necessary)
    if (pos < sequence_length)
        regions->push_back(new SeqRegion(TYPE_R, pos));
    
    return regions;
}
