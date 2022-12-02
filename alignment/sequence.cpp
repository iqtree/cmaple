//
//  sequence.cpp
//  alignment
//
//  Created by NhanLT on 31/3/2022.
//

#include "sequence.h"
using namespace std;
Sequence::Sequence(string n_seq_name)
{
    seq_name = std::move(n_seq_name);
    resize(0);
}

Sequence::Sequence(string n_seq_name, vector<Mutation*> n_mutations):vector<Mutation*>(std::move(n_mutations))
{
    seq_name = std::move(n_seq_name);
}

Sequence::~Sequence()
{
    for (Mutation *mutation : (*this))
        delete mutation;
}

SeqRegions* Sequence::getLowerLhVector(PositionType sequence_length, StateType num_states, SeqType seq_type)
{
    SeqRegions* regions = new SeqRegions();
    regions->reserve(this->size() * 2); // avoid realloc of vector data (alternatively, scan through vector first to determine final # of elements)
    PositionType pos = 0;
    
    for (Mutation* mutation: (*this))
    {
        // insert Region of type R (if necessary)
        if (mutation->position > pos)
            regions->emplace_back(TYPE_R, mutation->position - 1);
        
        // convert the current mutation
        pos = mutation->position + mutation->getLength();
        regions->emplace_back(mutation, seq_type, num_states);
    }
    
    // insert the last Region of type R (if necessary)
    if (pos < sequence_length)
        regions->emplace_back(TYPE_R, sequence_length - 1);
    
    return regions;
}
