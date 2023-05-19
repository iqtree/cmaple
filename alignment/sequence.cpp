//
//  sequence.cpp
//  alignment
//
//  Created by NhanLT on 31/3/2022.
//

#include "sequence.h"
using namespace std;
using namespace cmaple;

cmaple::Sequence::Sequence(string n_seq_name)
  : seq_name(std::move(n_seq_name))
{
}

cmaple::Sequence::Sequence(string n_seq_name, vector<Mutation> n_mutations):
  vector<Mutation>(std::move(n_mutations)),
  seq_name(std::move(n_seq_name))
{    
}

std::unique_ptr<SeqRegions> cmaple::Sequence::getLowerLhVector(const PositionType sequence_length, const StateType num_states, const SeqType seq_type)
{
    std::unique_ptr<SeqRegions> regions = make_unique<SeqRegions>(SeqRegions());
    regions->reserve(size() * 2); // avoid realloc of vector data (we could also count explicitly - see below), but that seems a bit slower (yet potentially more memory efficient)

    // count number of items we will need
    /**
    size_t count{this->size()};
    PositionType pos = 0;
    for (auto& mutation : (*this))
    { // insert Region of type R (if necessary)
      if (mutation.position > pos) ++count;
      pos = mutation.position + mutation.getLength();
    }
    regions->reserve(++count); // avoid realloc of vector data 
    
    pos = 0;
    */

    PositionType pos = 0;
    for (auto& mutation: (*this))
    {
        // insert Region of type R (if necessary)
        if (mutation.position > pos)
            regions->emplace_back(TYPE_R, mutation.position - 1);
        
        // convert the current mutation
        pos = mutation.position + mutation.getLength();
        regions->emplace_back(&mutation, seq_type, num_states);
    }
    
    // insert the last Region of type R (if necessary)
    if (pos < sequence_length)
        regions->emplace_back(TYPE_R, sequence_length - 1);
    
    return regions;
}
