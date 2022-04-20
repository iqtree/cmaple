//
//  sequence.h
//  alignment
//
//  Created by Nhan Ly-Trong on 24/01/2022.
//

#include "mutationindel.h"
#include "region.h"

#ifndef SEQUENCE_H
#define SEQUENCE_H

class Sequence {
public:
    string seq_name;
    vector<Mutation*> mutations;
    
    /**
    *  Sequence constructor
    */
    Sequence(){};
    
    /**
    *  Sequence constructor
    */
    Sequence(string n_seq_name);
    
    /**
    *  Sequence deconstructor
    */
    ~Sequence();
    
    /**
    *  Sequence constructor
    */
    Sequence(string n_seq_name, vector<Mutation*> n_mutations);
    
    /**
    *  get a vector of Regions from a vector of Mutations
    */
    vector<Region*> getRegionsFromMutations(PositionType sequence_length, SeqType seq_type, int max_num_states);
};
#endif
