//
//  sequence.h
//  alignment
//
//  Created by Nhan Ly-Trong on 24/01/2022.
//

#include "mutationindel.h"
#include "regions.h"

#ifndef SEQUENCE_H
#define SEQUENCE_H

class Alignment;
class Regions;

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
    *  Sequence constructor
    */
    Sequence(string n_seq_name, vector<Mutation*> n_mutations);
    
    /**
    *  Sequence deconstructor
    */
    ~Sequence();
    
    /**
    *  get lower likelihood vector (converting a vector of Mutations into a vector of Regions)
    */
    Regions* getLowerLhVector(PositionType sequence_length, StateType num_states, SeqType seq_type);
};
#endif
