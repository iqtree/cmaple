//
//  sequence.h
//  alignment
//
//  Created by Nhan Ly-Trong on 24/01/2022.
//

#include "mutationindel.h"
#include "seqregions.h"

#ifndef SEQUENCE_H
#define SEQUENCE_H

class Alignment;
class SeqRegions;

/** Class present a sequence */
class Sequence: public vector<Mutation*> {
public:
    /**
        Name of the sequence
     */
    string seq_name;
    
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
        Extract the lower likelihood vector (converting a vector of Mutations into a vector of Regions)
    */
    SeqRegions* getLowerLhVector(PositionType sequence_length, StateType num_states, SeqType seq_type);
};
#endif
