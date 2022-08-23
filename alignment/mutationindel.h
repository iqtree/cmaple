//
//  mutationindel.h
//  alignment
//
//  Created by Nhan Ly-Trong on 24/01/2022.
//

#include "mutation.h"

#ifndef MUTATIONINDEL_H
#define MUTATIONINDEL_H

/** An indel mutation contains gaps */
class MutationIndel: public Mutation {
private:
    /**
        Length of mutation (continuous sites)
     */
    PositionType length;
    
public:
    
    /**
        MutationIndel constructor
    */
    MutationIndel();
    
    /**
        MutationIndel constructor
    */
    MutationIndel(StateType n_type, PositionType n_position, PositionType n_length);
    
    /**
        Return Length of the current mutation
    */
    virtual PositionType getLength() { return length;};
};
#endif
