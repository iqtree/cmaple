//
//  seqregion.h
//  alignment
//
//  Created by Nhan Ly-Trong on 24/01/2022.
//

#include "mutation.h"

#ifndef SEQREGION_H
#define SEQREGION_H

/** A region in a sequence */
class SeqRegion: public Mutation {
private:
    /**
    *  Compute the relative likelihood for an ambiguious state
    */
    void computeLhAmbiguity(IntVector entries);
    
    /**
    *  Convert Ambiguious state into typical states: nucleotides/amino-acids...; N; O; R
    */
    void convertAmbiguiousState(SeqType seq_type, int max_num_states);
    
    /**
    *  Convert Ambiguious state (of DNA data) into typical states: A/C/G/T; N; O; R
    */
    void convertAmbiguiousStateDNA(int max_num_states);
    
public:
    /**
        Length of the path between the current phylo node and the node where the likelihood is calculated
     */
    RealNumType plength_observation2node;
    
    /**
        Distance separates the observation at root from the observation at the current node
        To take into account that the observation might have occurred on the other side of the phylogeny with respect to the root;
        for example if the entry is (type = 1, pos = 234, plength_observation = 0.0001, plength_from_root = 0.0002) then this means that a "C" was observed at genome position 234, and the observation is separated from the root by a distance of 0.0001, while a distance of 0.0002 separates the root from the current node (or position along a branch) considered.
     */
    RealNumType plength_observation2root;
    
    /**
        The relative partial likelihood
     */
    RealNumType* likelihood;
    
    /**
    *  Region constructor
    */
    SeqRegion();
    
    /**
    *  Region constructor
    */
    SeqRegion(StateType n_type, PositionType n_position, RealNumType n_plength_observation = -1, RealNumType n_plength_from_root = -1, RealNumType* n_likelihood = NULL);
    
    /**
    *  Region constructor
    */
    SeqRegion(StateType n_type, PositionType n_position, SeqType seq_type, int max_num_states);
    
    /**
    *  Region constructor
    */
    SeqRegion(Mutation* n_mutation, SeqType seq_type, int max_num_states);
    
    /**
    *  Region constructor
    */
    SeqRegion(SeqRegion* region, StateType num_states, bool copy_likelihood = true);
    
    /**
    *  Region deconstructor
    */
    ~SeqRegion();
};
#endif
