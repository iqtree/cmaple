//
//  regions.h
//  alignment
//
//  Created by Nhan Ly-Trong on 24/01/2022.
//

#include "region.h"

#ifndef REGIONS_H
#define REGIONS_H

class Regions: public vector<Region*> {
private:
    
public:
    
    /**
    *  Regions constructor
    */
    Regions();
    
    /**
    *  Regions deconstructor
    */
    ~Regions();
    
    /**
    *  copy regions
    */
    void copyRegions(Regions* n_regions, StateType num_states);
    
    /**
    *  get a Region
    */
    Region* getRegion(PositionType id);
    
    /**
    *  Delete all region items one by one
    */
    void deleteRegions();
    
    /**
        merge consecutive R regions
     */
    void mergeRegionR(StateType num_states, double threshold);
    
    /**
        move to the next region in the vector of regions
        @param sequence: a vector of regions; region_index: the index of the next region
        @return region: a region; current_pos: the current position; end_pos: the end position
     */
    static void move2NextRegion(Regions* sequence, PositionType region_index, Region* &region, PositionType &current_pos, PositionType &end_pos, PositionType seq_length);
    
    /**
        get the shared segment between the next regions of two sequences
        @param current_pos: current site posisition; sequence1, sequence2: vectors of regions; seq1_index, seq2_index: the indexes of the current regions; seq1_end_pos, seq2_end_pos: the end positions of the current regions
        @return seq1_region, seq2_region: the regions contains the shared segment; length: length of the shared segment
     */
    static void getNextSharedSegment(PositionType current_pos, PositionType seq_length, Regions* sequence1, Regions* sequence2, PositionType &seq1_index, PositionType &seq2_index, Region* &seq1_region, Region* &seq2_region, PositionType &seq1_end_pos, PositionType &seq2_end_pos, PositionType &length);
    
    /**
        compare the current sequence with another sequence regarding the amount of information
        @param sequence2 the sequence to compare
        @param seq_length the sequence length
        @param num_states the number of states
        @return 0: if the two sequences are incomparable; 1: if the current sequence is more or equally informative than/to sequence2; -1; if the current sequence  is less informative than sequence2
     */
    int compareWithSample(Regions* sequence2, PositionType seq_length, StateType num_states);
    
    /**
        Check if the current regions and regions2 represent the same partial likelihoods or not -> be used to stop traversing the tree further for updating partial likelihoods
     */
    bool areDiffFrom(Regions* regions2, PositionType seq_length, StateType num_states, Params* params);
};
#endif
