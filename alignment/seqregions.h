//
//  seqregions.h
//  alignment
//
//  Created by Nhan Ly-Trong on 24/01/2022.
//

#ifndef SEQREGIONS_H
#define SEQREGIONS_H

#include "seqregion.h"
#include "alignment.h"
#include "model/model.h"

class Alignment;

/** Vector of sequence regions, used to represent/compute partial/total likelihood */
class SeqRegions: public vector<SeqRegion*> {
private:
    
public:
    
    /**
    *  Regions constructor
    */
    SeqRegions();
    
    /**
    *  Regions deconstructor
    */
    ~SeqRegions();
    
    /**
    *  Copy regions
    */
    void copyRegions(SeqRegions* n_regions, StateType num_states);
    
    /**
    *  Delete all region items one by one
    */
    void deleteRegions();
    
    /**
        Merge consecutive R regions
     */
    void mergeRegionR(StateType num_states, RealNumType threshold);
    
    /**
        Move to the next region in the vector of regions
        @param sequence: a vector of regions; region_index: the index of the next region
        @return region: a region; current_pos: the current position; end_pos: the end position
     */
    static void move2NextRegion(SeqRegions* sequence, PositionType region_index, SeqRegion* &region, PositionType &current_pos, PositionType &end_pos, PositionType seq_length);
    
    /**
        Get the shared segment between the next regions of two sequences
        @param current_pos: current site posisition; sequence1, sequence2: vectors of regions; seq1_index, seq2_index: the indexes of the current regions; seq1_end_pos, seq2_end_pos: the end positions of the current regions
        @return seq1_region, seq2_region: the regions contains the shared segment; length: length of the shared segment
     */
    static void getNextSharedSegment(PositionType current_pos, PositionType seq_length, SeqRegions* sequence1, SeqRegions* sequence2, PositionType &seq1_index, PositionType &seq2_index, SeqRegion* &seq1_region, SeqRegion* &seq2_region, PositionType &seq1_end_pos, PositionType &seq2_end_pos, PositionType &length);
    
    /**
        Compare the current sequence with another sequence regarding the amount of information
        @param sequence2 the sequence to compare
        @param seq_length the sequence length
        @param num_states the number of states
        @return 0: if the two sequences are incomparable; 1: if the current sequence is more or equally informative than/to sequence2; -1; if the current sequence  is less informative than sequence2
     */
    int compareWithSample(SeqRegions* sequence2, PositionType seq_length, StateType num_states);
    
    /**
        Check if the current regions and regions2 represent the same partial likelihoods or not -> be used to stop traversing the tree further for updating partial likelihoods
     */
    bool areDiffFrom(SeqRegions* regions2, PositionType seq_length, StateType num_states, Params* params);
    
    /**
        Merge two likelihood vectors, one from above and one from below
        This regions is the upper regions
        @param lower_regions the lower regions
        @param upper_plength the length of the upper branch
        @param lower_plength the length of the lower branch
        @param merged_regions the output regions by merging the upper and lower regions
        @param aln the alignment
        @param model the model of evolution
        @param threshold the threshold for approximation
     */
    void mergeUpperLower(SeqRegions* &merged_regions, RealNumType upper_plength, SeqRegions* lower_regions, RealNumType lower_plength, Alignment* aln, Model* model, RealNumType threshold);
    
    /**
        Merge two lower likelihood vectors
        This regions is one of the two lower regions
        @param regions2 the other lower regions
        @param plength1, plength2 the lengths of the lower branches
        @param merged_regions the output regions by merging the upper and lower regions
        @param aln the alignment
        @param model the model of evolution
        @param threshold the threshold for approximation
        @param cumulative_rate the cumulative rates of the reference sequence
        @param return_log_lh TRUE to return the log likelihood
     */
    RealNumType mergeTwoLowers(SeqRegions* &merged_regions, RealNumType plength1, SeqRegions* regions2, RealNumType plength2, Alignment* aln, Model* model, RealNumType threshold, RealNumType* cumulative_rate, bool return_log_lh = false);
    
    /**
        Compute total lh/upper left_right for root node
        @param blength the branch length; (-1 by default).
     */
    SeqRegions* computeTotalLhAtRoot(StateType num_states, Model* model, RealNumType blength = -1);
    
    /**
        Compute the likelihood by merging the lower lh with root frequencies
     */
    RealNumType computeAbsoluteLhAtRoot(Alignment* aln, Model* model, vector< vector<PositionType> > &cumulative_base);
    
    /**
        Convert an entry 'O' into a normal nucleotide if its probability dominated others
     */
    StateType simplifyO(RealNumType* &partial_lh, StateType ref_state, StateType num_states, RealNumType threshold);
};
#endif
