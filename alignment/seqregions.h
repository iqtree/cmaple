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
#include <algorithm>

class Alignment;

inline PositionType minFast(PositionType a, PositionType b)
{ // some compilers cannot optimize std::min as well as the ternary op
  return a < b ? a : b;
}

/** Vector of sequence regions, used to represent/compute partial/total likelihood */
class SeqRegions: public std::vector<SeqRegion> {
private:
    
public:
    
    /**
    *  Regions constructor
    */
    SeqRegions() = default;
    
    /**
    *  Regions constructor
    */
    explicit SeqRegions(SeqRegions* n_regions);
    
    /**
    *  Regions destructor
    */
    ~SeqRegions();
    
    /**
    *  Delete all region items one by one
    */
    void deleteRegions();
    
    /**
        Merge consecutive R regions
     */
    void mergeRegionR(StateType num_states, RealNumType threshold);
    
    /**
        Get the shared segment between the next regions of two sequences
        @param current_pos: current site position; 
        @return seq1_region, seq2_region: the regions contains the shared segment; end_pos: ending position of the shared segment
     */
    inline static void getNextSharedSegment(PositionType current_pos, const SeqRegions& seq1_region, const SeqRegions& seq2_region, size_t& i1, size_t& i2, PositionType &end_pos)
    {  // 14% of runtime, invoked from Tree::calculateSubTreePlacementCostTemplate
      if (current_pos > seq1_region[i1].position) ++i1;
      if (current_pos > seq2_region[i2].position) ++i2;
        
      // compute the end_pos for the shared segment
      end_pos = minFast(seq1_region[i1].position, seq2_region[i2].position);
    }
    
    /**
        Compare the current sequence with another sequence regarding the amount of information
        @param sequence2 the sequence to compare
        @param seq_length the sequence length
        @param num_states the number of states
        @return 0: if the two sequences are incomparable; 1: if the current sequence is more or equally informative than/to sequence2; -1; if the current sequence  is less informative than sequence2
     */
    int compareWithSample(const SeqRegions& sequence2, PositionType seq_length, StateType num_states) const;
    
    /**
        Check if the current regions and regions2 represent the same partial likelihoods or not -> be used to stop traversing the tree further for updating partial likelihoods
     */
    bool areDiffFrom(const SeqRegions& regions2, PositionType seq_length, StateType num_states, const Params* params) const;
    
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
    void mergeUpperLower(SeqRegions* &merged_regions, RealNumType upper_plength, const SeqRegions& lower_regions, RealNumType lower_plength, const
                         Alignment& aln, const Model& model, RealNumType threshold) const;
    
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
    RealNumType mergeTwoLowers(SeqRegions* &merged_regions, RealNumType plength1, const SeqRegions* const regions2, RealNumType plength2, const Alignment& aln, const
                               Model& model, RealNumType threshold, RealNumType* cumulative_rate, bool return_log_lh = false);
    
    /**
        Compute total lh/upper left_right for root node
        @param blength the branch length; (-1 by default).
     */
    SeqRegions* computeTotalLhAtRoot(StateType num_states, const Model& model, RealNumType blength = -1);
    
    /**
        Compute the likelihood by merging the lower lh with root frequencies
     */
    RealNumType computeAbsoluteLhAtRoot(const Alignment& aln, const Model& model, std::vector< std::vector<PositionType> > &cumulative_base);
    
    /**
        Convert an entry 'O' into a normal nucleotide if its probability dominated others
     */
    StateType simplifyO(RealNumType* const partial_lh, StateType ref_state, StateType num_states, RealNumType threshold) const;
};
#endif
