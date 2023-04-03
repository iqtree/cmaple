//
//  seqregions.h
//  alignment
//
//  Created by Nhan Ly-Trong on 24/01/2022.
//

#pragma once

#include "seqregion.h"
#include "alignment.h"
#include "model/model.h"
#include <algorithm>

class Alignment;
class Model;

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
    explicit SeqRegions(const std::unique_ptr<SeqRegions>& n_regions);
    
    /// Move CTor
    SeqRegions(SeqRegions&& regions) = default;
    /// Move Assignment
    SeqRegions& operator=(SeqRegions&& regions) = default;
    
    /**
        Add a new region and automatically merged consecutive R regions
     */
    static void addNonConsecutiveRRegion(SeqRegions& regions, const StateType new_region_type, const RealNumType plength_observation2node, const RealNumType plength_observation2root, const PositionType end_pos, const RealNumType threshold_prob)
    {
        // cannot merge consecutive R regions if no region exists in regions
        if (!regions.empty())
        {
            // try to merge consecutive R regions
            SeqRegion& last_region = regions.back();
            if (new_region_type == TYPE_R
                && last_region.type == TYPE_R
                && fabs(last_region.plength_observation2node - plength_observation2node) < threshold_prob
                && fabs(last_region.plength_observation2root - plength_observation2root) < threshold_prob)
            {
                last_region.position = end_pos;
                last_region.plength_observation2node = plength_observation2node;
                last_region.plength_observation2root = plength_observation2root;
                return;
            }
        }
        
        // if we cannot merge new region into existing R region => just add a new one
        regions.emplace_back(new_region_type, end_pos, plength_observation2node, plength_observation2root);
    }
    
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
        Count the number of shared segments
     */
    size_t countSharedSegments(const SeqRegions& seq2_regions, const size_t seq_length) const;
    
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
    void mergeUpperLower(std::unique_ptr<SeqRegions>& merged_regions, RealNumType upper_plength, const SeqRegions& lower_regions, RealNumType lower_plength, const
                         Alignment& aln, const Model& model, const RealNumType threshold) const;
    
    /**
        Merge two lower likelihood vectors
        This regions is one of the two lower regions
        @param regions2 the other lower regions
        @param plength1, plength2 the lengths of the lower branches
        @param merged_regions the output regions by merging the upper and lower regions
        @param aln the alignment
        @param model the model of evolution
        @param threshold the threshold for approximation
        @param return_log_lh TRUE to return the log likelihood
     */
    RealNumType mergeTwoLowers(std::unique_ptr<SeqRegions>& merged_regions, const RealNumType plength1, const SeqRegions& regions2, const RealNumType plength2, const Alignment& aln, const Model& model, const RealNumType threshold_prob, const bool return_log_lh = false) const;
    
    /**
        Compute total lh/upper left_right for root node
        @param blength the branch length; (-1 by default).
     */
    void computeTotalLhAtRoot(std::unique_ptr<SeqRegions>& total_lh, const StateType num_states, const Model& model, RealNumType blength = -1) const;
    
    /**
        Compute the likelihood by merging the lower lh with root frequencies
     */
    RealNumType computeAbsoluteLhAtRoot(const StateType num_states, const Model& model);
    
    /**
        Convert an entry 'O' into a normal nucleotide if its probability dominated others
     */
    static StateType simplifyO(RealNumType* const partial_lh, StateType ref_state, StateType num_states, RealNumType threshold)
    {
        // dummy variables
        ASSERT(partial_lh);
        RealNumType max_prob = 0;
        StateType max_index = 0;
        StateType high_prob_count = 0;
        
        // Check all states one by one
        for (StateType i = 0; i < num_states; ++i)
        {
            // record the state with the highest likelihood
            if (partial_lh[i] > max_prob)
            {
                max_prob = partial_lh[i];
                max_index = i;
            }
            
            // count the number of states that have the likelihood greater than a threshold
            if (partial_lh[i] > threshold)
                ++high_prob_count;
        }
        
        // if the partial lh concentrates at a specific state -> return new state
        if (high_prob_count == 1)
        {
            // new state matches with the reference state
            if (max_index == ref_state)
                return TYPE_R;
            // return a nucleotide
            else
                return max_index;
        }
        // otherwise, cannot simplify
        else
            return TYPE_O;
    }
    
    static void addSimplifiedO(const PositionType end_pos, SeqRegion::LHType& new_lh, const Alignment& aln, const RealNumType threshold_prob, SeqRegions& merged_regions)
    {
        StateType new_state = SeqRegions::simplifyO(new_lh.data(), aln.ref_seq[end_pos], aln.num_states, threshold_prob);

        if (new_state == TYPE_O)
            merged_regions.emplace_back(TYPE_O, end_pos, 0, 0, std::move(new_lh));
        else
        {
            // add a new region and try to merge consecutive R regions together
            SeqRegions::addNonConsecutiveRRegion(merged_regions, new_state, -1, -1, end_pos, threshold_prob);
        }
    }
    
    /**
        For testing only, export codes to re-contruct this seqregions
     */
    void writeConstructionCodes(const std::string regions_name, std::ofstream& out, const StateType num_states) const;
    
    /**
        Compare two regions
     */
    bool operator==(const SeqRegions& seqregions_1) const;
};

/**
    MergeUpperLower case N with O
 */
void merge_N_O(const RealNumType lower_plength, const SeqRegion& reg_o, const Model& model,
               const PositionType end_pos, const StateType num_states, SeqRegions& merged_target);

/**
    MergeUpperLower case N_RACGT
 */
void merge_N_RACGT(const SeqRegion& reg_racgt, const RealNumType lower_plength, const PositionType end_pos,
                   const RealNumType threshold_prob, SeqRegions& merged_regions);

/**
    MergeUpperLower case O_N
 */
void merge_O_N(const SeqRegion& reg_o, const RealNumType upper_plength, const PositionType end_pos, const Model& model, const StateType num_states, SeqRegions& merged_regions);

/**
    MergeUpperLower case RACGT_N
 */
void merge_RACGT_N(const SeqRegion& reg_n, const RealNumType upper_plength, const PositionType end_pos,
                   const RealNumType threshold_prob, SeqRegions& merged_regions);

/**
    MergeUpperLower case Zero_Distance
 */
bool merge_Zero_Distance(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType total_blength_1, const RealNumType total_blength_2, const PositionType end_pos, const RealNumType threshold_prob, const StateType num_states, std::unique_ptr<SeqRegions>& merged_regions);

/**
    MergeUpperLower case O_ORACGT
 */
void merge_O_ORACGT(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType total_blength_1, const RealNumType total_blength_2, const PositionType end_pos, const RealNumType threshold_prob, const Model& model, const Alignment& aln, SeqRegions& merged_regions);

/**
    MergeUpperLower case RACGT_O
 */
void merge_RACGT_O(const SeqRegion& seq2_region, const RealNumType total_blength_2, const PositionType end_pos, SeqRegion::LHType& new_lh, const RealNumType threshold_prob, const Model& model, const Alignment& aln, SeqRegions& merged_regions);

/**
    MergeUpperLower case RACGT_RACGT
 */
void merge_RACGT_RACGT(const SeqRegion& seq2_region, const RealNumType total_blength_2, const PositionType end_pos, SeqRegion::LHType& new_lh, const Model& model, const Alignment& aln, SeqRegions& merged_regions);

/**
    MergeUpperLower case RACGT_ORACGT
 */
void merge_RACGT_ORACGT(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType total_blength_1, const RealNumType total_blength_2, const RealNumType upper_plength, const PositionType end_pos, const RealNumType threshold_prob, const Model& model, const Alignment& aln, SeqRegions& merged_regions);

/**
    MergeTwoLowers case N with O
 */
void merge_N_O_TwoLowers(const SeqRegion& seq2_region, const PositionType end_pos, const RealNumType plength2, SeqRegions& merged_regions);

/**
    MergeTwoLowers case N_RACGT
 */
void merge_N_RACGT_TwoLowers(const SeqRegion& seq2_region, const PositionType end_pos, const RealNumType plength2, const RealNumType threshold_prob, SeqRegions& merged_regions);

/**
    MergeTwoLowers case identicalRACGT
 */
void merge_identicalRACGT_TwoLowers(const SeqRegion& seq1_region, const PositionType end_pos, RealNumType total_blength_1, RealNumType total_blength_2, const PositionType pos, const RealNumType threshold_prob, const Model& model, RealNumType &log_lh, SeqRegions& merged_regions, const bool return_log_lh);

/**
    MergeTwoLowers case O_O
 */
bool merge_O_O_TwoLowers(const SeqRegion& seq2_region, RealNumType total_blength_2, const PositionType end_pos, const Alignment& aln, const Model& model, const RealNumType threshold_prob, RealNumType &log_lh, SeqRegion::LHType& new_lh, std::unique_ptr<SeqRegions>& merged_regions, const bool return_log_lh);

/**
    MergeTwoLowers case O_RACGT
 */
bool merge_O_RACGT_TwoLowers(const SeqRegion& seq2_region, RealNumType total_blength_2, const PositionType end_pos, const Alignment& aln, const Model& model, const RealNumType threshold_prob, RealNumType &log_lh, SeqRegion::LHType& new_lh, RealNumType& sum_lh, std::unique_ptr<SeqRegions>& merged_regions, const bool return_log_lh);

/**
    MergeTwoLowers case O_ORACGT
 */
bool merge_O_ORACGT_TwoLowers(const SeqRegion& seq1_region, const SeqRegion& seq2_region, RealNumType total_blength_1, RealNumType total_blength_2, const PositionType end_pos, const Alignment& aln, const Model& model, const RealNumType threshold_prob, RealNumType &log_lh, std::unique_ptr<SeqRegions>& merged_regions, const bool return_log_lh);

/**
    MergeTwoLowers case RACGT_O
 */
bool merge_RACGT_O_TwoLowers(const SeqRegion& seq2_region, RealNumType total_blength_2, const PositionType end_pos, const Alignment& aln, const Model& model, const RealNumType threshold_prob, SeqRegion::LHType& new_lh, RealNumType &log_lh, std::unique_ptr<SeqRegions>& merged_regions, const bool return_log_lh);

/**
    MergeTwoLowers case RACGT_RACGT
 */
bool merge_RACGT_RACGT_TwoLowers(const SeqRegion& seq2_region, RealNumType total_blength_2, const PositionType end_pos, const Alignment& aln, const Model& model, const RealNumType threshold_prob, SeqRegion::LHType& new_lh, RealNumType& sum_lh, RealNumType &log_lh, std::unique_ptr<SeqRegions>& merged_regions, const bool return_log_lh);

/**
    MergeTwoLowers case RACGT_ORACGT
 */
bool merge_RACGT_ORACGT_TwoLowers(const SeqRegion& seq1_region, const SeqRegion& seq2_region, RealNumType total_blength_1, RealNumType total_blength_2, const PositionType end_pos, const Alignment& aln, const Model& model, const RealNumType threshold_prob, RealNumType &log_lh, std::unique_ptr<SeqRegions>& merged_regions, const bool return_log_lh);

/**
    MergeTwoLowers case notN_notN
 */
bool merge_notN_notN_TwoLowers(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType plength1, const RealNumType plength2, const PositionType end_pos, const PositionType pos, const Alignment& aln, const Model& model, const RealNumType threshold_prob, RealNumType &log_lh, std::unique_ptr<SeqRegions>& merged_regions, const bool return_log_lh);

