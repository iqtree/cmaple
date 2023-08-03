#pragma once

#include "seqregion.h"
#include "alignment.h"
#include "../model/modelbase.h"
#include <algorithm>

namespace cmaple
{
    class Alignment;
    class ModelBase;

    inline cmaple::PositionType minFast(cmaple::PositionType a, cmaple::PositionType b)
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
         *  @throw std::invalid\_argument if n\_regions is null
         */
        explicit SeqRegions(const std::unique_ptr<SeqRegions>& n_regions);
        
        /// Move CTor
        SeqRegions(SeqRegions&& regions) = default;
        /// Move Assignment
        SeqRegions& operator=(SeqRegions&& regions) = default;
        
        /**
         Add a new region and automatically merged consecutive R regions
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        static void addNonConsecutiveRRegion(SeqRegions& regions, const cmaple::StateType new_region_type, const cmaple::RealNumType plength_observation2node, const cmaple::RealNumType plength_observation2root, const cmaple::PositionType end_pos, const cmaple::RealNumType threshold_prob)
        {
            // cannot merge consecutive R regions if no region exists in regions
            if (!regions.empty())
            {
                // try to merge consecutive R regions
                SeqRegion& last_region = regions.back();
                if (new_region_type == cmaple::TYPE_R
                    && last_region.type == cmaple::TYPE_R
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
        inline static void getNextSharedSegment(cmaple::PositionType current_pos, const SeqRegions& seq1_region, const SeqRegions& seq2_region, size_t& i1, size_t& i2, cmaple::PositionType &end_pos)
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
        int compareWithSample(const SeqRegions& sequence2, cmaple::PositionType seq_length, cmaple::StateType num_states) const;
        
        /**
         Check if the current regions and regions2 represent the same partial likelihoods or not -> be used to stop traversing the tree further for updating partial likelihoods
         */
        bool areDiffFrom(const SeqRegions& regions2, cmaple::PositionType seq_length, cmaple::StateType num_states, const cmaple::Params& params) const;
        
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
         
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType num_states>
        void mergeUpperLower(std::unique_ptr<SeqRegions>& merged_regions, cmaple::RealNumType upper_plength, const SeqRegions& lower_regions, cmaple::RealNumType lower_plength, const
                             Alignment* aln, const ModelBase* model, const cmaple::RealNumType threshold) const;
        
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
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType num_states>
        cmaple::RealNumType mergeTwoLowers(std::unique_ptr<SeqRegions>& merged_regions, const cmaple::RealNumType plength1, const SeqRegions& regions2, const cmaple::RealNumType plength2, const Alignment* aln, const ModelBase* model, const cmaple::RealNumType threshold_prob, const bool return_log_lh = false) const;
        
        /**
         Calculate site-lh contributions by merging two lower likelihood vectors
         This regions is one of the two lower regions
         @param regions2 the other lower regions
         @param plength1, plength2 the lengths of the lower branches
         @param merged_regions the output regions by merging the upper and lower regions
         @param aln the alignment
         @param model the model of evolution
         @param threshold the threshold for approximation
         @param return_log_lh TRUE to return the log likelihood
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType num_states>
        cmaple::RealNumType calculateSiteLhContributions(std::vector<cmaple::RealNumType>& site_lh_contributions, std::unique_ptr<SeqRegions>& merged_regions, const cmaple::RealNumType plength1, const SeqRegions& regions2, const cmaple::RealNumType plength2, const Alignment* aln, const ModelBase* model, const cmaple::RealNumType threshold_prob) const;
        
        /**
         Compute total lh/upper left_right for root node
         @param blength the branch length; (-1 by default).
         
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType num_states>
        void computeTotalLhAtRoot(std::unique_ptr<SeqRegions>& total_lh, const ModelBase* model, cmaple::RealNumType blength = -1) const;
        
        /**
         Compute the likelihood by merging the lower lh with root frequencies
         */
        template <const cmaple::StateType num_states>
        cmaple::RealNumType computeAbsoluteLhAtRoot(const ModelBase* model);
        
        /**
         Compute the site likelihood at root by merging the lower lh with root frequencies
         */
        template <const cmaple::StateType num_states>
        cmaple::RealNumType computeSiteLhAtRoot(std::vector<cmaple::RealNumType>& site_lh_contributions, const ModelBase* model);
        
        /**
         Convert an entry 'O' into a normal nucleotide if its probability dominated others
         */
        static cmaple::StateType simplifyO(cmaple::RealNumType* const partial_lh, cmaple::StateType ref_state, cmaple::StateType num_states, cmaple::RealNumType threshold)
        {
            // dummy variables
            ASSERT(partial_lh);
            cmaple::RealNumType max_prob = 0;
            cmaple::StateType max_index = 0;
            cmaple::StateType high_prob_count = 0;
            
            // Check all states one by one
            for (cmaple::StateType i = 0; i < num_states; ++i)
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
                    return cmaple::TYPE_R;
                // return a nucleotide
                else
                    return max_index;
            }
            // otherwise, cannot simplify
            else
                return cmaple::TYPE_O;
        }
        
        /*
         Add an simplified O seqregion to the current vectors of seqregions
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        static void addSimplifiedO(const cmaple::PositionType end_pos, SeqRegion::LHType& new_lh, const Alignment* aln, const cmaple::RealNumType threshold_prob, SeqRegions& merged_regions)
        {
            cmaple::StateType new_state = SeqRegions::simplifyO(new_lh.data(), aln->ref_seq[end_pos], aln->num_states, threshold_prob);
            
            if (new_state == cmaple::TYPE_O)
                merged_regions.emplace_back(cmaple::TYPE_O, end_pos, 0, 0, std::move(new_lh));
            else
            {
                // add a new region and try to merge consecutive R regions together
                SeqRegions::addNonConsecutiveRRegion(merged_regions, new_state, -1, -1, end_pos, threshold_prob);
            }
        }
        
        /**
         For testing only, export codes to re-contruct this seqregions
         */
        void writeConstructionCodes(const std::string regions_name, std::ofstream& out, const cmaple::StateType num_states) const;
        
        /**
         Compare two regions
         */
        bool operator==(const SeqRegions& seqregions_1) const;
    };

    /**
     MergeUpperLower case N with O
     @throw std::logic\_error if unexpected values/behaviors found during the operations
     */
    template <const cmaple::StateType num_states>
    void merge_N_O(const cmaple::RealNumType lower_plength, const SeqRegion& reg_o, const ModelBase* model,
                   const cmaple::PositionType end_pos, SeqRegions& merged_target);

    /**
     MergeUpperLower case N_RACGT
     
     @throw std::logic\_error if unexpected values/behaviors found during the operations
     */
    void merge_N_RACGT(const SeqRegion& reg_racgt, const cmaple::RealNumType lower_plength, const cmaple::PositionType end_pos,
                       const cmaple::RealNumType threshold_prob, SeqRegions& merged_regions);

    /**
     MergeUpperLower case O_N
     
     @throw std::logic\_error if unexpected values/behaviors found during the operations
     */
    template <const cmaple::StateType num_states>
    void merge_O_N(const SeqRegion& reg_o, const cmaple::RealNumType upper_plength, const cmaple::PositionType end_pos, const ModelBase* model, SeqRegions& merged_regions);

    /**
     MergeUpperLower case RACGT_N
     
     @throw std::logic\_error if unexpected values/behaviors found during the operations
     */
    void merge_RACGT_N(const SeqRegion& reg_n, const cmaple::RealNumType upper_plength, const cmaple::PositionType end_pos,
                       const cmaple::RealNumType threshold_prob, SeqRegions& merged_regions);

    /**
     MergeUpperLower case Zero_Distance
     
     @throw std::logic\_error if unexpected values/behaviors found during the operations
     */
    bool merge_Zero_Distance(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const cmaple::RealNumType total_blength_1, const cmaple::RealNumType total_blength_2, const cmaple::PositionType end_pos, const cmaple::RealNumType threshold_prob, const cmaple::StateType num_states, std::unique_ptr<SeqRegions>& merged_regions);

    /**
     MergeUpperLower case O_ORACGT
     
     @throw std::logic\_error if unexpected values/behaviors found during the operations
     */
    template <const cmaple::StateType num_states>
    void merge_O_ORACGT(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const cmaple::RealNumType total_blength_1, const cmaple::RealNumType total_blength_2, const cmaple::PositionType end_pos, const cmaple::RealNumType threshold_prob, const ModelBase* model, const Alignment* aln, SeqRegions& merged_regions);

    /**
     MergeUpperLower case RACGT_O
     
     @throw std::logic\_error if unexpected values/behaviors found during the operations
     */
    template <const cmaple::StateType num_states>
    void merge_RACGT_O(const SeqRegion& seq2_region, const cmaple::RealNumType total_blength_2, const cmaple::PositionType end_pos, SeqRegion::LHType& new_lh, const cmaple::RealNumType threshold_prob, const ModelBase* model, const Alignment* aln, SeqRegions& merged_regions);

    /**
     MergeUpperLower case RACGT_RACGT
     
     @throw std::logic\_error if unexpected values/behaviors found during the operations
     */
    template <const cmaple::StateType num_states>
    void merge_RACGT_RACGT(const SeqRegion& seq2_region, const cmaple::RealNumType total_blength_2, const cmaple::PositionType end_pos, SeqRegion::LHType& new_lh, const ModelBase* model, const Alignment* aln, SeqRegions& merged_regions);

    /**
     MergeUpperLower case RACGT_ORACGT
     
     @throw std::logic\_error if unexpected values/behaviors found during the operations
     */
    template <const cmaple::StateType num_states>
    void merge_RACGT_ORACGT(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const cmaple::RealNumType total_blength_1, const cmaple::RealNumType total_blength_2, const cmaple::RealNumType upper_plength, const cmaple::PositionType end_pos, const cmaple::RealNumType threshold_prob, const ModelBase* model, const Alignment* aln, SeqRegions& merged_regions);

    /**
     MergeTwoLowers case N with O
     */
    void merge_N_O_TwoLowers(const SeqRegion& seq2_region, const cmaple::PositionType end_pos, const cmaple::RealNumType plength2, SeqRegions& merged_regions);

    /**
     MergeTwoLowers case N_RACGT
     
     @throw std::logic\_error if unexpected values/behaviors found during the operations
     */
    void merge_N_RACGT_TwoLowers(const SeqRegion& seq2_region, const cmaple::PositionType end_pos, const cmaple::RealNumType plength2, const cmaple::RealNumType threshold_prob, SeqRegions& merged_regions);

    /**
     MergeTwoLowers case identicalRACGT
     @throw std::logic\_error if unexpected values/behaviors found during the operations
     */
    void merge_identicalRACGT_TwoLowers(const SeqRegion& seq1_region, const cmaple::PositionType end_pos, cmaple::RealNumType total_blength_1, cmaple::RealNumType total_blength_2, const cmaple::PositionType pos, const cmaple::RealNumType threshold_prob, const ModelBase* model, cmaple::RealNumType &log_lh, SeqRegions& merged_regions, const bool return_log_lh);

    /**
     MergeTwoLowers case O_O
     
     @throw std::logic\_error if unexpected values/behaviors found during the operations
     */
    template <const cmaple::StateType num_states>
    bool merge_O_O_TwoLowers(const SeqRegion& seq2_region, cmaple::RealNumType total_blength_2, const cmaple::PositionType end_pos, const Alignment* aln, const ModelBase* model, const cmaple::RealNumType threshold_prob, cmaple::RealNumType &log_lh, SeqRegion::LHType& new_lh, std::unique_ptr<SeqRegions>& merged_regions, const bool return_log_lh);

    /**
     MergeTwoLowers case O_RACGT
     
     @throw std::logic\_error if unexpected values/behaviors found during the operations
     */
    template <const cmaple::StateType num_states>
    bool merge_O_RACGT_TwoLowers(const SeqRegion& seq2_region, cmaple::RealNumType total_blength_2, const cmaple::PositionType end_pos, const Alignment* aln, const ModelBase* model, const cmaple::RealNumType threshold_prob, cmaple::RealNumType &log_lh, SeqRegion::LHType& new_lh, cmaple::RealNumType& sum_lh, std::unique_ptr<SeqRegions>& merged_regions, const bool return_log_lh);

    /**
     MergeTwoLowers case O_ORACGT
     
     @throw std::logic\_error if unexpected values/behaviors found during the operations
     */
    template <const cmaple::StateType num_states>
    bool merge_O_ORACGT_TwoLowers(const SeqRegion& seq1_region, const SeqRegion& seq2_region, cmaple::RealNumType total_blength_1, cmaple::RealNumType total_blength_2, const cmaple::PositionType end_pos, const Alignment* aln, const ModelBase* model, const cmaple::RealNumType threshold_prob, cmaple::RealNumType &log_lh, std::unique_ptr<SeqRegions>& merged_regions, const bool return_log_lh);

    /**
     MergeTwoLowers case RACGT_O
     
     @throw std::logic\_error if unexpected values/behaviors found during the operations
     */
    template <const cmaple::StateType num_states>
    bool merge_RACGT_O_TwoLowers(const SeqRegion& seq2_region, cmaple::RealNumType total_blength_2, const cmaple::PositionType end_pos, const Alignment* aln, const ModelBase* model, const cmaple::RealNumType threshold_prob, SeqRegion::LHType& new_lh, cmaple::RealNumType &log_lh, std::unique_ptr<SeqRegions>& merged_regions, const bool return_log_lh);

    /**
     MergeTwoLowers case RACGT_RACGT
     
     @throw std::logic\_error if unexpected values/behaviors found during the operations
     */
    template <const cmaple::StateType num_states>
    bool merge_RACGT_RACGT_TwoLowers(const SeqRegion& seq2_region, cmaple::RealNumType total_blength_2, const cmaple::PositionType end_pos, const Alignment* aln, const ModelBase* model, const cmaple::RealNumType threshold_prob, SeqRegion::LHType& new_lh, cmaple::RealNumType& sum_lh, cmaple::RealNumType &log_lh, std::unique_ptr<SeqRegions>& merged_regions, const bool return_log_lh);

    /**
     MergeTwoLowers case RACGT_ORACGT
     
     @throw std::logic\_error if unexpected values/behaviors found during the operations
     */
    template <const cmaple::StateType num_states>
    bool merge_RACGT_ORACGT_TwoLowers(const SeqRegion& seq1_region, const SeqRegion& seq2_region, cmaple::RealNumType total_blength_1, cmaple::RealNumType total_blength_2, const cmaple::PositionType end_pos, const Alignment* aln, const ModelBase* model, const cmaple::RealNumType threshold_prob, cmaple::RealNumType &log_lh, std::unique_ptr<SeqRegions>& merged_regions, const bool return_log_lh);

    /**
     MergeTwoLowers case notN_notN
     
     @throw std::logic\_error if unexpected values/behaviors found during the operations
     */
    template <const cmaple::StateType num_states>
    bool merge_notN_notN_TwoLowers(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const cmaple::RealNumType plength1, const cmaple::RealNumType plength2, const cmaple::PositionType end_pos, const cmaple::PositionType pos, const Alignment* aln, const ModelBase* model, const cmaple::RealNumType threshold_prob, cmaple::RealNumType &log_lh, std::unique_ptr<SeqRegions>& merged_regions, const bool return_log_lh);

}
