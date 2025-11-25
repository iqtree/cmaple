#pragma once

#include <algorithm>
#include "../model/modelbase.h"
#include "alignment.h"
#include "seqregion.h"
#include "../utils/tools.h"

namespace cmaple {

class ModelBase;

inline cmaple::PositionType minFast(
    cmaple::PositionType a,
    cmaple::PositionType b) {  // some compilers cannot optimize std::min as
                               // well as the ternary op
  return a < b ? a : b;
}

using DoubleState = uint16_t;
static constexpr DoubleState NN = (DoubleState(TYPE_N) << 8) | TYPE_N;
static constexpr DoubleState NO = (DoubleState(TYPE_N) << 8) | TYPE_O;
static constexpr DoubleState RR = (DoubleState(TYPE_R) << 8) | TYPE_R;
static constexpr DoubleState RO = (DoubleState(TYPE_R) << 8) | TYPE_O;
static constexpr DoubleState OO = (DoubleState(TYPE_O) << 8) | TYPE_O;
static constexpr DoubleState ON = (DoubleState(TYPE_O) << 8) | TYPE_N;

/** Vector of sequence regions, used to represent/compute partial/total
 *  likelihood
 */
class SeqRegions : public std::vector<SeqRegion> {
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
     Clone a SeqRegions instance
     */
    static inline SeqRegions clone(const SeqRegions& from) {
        // clone regions one by one
        SeqRegions new_seq_regions;
        new_seq_regions.reserve(from.size());
        for (const auto& region : from) {
            new_seq_regions.push_back(SeqRegion::clone(region));
        }
        return new_seq_regions;
    }

  /**
   Add a new region and automatically merged consecutive R regions
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  static void addNonConsecutiveRRegion(
      SeqRegions& regions,
      const cmaple::StateType new_region_type,
      const cmaple::StateType new_region_prev_state,
      const cmaple::RealNumType plength_observation2node,
      const cmaple::RealNumType plength_observation2root,
      const cmaple::PositionType end_pos,
      const cmaple::RealNumType threshold_prob);

  /**
   Get the shared segment between the next regions of two sequences
   @param current_pos: current site position;
   @return seq1_region, seq2_region: the regions contains the shared segment;
   end_pos: ending position of the shared segment
   */
  inline static void getNextSharedSegment(
      cmaple::PositionType current_pos,
      const SeqRegions& seq1_region,
      const SeqRegions& seq2_region,
      size_t& i1,
      size_t& i2,
      cmaple::PositionType&
          end_pos) {  // 14% of runtime, invoked from
                      // Tree::calculateSubTreePlacementCostTemplate
    assert(seq1_region.size() > i1);
    assert(seq2_region.size() > i2);
      
    if (current_pos > seq1_region[i1].position)
      ++i1;
    if (current_pos > seq2_region[i2].position)
      ++i2;

    // compute the end_pos for the shared segment
    end_pos = minFast(seq1_region[i1].position, seq2_region[i2].position);
  }

  /**
   Count the number of shared segments
   */
  size_t countSharedSegments(const SeqRegions& seq2_regions,
                             const size_t seq_length) const;

  /**
   Compare the current sequence with another sequence regarding the amount of
   information
   @param sequence2 the sequence to compare
   @param seq_length the sequence length
   @param Alignment the alignment
   @param check_ident_only TRUE if only checking whether two samples are identical or not
   @return 0: if the two sequences are incomparable; 1: if the current sequence
   is more or equally informative than/to sequence2; -1; if the current sequence
   is less informative than sequence2; if check\_ident\_only is TRUE,
   return 1 if two samples are identical, otherwise, return 0
   */
  int compareWithSample(const SeqRegions& sequence2,
                        cmaple::PositionType seq_length,
                        const Alignment* aln,
                        const bool check_ident_only = false) const;

  /**
   Check if the current regions and regions2 represent the same partial
   likelihoods or not -> be used to stop traversing the tree further for
   updating partial likelihoods
   */
  bool areDiffFrom(const std::unique_ptr<SeqRegions>& regions2,
                   cmaple::PositionType seq_length,
                   cmaple::StateType num_states,
                   const cmaple::Params& params) const;
    
    /**
     Integrate mutations at a branch to a likelihood vector
     This regions is the input  likelihood vector
     @param mutations the vector of mutations
     @param aln the alignment
     @param inverse True to inverse the mutations
     @return the output regions
     @throw std::logic\_error if unexpected values/behaviors found during the
     operations
     */
    template <const cmaple::StateType num_states>
    auto integrateMutations(
                         std::unique_ptr<SeqRegions>& mutations,
                         const Alignment* aln,
                         const bool inverse = false) const
        -> std::unique_ptr<SeqRegions>;
    
    /**
     Merge two local references
     This regions is the top/new one
     @param mutations_2 the second local reference
     @param aln the alignment
     @param downward True if merging different sides
     @return the output regions
     @throw std::logic\_error if unexpected values/behaviors found during the
     operations
     */
    template <const cmaple::StateType num_states>
    auto mergeTwoRefs(std::unique_ptr<SeqRegions>& mutations_2,
                      const Alignment* aln, const cmaple::RealNumType threshold_prob,
                      const bool downward = false) const
        -> std::unique_ptr<SeqRegions>;
    
    /**
     Flip all mutations
     */
    template <const cmaple::StateType num_states>
    void flipMutations();
    
    /**
     Check if this likelihood vector contains at least N mutations
     @param min_mut the number of mutations
     @return TRUE if this likelihood vector contains at least N mutations
     */
    template <const cmaple::StateType num_states>
    auto containAtLeastNMuts(const int min_mut) const
        -> bool;

  /**
   Merge two likelihood vectors, one from above and one from below
   This regions is the upper regions
   @param lower_regions the lower regions
   @param upper_plength the length of the upper branch
   @param lower_plength the length of the lower branch
   @param merged_regions the output regions by merging the upper and lower
   regions
   @param aln the alignment
   @param model the model of evolution
   @param threshold the threshold for approximation

   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void mergeUpperLower(std::unique_ptr<SeqRegions>& merged_regions,
                       cmaple::RealNumType upper_plength,
                       const SeqRegions& lower_regions,
                       cmaple::RealNumType lower_plength,
                       const Alignment* aln,
                       const ModelBase* model,
                       const cmaple::RealNumType threshold) const;

  /**
   Merge two lower likelihood vectors
   This regions is one of the two lower regions
   @param regions2 the other lower regions
   @param plength1, plength2 the lengths of the lower branches
   @param merged_regions the output regions by merging the upper and lower
   regions
   @param aln the alignment
   @param model the model of evolution
   @param threshold the threshold for approximation
   @param return_log_lh TRUE to return the log likelihood
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  cmaple::RealNumType mergeTwoLowers(
      std::unique_ptr<SeqRegions>& merged_regions,
      const cmaple::RealNumType plength1,
      const SeqRegions& regions2,
      const cmaple::RealNumType plength2,
      const Alignment* aln,
      const ModelBase* model,
      const RealNumType* const cumulative_rate,
      const cmaple::RealNumType threshold_prob,
      const bool return_log_lh = false) const;

  /**
   Calculate site-lh contributions by merging two lower likelihood vectors
   This regions is one of the two lower regions
   @param regions2 the other lower regions
   @param plength1, plength2 the lengths of the lower branches
   @param merged_regions the output regions by merging the upper and lower
   regions
   @param aln the alignment
   @param model the model of evolution
   @param threshold the threshold for approximation
   @param return_log_lh TRUE to return the log likelihood
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  cmaple::RealNumType calculateSiteLhContributions(
      std::vector<cmaple::RealNumType>& site_lh_contributions,
      std::unique_ptr<SeqRegions>& merged_regions,
      const cmaple::RealNumType plength1,
      const SeqRegions& regions2,
      const cmaple::RealNumType plength2,
      const Alignment* aln,
      const ModelBase* model,
      const RealNumType* const cumulative_rate,
      const cmaple::RealNumType threshold_prob) const;

  /**
   Compute total lh/upper left_right for root node
   @param blength the branch length; (-1 by default).

   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void computeTotalLhAtRoot(std::unique_ptr<SeqRegions>& total_lh,
                            const ModelBase* model,
                            cmaple::RealNumType blength = -1) const;

  /**
   Compute the likelihood by merging the lower lh with root frequencies
   */
  template <const cmaple::StateType num_states>
  cmaple::RealNumType computeAbsoluteLhAtRoot(
      Alignment* aln,
      const ModelBase* model,
      const std::vector<std::vector<PositionType>>& cumulative_base);

  /**
   Compute the site likelihood at root by merging the lower lh with root
   frequencies
   */
  template <const cmaple::StateType num_states>
  cmaple::RealNumType computeSiteLhAtRoot(
      std::vector<cmaple::RealNumType>& site_lh_contributions,
      const ModelBase* model,
      const std::vector<std::vector<PositionType>>& cumulative_base);

  /**
   Convert an entry 'O' into a normal nucleotide if its probability dominated
   others
   */
  static cmaple::StateType simplifyO(cmaple::RealNumType* const partial_lh,
                                     cmaple::StateType ref_state,
                                     cmaple::StateType num_states,
                                     cmaple::RealNumType threshold);

  /*
   Add an simplified O seqregion to the current vectors of seqregions
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  static void addSimplifiedO(const StateType ref_state,
                             const cmaple::PositionType end_pos,
                             SeqRegion::LHType& new_lh,
                             const Alignment* aln,
                             const cmaple::RealNumType threshold_prob,
                             SeqRegions& merged_regions);

  /**
   Compare two regions
   */
  bool operator==(const SeqRegions& seqregions_1) const;
};

/**
 MergeUpperLower case N with O
 @throw std::logic\_error if unexpected values/behaviors found during the
 operations
 */
template <const cmaple::StateType num_states>
void merge_N_O(const cmaple::RealNumType lower_plength,
               const SeqRegion& reg_o,
               const ModelBase* model,
               const cmaple::PositionType end_pos,
               SeqRegions& merged_target);

/**
 MergeUpperLower case N_RACGT

 @throw std::logic\_error if unexpected values/behaviors found during the
 operations
 */
void merge_N_RACGT(const SeqRegion& reg_racgt,
                   const cmaple::RealNumType lower_plength,
                   const cmaple::PositionType end_pos,
                   const cmaple::RealNumType threshold_prob,
                   SeqRegions& merged_regions);

/**
 MergeUpperLower case O_N

 @throw std::logic\_error if unexpected values/behaviors found during the
 operations
 */
template <const cmaple::StateType num_states>
void merge_O_N(const SeqRegion& reg_o,
               const cmaple::RealNumType upper_plength,
               const cmaple::PositionType end_pos,
               const ModelBase* model,
               SeqRegions& merged_regions);

/**
 MergeUpperLower case RACGT_N

 @throw std::logic\_error if unexpected values/behaviors found during the
 operations
 */
void merge_RACGT_N(const SeqRegion& reg_n,
                   const cmaple::RealNumType upper_plength,
                   const cmaple::PositionType end_pos,
                   const cmaple::RealNumType threshold_prob,
                   SeqRegions& merged_regions);

/**
 MergeUpperLower case Zero_Distance

 @throw std::logic\_error if unexpected values/behaviors found during the
 operations
 */
bool merge_Zero_Distance(const SeqRegion& seq1_region,
                         const SeqRegion& seq2_region,
                         const cmaple::RealNumType total_blength_1,
                         const cmaple::RealNumType total_blength_2,
                         const cmaple::PositionType end_pos,
                         const cmaple::RealNumType threshold_prob,
                         const cmaple::StateType num_states,
                         std::unique_ptr<SeqRegions>& merged_regions);

/**
 MergeUpperLower case O_ORACGT

 @throw std::logic\_error if unexpected values/behaviors found during the
 operations
 */
/*template <const cmaple::StateType num_states>
void merge_O_ORACGT(const SeqRegion& seq1_region,
                    const SeqRegion& seq2_region,
                    const cmaple::RealNumType total_blength_1,
                    const cmaple::RealNumType total_blength_2,
                    const cmaple::PositionType end_pos,
                    const cmaple::RealNumType threshold_prob,
                    const ModelBase* model,
                    const Alignment* aln,
                    SeqRegions& merged_regions);*/

/**
 MergeUpperLower case RACGT_O

 @throw std::logic\_error if unexpected values/behaviors found during the
 operations
 */
/*template <const cmaple::StateType num_states>
void merge_RACGT_O(const SeqRegion& seq2_region,
                   const cmaple::RealNumType total_blength_2,
                   const cmaple::PositionType end_pos,
                   SeqRegion::LHType& new_lh,
                   const cmaple::RealNumType threshold_prob,
                   const ModelBase* model,
                   const Alignment* aln,
                   SeqRegions& merged_regions);*/

/**
 MergeUpperLower case RACGT_RACGT

 @throw std::logic\_error if unexpected values/behaviors found during the
 operations
 */
/*template <const cmaple::StateType num_states>
void merge_RACGT_RACGT(const SeqRegion& seq2_region,
                       const cmaple::RealNumType total_blength_2,
                       const cmaple::PositionType end_pos,
                       SeqRegion::LHType& new_lh,
                       const ModelBase* model,
                       const Alignment* aln,
                       SeqRegions& merged_regions);*/

/**
 MergeUpperLower case RACGT_ORACGT

 @throw std::logic\_error if unexpected values/behaviors found during the
 operations
 */
/*template <const cmaple::StateType num_states>
void merge_RACGT_ORACGT(const SeqRegion& seq1_region,
                        const SeqRegion& seq2_region,
                        const cmaple::RealNumType total_blength_1,
                        const cmaple::RealNumType total_blength_2,
                        const cmaple::RealNumType upper_plength,
                        const cmaple::PositionType end_pos,
                        const cmaple::RealNumType threshold_prob,
                        const ModelBase* model,
                        const Alignment* aln,
                        SeqRegions& merged_regions);*/

/**
 MergeTwoLowers case N with O
 */
void merge_N_O_TwoLowers(const SeqRegion& seq2_region,
                         const cmaple::PositionType end_pos,
                         const cmaple::RealNumType plength2,
                         SeqRegions& merged_regions);

/**
 MergeTwoLowers case N_RACGT

 @throw std::logic\_error if unexpected values/behaviors found during the
 operations
 */
void merge_N_RACGT_TwoLowers(const SeqRegion& seq2_region,
                             const cmaple::PositionType end_pos,
                             const cmaple::RealNumType plength2,
                             const cmaple::RealNumType threshold_prob,
                             SeqRegions& merged_regions);

/**
 MergeTwoLowers case identicalRACGT
 @throw std::logic\_error if unexpected values/behaviors found during the
 operations
 */
void merge_identicalRACGT_TwoLowers(const SeqRegion& seq1_region,
                                    const cmaple::PositionType end_pos,
                                    cmaple::RealNumType total_blength_1,
                                    cmaple::RealNumType total_blength_2,
                                    const cmaple::RealNumType blength_1,
                                    const cmaple::RealNumType blength_2,
                                    const cmaple::PositionType pos,
                                    const cmaple::RealNumType threshold_prob,
                                    const ModelBase* model,
                                    const RealNumType* const cumulative_rate,
                                    cmaple::RealNumType& log_lh,
                                    SeqRegions& merged_regions,
                                    const bool return_log_lh);

/**
 MergeTwoLowers case O_O

 @throw std::logic\_error if unexpected values/behaviors found during the
 operations
 */
/*template <const cmaple::StateType num_states>
bool merge_O_O_TwoLowers(const SeqRegion& seq2_region,
                         cmaple::RealNumType total_blength_2,
                         const cmaple::PositionType end_pos,
                         const Alignment* aln,
                         const ModelBase* model,
                         const cmaple::RealNumType threshold_prob,
                         cmaple::RealNumType& log_lh,
                         SeqRegion::LHType& new_lh,
                         std::unique_ptr<SeqRegions>& merged_regions,
                         const bool return_log_lh);*/

/**
 MergeTwoLowers case O_RACGT

 @throw std::logic\_error if unexpected values/behaviors found during the
 operations
 */
/*template <const cmaple::StateType num_states>
bool merge_O_RACGT_TwoLowers(const SeqRegion& seq2_region,
                             cmaple::RealNumType total_blength_2,
                             const cmaple::PositionType end_pos,
                             const Alignment* aln,
                             const ModelBase* model,
                             const cmaple::RealNumType threshold_prob,
                             cmaple::RealNumType& log_lh,
                             SeqRegion::LHType& new_lh,
                             cmaple::RealNumType& sum_lh,
                             std::unique_ptr<SeqRegions>& merged_regions,
                             const bool return_log_lh);*/

/**
 * Initialize partial lh vector when merging two regions
 * @param ori_type the original state
 * @param total_blength total branch length
 * @param model the model
 * @param input_lh_vec the input likelihood vector
 * @param inverse TRUE if going upward a branch
 * @return a likelihood vector
*/
template <const cmaple::StateType num_states>
SeqRegion::LHType initLhVecMerge(
    const StateType ori_type, const RealNumType total_blength,
    const ModelBase* model, const bool inverse,
    const SeqRegion::LHPtrType& input_lh_vec = nullptr);

/**
 MergeTwoLowers case O_ORACGT

 @throw std::logic\_error if unexpected values/behaviors found during the
 operations
 */
/*template <const cmaple::StateType num_states>
bool merge_O_ORACGT_TwoLowers(const SeqRegion& seq1_region,
                              const SeqRegion& seq2_region,
                              cmaple::RealNumType total_blength_1,
                              cmaple::RealNumType total_blength_2,
                              const cmaple::PositionType end_pos,
                              const Alignment* aln,
                              const ModelBase* model,
                              const cmaple::RealNumType threshold_prob,
                              cmaple::RealNumType& log_lh,
                              std::unique_ptr<SeqRegions>& merged_regions,
                              const bool return_log_lh);*/

/**
 MergeTwoLowers case RACGT_O

 @throw std::logic\_error if unexpected values/behaviors found during the
 operations
 */
/*template <const cmaple::StateType num_states>
bool merge_RACGT_O_TwoLowers(const SeqRegion& seq2_region,
                             cmaple::RealNumType total_blength_2,
                             const cmaple::PositionType end_pos,
                             const Alignment* aln,
                             const ModelBase* model,
                             const cmaple::RealNumType threshold_prob,
                             SeqRegion::LHType& new_lh,
                             cmaple::RealNumType& log_lh,
                             std::unique_ptr<SeqRegions>& merged_regions,
                             const bool return_log_lh);*/

/**
 MergeTwoLowers case RACGT_RACGT

 @throw std::logic\_error if unexpected values/behaviors found during the
 operations
 */
/*template <const cmaple::StateType num_states>
bool merge_RACGT_RACGT_TwoLowers(const SeqRegion& seq2_region,
                                 cmaple::RealNumType total_blength_2,
                                 const cmaple::PositionType end_pos,
                                 const Alignment* aln,
                                 const ModelBase* model,
                                 const cmaple::RealNumType threshold_prob,
                                 SeqRegion::LHType& new_lh,
                                 cmaple::RealNumType& sum_lh,
                                 cmaple::RealNumType& log_lh,
                                 std::unique_ptr<SeqRegions>& merged_regions,
                                 const bool return_log_lh);*/

/**
 MergeTwoLowers case RACGT_ORACGT

 @throw std::logic\_error if unexpected values/behaviors found during the
 operations
 */
/*template <const cmaple::StateType num_states>
bool merge_RACGT_ORACGT_TwoLowers(const SeqRegion& seq1_region,
                                  const SeqRegion& seq2_region,
                                  cmaple::RealNumType total_blength_1,
                                  cmaple::RealNumType total_blength_2,
                                  const cmaple::PositionType end_pos,
                                  const Alignment* aln,
                                  const ModelBase* model,
                                  const cmaple::RealNumType threshold_prob,
                                  cmaple::RealNumType& log_lh,
                                  std::unique_ptr<SeqRegions>& merged_regions,
                                  const bool return_log_lh);*/

/**
 MergeTwoLowers case ORACGT_ORACGT

 @throw std::logic\_error if unexpected values/behaviors found during the
 operations
 */
template <const cmaple::StateType num_states>
auto do_merge_ORACGT_ORACGT(
                              const SeqRegion& seq1_region,
                              const SeqRegion& seq2_region,
                              const cmaple::RealNumType total_blength_1,
                              const cmaple::RealNumType total_blength_2,
                              const cmaple::RealNumType blength_1,
                              const cmaple::PositionType end_pos,
                              const Alignment* aln,
                              const ModelBase* model,
                              const cmaple::RealNumType threshold_prob,
                              std::unique_ptr<SeqRegions>& merged_regions,
                              const bool inverse,
                              cmaple::RealNumType& total_prob) -> bool;

/**
 MergeTwoLowers case ORACGT_ORACGT

 @throw std::logic\_error if unexpected values/behaviors found during the
 operations
 */
template <const cmaple::StateType num_states>
bool merge_ORACGT_ORACGT(const SeqRegion& seq1_region,
                                   const SeqRegion& seq2_region,
                                   const cmaple::RealNumType total_blength_1,
                                   const cmaple::RealNumType total_blength_2,
                                   const cmaple::RealNumType blength_1,
                                   const cmaple::PositionType end_pos,
                                   const Alignment* aln,
                                   const ModelBase* model,
                                   const cmaple::RealNumType threshold_prob,
                                   std::unique_ptr<SeqRegions>& merged_regions,
                                   const bool inverse,
                                   const bool return_log_lh,
                                   cmaple::RealNumType& total_factor);

/**
 MergeTwoLowers case notN_notN

 @throw std::logic\_error if unexpected values/behaviors found during the
 operations
 */
template <const cmaple::StateType num_states>
bool merge_notN_notN_TwoLowers(const SeqRegion& seq1_region,
                               const SeqRegion& seq2_region,
                               const cmaple::RealNumType plength1,
                               const cmaple::RealNumType plength2,
                               const cmaple::PositionType end_pos,
                               const cmaple::PositionType pos,
                               const Alignment* aln,
                               const ModelBase* model,
                               const RealNumType* const cumulative_rate,
                               const cmaple::RealNumType threshold_prob,
                               cmaple::RealNumType& log_lh,
                               cmaple::RealNumType& total_factor,
                               std::unique_ptr<SeqRegions>& merged_regions,
                               const bool return_log_lh);

template <const StateType num_states>
auto updateLHwithModel(const ModelBase* model,
                       const SeqRegion::LHType& prior,
                       SeqRegion::LHType& posterior,
                       const RealNumType total_blength,
                       const PositionType pos) -> RealNumType {
  assert(model);
  
  bool negative_total = false;
  const RealNumType equal_prob = 1.0 / num_states;
  RealNumType sum_lh = 0;
  for (StateType i = 0; i < num_states; ++i) {
    RealNumType tot = 0;
    if (total_blength > 0)  // TODO: avoid
    {

      const RealNumType* mutation_mat_row = model->getMutationMatrixRow(i, pos);
      tot += dotProduct<num_states>(&(prior)[0], mutation_mat_row);

      tot *= total_blength;
    }

    tot += prior[i];
      
      // special treatment, if tot is negative
      // then return equal probability
      if (tot < 0)
          negative_total = true;
      if (negative_total)
          tot = equal_prob;
      
      posterior[i] = tot * model->getRootFreq(i);
    sum_lh += posterior[i];
  }
  return sum_lh;
}

template <const StateType num_states>
auto updateLHwithMat(const RealNumType* mat_row,
                     const SeqRegion::LHType& prior,
                     SeqRegion::LHType& posterior,
                     const RealNumType total_blength) -> RealNumType {
  assert(mat_row);
  RealNumType sum_lh = 0;
  bool negative_tot = false;
  const RealNumType equal_prob = 1.0 / num_states;
  for (StateType i = 0; i < num_states; ++i, mat_row += num_states) {
    RealNumType tot = 0;
    tot += dotProduct<num_states>(&(prior)[0], mat_row);
    tot *= total_blength;
    tot += prior[i];
    
    // record negative tot
    if (tot < 0)
        negative_tot = true;
      
    posterior[i] = tot;
    sum_lh += tot;
  }
    
    // if negative tot found -> return a vector of equal probabilities
    if (negative_tot)
    {
        for (auto i = 0; i < num_states; ++i)
            posterior[i] = equal_prob;
    }
    
  return sum_lh;
}

template <const StateType num_states>
auto updateMultLHwithMat(const RealNumType* mat_row,
                         const SeqRegion::LHType& prior,
                         SeqRegion::LHType& posterior,
                         const RealNumType total_blength) -> RealNumType {
  assert(mat_row);
  RealNumType sum_lh = 0;
  for (StateType i = 0; i < num_states; ++i, mat_row += num_states) {
    RealNumType tot = 0;
    if (total_blength > 0)  // TODO: avoid
    {
      tot += dotProduct<num_states>(&(prior)[0], mat_row);
      tot *= total_blength;
    }
    tot += prior[i];
    posterior[i] *= tot;
    sum_lh += posterior[i];
  }
  return sum_lh;
}

template <const StateType num_states>
auto SeqRegions::integrateMutations(
                                    std::unique_ptr<SeqRegions>& mutations,
                                    const Alignment* aln,
                                    const bool inverse) const
        -> std::unique_ptr<SeqRegions>
{
    // if there is no mutations clone and return the original likelihood vector
    if (!mutations || !mutations->size())
    {
        return cmaple::make_unique<SeqRegions>(SeqRegions::clone(*this));
    }
            
    // init merged_regions
    std::unique_ptr<SeqRegions> output_regions = cmaple::make_unique<SeqRegions>();
    
    assert(mutations->size() > 0);
    assert(aln);
    
  // init variables
  PositionType pos = 0;
  const SeqRegions& seq_regions = *this;
  const SeqRegions& mutation_regions = *mutations;
  size_t iseq = 0;
  size_t imut = 0;
  const PositionType seq_length = static_cast<PositionType>(aln->ref_seq.size());
  const RealNumType threshold_prob = 1e-8;

  // avoid realloc of vector data (minimize memory footprint)
    output_regions->reserve(countSharedSegments(
        mutation_regions, static_cast<size_t>(seq_length)));
    
#ifdef DEBUG
  // remember capacity (may be more than we 'reserved')
  const size_t max_elements = output_regions->capacity();
#endif

  while (pos < seq_length) {
    PositionType end_pos;

    // get the next shared segment in the two sequences
    cmaple::SeqRegions::getNextSharedSegment(pos, seq_regions, mutation_regions,
                                             iseq, imut, end_pos);
    const auto* const seq_region = &seq_regions[iseq];
    const auto* const mutation = &mutation_regions[imut];
      
      // extract the new state of the mutation
      // according to the direction
      cmaple::StateType mut_new_state = inverse ?
        mutation->prev_state : mutation->type;
      
      // extract states of seq_region
      StateType region_type = seq_region->type;
      StateType region_prev_state = seq_region->prev_state;
      
      // case 1: seq_region = N
      if (region_type == TYPE_N)
      {
          // force moving to the end of this region
          while (end_pos < seq_region->position)
          {
              // update pos
              pos = end_pos + 1;
              
              // move to the next share segment
              cmaple::SeqRegions::getNextSharedSegment(
                    pos, seq_regions, mutation_regions,
                    iseq, imut, end_pos);
          }
          assert(end_pos == seq_region->position);
      }
      // case 2: seq_region = A/C/G/T
      else if (region_type < num_states)
      {
          // if a mutation occurs at this position,
          // modify the type (i.e., current state) and the previous state
          if (mutation->type < num_states)
          {
              // if the current state is identical to the new state of the mutation
              // change it to R
              if (region_type == mut_new_state)
              {
                  region_type = TYPE_R;
                  region_prev_state = TYPE_N;
              }
              // otherwise, keep the current state but update the previous state
              else
              {
                  region_prev_state = mut_new_state;
              }
          }
      }
      // case 3: seq_region = R
      else if (region_type == TYPE_R)
      {
          // if a mutation occurs at this position,
          // modify the type (i.e., current state) and the previous state
          if (mutation->type < num_states)
          {
              // set the previous state of the newly-added region
              // as the new state of the mutation
              region_prev_state = mut_new_state;
              
              // extract the previous state of the mutation
              // according to the direction
              cmaple::StateType mut_prev_state = inverse ?
                mutation->type : mutation->prev_state;
              
              // the new state of the newly-added region
              // is R which is also the previous state of the mutation
              region_type = mut_prev_state;
          }
      }
      // case 4: seq_region = O
      else
      {
          // if a mutation occurs at this position,
          // modify the previous state
          if (mutation->type < num_states)
          {
              region_prev_state = mut_new_state;
          }
      }
      
      // add the region the the output
      if (region_type == cmaple::TYPE_O) {
          // first add seq_region to the output lh vector
          output_regions->push_back(SeqRegion::clone(*seq_region));
          // extract the newly added seq_region
          SeqRegion& last_region = output_regions->back();
          last_region.prev_state = region_prev_state;
      } else {
        // add a new region and try to merge consecutive R regions together
          cmaple::SeqRegions::addNonConsecutiveRRegion(*output_regions, region_type,
                region_prev_state, seq_region->plength_observation2node,
                seq_region->plength_observation2root, end_pos, threshold_prob);
      }

    // update pos
    pos = end_pos + 1;
  }

#ifdef DEBUG
  // ensure we did the correct reserve, otherwise it was
  // a pessimization
  assert(output_regions->capacity() == max_elements);
#endif
    return output_regions;
}

template <const StateType num_states>
auto SeqRegions::mergeTwoRefs(std::unique_ptr<SeqRegions>& mutations_2,
                            const Alignment* aln, const cmaple::RealNumType threshold_prob,
                            const bool downward) const -> std::unique_ptr<SeqRegions>
{
    assert(aln);
    assert(size());
    assert(mutations_2 && mutations_2->size());
            
    // init merged_regions
    std::unique_ptr<SeqRegions> output_regions = cmaple::make_unique<SeqRegions>();
    
  // init variables
  PositionType pos = 0;
  const SeqRegions& seq_regions_1 = *this;
  const SeqRegions& seq_regions_2 = *mutations_2;
  size_t iseq = 0;
  size_t imut = 0;
  const PositionType seq_length = static_cast<PositionType>(aln->ref_seq.size());

  // avoid realloc of vector data (minimize memory footprint)
    output_regions->reserve(countSharedSegments(
        seq_regions_2, static_cast<size_t>(seq_length)));
    
#ifdef DEBUG
  // remember capacity (may be more than we 'reserved')
  const size_t max_elements = output_regions->capacity();
#endif

  while (pos < seq_length) {
    PositionType end_pos;

    // get the next shared segment in the two sequences
    cmaple::SeqRegions::getNextSharedSegment(pos, seq_regions_1, seq_regions_2,
                                             iseq, imut, end_pos);
    const auto* const seq_region_1 = &seq_regions_1[iseq];
    const auto* const seq_region_2 = &seq_regions_2[imut];
      
    // if seq_region_2 is a mutation
    if (seq_region_2->type < num_states)
    {
        // if seq_region_1 is also a mutation
        if (seq_region_1->type < num_states)
        {
            StateType from_state = seq_region_1->prev_state;
            StateType to_state = seq_region_1->type;
            if (downward)
            {
                from_state = seq_region_1->type;
                to_state = seq_region_1->prev_state;
            }
            
            // validate the two mutations
            if (to_state != seq_region_2->prev_state)
            {
                std::cout << "WARNING: inconsistent mutations " << std::endl;
            }
            
            // if the two mutations don't cancel out each other
            // record them
            if (from_state != seq_region_2->type)
            {
                // clone the mutation from the ref 2
                output_regions->push_back(SeqRegion::clone(*seq_region_2));
                
                // update the previous state of the newly added mutation
                SeqRegion& last_region = output_regions->back();
                last_region.prev_state = from_state;
            }
            // otherwise, add an R to make sure continue regions
            else
            {
                cmaple::SeqRegions::addNonConsecutiveRRegion(*output_regions, TYPE_R,
                                        TYPE_N, -1, -1, end_pos, threshold_prob);
            }
        }
        // otherwise, simply add seq_region_2
        else
        {
            output_regions->push_back(SeqRegion::clone(*seq_region_2));
        }
    }
    else
    {
        // if seq_region_1 is a mutation
        if (seq_region_1->type < num_states)
        {
            // add mutation from the reference 1
            output_regions->push_back(SeqRegion::clone(*seq_region_1));
            // swap the direction if needed
            if (downward)
            {
                SeqRegion& last_region = output_regions->back();
                const StateType ori_state = last_region.type;
                last_region.type = last_region.prev_state;
                last_region.prev_state = ori_state;
            }
        }
        // otherwise, both are not mutations
        // add an R region
        else
        {
            cmaple::SeqRegions::addNonConsecutiveRRegion(*output_regions, TYPE_R,
                                    TYPE_N, -1, -1, end_pos, threshold_prob);
        }
    }

    // update pos
    pos = end_pos + 1;
  }

#ifdef DEBUG
  // ensure we did the correct reserve, otherwise it was
  // a pessimization
  assert(output_regions->capacity() == max_elements);
#endif
    return output_regions;
}

template <const StateType num_states>
void merge_N_O(const RealNumType lower_plength,
               const SeqRegion& reg_o,
               const ModelBase* model,
               const PositionType end_pos,
               SeqRegions& merged_target) {
  assert(reg_o.type == TYPE_O);
  assert(model);
    
  RealNumType total_blength = lower_plength;
  if (reg_o.plength_observation2node >= 0) {
    total_blength = reg_o.plength_observation2node +
                    (lower_plength > 0 ? lower_plength : 0);
  }
  auto new_lh =
      cmaple::make_unique<SeqRegion::LHType>();  // = new RealNumType[num_states];
  RealNumType sum_lh = updateLHwithModel<num_states>(model, *reg_o.likelihood,
                                                     (*new_lh), total_blength,
                                                     reg_o.position);
  // normalize the new partial likelihood
  normalize_arr(new_lh->data(), num_states, sum_lh);

  // add merged region into merged_regions
  merged_target.emplace_back(TYPE_O, end_pos, reg_o.prev_state , 0, 0, std::move(new_lh));
}

template <const StateType num_states>
void merge_O_N(const SeqRegion& reg_o,
               const RealNumType upper_plength,
               const PositionType end_pos,
               const ModelBase* model,
               SeqRegions& merged_regions) {
  /*
   NhanLT: total_blength may be wrong if upper_plength == -1 and
  reg_o.plength_observation2node > 0 RealNumType total_blength = upper_plength;
  if (reg_o.plength_observation2node > 0)
  {
    total_blength += reg_o.plength_observation2node;
  }*/
  assert(reg_o.type == TYPE_O);
  assert(model);

  RealNumType total_blength = -1;

  if (reg_o.plength_observation2node >= 0) {
    total_blength = reg_o.plength_observation2node;
    if (upper_plength > 0) {
      total_blength += upper_plength;
    }
  } else if (upper_plength > 0) {
    total_blength = upper_plength;
  }

  if (total_blength > 0) {
    auto new_lh = cmaple::make_unique<SeqRegion::LHType>();  // = new
    // RealNumType[num_states];
    RealNumType sum_lh = updateLHwithMat<num_states>(
        model->getTransposedMutationMatrix(end_pos), *(reg_o.likelihood), *new_lh, total_blength);

    // normalize the new partial likelihood
    normalize_arr(new_lh->data(), num_states, sum_lh);

    // add merged region into merged_regions
    merged_regions.emplace_back(TYPE_O, end_pos, reg_o.prev_state, 0, 0, std::move(new_lh));
  } else {
    // add merged region into merged_regions
    merged_regions.emplace_back(TYPE_O, end_pos, reg_o.prev_state, 0, 0, *(reg_o.likelihood));
  }
}

/*template <const StateType num_states>
void merge_O_ORACGT(const SeqRegion& seq1_region,
                    const SeqRegion& seq2_region,
                    const RealNumType total_blength_1,
                    const RealNumType total_blength_2,
                    const PositionType end_pos,
                    const RealNumType threshold_prob,
                    const ModelBase* model,
                    const Alignment* aln,
                    SeqRegions& merged_regions) {
  assert(seq1_region.type == TYPE_O);
  assert(seq2_region.type != TYPE_N);
  assert(model);
  assert(aln);
    
  auto new_lh =
      cmaple::make_unique<SeqRegion::LHType>();  // = new RealNumType[num_states];
  auto& new_lh_value = *new_lh;

  // if total_blength_1 > 0 => compute new partial likelihood
  if (total_blength_1 > 0) {
    updateLHwithMat<num_states>(model->getTransposedMutationMatrix(end_pos),
                                *(seq1_region.likelihood), *new_lh,
                                total_blength_1);
    // otherwise, clone the partial likelihood from seq1
  } else {
    *new_lh = *seq1_region.likelihood;
  }

  RealNumType sum_new_lh = 0;

  // seq1 = seq2 = O
  if (seq2_region.type == TYPE_O) {
    sum_new_lh = updateMultLHwithMat<num_states>(model->getMutationMatrix(end_pos),
                                                 *(seq2_region.likelihood),
                                                 *new_lh, total_blength_2);
  }
  // seq1 = "O" and seq2 = R or ACGT
  else {
    StateType seq2_state = seq2_region.type;
    if (seq2_state == TYPE_R) {
      seq2_state = aln->ref_seq[static_cast<std::vector<cmaple::StateType>::size_type>(end_pos)];
    }

    if (total_blength_2 > 0) {
      sum_new_lh += updateVecWithState<num_states>(
          new_lh_value.data(), seq2_state, 
          model->getTransposedMutationMatrixRow(seq2_state, end_pos),
          total_blength_2);
    } else {
      sum_new_lh += resetLhVecExceptState<num_states>(
          new_lh_value.data(), seq2_state, new_lh_value[seq2_state]);
    }
  }

  // normalize the new partial lh
  if (sum_new_lh == 0) {
    throw std::logic_error("Sum of new partital lh is zero.");
  }

  // normalize the new partial likelihood
  normalize_arr(new_lh->data(), num_states, sum_new_lh);
  cmaple::SeqRegions::addSimplifiedO(TYPE_N, end_pos, new_lh_value, aln, threshold_prob,
                                     merged_regions);
}*/

/*template <const StateType num_states>
void merge_RACGT_O(const SeqRegion& seq2_region,
                   const RealNumType total_blength_2,
                   const PositionType end_pos,
                   SeqRegion::LHType& new_lh,
                   const RealNumType threshold_prob,
                   const ModelBase* model,
                   const Alignment* aln,
                   SeqRegions& merged_regions) {
  assert(seq2_region.type == TYPE_O);
  assert(model);
  assert(aln);

  RealNumType sum_new_lh = updateMultLHwithMat<num_states>(
      model->getMutationMatrix(end_pos), *(seq2_region.likelihood), new_lh, total_blength_2);

  // normalize the new partial likelihood
  normalize_arr(new_lh.data(), num_states, sum_new_lh);
  cmaple::SeqRegions::addSimplifiedO(TYPE_N, end_pos, new_lh, aln, threshold_prob,
                                     merged_regions);
}

template <const StateType num_states>
void merge_RACGT_RACGT(const SeqRegion& seq2_region,
                       const RealNumType total_blength_2,
                       const PositionType end_pos,
                       SeqRegion::LHType& new_lh,
                       const ModelBase* model,
                       const Alignment* aln,
                       SeqRegions& merged_regions) {
  assert(seq2_region.type != TYPE_N && seq2_region.type != TYPE_O);
  assert(model);
  assert(aln);
    
  RealNumType sum_new_lh = 0;
  StateType seq2_state = seq2_region.type;

  if (seq2_state == TYPE_R) {
    seq2_state = aln->ref_seq[static_cast<std::vector<cmaple::StateType>::size_type>(end_pos)];
  }

  // TODO: this seems a weird operation on `new_lh_value` (since it was just
  // created anew and is all 00000)! please check this makes sense! CHECKED:
  // NHANLT: new_lh_value or new_lh is already initialized or updated in the
  // section "init or update new_lh/new_lh_value" above
  if (total_blength_2 > 0) {
    sum_new_lh = updateVecWithState<num_states>(
        new_lh.data(), seq2_state, 
        model->getTransposedMutationMatrixRow(seq2_state, end_pos), 
        total_blength_2);
  } else {
    sum_new_lh = resetLhVecExceptState<num_states>(new_lh.data(), seq2_state,
                                                   new_lh[seq2_state]);
  }

  // normalize the new partial likelihood
  normalize_arr(new_lh.data(), num_states, sum_new_lh);

  // add new region into the merged regions
  merged_regions.emplace_back(TYPE_O, end_pos, TYPE_N, 0, 0, std::move(new_lh));
}*/

/*template <const StateType num_states>
void merge_RACGT_ORACGT(const SeqRegion& seq1_region,
                        const SeqRegion& seq2_region,
                        const RealNumType total_blength_1,
                        const RealNumType total_blength_2,
                        const RealNumType upper_plength,
                        const PositionType end_pos,
                        const RealNumType threshold_prob,
                        const ModelBase* model,
                        const Alignment* aln,
                        SeqRegions& merged_regions) {
  assert(seq1_region.type != TYPE_N && seq1_region.type != TYPE_O);
  assert(seq2_region.type != TYPE_N);
  assert(model);
  assert(aln);
    
  StateType seq1_state = seq1_region.type;
  if (seq1_state == TYPE_R) {
    seq1_state = aln->ref_seq[static_cast<std::vector<cmaple::StateType>::size_type>(end_pos)];
  }

  auto new_lh =
      cmaple::make_unique<SeqRegion::LHType>();  // = new RealNumType[num_states];
  auto& new_lh_value = *new_lh;

  // init or update new_lh/new_lh_value
  if (seq1_region.plength_observation2root >= 0) {
    RealNumType length_to_root = seq1_region.plength_observation2root;
    if (upper_plength > 0) {
      length_to_root += upper_plength;
    }
    SeqRegion::LHType root_vec;
    
    /*memcpy(root_vec.data(), model->root_freqs,
           sizeof(RealNumType) * num_states);*/
    /*memcpy(root_vec.data(), model->root_freqs, sizeof(SeqRegion::LHType));


    updateVecWithState<num_states>(root_vec.data(), seq1_state,
                                    model->getTransposedMutationMatrixRow(seq1_state, end_pos),
                                   seq1_region.plength_observation2node);

    updateLHwithMat<num_states>(model->getTransposedMutationMatrix(end_pos), 
                                root_vec, *new_lh, length_to_root);
  } else {
    if (total_blength_1 > 0) {
      const RealNumType* mutation_mat_row = model->getMutationMatrixRow(seq1_state, end_pos);
      setVecWithState<num_states>(new_lh_value.data(), seq1_state,
                                  mutation_mat_row, total_blength_1);
    } else {
      for (auto& v : new_lh_value) {
        v = 0;
      }
      new_lh_value[seq1_state] = 1;
    }
  }

  // seq1 = R/ACGT and seq2 = "O"
  if (seq2_region.type == TYPE_O) {
    merge_RACGT_O<num_states>(seq2_region, total_blength_2, end_pos, *new_lh,
                              threshold_prob, model, aln, merged_regions);
  }
  // seq1 = R/ACGT and different from seq2 = R/ACGT
  else {
    merge_RACGT_RACGT<num_states>(seq2_region, total_blength_2, end_pos,
                                  *new_lh, model, aln, merged_regions);
  }
}*/

template <const StateType num_states>
void SeqRegions::mergeUpperLower(std::unique_ptr<SeqRegions>& merged_regions,
                                 RealNumType upper_plength,
                                 const SeqRegions& lower_regions,
                                 RealNumType lower_plength,
                                 const Alignment* aln,
                                 const ModelBase* model,
                                 const RealNumType threshold_prob) const {
  assert(lower_regions.size() > 0);
  assert(model);
  assert(aln);
    
  // init variables
  PositionType pos = 0;
  const SeqRegions& seq1_regions = *this;
  const SeqRegions& seq2_regions = lower_regions;
  size_t iseq1 = 0;
  size_t iseq2 = 0;
  const PositionType seq_length = static_cast<PositionType>(aln->ref_seq.size());

  // init merged_regions
  if (merged_regions) {
    merged_regions->clear();
  } else {
    merged_regions = cmaple::make_unique<SeqRegions>();
  }

  // avoid realloc of vector data (minimize memory footprint)
  merged_regions->reserve(countSharedSegments(
      seq2_regions, static_cast<size_t>(seq_length)));  // avoid realloc of vector data
#ifdef DEBUG
  const size_t max_elements =
      merged_regions
          ->capacity();  // remember capacity (may be more than we 'reserved')
#endif

  while (pos < seq_length) {
    PositionType end_pos;

    // get the next shared segment in the two sequences
    cmaple::SeqRegions::getNextSharedSegment(pos, seq1_regions, seq2_regions,
                                             iseq1, iseq2, end_pos);
    const auto* const seq1_region = &seq1_regions[iseq1];
    const auto* const seq2_region = &seq2_regions[iseq2];
    const DoubleState s1s2 =
        (DoubleState(seq1_region->type) << 8) | seq2_region->type;

    // seq1_entry = 'N'
    // seq1_entry = 'N' and seq2_entry = 'N'
    if (s1s2 == NN) {
      merged_regions->emplace_back(TYPE_N, end_pos, TYPE_N);
      // seq1_entry = 'N' and seq2_entry = O/R/ACGT
      // seq1_entry = 'N' and seq2_entry = O
    } else if (s1s2 == NO) {
      merge_N_O<num_states>(lower_plength, *seq2_region, model, end_pos,
                            *merged_regions);
    }
    // seq1_entry = 'N' and seq2_entry = R/ACGT
    else if (seq1_region->type == TYPE_N) {
      merge_N_RACGT(*seq2_region, lower_plength, end_pos, threshold_prob,
                    *merged_regions);
    }
    // seq2_entry = 'N'
    // seq1_entry = 'O' and seq2_entry = N
    else if (s1s2 == ON) {
      merge_O_N<num_states>(*seq1_region, upper_plength, end_pos, model,
                            *merged_regions);
    }
    // seq2_entry = 'N' and seq1_entry = R/ACGT
    else if (seq2_region->type == TYPE_N) {
      merge_RACGT_N(*seq1_region, upper_plength, end_pos, threshold_prob,
                    *merged_regions);
    }
    // seq1_entry = seq2_entry = R/ACGT
    // todo: improve condition: a == b & a == RACGT
    else if (seq1_region->type == seq2_region->type &&
             (seq1_region->type < num_states || seq1_region->type == TYPE_R)) {
      // add a new region and try to merge consecutive R regions together
      addNonConsecutiveRRegion(*merged_regions, seq1_region->type, seq1_region->prev_state,
                               -1, -1, end_pos, threshold_prob);
    }
    // cases where the new genome list entry will likely be of type "O"
    else {
      RealNumType total_blength_1 = upper_plength;
      if (seq1_region->plength_observation2node >= 0) {
        total_blength_1 = seq1_region->plength_observation2node;
        if (upper_plength > 0) {
          total_blength_1 += upper_plength;
        }

        if (seq1_region->type != TYPE_O &&
            seq1_region->plength_observation2root >= 0) {
          total_blength_1 += seq1_region->plength_observation2root;
        }
      }

      RealNumType total_blength_2 = lower_plength;
      if (seq2_region->plength_observation2node >= 0) {
        total_blength_2 = seq2_region->plength_observation2node;
        if (lower_plength > 0) {
          total_blength_2 += lower_plength;
        }
      }

      // special cases when either total_blength_1 or total_blength_2 is zero
      if (merge_Zero_Distance(*seq1_region, *seq2_region, total_blength_1,
                              total_blength_2, end_pos, threshold_prob,
                              num_states, merged_regions)) {
        if (!merged_regions) {
          return;
        }
      }
      // seq1 and seq2 are ORACGT
      else
      {
          // Dummy variables
          RealNumType total_factor = 1;
          merge_ORACGT_ORACGT<num_states>(*seq1_region, *seq2_region,
               total_blength_1, total_blength_2, upper_plength, end_pos, aln, model,
               threshold_prob, merged_regions, true, false, total_factor);
      }
      /* ------ OLD-CODE REMOVED AFTER UPDATING TO MAPLE 0.7.5 --------------
      // seq1_entry = O and seq2_entry = O/R/ACGT
      else if (seq1_region->type == TYPE_O) {
        merge_O_ORACGT<num_states>(*seq1_region, *seq2_region, total_blength_1,
                                   total_blength_2, end_pos, threshold_prob,
                                   model, aln, *merged_regions);
      }
      // seq1_entry = R/ACGT and seq2_entry = O/R/ACGT
      else {
        merge_RACGT_ORACGT<num_states>(*seq1_region, *seq2_region,
                                       total_blength_1, total_blength_2,
                                       upper_plength, end_pos, threshold_prob,
                                       model, aln, *merged_regions);
      }*/

      // NHANLT: LOGS FOR DEBUGGING
      /*if (Params::getInstance().debug &&
      merged_regions->at(merged_regions->size()-1).type == TYPE_O)
      {
          SeqRegion::LHType lh =
      *merged_regions->at(merged_regions->size()-1).likelihood; std::cout <<
      "mergeUpLow " << pos << " " <<  std::setprecision(20) << lh[0] << " " <<
      lh[1] << " " << lh[2] << " " << lh[3] << " " << std::endl;
      }*/
    }

    // update pos
    pos = end_pos + 1;
  }

#ifdef DEBUG
  assert(merged_regions->capacity() ==
         max_elements);  // ensure we did the correct reserve, otherwise it was
  // a pessimization
#endif
}

/*template <const StateType num_states>
auto merge_O_O_TwoLowers(const SeqRegion& seq2_region,
                         RealNumType total_blength_2,
                         const PositionType end_pos,
                         const Alignment* aln,
                         const ModelBase* model,
                         const RealNumType threshold_prob,
                         RealNumType& log_lh,
                         SeqRegion::LHType& new_lh,
                         std::unique_ptr<SeqRegions>& merged_regions,
                         const bool return_log_lh) -> bool {
  assert(seq2_region.type == TYPE_O);
  assert(model);
  assert(aln);

  RealNumType sum_lh = updateMultLHwithMat<num_states>(
      model->getMutationMatrix(end_pos), *seq2_region.likelihood, new_lh, total_blength_2);

  if (sum_lh == 0) {
    merged_regions = nullptr;
    return false;
  }

  // normalize the new partial likelihood
  normalize_arr(new_lh.data(), num_states, sum_lh);
  cmaple::SeqRegions::addSimplifiedO(TYPE_N, end_pos, new_lh, aln, threshold_prob,
                                     *merged_regions);

  if (return_log_lh) {
    log_lh += log(sum_lh);
  }

  // no error
  return true;
}*/

/*template <const StateType num_states>
auto merge_O_RACGT_TwoLowers(const SeqRegion& seq2_region,
                             RealNumType total_blength_2,
                             const PositionType end_pos,
                             const Alignment* aln,
                             const ModelBase* model,
                             const RealNumType threshold_prob,
                             RealNumType& log_lh,
                             SeqRegion::LHType& new_lh,
                             RealNumType& sum_lh,
                             std::unique_ptr<SeqRegions>& merged_regions,
                             const bool return_log_lh) -> bool {
  assert(seq2_region.type != TYPE_N && seq2_region.type != TYPE_O);
  assert(model);
  assert(aln);
  StateType seq2_state = seq2_region.type;
  if (seq2_state == TYPE_R) {
    seq2_state = aln->ref_seq[static_cast<std::vector<cmaple::StateType>::size_type>(end_pos)];
  }

  if (total_blength_2 > 0) {
    sum_lh = updateVecWithState<num_states>(
        new_lh.data(), seq2_state,  
        model->getTransposedMutationMatrixRow(seq2_state, end_pos), 
        total_blength_2);

    // normalize new partial lh
    // normalize the new partial likelihood
    normalize_arr(new_lh.data(), num_states, sum_lh);
    cmaple::SeqRegions::addSimplifiedO(TYPE_N, end_pos, new_lh, aln, threshold_prob,
                                       *merged_regions);

    if (return_log_lh) {
      log_lh += log(sum_lh);
    }
  } else {
    if (new_lh[seq2_state] == 0) {
      merged_regions = nullptr;
      return false;
    }

    // add a new region and try to merge consecutive R regions together
    cmaple::SeqRegions::addNonConsecutiveRRegion(
        *merged_regions, seq2_region.type, seq2_region.prev_state,
        -1, -1, end_pos, threshold_prob);

    if (return_log_lh) {
      log_lh += log(new_lh[seq2_state]);
    }
  }
  // no error
  return true;
}*/

template <const cmaple::StateType num_states>
auto initLhVecMerge(
    const StateType ori_type, const RealNumType total_blength,
    const ModelBase* model, const bool inverse, const SeqRegion::LHPtrType& input_lh_vec) -> SeqRegion::LHType
{
    // dummy variables
    SeqRegion::LHType output_lh_vec;
    
    // case 1: ori_type is O
    if (ori_type == TYPE_O)
    {
        // special case: return the existed vector
        if (total_blength <= 0)
        {
            output_lh_vec = *input_lh_vec;
        }
        // normal case
        else
        {
            // used when merging upper and lower regions
            if (inverse)
            {
                updateLHwithMat<num_states>(model->transposed_mut_mat, *input_lh_vec, output_lh_vec,
                                            total_blength);
            }
            // used when mergin two lower regions
            else
            {
                updateLHwithMat<num_states>(model->mutation_mat, *input_lh_vec,
                                            output_lh_vec, total_blength);
            }
        }
    }
    // case 2: ori_type is RACGT
    else
    {
        // special case: return a default vector concentrated at the original state
        if (total_blength <= 0)
        {
            resetLhVecExceptState<num_states>(output_lh_vec.data(), ori_type, 1);
        }
        // normal case
        else
        {
            // used when merging upper and lower regions
            if (inverse)
            {
                RealNumType* mutation_mat_row =
                    model->mutation_mat + model->row_index[ori_type];
                setVecWithState<num_states>(output_lh_vec.data(), ori_type,
                                            mutation_mat_row, total_blength);
            }
            // used when mergin two lower regions
            else
            {
                RealNumType* transposed_mut_mat_row =
                model->transposed_mut_mat + model->row_index[ori_type];
                setVecWithState<num_states>(output_lh_vec.data(), ori_type,
                                            transposed_mut_mat_row, total_blength);
            }
            
            // special treatment in abnormal case
            // output_lh_vec[ori_type] is negative
            if (output_lh_vec[ori_type] < 0)
            {
                const RealNumType equal_prob = 1.0 / num_states;
                for (auto i = 0; i < num_states; ++i)
                    output_lh_vec[i] = equal_prob;
            }
        }
    }
    return output_lh_vec;
}


/*template <const StateType num_states>
auto merge_O_ORACGT_TwoLowers(const SeqRegion& seq1_region,
                              const SeqRegion& seq2_region,
                              RealNumType total_blength_1,
                              RealNumType total_blength_2,
                              const PositionType end_pos,
                              const Alignment* aln,
                              const ModelBase* model,
                              const RealNumType threshold_prob,
                              RealNumType& log_lh,
                              std::unique_ptr<SeqRegions>& merged_regions,
                              const bool return_log_lh) -> bool {
  assert(seq1_region.type == TYPE_O);
  assert(seq2_region.type != TYPE_N);
  assert(model);
  assert(aln);
    
  auto new_lh =
      cmaple::make_unique<SeqRegion::LHType>();  // = new RealNumType[num_states];
  RealNumType sum_lh = 0;
  if (total_blength_1 > 0) {
    updateLHwithMat<num_states>(model->getMutationMatrix(end_pos), *(seq1_region.likelihood),
                                *new_lh, total_blength_1);
    // otherwise, clone the partial likelihood from seq1
  } else {
    *new_lh = *seq1_region.likelihood;
  }

  // seq1_entry = O and seq2_entry = O
  if (seq2_region.type == TYPE_O) {
    auto ret = merge_O_O_TwoLowers<num_states>(
        seq2_region, total_blength_2, end_pos, aln, model, threshold_prob,
        log_lh, *new_lh, merged_regions, return_log_lh);
    return ret;
  }
  // seq1_entry = O and seq2_entry = R/ACGT
  else {
    auto ret = merge_O_RACGT_TwoLowers<num_states>(
        seq2_region, total_blength_2, end_pos, aln, model, threshold_prob,
        log_lh, *new_lh, sum_lh, merged_regions, return_log_lh);
    return ret;
  }

  // no error
  return true;
}*/

/*template <const StateType num_states>
auto merge_RACGT_O_TwoLowers(const SeqRegion& seq2_region,
                             RealNumType total_blength_2,
                             const PositionType end_pos,
                             const Alignment* aln,
                             const ModelBase* model,
                             const RealNumType threshold_prob,
                             SeqRegion::LHType& new_lh,
                             RealNumType& log_lh,
                             std::unique_ptr<SeqRegions>& merged_regions,
                             const bool return_log_lh) -> bool {
  assert(seq2_region.type == TYPE_O);
  assert(model);
  assert(aln);

  RealNumType sum_lh = updateMultLHwithMat<num_states>(
      model->getMutationMatrix(end_pos), *(seq2_region.likelihood), new_lh, total_blength_2);

  if (sum_lh == 0) {
    merged_regions = nullptr;
    return false;
  }

  // normalize the new partial likelihood
  normalize_arr(new_lh.data(), num_states, sum_lh);
  cmaple::SeqRegions::addSimplifiedO(TYPE_N, end_pos, new_lh, aln, threshold_prob,
                                     *merged_regions);

  if (return_log_lh) {
    log_lh += log(sum_lh);
  }

  // no error
  return true;
}*/

/*template <const StateType num_states>
auto merge_RACGT_RACGT_TwoLowers(const SeqRegion& seq2_region,
                                 RealNumType total_blength_2,
                                 const PositionType end_pos,
                                 const Alignment* aln,
                                 const ModelBase* model,
                                 const RealNumType threshold_prob,
                                 SeqRegion::LHType& new_lh,
                                 RealNumType& sum_lh,
                                 RealNumType& log_lh,
                                 std::unique_ptr<SeqRegions>& merged_regions,
                                 const bool return_log_lh) -> bool {
  assert(seq2_region.type != TYPE_N && seq2_region.type != TYPE_O);
  assert(model);
  assert(aln);
    
  StateType seq2_state = seq2_region.type;
  if (seq2_state == TYPE_R) {
    seq2_state = aln->ref_seq[static_cast<std::vector<cmaple::StateType>::size_type>(end_pos)];
  }

  if (total_blength_2 > 0) {
    sum_lh += updateVecWithState<num_states>(
        new_lh.data(), seq2_state,  
        model->getTransposedMutationMatrixRow(seq2_state, end_pos), 
        total_blength_2);

    // normalize the new partial likelihood
    normalize_arr(new_lh.data(), num_states, sum_lh);
    cmaple::SeqRegions::addSimplifiedO(TYPE_N, end_pos, new_lh, aln, threshold_prob,
                                       *merged_regions);

    if (return_log_lh) {
      log_lh += log(sum_lh);
    }
  } else {
    // add a new region and try to merge consecutive R regions together
    cmaple::SeqRegions::addNonConsecutiveRRegion(
        *merged_regions, seq2_region.type, seq2_region.prev_state,
        -1, -1, end_pos, threshold_prob);

    if (return_log_lh) {
      log_lh += log(new_lh[seq2_state]);
    }
  }

  // no error
  return true;
}*/

/*template <const StateType num_states>
auto merge_RACGT_ORACGT_TwoLowers(const SeqRegion& seq1_region,
                                  const SeqRegion& seq2_region,
                                  RealNumType total_blength_1,
                                  RealNumType total_blength_2,
                                  const PositionType end_pos,
                                  const Alignment* aln,
                                  const ModelBase* model,
                                  const RealNumType threshold_prob,
                                  RealNumType& log_lh,
                                  std::unique_ptr<SeqRegions>& merged_regions,
                                  const bool return_log_lh) -> bool {
  assert(seq1_region.type != TYPE_O && seq1_region.type != TYPE_N);
  assert(seq2_region.type != TYPE_N);
  assert(model);
  assert(aln);
    
  StateType seq1_state = seq1_region.type;
  if (seq1_state == TYPE_R) {
    seq1_state = aln->ref_seq[static_cast<std::vector<cmaple::StateType>::size_type>(end_pos)];
  }

  auto new_lh =
      cmaple::make_unique<SeqRegion::LHType>();  // = new RealNumType[num_states];
  // auto& new_lh_value = *new_lh;
  RealNumType sum_lh = 0;

  if (total_blength_1 > 0) {
    setVecWithState<num_states>(new_lh->data(), seq1_state,
                                 model->getTransposedMutationMatrixRow(seq1_state, end_pos), 
                                 total_blength_1);
  } else {
    resetLhVecExceptState<num_states>(new_lh->data(), seq1_state, 1);
  }

  // seq1_entry = R/ACGT and seq2_entry = O
  if (seq2_region.type == TYPE_O) {
    return merge_RACGT_O_TwoLowers<num_states>(
        seq2_region, total_blength_2, end_pos, aln, model, threshold_prob,
        *new_lh, log_lh, merged_regions, return_log_lh);
  }

  // otherwise, seq1_entry = R/ACGT and seq2_entry = R/ACGT
  return merge_RACGT_RACGT_TwoLowers<num_states>(
      seq2_region, total_blength_2, end_pos, aln, model, threshold_prob,
      *new_lh, sum_lh, log_lh, merged_regions, return_log_lh);
}*/

template <const cmaple::StateType num_states>
auto do_merge_ORACGT_ORACGT(
                              const SeqRegion& seq1_region,
                              const SeqRegion& seq2_region,
                              const cmaple::RealNumType total_blength_1,
                              const cmaple::RealNumType total_blength_2,
                              const cmaple::RealNumType blength_1,
                              const cmaple::PositionType end_pos,
                              const Alignment* aln,
                              const ModelBase* model,
                              const cmaple::RealNumType threshold_prob,
                              std::unique_ptr<SeqRegions>& merged_regions,
                              const bool inverse,
                              cmaple::RealNumType& total_prob) -> bool
{
    // dummy variables
    SeqRegion::LHType seq1_lh_vec, seq2_lh_vec;
    
    // extract states of seq1 and the reference
    StateType seq1_state, ref_state;
    if (seq1_region.type == TYPE_R)
    {
        ref_state = seq2_region.prev_state;
        seq1_state = ref_state;
    }
    else
    {
       ref_state = seq1_region.prev_state;
       seq1_state = seq1_region.type;
    }
    
    // extract the lh vector for seq1,
    // taking into account the evolution along the branch
    // seq1_state is O
    if (seq1_state == TYPE_O)
    {
       if (total_blength_1 > 0)
       {
           assert(seq1_state == TYPE_O);
           seq1_lh_vec = initLhVecMerge<num_states>(TYPE_O,
                           total_blength_1, model, inverse, seq1_region.likelihood);
       }
       else
       {
           seq1_lh_vec = *seq1_region.likelihood;
       }
    }
    // seq1_state is RACGT
    else
    {
       if (total_blength_1 > 0)
       {
           if (inverse && seq1_region.plength_observation2root > 0)
           {
               // init lh vector for seq1, taking into account
               // the evolution along plength_observation2node
               seq1_lh_vec = initLhVecMerge<num_states>(seq1_state,
                                seq1_region.plength_observation2node, model, false);
               
               // combine with the root freqs
               for (auto i = 0; i < num_states; ++i)
               {
                   seq1_lh_vec[i] *= model->root_freqs[i];
               }
               
               // account for the evolution along plength_observation2root
               // and blength_1
               RealNumType seq1_total_observation2root =
                    seq1_region.plength_observation2root > 0 ?
                    seq1_region.plength_observation2root : 0;
               seq1_total_observation2root +=
                    blength_1 > 0 ? blength_1 : 0;
               if (seq1_total_observation2root)
               {
                   seq1_lh_vec = initLhVecMerge<num_states>(TYPE_O,
                        seq1_total_observation2root, model, true,
                        cmaple::make_unique<SeqRegion::LHType>(std::move(seq1_lh_vec)));
               }
               
           }
           else
           {
               seq1_lh_vec = initLhVecMerge<num_states>(seq1_state,
                                total_blength_1, model, inverse);
           }
       }
       else
       {
           resetLhVecExceptState<num_states>(seq1_lh_vec.data(), seq1_state, 1);
       }
    }

    // extract seq2_state
    const StateType seq2_state = seq2_region.type == TYPE_R ? ref_state : seq2_region.type;

    // extract the lh vector for seq2,
    // taking into account the evolution along the branch
    // seq2 is O
    if (seq2_state == TYPE_O)
    {
        if (total_blength_2 > 0)
        {
            assert(seq2_state == TYPE_O);
            seq2_lh_vec = initLhVecMerge<num_states>(TYPE_O,
                            total_blength_2, model, false, seq2_region.likelihood);
        }
        else
        {
            seq2_lh_vec = *seq2_region.likelihood;
        }
   }
   // seq2 is RACGT
   else
   {
       if (total_blength_2 > 0)
       {
           seq2_lh_vec = initLhVecMerge<num_states>(seq2_state,
                               total_blength_2, model, false);
       }
       else
       {
           resetLhVecExceptState<num_states>(seq2_lh_vec.data(), seq2_state, 1);
       }
   }

    // merge two lh vectors
    SeqRegion::LHType merged_lh_vec;
    setVecByProduct<num_states>(merged_lh_vec.data(), seq1_lh_vec.data(), seq2_lh_vec.data());
    total_prob = 0;
    for (auto i = 0; i < num_states; ++i)
        total_prob += merged_lh_vec[i];
    
    // handle abnormal case
    if (total_prob == 0)
    {
        merged_regions = nullptr;
        return false;
    }
    
    // normalize the new partial likelihood
    normalize_arr(merged_lh_vec.data(), num_states, total_prob);
    cmaple::SeqRegions::addSimplifiedO(ref_state, end_pos, merged_lh_vec, aln, threshold_prob,
                                         *merged_regions);
    
    return true;
}

template <const cmaple::StateType num_states>
auto merge_ORACGT_ORACGT(const SeqRegion& seq1_region,
                       const SeqRegion& seq2_region,
                       const cmaple::RealNumType total_blength_1,
                       const cmaple::RealNumType total_blength_2,
                       const cmaple::RealNumType blength_1,
                       const cmaple::PositionType end_pos,
                       const Alignment* aln,
                       const ModelBase* model,
                       const cmaple::RealNumType threshold_prob,
                       std::unique_ptr<SeqRegions>& merged_regions,
                       const bool inverse,
                       const bool return_log_lh,
                       cmaple::RealNumType& total_factor) -> bool
{
    RealNumType total_prob = 0;
    if (!do_merge_ORACGT_ORACGT<num_states>(seq1_region, seq2_region,
                                total_blength_1, total_blength_2,
                                blength_1, end_pos, aln, model,
                                threshold_prob, merged_regions,
                                inverse, total_prob))
    {
        return false;
    }

    if (return_log_lh)
    {
       total_factor *= total_prob;
    }
    
    return true;
}

template <const StateType num_states>
auto merge_notN_notN_TwoLowers(const SeqRegion& seq1_region,
                               const SeqRegion& seq2_region,
                               const RealNumType plength1,
                               const RealNumType plength2,
                               const PositionType end_pos,
                               const PositionType pos,
                               const Alignment* aln,
                               const ModelBase* model,
                               const RealNumType* const cumulative_rate,
                               const RealNumType threshold_prob,
                               RealNumType& log_lh,
                               cmaple::RealNumType& total_factor,
                               std::unique_ptr<SeqRegions>& merged_regions,
                               const bool return_log_lh) -> bool {
  assert(seq1_region.type != TYPE_N);
  assert(seq2_region.type != TYPE_N);
  assert(model);
  assert(aln);
  assert(cumulative_rate);

  RealNumType total_blength_1 = plength1;
  if (seq1_region.plength_observation2node >= 0) {
    total_blength_1 = seq1_region.plength_observation2node;
    if (plength1 > 0) {
      total_blength_1 += plength1;
    }
  }

  RealNumType total_blength_2 = plength2;
  if (seq2_region.plength_observation2node >= 0) {
    total_blength_2 = seq2_region.plength_observation2node;
    if (plength2 > 0) {
      total_blength_2 += plength2;
    }
  }

  // seq1_entry and seq2_entry are identical seq1_entry = R/ACGT
  if (seq1_region.type == seq2_region.type &&
      (seq1_region.type == TYPE_R || seq1_region.type < num_states)) {
    merge_identicalRACGT_TwoLowers(seq1_region, end_pos, total_blength_1,
                                   total_blength_2, plength1, plength2,
                                   pos, threshold_prob, model,
                                   cumulative_rate, log_lh, *merged_regions,
                                   return_log_lh);
  }
  // #0 distance between different nucleotides: merge is not possible
  else if (total_blength_1 <= 0 && total_blength_2 <= 0 &&
           (seq1_region.type == TYPE_R || seq1_region.type < num_states) &&
           (seq2_region.type == TYPE_R || seq2_region.type < num_states)) {
    merged_regions = nullptr;
    return false;
  }
 // seq1 and seq2 are ORACGT
 else
 {
     return merge_ORACGT_ORACGT<num_states>(seq1_region, seq2_region,
        total_blength_1, total_blength_2, plength1, end_pos, aln, model,
        threshold_prob, merged_regions, false, return_log_lh, total_factor);
 }
/* ------ OLD-CODE REMOVED AFTER UPDATING TO MAPLE 0.7.5 --------------
  // seq1_entry = O
  else if (seq1_region.type == TYPE_O) {
    auto ret = merge_O_ORACGT_TwoLowers<num_states>(
        seq1_region, seq2_region, total_blength_1, total_blength_2, end_pos,
        aln, model, threshold_prob, log_lh, merged_regions, return_log_lh);
    return ret;
  }
  // seq1_entry = R/ACGT
  else {
    auto ret = merge_RACGT_ORACGT_TwoLowers<num_states>(
        seq1_region, seq2_region, total_blength_1, total_blength_2, end_pos,
        aln, model, threshold_prob, log_lh, merged_regions, return_log_lh);
    return ret;
  }
*/
  // no error
  return true;
}

template <const StateType num_states>
RealNumType SeqRegions::mergeTwoLowers(
    std::unique_ptr<SeqRegions>& merged_regions,
    const RealNumType plength1,
    const SeqRegions& regions2,
    const RealNumType plength2,
    const Alignment* aln,
    const ModelBase* model,
    const RealNumType* const cumulative_rate,
    const RealNumType threshold_prob,
    const bool return_log_lh) const {
  assert(regions2.size() > 0);
  assert(model);
  assert(aln);
  assert(cumulative_rate);

  // init variables
  RealNumType total_factor = 1.0;
  RealNumType log_lh = 0;
  PositionType pos = 0;
  const SeqRegions& seq1_regions = *this;
  const SeqRegions& seq2_regions = regions2;
  size_t iseq1 = 0;
  size_t iseq2 = 0;
  const PositionType seq_length = static_cast<PositionType>(aln->ref_seq.size());
  // contribution to non-mutation for the whole genome
  const RealNumType global_total_rate = -seq_length;
  const RealNumType total_blength = ((plength1 > 0 ? plength1 : 0) + (plength2 > 0 ? plength2 : 0));
  if (return_log_lh)
      log_lh = total_blength * global_total_rate;

  // init merged_regions
  if (merged_regions) {
    merged_regions->clear();
  } else {
    merged_regions = cmaple::make_unique<SeqRegions>();
  }

  // avoid realloc of vector data (minimize memory footprint)
  merged_regions->reserve(countSharedSegments(
      seq2_regions, static_cast<size_t>(seq_length)));  // avoid realloc of vector data
#ifdef DEBUG
  const size_t max_elements =
      merged_regions
          ->capacity();  // remember capacity (may be more than we 'reserved')
#endif

  while (pos < seq_length) {
    PositionType end_pos;
    // get the next shared segment in the two sequences
    cmaple::SeqRegions::getNextSharedSegment(pos, seq1_regions, seq2_regions,
                                             iseq1, iseq2, end_pos);
    const auto* const seq1_region = &seq1_regions[iseq1];
    const auto* const seq2_region = &seq2_regions[iseq2];
    const DoubleState s1s2 =
        (DoubleState(seq1_region->type) << 8) | seq2_region->type;

    // seq1_entry = 'N'
    // seq1_entry = 'N' and seq2_entry = 'N'
    if (s1s2 == NN) {
      merged_regions->emplace_back(TYPE_N, end_pos, TYPE_N);
      // seq1_entry = 'N' and seq2_entry = O/R/ACGT
      // seq1_entry = 'N' and seq2_entry = O
    } else if (s1s2 == NO) {
      merge_N_O_TwoLowers(*seq2_region, end_pos, plength2, *merged_regions);
    }
    // seq1_entry = 'N' and seq2_entry = R/ACGT
    else if (seq1_region->type == TYPE_N) {
      merge_N_RACGT_TwoLowers(*seq2_region, end_pos, plength2, threshold_prob,
                              *merged_regions);
    }
    // seq2_entry = 'N'
    // seq1_entry = 'O' and seq2_entry = N
    else if (s1s2 == ON) {
      // NOTE: merge_N_O_TwoLowers can be used to merge O_N
      merge_N_O_TwoLowers(*seq1_region, end_pos, plength1, *merged_regions);
    }
    // seq1_entry = R/ACGT and seq2_entry = 'N'
    else if (seq2_region->type == TYPE_N) {
      // NOTE: merge_N_RACGT can be used to merge RACGT_N
      merge_N_RACGT_TwoLowers(*seq1_region, end_pos, plength1, threshold_prob,
                              *merged_regions);
    }
    // neither seq1_entry nor seq2_entry = N
    else {
        // compute lh, if needed
        // removing pre-calculated overall contributions to the likelihood
        // from this position if they're not both R
        if (return_log_lh
            && !(seq1_region->type == TYPE_R && seq2_region->type == TYPE_R))
        {
            const StateType ref_state = seq1_region->type != TYPE_R ?
                                    seq1_region->prev_state : seq2_region->prev_state;
            log_lh -= model->diagonal_mut_mat[ref_state] * total_blength;
        }
        
      if (!merge_notN_notN_TwoLowers<num_states>(
              *seq1_region, *seq2_region, plength1, plength2, end_pos, pos, aln,
              model, cumulative_rate, threshold_prob, log_lh, total_factor, merged_regions,
              return_log_lh)) {
        return MIN_NEGATIVE;
      }
    }
      
    // calculate lh if needed
    if (return_log_lh
        && (seq1_region->type == TYPE_N || seq2_region->type == TYPE_N))
    {
        // removing pre-calculated overall contributions to the likelihood from these positions?
        log_lh += total_blength * (cumulative_rate[pos] - cumulative_rate[end_pos + 1]);
    }

    // NHANLT: LOGS FOR DEBUGGING
    /*if (Params::getInstance().debug &&
    merged_regions->at(merged_regions->size()-1).type == TYPE_O)
    {
        SeqRegion::LHType lh =
    *merged_regions->at(merged_regions->size()-1).likelihood; std::cout <<
    "merge2Low " << pos << " " << std::setprecision(20) << lh[0] << " " << lh[1]
    << " " << lh[2] << " " << lh[3] << " " << std::endl;
    }*/
      
    // avoid underflow on total_factor
    // approximately update lh_cost and total_factor
    if (return_log_lh && total_factor <= MIN_CARRY_OVER) {
        if (total_factor < MIN_POSITIVE) {
           return MIN_NEGATIVE;
        }

        // lh_cost += log(total_factor);
        // total_factor = 1.0;
        total_factor *= MAX_POSITIVE;
        log_lh -= LOG_MAX_POSITIVE;
    }

    // update pos
    pos = end_pos + 1;
  }

#ifdef DEBUG
  assert(merged_regions->capacity() ==
         max_elements);  // ensure we did the correct reserve, otherwise it was
  // a pessimization
#endif
        
  // update log_lh
  log_lh += log(total_factor);

  return log_lh;
}

template <const StateType num_states>
auto SeqRegions::computeAbsoluteLhAtRoot(
    Alignment* aln,
    const ModelBase* model,
    const std::vector<std::vector<PositionType>>& cumulative_base)
    -> RealNumType {
  assert(model);
  assert(size() > 0);
        
  // dummy variables
  RealNumType log_lh = 0;
  RealNumType log_factor = 1;
  PositionType start_pos = 0;
  const SeqRegions& regions = *this;

  // browse regions one by one to compute the likelihood of each region
  for (const SeqRegion& region : regions) {
    // type R
    if (region.type == TYPE_R) {
      for (StateType i = 0; i < num_states; ++i) {
        log_lh += model->getRootLogFreq(i) *
                  (cumulative_base[static_cast<size_t>(region.position) + 1][i] -
                   cumulative_base[static_cast<size_t>(start_pos)][i]); 
      }
    }
    // type ACGT
    else if (region.type < num_states) {
      log_lh += model->getRootLogFreq(region.type);
      // type O
    } else if (region.type == TYPE_O) {
      RealNumType tot = 0;
      tot +=
          dotProduct<num_states>(&(*region.likelihood)[0], model->getRootFreqs());
      log_factor *= tot;
    }

    // maintain start_pos
    start_pos = region.position + 1;

    // NhanLT: avoid underflow on log_factor
    if (log_factor <= MIN_CARRY_OVER) {
      // if (log_factor < MIN_POSITIVE)
      //     return MIN_NEGATIVE;

      // lh_cost += log(total_factor);
      // total_factor = 1.0;
      log_factor *= MAX_POSITIVE;
      log_lh -= LOG_MAX_POSITIVE;
    }
  }

  // update log_lh
  log_lh += log(log_factor);

  // return the absolute likelihood
  return log_lh;
}

template <const StateType num_states>
RealNumType SeqRegions::computeSiteLhAtRoot(
    std::vector<RealNumType>& site_lh_contributions,
    const ModelBase* model,
    const std::vector<std::vector<PositionType>>& cumulative_base) {
  assert(model);
  assert(size() > 0);
    
  // dummy variables
  RealNumType log_lh = 0;
  RealNumType log_factor = 1;
  PositionType start_pos = 0;
  const SeqRegions& regions = *this;

  // browse regions one by one to compute the likelihood of each region
  for (const SeqRegion& region : regions) {
    // type R
    if (region.type == TYPE_R) {
      for (StateType i = 0; i < num_states; ++i) {
        log_lh += model->getRootLogFreq(i) *
                  (cumulative_base[static_cast<size_t>(region.position) + 1][i] -
                   cumulative_base[static_cast<size_t>(start_pos)][i]);
      }

      // calculate site lhs
      for (PositionType pos = start_pos; pos < region.position + 1; ++pos) {
        for (StateType i = 0; i < num_states; ++i) {
          site_lh_contributions[static_cast<std::vector<RealNumType>::size_type>(pos)] +=
              model->getRootLogFreq(i) *
              (cumulative_base[static_cast<size_t>(pos) + 1][i] -
               cumulative_base[static_cast<size_t>(pos)][i]);
        }
      }
    }
    // type ACGT
    else if (region.type < num_states) {
      RealNumType lh_contribution = model->getRootLogFreq(region.type);
      log_lh += lh_contribution;

      // calculate site lhs
      site_lh_contributions[static_cast<std::vector<RealNumType>
                            ::size_type>(start_pos)] += lh_contribution;
    }
    // type O
    else if (region.type == TYPE_O) {
      RealNumType tot = 0;
      tot +=
          dotProduct<num_states>(&(*region.likelihood)[0], model->getRootFreqs());
      log_factor *= tot;

      // calculate site lhs
      site_lh_contributions[static_cast<std::vector<RealNumType>
                            ::size_type>(start_pos)] += log(tot);
    }

    // maintain start_pos
    start_pos = region.position + 1;

    // NhanLT: avoid underflow on log_factor
    if (log_factor <= MIN_CARRY_OVER) {
      // if (log_factor < MIN_POSITIVE)
      //     return MIN_NEGATIVE;

      // lh_cost += log(total_factor);
      // total_factor = 1.0;
      log_factor *= MAX_POSITIVE;
      log_lh -= LOG_MAX_POSITIVE;
    }
  }

  // update log_lh
  log_lh += log(log_factor);

  // return the absolute likelihood
  return log_lh;
}

template <const StateType num_states>
void SeqRegions::computeTotalLhAtRoot(std::unique_ptr<SeqRegions>& total_lh,
                                      const ModelBase* model,
                                      RealNumType blength) const {
  assert(model);
  assert(size() > 0);
    
  if (total_lh) {
    total_lh->clear();
  } else {
    total_lh = cmaple::make_unique<SeqRegions>();
  }

  total_lh->reserve(size());  // avoid realloc of vector data
  for (const SeqRegion& elem : (*this)) {
    const auto* const region = &elem;
    // type N
    if (region->type == TYPE_N) {
        // need double-check region->prev_state
      total_lh->emplace_back(region->type, region->position, region->prev_state,
                             region->plength_observation2node,
                             region->plength_observation2root);
    } else {
      // type O
      if (region->type == TYPE_O) {
        // compute total blength
        RealNumType total_blength = blength;
        if (region->plength_observation2node >= 0) {
          total_blength = region->plength_observation2node;
          if (blength > 0) {
            total_blength += blength;
          }
        }

        // init new likelihood
        auto new_lh = cmaple::make_unique<SeqRegion::LHType>();  // = new
        // RealNumType[num_states];
        RealNumType sum_lh = updateLHwithModel<num_states>(
            model, *region->likelihood, (*new_lh), total_blength, region->position);
        // normalize the new partial likelihood
        normalize_arr(new_lh->data(), num_states, sum_lh);

        // add new region to the total_lh_regions
        // need double-check region->prev_state
        total_lh->emplace_back(
            region->type, region->position, region->prev_state, region->plength_observation2node,
            region->plength_observation2root, std::move(new_lh));
      }
      // other types: R or A/C/G/T
      else {
        // add new region to the total_lh_regions
#if __cplusplus >= 201703L
        // need double-check region->prev_state
        SeqRegion& new_region = total_lh->emplace_back(
            region->type, region->position, region->prev_state, region->plength_observation2node,
            region->plength_observation2root);
#else
          // need double-check region->prev_state
          total_lh->emplace_back(
              region->type, region->position, region->prev_state, region->plength_observation2node,
              region->plength_observation2root);
          SeqRegion& new_region = total_lh->at(total_lh->size() - 1);
#endif

        if (new_region.plength_observation2node >= 0) {
          if (blength > 0) {
            new_region.plength_observation2node += blength;
          }

          new_region.plength_observation2root = 0;
        } else if (blength > 0) {
          new_region.plength_observation2node = blength;
          new_region.plength_observation2root = 0;
        }
      }
    }
  }
}

void calSiteLhs_identicalRACGT(std::vector<RealNumType>& site_lh_contributions,
                               const SeqRegion& seq1_region,
                               const PositionType end_pos,
                               RealNumType total_blength_1,
                               RealNumType total_blength_2,
                               const cmaple::RealNumType blength_1,
                               const cmaple::RealNumType blength_2,
                               const PositionType pos,
                               const RealNumType threshold_prob,
                               const ModelBase* model,
                               const RealNumType* const cumulative_rate,
                               RealNumType& log_lh,
                               SeqRegions& merged_regions);

template <const StateType num_states>
inline void addSimplifyOAndCalSiteLh(std::vector<RealNumType>& site_lh_contributions,
    RealNumType& log_lh, SeqRegion::LHType& new_lh, RealNumType& sum_lh,
    const PositionType end_pos, const Alignment* aln, const RealNumType threshold_prob,
    std::unique_ptr<SeqRegions>& merged_regions)
{
    // normalize the new partial likelihood
    normalize_arr(new_lh.data(), num_states, sum_lh);
    cmaple::SeqRegions::addSimplifiedO(TYPE_N, end_pos, new_lh, aln, threshold_prob,
                                       *merged_regions);

    // compute (site) lh contributions
    RealNumType lh_contribution = log(sum_lh);
    log_lh += lh_contribution;
    site_lh_contributions[static_cast<std::vector<RealNumType>
                        ::size_type>(end_pos)] += lh_contribution;
}

/*template <const StateType num_states>
bool calSiteLhs_O_O(std::vector<RealNumType>& site_lh_contributions,
                    const SeqRegion& seq2_region,
                    RealNumType total_blength_2,
                    const PositionType end_pos,
                    const Alignment* aln,
                    const ModelBase* model,
                    const RealNumType threshold_prob,
                    RealNumType& log_lh,
                    SeqRegion::LHType& new_lh,
                    std::unique_ptr<SeqRegions>& merged_regions) {
  assert(seq2_region.type == TYPE_O);
  assert(aln);
  assert(model);

  RealNumType sum_lh = updateMultLHwithMat<num_states>(
      model->getMutationMatrix(end_pos), *seq2_region.likelihood, new_lh, total_blength_2);

  if (sum_lh == 0) {
    merged_regions = nullptr;
    return false;
  }
    
  // normalize the new partial likelihood
  // add simplify O
  // compute (site) lh contributions
  addSimplifyOAndCalSiteLh<num_states>(site_lh_contributions, log_lh,
        new_lh, sum_lh, end_pos, aln, threshold_prob, merged_regions);

  // no error
  return true;
}*/

/*template <const StateType num_states>
bool calSiteLhs_O_RACGT(std::vector<RealNumType>& site_lh_contributions,
                        const SeqRegion& seq2_region,
                        const StateType& ref_state,
                        RealNumType total_blength_2,
                        const PositionType end_pos,
                        const Alignment* aln,
                        const ModelBase* model,
                        const RealNumType threshold_prob,
                        RealNumType& log_lh,
                        SeqRegion::LHType& new_lh,
                        RealNumType& sum_lh,
                        std::unique_ptr<SeqRegions>& merged_regions) {
  assert(seq2_region.type != TYPE_N && seq2_region.type != TYPE_O);
  assert(aln);
  assert(model);
    
  StateType seq2_state = seq2_region.type;
  if (seq2_state == TYPE_R) {
    seq2_state = ref_state;
  }

  if (total_blength_2 > 0) {
    sum_lh = updateVecWithState<num_states>(
        new_lh.data(), seq2_state,  
        model->getTransposedMutationMatrixRow(seq2_state, end_pos), 
        total_blength_2);

  // normalize the new partial likelihood
  // add simplify O
  // compute (site) lh contributions
  addSimplifyOAndCalSiteLh<num_states>(site_lh_contributions, log_lh, new_lh,
            sum_lh, end_pos, aln, threshold_prob, merged_regions);
  } else {
    if (new_lh[seq2_state] == 0) {
      merged_regions = nullptr;
      return false;
    }

    // add a new region and try to merge consecutive R regions together
    cmaple::SeqRegions::addNonConsecutiveRRegion(
        *merged_regions, seq2_region.type, seq2_region.prev_state,
        -1, -1, end_pos, threshold_prob);

    // compute (site) lh contributions
    RealNumType lh_contribution = log(new_lh[seq2_state]);
    log_lh += lh_contribution;
    site_lh_contributions[static_cast<std::vector<RealNumType>
                        ::size_type>(end_pos)] += lh_contribution;
  }

  // no error
  return true;
}*/

/*template <const StateType num_states>
bool calSiteLhs_O_ORACGT(std::vector<RealNumType>& site_lh_contributions,
                         const SeqRegion& seq1_region,
                         const SeqRegion& seq2_region,
                         RealNumType total_blength_1,
                         RealNumType total_blength_2,
                         const PositionType end_pos,
                         const Alignment* aln,
                         const ModelBase* model,
                         const RealNumType threshold_prob,
                         RealNumType& log_lh,
                         std::unique_ptr<SeqRegions>& merged_regions) {
  assert(seq1_region.type == TYPE_O);
  assert(seq2_region.type != TYPE_N);
  assert(aln);
  assert(model);
    
  auto new_lh =
      cmaple::make_unique<SeqRegion::LHType>();  // = new RealNumType[num_states];
  RealNumType sum_lh = 0;

  if (total_blength_1 > 0) {
    updateLHwithMat<num_states>(model->getMutationMatrix(end_pos), *(seq1_region.likelihood),
                                *new_lh, total_blength_1);
    // otherwise, clone the partial likelihood from seq1
  } else {
    *new_lh = *seq1_region.likelihood;
  }

  // seq1_entry = O and seq2_entry = O
  if (seq2_region.type == TYPE_O) {
    return calSiteLhs_O_O<num_states>(
        site_lh_contributions, seq2_region, total_blength_2, end_pos, aln,
        model, threshold_prob, log_lh, *new_lh, merged_regions);
  }
  // seq1_entry = O and seq2_entry = R/ACGT
  else {
    return calSiteLhs_O_RACGT<num_states>(
        site_lh_contributions, seq2_region, seq1_region.prev_state, total_blength_2,
        end_pos, aln, model, threshold_prob, log_lh, *new_lh, sum_lh, merged_regions);
  }

  // no error
  return true;
}*/

/*template <const StateType num_states>
bool calSiteLhs_RACGT_O(std::vector<RealNumType>& site_lh_contributions,
                        const SeqRegion& seq2_region,
                        RealNumType total_blength_2,
                        const PositionType end_pos,
                        const Alignment* aln,
                        const ModelBase* model,
                        const RealNumType threshold_prob,
                        SeqRegion::LHType& new_lh,
                        RealNumType& log_lh,
                        std::unique_ptr<SeqRegions>& merged_regions) {
  assert(seq2_region.type == TYPE_O);
  assert(aln);
  assert(model);
  RealNumType sum_lh = updateMultLHwithMat<num_states>(
      model->getMutationMatrix(end_pos), *(seq2_region.likelihood), new_lh, total_blength_2);

  if (sum_lh == 0) {
    merged_regions = nullptr;
    return false;
  }

  // normalize the new partial likelihood
  // add simplify O
  // compute (site) lh contributions
  addSimplifyOAndCalSiteLh<num_states>(site_lh_contributions, log_lh,
        new_lh, sum_lh, end_pos, aln, threshold_prob, merged_regions);

  // no error
  return true;
}*/

/*template <const StateType num_states>
bool calSiteLhs_RACGT_RACGT(std::vector<RealNumType>& site_lh_contributions,
                            const SeqRegion& seq2_region,
                            const StateType& ref_state,
                            RealNumType total_blength_2,
                            const PositionType end_pos,
                            const Alignment* aln,
                            const ModelBase* model,
                            const RealNumType threshold_prob,
                            SeqRegion::LHType& new_lh,
                            RealNumType& sum_lh,
                            RealNumType& log_lh,
                            std::unique_ptr<SeqRegions>& merged_regions) {
  assert(seq2_region.type != TYPE_N && seq2_region.type != TYPE_O);
  assert(aln);
  assert(model);
    
  StateType seq2_state = seq2_region.type;
  if (seq2_state == TYPE_R) {
    seq2_state = ref_state;
  }

  if (total_blength_2 > 0) {
    sum_lh += updateVecWithState<num_states>(
        new_lh.data(), seq2_state,  
        model->getTransposedMutationMatrixRow(seq2_state, end_pos), 
        total_blength_2);

  // normalize the new partial likelihood
  // add simplify O
  // compute (site) lh contributions
  addSimplifyOAndCalSiteLh<num_states>(site_lh_contributions, log_lh,
                new_lh, sum_lh, end_pos, aln, threshold_prob, merged_regions);
  } else {
    // add a new region and try to merge consecutive R regions together
    cmaple::SeqRegions::addNonConsecutiveRRegion(
        *merged_regions, seq2_region.type, seq2_region.prev_state,
        -1, -1, end_pos, threshold_prob);

    // compute (site) lh contributions
    RealNumType lh_contribution = log(new_lh[seq2_state]);
    log_lh += lh_contribution;
    site_lh_contributions[static_cast<std::vector<RealNumType>
                            ::size_type>(end_pos)] += lh_contribution;
  }

  // no error
  return true;
}*/

/*template <const StateType num_states>
bool calSiteLhs_RACGT_ORACGT(std::vector<RealNumType>& site_lh_contributions,
                             const SeqRegion& seq1_region,
                             const SeqRegion& seq2_region,
                             RealNumType total_blength_1,
                             RealNumType total_blength_2,
                             const PositionType end_pos,
                             const Alignment* aln,
                             const ModelBase* model,
                             const RealNumType threshold_prob,
                             RealNumType& log_lh,
                             std::unique_ptr<SeqRegions>& merged_regions) {
  assert(seq1_region.type != TYPE_N && seq1_region.type != TYPE_O);
  assert(seq2_region.type != TYPE_N);
  assert(aln);
  assert(model);
    
  StateType seq1_state = seq1_region.type;
  if (seq1_state == TYPE_R) {
      seq1_state = seq2_region.prev_state;
  }

  auto new_lh =
      cmaple::make_unique<SeqRegion::LHType>();  // = new RealNumType[num_states];
  // auto& new_lh_value = *new_lh;
  RealNumType sum_lh = 0;

  if (total_blength_1 > 0) {
    setVecWithState<num_states>(new_lh->data(), seq1_state,
                                 model->getTransposedMutationMatrixRow(seq1_state, end_pos), total_blength_1);
  } else {
    resetLhVecExceptState<num_states>(new_lh->data(), seq1_state, 1);
  }

  // seq1_entry = R/ACGT and seq2_entry = O
  if (seq2_region.type == TYPE_O) {
    return calSiteLhs_RACGT_O<num_states>(
        site_lh_contributions, seq2_region, total_blength_2, end_pos, aln,
        model, threshold_prob, *new_lh, log_lh, merged_regions);
  }

  // otherwise, seq1_entry = R/ACGT and seq2_entry = R/ACGT
  return calSiteLhs_RACGT_RACGT<num_states>(
      site_lh_contributions, seq2_region, seq1_region.prev_state, total_blength_2,
      end_pos, aln, model, threshold_prob, *new_lh, sum_lh, log_lh, merged_regions);
}*/

template <const cmaple::StateType num_states>
auto calSiteLhs_ORACGT_ORACGT(std::vector<RealNumType>& site_lh_contributions,
                              const SeqRegion& seq1_region,
                              const SeqRegion& seq2_region,
                              const cmaple::RealNumType total_blength_1,
                              const cmaple::RealNumType total_blength_2,
                              const cmaple::RealNumType blength_1,
                              const cmaple::PositionType end_pos,
                              const Alignment* aln,
                              const ModelBase* model,
                              const cmaple::RealNumType threshold_prob,
                              std::unique_ptr<SeqRegions>& merged_regions,
                              cmaple::RealNumType& total_factor) -> bool
{
    RealNumType total_prob = 0;
    if (!do_merge_ORACGT_ORACGT<num_states>(seq1_region, seq2_region,
                                total_blength_1, total_blength_2,
                                blength_1, end_pos, aln, model,
                                threshold_prob, merged_regions,
                                false, total_prob))
    {
        return false;
    }

    total_factor *= total_prob;
    
    // update site_lh_contributions
    site_lh_contributions[end_pos] += log(total_prob);
    
    return true;
}

template <const StateType num_states>
bool calSiteLhs_notN_notN(std::vector<RealNumType>& site_lh_contributions,
                          const SeqRegion& seq1_region,
                          const SeqRegion& seq2_region,
                          const RealNumType plength1,
                          const RealNumType plength2,
                          const PositionType end_pos,
                          const PositionType pos,
                          const Alignment* aln,
                          const ModelBase* model,
                          const RealNumType* const cumulative_rate,
                          const RealNumType threshold_prob,
                          RealNumType& log_lh,
                          cmaple::RealNumType& total_factor,
                          std::unique_ptr<SeqRegions>& merged_regions) {
  assert(seq1_region.type != TYPE_N);
  assert(seq2_region.type != TYPE_N);
  assert(aln);
  assert(model);
  assert(cumulative_rate);
    
  RealNumType total_blength_1 = plength1;
  if (seq1_region.plength_observation2node >= 0) {
    total_blength_1 = seq1_region.plength_observation2node;
    if (plength1 > 0) {
      total_blength_1 += plength1;
    }
  }

  RealNumType total_blength_2 = plength2;
  if (seq2_region.plength_observation2node >= 0) {
    total_blength_2 = seq2_region.plength_observation2node;
    if (plength2 > 0) {
      total_blength_2 += plength2;
    }
  }

  // seq1_entry and seq2_entry are identical seq1_entry = R/ACGT
  if (seq1_region.type == seq2_region.type &&
      (seq1_region.type == TYPE_R || seq1_region.type < num_states)) {
    calSiteLhs_identicalRACGT(site_lh_contributions, seq1_region, end_pos,
                              total_blength_1, total_blength_2, plength1,
                              plength2, pos, threshold_prob, model,
                              cumulative_rate, log_lh, *merged_regions);
  }
  // #0 distance between different nucleotides: merge is not possible
  else if (total_blength_1 == 0 && total_blength_2 == 0 &&
           (seq1_region.type == TYPE_R || seq1_region.type < num_states) &&
           (seq2_region.type == TYPE_R || seq2_region.type < num_states)) {
    merged_regions = nullptr;
    return false;
  }
  // seq1 and seq2 are ORACGT
  else
  {
      return calSiteLhs_ORACGT_ORACGT<num_states>(site_lh_contributions, seq1_region,
        seq2_region, total_blength_1, total_blength_2, plength1, end_pos, aln, model,
        threshold_prob, merged_regions, total_factor);
  }
  /* ------ OLD-CODE REMOVED AFTER UPDATING TO MAPLE 0.7.5 --------------
  // seq1_entry = O
  else if (seq1_region.type == TYPE_O) {
    return calSiteLhs_O_ORACGT<num_states>(
        site_lh_contributions, seq1_region, seq2_region, total_blength_1,
        total_blength_2, end_pos, aln, model, threshold_prob, log_lh,
        merged_regions);
  }
  // seq1_entry = R/ACGT
  else {
    return calSiteLhs_RACGT_ORACGT<num_states>(
        site_lh_contributions, seq1_region, seq2_region, total_blength_1,
        total_blength_2, end_pos, aln, model, threshold_prob, log_lh,
        merged_regions);
  }*/

  // no error
  return true;
}

template <const StateType num_states>
RealNumType SeqRegions::calculateSiteLhContributions(
    std::vector<RealNumType>& site_lh_contributions,
    std::unique_ptr<SeqRegions>& merged_regions,
    const RealNumType plength1,
    const SeqRegions& regions2,
    const RealNumType plength2,
    const Alignment* aln,
    const ModelBase* model,
    const RealNumType* const cumulative_rate,
    const RealNumType threshold_prob) const {
  assert(aln);
  assert(model);
  assert(cumulative_rate);
        
  // init variables
  RealNumType total_factor = 1.0;
  RealNumType log_lh = 0;
  PositionType pos = 0;
  const SeqRegions& seq1_regions = *this;
  const SeqRegions& seq2_regions = regions2;
  size_t iseq1 = 0;
  size_t iseq2 = 0;
  const PositionType seq_length = static_cast<PositionType>(aln->ref_seq.size());
  // contribution to non-mutation for the whole genome
  const RealNumType global_total_rate = -seq_length;
  const RealNumType total_blength = ((plength1 > 0 ? plength1 : 0) + (plength2 > 0 ? plength2 : 0));
  assert(site_lh_contributions.size() == seq_length);
  log_lh = total_blength * global_total_rate;
  for (RealNumType &site_lh_contribution : site_lh_contributions)
  {
      site_lh_contribution -= total_blength;
  }

  // init merged_regions
  if (merged_regions) {
    merged_regions->clear();
  } else {
    merged_regions = cmaple::make_unique<SeqRegions>();
  }

  // avoid realloc of vector data (minimize memory footprint)
  merged_regions->reserve(countSharedSegments(
      seq2_regions, static_cast<size_t>(seq_length)));  // avoid realloc of vector data
#ifdef DEBUG
  const size_t max_elements =
      merged_regions
          ->capacity();  // remember capacity (may be more than we 'reserved')
#endif

  while (pos < seq_length) {
    PositionType end_pos;

    // get the next shared segment in the two sequences
    cmaple::SeqRegions::getNextSharedSegment(pos, seq1_regions, seq2_regions,
                                             iseq1, iseq2, end_pos);
    const auto* const seq1_region = &seq1_regions[iseq1];
    const auto* const seq2_region = &seq2_regions[iseq2];
    const DoubleState s1s2 =
        (DoubleState(seq1_region->type) << 8) | seq2_region->type;

    // seq1_entry = 'N'
    // seq1_entry = 'N' and seq2_entry = 'N'
    if (s1s2 == NN) {
      merged_regions->emplace_back(TYPE_N, end_pos, TYPE_N);
      // seq1_entry = 'N' and seq2_entry = O/R/ACGT
      // seq1_entry = 'N' and seq2_entry = O
    } else if (s1s2 == NO) {
      merge_N_O_TwoLowers(*seq2_region, end_pos, plength2, *merged_regions);
    }
    // seq1_entry = 'N' and seq2_entry = R/ACGT
    else if (seq1_region->type == TYPE_N) {
      merge_N_RACGT_TwoLowers(*seq2_region, end_pos, plength2, threshold_prob,
                              *merged_regions);
    }
    // seq2_entry = 'N'
    // seq1_entry = 'O' and seq2_entry = N
    else if (s1s2 == ON) {
      // NOTE: merge_N_O_TwoLowers can be used to merge O_N
      merge_N_O_TwoLowers(*seq1_region, end_pos, plength1, *merged_regions);
    }
    // seq1_entry = R/ACGT and seq2_entry = 'N'
    else if (seq2_region->type == TYPE_N) {
      // NOTE: merge_N_RACGT can be used to merge RACGT_N
      merge_N_RACGT_TwoLowers(*seq1_region, end_pos, plength1, threshold_prob,
                              *merged_regions);
    }
    // neither seq1_entry nor seq2_entry = N
    else {
        // compute lh
        // removing pre-calculated overall contributions to the likelihood
        // from this position if they're not both R
        if (!(seq1_region->type == TYPE_R && seq2_region->type == TYPE_R))
        {
            const StateType ref_state = seq1_region->type != TYPE_R ?
                                    seq1_region->prev_state : seq2_region->prev_state;
            const RealNumType lh_deduction = model->diagonal_mut_mat[ref_state] * total_blength;
            log_lh -= lh_deduction;
            // update site_lh_contributions
            site_lh_contributions[pos] -= lh_deduction;
        }
        
      if (!calSiteLhs_notN_notN<num_states>(
              site_lh_contributions, *seq1_region, *seq2_region, plength1,
              plength2, end_pos, pos, aln, model, cumulative_rate,
              threshold_prob, log_lh, total_factor, merged_regions)) {
        if (cmaple::verbose_mode >= cmaple::VB_DEBUG) {
          outWarning("calculateSiteLhContributions() returns MIN_NEGATIVE!");
        }
        return MIN_NEGATIVE;
      }
    }
      
      // calculate lh
      if (seq1_region->type == TYPE_N || seq2_region->type == TYPE_N)
      {
          // removing pre-calculated overall contributions to the likelihood from these positions?
          log_lh += total_blength * (cumulative_rate[pos] - cumulative_rate[end_pos + 1]);
          
          // update site_lh_contributions
          for (PositionType site_index = pos; site_index <= end_pos; ++site_index)
          {
              site_lh_contributions[site_index] += total_blength
                * (cumulative_rate[site_index] - cumulative_rate[site_index + 1]);
          }
      }

    // NHANLT: LOGS FOR DEBUGGING
    /*if (Params::getInstance().debug &&
    merged_regions->at(merged_regions->size()-1).type == TYPE_O)
    {
        SeqRegion::LHType lh =
    *merged_regions->at(merged_regions->size()-1).likelihood; std::cout <<
    "merge2Low " << pos << " " << std::setprecision(20) << lh[0] << " " << lh[1]
    << " " << lh[2] << " " << lh[3] << " " << std::endl;
    }*/
      
      // avoid underflow on total_factor
      // approximately update lh_cost and total_factor
      if (total_factor <= MIN_CARRY_OVER) {
          if (total_factor < MIN_POSITIVE) {
             return MIN_NEGATIVE;
          }

          // lh_cost += log(total_factor);
          // total_factor = 1.0;
          total_factor *= MAX_POSITIVE;
          log_lh -= LOG_MAX_POSITIVE;
      }

    // update pos
    pos = end_pos + 1;
  }

#ifdef DEBUG
  assert(merged_regions->capacity() ==
         max_elements);  // ensure we did the correct reserve, otherwise it was
#endif
  // a pessimization
        
  // update log_lh
  log_lh += log(total_factor);

  return log_lh;
}

template <const StateType num_states>
auto cmaple::SeqRegions::containAtLeastNMuts(const int min_mut) const
    -> bool
{
    int count_mutations = 0;
    
    // loop over the vector of regions
    for (auto i = 0; i < size(); ++i)
    {
        const SeqRegion& seq_region = at(i);
        
        if (seq_region.type < num_states)
        {
            ++count_mutations;
            
            // check if it meets the requirement
            if (count_mutations >= min_mut)
                return true;
        }
    }
    
    return false;
}

template <const cmaple::StateType num_states>
void cmaple::SeqRegions::flipMutations()
{
    for (SeqRegion& region : *this)
    {
        // only flip mutations
        if (region.type < num_states)
        {
            const StateType prev_state = region.prev_state;
            region.prev_state = region.type;
            region.type = prev_state;
        }
    }
}

}  // namespace cmaple
