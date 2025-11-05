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
   Add a new region and automatically merged consecutive R regions
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  static void addNonConsecutiveRRegion(
      SeqRegions& regions,
      const cmaple::StateType new_region_type,
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
     @param output_regions the output regions
     @param mutations the vector of mutations
     @param aln the alignment
     @param inverse True to inverse the mutations
     @throw std::logic\_error if unexpected values/behaviors found during the
     operations
     */
    template <const cmaple::StateType num_states>
    void integrateMutations(std::unique_ptr<SeqRegions>& output_regions,
                         const SeqRegions& mutations,
                         const Alignment* aln,
                         const bool inverse = false) const;

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
  static void addSimplifiedO(const cmaple::PositionType end_pos,
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
template <const cmaple::StateType num_states>
void merge_O_ORACGT(const SeqRegion& seq1_region,
                    const SeqRegion& seq2_region,
                    const cmaple::RealNumType total_blength_1,
                    const cmaple::RealNumType total_blength_2,
                    const cmaple::PositionType end_pos,
                    const cmaple::RealNumType threshold_prob,
                    const ModelBase* model,
                    const Alignment* aln,
                    SeqRegions& merged_regions);

/**
 MergeUpperLower case RACGT_O

 @throw std::logic\_error if unexpected values/behaviors found during the
 operations
 */
template <const cmaple::StateType num_states>
void merge_RACGT_O(const SeqRegion& seq2_region,
                   const cmaple::RealNumType total_blength_2,
                   const cmaple::PositionType end_pos,
                   SeqRegion::LHType& new_lh,
                   const cmaple::RealNumType threshold_prob,
                   const ModelBase* model,
                   const Alignment* aln,
                   SeqRegions& merged_regions);

/**
 MergeUpperLower case RACGT_RACGT

 @throw std::logic\_error if unexpected values/behaviors found during the
 operations
 */
template <const cmaple::StateType num_states>
void merge_RACGT_RACGT(const SeqRegion& seq2_region,
                       const cmaple::RealNumType total_blength_2,
                       const cmaple::PositionType end_pos,
                       SeqRegion::LHType& new_lh,
                       const ModelBase* model,
                       const Alignment* aln,
                       SeqRegions& merged_regions);

/**
 MergeUpperLower case RACGT_ORACGT

 @throw std::logic\_error if unexpected values/behaviors found during the
 operations
 */
template <const cmaple::StateType num_states>
void merge_RACGT_ORACGT(const SeqRegion& seq1_region,
                        const SeqRegion& seq2_region,
                        const cmaple::RealNumType total_blength_1,
                        const cmaple::RealNumType total_blength_2,
                        const cmaple::RealNumType upper_plength,
                        const cmaple::PositionType end_pos,
                        const cmaple::RealNumType threshold_prob,
                        const ModelBase* model,
                        const Alignment* aln,
                        SeqRegions& merged_regions);

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
template <const cmaple::StateType num_states>
bool merge_O_O_TwoLowers(const SeqRegion& seq2_region,
                         cmaple::RealNumType total_blength_2,
                         const cmaple::PositionType end_pos,
                         const Alignment* aln,
                         const ModelBase* model,
                         const cmaple::RealNumType threshold_prob,
                         cmaple::RealNumType& log_lh,
                         SeqRegion::LHType& new_lh,
                         std::unique_ptr<SeqRegions>& merged_regions,
                         const bool return_log_lh);

/**
 MergeTwoLowers case O_RACGT

 @throw std::logic\_error if unexpected values/behaviors found during the
 operations
 */
template <const cmaple::StateType num_states>
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
                             const bool return_log_lh);

/**
 MergeTwoLowers case O_ORACGT

 @throw std::logic\_error if unexpected values/behaviors found during the
 operations
 */
template <const cmaple::StateType num_states>
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
                              const bool return_log_lh);

/**
 MergeTwoLowers case RACGT_O

 @throw std::logic\_error if unexpected values/behaviors found during the
 operations
 */
template <const cmaple::StateType num_states>
bool merge_RACGT_O_TwoLowers(const SeqRegion& seq2_region,
                             cmaple::RealNumType total_blength_2,
                             const cmaple::PositionType end_pos,
                             const Alignment* aln,
                             const ModelBase* model,
                             const cmaple::RealNumType threshold_prob,
                             SeqRegion::LHType& new_lh,
                             cmaple::RealNumType& log_lh,
                             std::unique_ptr<SeqRegions>& merged_regions,
                             const bool return_log_lh);

/**
 MergeTwoLowers case RACGT_RACGT

 @throw std::logic\_error if unexpected values/behaviors found during the
 operations
 */
template <const cmaple::StateType num_states>
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
                                 const bool return_log_lh);

/**
 MergeTwoLowers case RACGT_ORACGT

 @throw std::logic\_error if unexpected values/behaviors found during the
 operations
 */
template <const cmaple::StateType num_states>
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
                                  const bool return_log_lh);

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
                               std::unique_ptr<SeqRegions>& merged_regions,
                               const bool return_log_lh);

template <const StateType num_states>
auto updateLHwithModel(const ModelBase* model,
                       const SeqRegion::LHType& prior,
                       SeqRegion::LHType& posterior,
                       const RealNumType total_blength,
                       const PositionType pos) -> RealNumType {
  assert(model);
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
  for (StateType i = 0; i < num_states; ++i, mat_row += num_states) {
    RealNumType tot = 0;
    tot += dotProduct<num_states>(&(prior)[0], mat_row);
    tot *= total_blength;
    tot += prior[i];
    posterior[i] = tot;
    sum_lh += tot;
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
void SeqRegions::integrateMutations(std::unique_ptr<SeqRegions>& output_regions,
                                    const SeqRegions& mutations,
                                    const Alignment* aln,
                                    const bool inverse) const {
    assert(mutations.size() > 0);
    assert(aln);
    
  // init variables
  PositionType pos = 0;
  const SeqRegions& seq_regions = *this;
  const SeqRegions& mutation_regions = mutations;
  size_t iseq = 0;
  size_t imut = 0;
  const PositionType seq_length = static_cast<PositionType>(aln->ref_seq.size());

  // init merged_regions
  if (output_regions) {
      output_regions->clear();
  } else {
      output_regions = cmaple::make_unique<SeqRegions>();
  }

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
      
      // first add seq_region to the output lh vector
      output_regions->push_back(SeqRegion::clone(*seq_region));
      // extract the newly added seq_region
      SeqRegion& last_region = output_regions->back();
      
      // case 1: seq_region = N
      if (seq_region->type == TYPE_N)
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
      else if (seq_region->type < num_states)
      {
          // if a mutation occurs at this position,
          // modify the type (i.e., current state) and the previous state
          if (mutation->type < num_states)
          {
              // if the current state is identical to the new state of the mutation
              // change it to R
              if (last_region.type == mut_new_state)
              {
                  last_region.type = TYPE_R;
                  last_region.prev_state = TYPE_N;
              }
              // otherwise, keep the current state but update the previous state
              else
              {
                  last_region.prev_state = mut_new_state;
              }
          }
      }
      // case 3: seq_region = R
      else if (seq_region->type == TYPE_R)
      {
          // if a mutation occurs at this position,
          // modify the type (i.e., current state) and the previous state
          if (mutation->type < num_states)
          {
              // set the previous state of the newly-added region
              // as the new state of the mutation
              last_region.prev_state = mut_new_state;
              
              // extract the previous state of the mutation
              // according to the direction
              cmaple::StateType mut_prev_state = inverse ?
                mutation->type : mutation->prev_state;
              
              // the new state of the newly-added region
              // is R which is also the previous state of the mutation
              last_region.type = mut_prev_state;
          }
          // otherwise, update the last position of the newly-added R region
          {
              last_region.position = end_pos;
          }
      }
      // case 4: seq_region = O
      else
      {
          // if a mutation occurs at this position,
          // modify the previous state
          if (mutation->type < num_states)
          {
              last_region.prev_state = mut_new_state;
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
  merged_target.emplace_back(TYPE_O, end_pos, 0, 0, std::move(new_lh));
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
    merged_regions.emplace_back(TYPE_O, end_pos, 0, 0, std::move(new_lh));
  } else {
    // add merged region into merged_regions
    merged_regions.emplace_back(TYPE_O, end_pos, 0, 0, *(reg_o.likelihood));
  }
}

template <const StateType num_states>
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
  cmaple::SeqRegions::addSimplifiedO(end_pos, new_lh_value, aln, threshold_prob,
                                     merged_regions);
}

template <const StateType num_states>
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
  cmaple::SeqRegions::addSimplifiedO(end_pos, new_lh, aln, threshold_prob,
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
  merged_regions.emplace_back(TYPE_O, end_pos, 0, 0, std::move(new_lh));
}

template <const StateType num_states>
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
    memcpy(root_vec.data(), model->root_freqs, sizeof(SeqRegion::LHType));


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
}

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
      merged_regions->emplace_back(TYPE_N, end_pos);
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
      addNonConsecutiveRRegion(*merged_regions, seq1_region->type, -1, -1,
                               end_pos, threshold_prob);
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
      }

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

template <const StateType num_states>
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
  cmaple::SeqRegions::addSimplifiedO(end_pos, new_lh, aln, threshold_prob,
                                     *merged_regions);

  if (return_log_lh) {
    log_lh += log(sum_lh);
  }

  // no error
  return true;
}

template <const StateType num_states>
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
    cmaple::SeqRegions::addSimplifiedO(end_pos, new_lh, aln, threshold_prob,
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
        *merged_regions, seq2_region.type, -1, -1, end_pos, threshold_prob);

    if (return_log_lh) {
      log_lh += log(new_lh[seq2_state]);
    }
  }
  // no error
  return true;
}

template <const StateType num_states>
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
}

template <const StateType num_states>
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
  cmaple::SeqRegions::addSimplifiedO(end_pos, new_lh, aln, threshold_prob,
                                     *merged_regions);

  if (return_log_lh) {
    log_lh += log(sum_lh);
  }

  // no error
  return true;
}

template <const StateType num_states>
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
    cmaple::SeqRegions::addSimplifiedO(end_pos, new_lh, aln, threshold_prob,
                                       *merged_regions);

    if (return_log_lh) {
      log_lh += log(sum_lh);
    }
  } else {
    // add a new region and try to merge consecutive R regions together
    cmaple::SeqRegions::addNonConsecutiveRRegion(
        *merged_regions, seq2_region.type, -1, -1, end_pos, threshold_prob);

    if (return_log_lh) {
      log_lh += log(new_lh[seq2_state]);
    }
  }

  // no error
  return true;
}

template <const StateType num_states>
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
                                   total_blength_2, pos, threshold_prob, model,
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
  RealNumType log_lh = 0;
  PositionType pos = 0;
  const SeqRegions& seq1_regions = *this;
  const SeqRegions& seq2_regions = regions2;
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
      merged_regions->emplace_back(TYPE_N, end_pos);
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
      if (!merge_notN_notN_TwoLowers<num_states>(
              *seq1_region, *seq2_region, plength1, plength2, end_pos, pos, aln,
              model, cumulative_rate, threshold_prob, log_lh, merged_regions,
              return_log_lh)) {
        return MIN_NEGATIVE;
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

    // update pos
    pos = end_pos + 1;
  }

#ifdef DEBUG
  assert(merged_regions->capacity() ==
         max_elements);  // ensure we did the correct reserve, otherwise it was
  // a pessimization
#endif
  return log_lh;
}

template <const StateType num_states>
auto SeqRegions::computeAbsoluteLhAtRoot(
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
      total_lh->emplace_back(region->type, region->position,
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
        total_lh->emplace_back(
            region->type, region->position, region->plength_observation2node,
            region->plength_observation2root, std::move(new_lh));
      }
      // other types: R or A/C/G/T
      else {
        // add new region to the total_lh_regions
#if __cplusplus >= 201703L
        SeqRegion& new_region = total_lh->emplace_back(
            region->type, region->position, region->plength_observation2node,
            region->plength_observation2root);
#else
          total_lh->emplace_back(
              region->type, region->position, region->plength_observation2node,
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
    cmaple::SeqRegions::addSimplifiedO(end_pos, new_lh, aln, threshold_prob,
                                       *merged_regions);

    // compute (site) lh contributions
    RealNumType lh_contribution = log(sum_lh);
    log_lh += lh_contribution;
    site_lh_contributions[static_cast<std::vector<RealNumType>
                        ::size_type>(end_pos)] += lh_contribution;
}

template <const StateType num_states>
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
}

template <const StateType num_states>
bool calSiteLhs_O_RACGT(std::vector<RealNumType>& site_lh_contributions,
                        const SeqRegion& seq2_region,
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
    seq2_state = aln->ref_seq[static_cast<std::vector<cmaple::StateType>::size_type>(end_pos)];
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
        *merged_regions, seq2_region.type, -1, -1, end_pos, threshold_prob);

    // compute (site) lh contributions
    RealNumType lh_contribution = log(new_lh[seq2_state]);
    log_lh += lh_contribution;
    site_lh_contributions[static_cast<std::vector<RealNumType>
                        ::size_type>(end_pos)] += lh_contribution;
  }

  // no error
  return true;
}

template <const StateType num_states>
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
        site_lh_contributions, seq2_region, total_blength_2, end_pos, aln,
        model, threshold_prob, log_lh, *new_lh, sum_lh, merged_regions);
  }

  // no error
  return true;
}

template <const StateType num_states>
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
}

template <const StateType num_states>
bool calSiteLhs_RACGT_RACGT(std::vector<RealNumType>& site_lh_contributions,
                            const SeqRegion& seq2_region,
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
    seq2_state = aln->ref_seq[static_cast<std::vector<cmaple::StateType>::size_type>(end_pos)];
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
        *merged_regions, seq2_region.type, -1, -1, end_pos, threshold_prob);

    // compute (site) lh contributions
    RealNumType lh_contribution = log(new_lh[seq2_state]);
    log_lh += lh_contribution;
    site_lh_contributions[static_cast<std::vector<RealNumType>
                            ::size_type>(end_pos)] += lh_contribution;
  }

  // no error
  return true;
}

template <const StateType num_states>
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
    seq1_state = aln->ref_seq[static_cast<std::vector<cmaple::StateType>::size_type>(end_pos)];
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
      site_lh_contributions, seq2_region, total_blength_2, end_pos, aln, model,
      threshold_prob, *new_lh, sum_lh, log_lh, merged_regions);
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
                              total_blength_1, total_blength_2, pos,
                              threshold_prob, model, cumulative_rate, log_lh,
                              *merged_regions);
  }
  // #0 distance between different nucleotides: merge is not possible
  else if (total_blength_1 == 0 && total_blength_2 == 0 &&
           (seq1_region.type == TYPE_R || seq1_region.type < num_states) &&
           (seq2_region.type == TYPE_R || seq2_region.type < num_states)) {
    merged_regions = nullptr;
    return false;
  }
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
  }

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
  RealNumType log_lh = 0;
  PositionType pos = 0;
  const SeqRegions& seq1_regions = *this;
  const SeqRegions& seq2_regions = regions2;
  size_t iseq1 = 0;
  size_t iseq2 = 0;
  const PositionType seq_length = static_cast<PositionType>(aln->ref_seq.size());
  assert(site_lh_contributions.size() == seq_length);

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
      merged_regions->emplace_back(TYPE_N, end_pos);
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
      if (!calSiteLhs_notN_notN<num_states>(
              site_lh_contributions, *seq1_region, *seq2_region, plength1,
              plength2, end_pos, pos, aln, model, cumulative_rate,
              threshold_prob, log_lh, merged_regions)) {
        if (cmaple::verbose_mode >= cmaple::VB_DEBUG) {
          outWarning("calculateSiteLhContributions() returns MIN_NEGATIVE!");
        }
        return MIN_NEGATIVE;
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

    // update pos
    pos = end_pos + 1;
  }

#ifdef DEBUG
  assert(merged_regions->capacity() ==
         max_elements);  // ensure we did the correct reserve, otherwise it was
#endif
  // a pessimization

  return log_lh;
}

}  // namespace cmaple
