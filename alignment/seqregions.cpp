//
//  regions.cpp
//  alignment
//
//  Created by NhanLT on 31/3/2022.
//

#include "seqregions.h"
#include <cassert>
#include "../utils/matrix.h"

using namespace cmaple;

cmaple::SeqRegions::SeqRegions(const std::unique_ptr<SeqRegions>& n_regions) {
  if (!n_regions) {
    throw std::invalid_argument("n_regions is null");
  }
    
  assert(n_regions);
    
  // clone regions one by one
  reserve(n_regions->size());
  for (const auto& region : *n_regions) {
    push_back(SeqRegion::clone(region));
  }
}

auto cmaple::SeqRegions::compareWithSample(const SeqRegions& sequence2,
                                           PositionType seq_length,
                                           const Alignment* aln) const -> int {
  assert(seq_length > 0);
  assert(sequence2.size() > 0);
  assert(size() > 0);
  assert(aln);
  assert(aln->ref_seq.size() == seq_length);

  // init dummy variables
  bool seq1_more_info = false;
  bool seq2_more_info = false;
  PositionType pos = 0;
  const SeqRegions& seq1_regions = *this;
  const SeqRegions& seq2_regions = sequence2;
  size_t iseq1 = 0;
  size_t iseq2 = 0;
  const StateType num_states = aln->num_states;

  while (pos < seq_length && (!seq1_more_info || !seq2_more_info)) {
    PositionType end_pos = 0;

    // get the next shared segment in the two sequences
    cmaple::SeqRegions::getNextSharedSegment(pos, seq1_regions, seq2_regions,
                                             iseq1, iseq2, end_pos);
    const auto* const seq1_region = &seq1_regions[iseq1];
    const auto* const seq2_region = &seq2_regions[iseq2];
    // The two regions have different types from each other
    if (seq1_region->type != seq2_region->type) {
      if (seq1_region->type == TYPE_N) {
        seq2_more_info = true;
      } else if (seq2_region->type == TYPE_N) {
        seq1_more_info = true;
      } else if (seq1_region->type == TYPE_O) {
        StateType seq2_state = seq2_region->type;
        if (seq2_state == TYPE_R)
            seq2_state = aln->ref_seq[static_cast<std::vector<cmaple::StateType>
                                        ::size_type>(end_pos)];
          
        if (seq1_region->getLH(seq2_state) > 0.1)
              seq2_more_info = true;
        else
              return 0;
        
      } else if (seq2_region->type == TYPE_O) {
        StateType seq1_state = seq1_region->type;
        if (seq1_state == TYPE_R)
            seq1_state = aln->ref_seq[static_cast<std::vector<cmaple::StateType>
                                        ::size_type>(end_pos)];
          
        if (seq2_region->getLH(seq1_state) > 0.1)
              seq1_more_info = true;
        else
              return 0;
      } else {
        return 0;
      }
    }
    // Both regions are type O
    else if (seq1_region->type == TYPE_O) {
      for (StateType i = 0; i < num_states; ++i) {
        if (seq2_region->getLH(i) > 0.1 && seq1_region->getLH(i) < 0.1) {
          seq1_more_info = true;
        } else if (seq1_region->getLH(i) > 0.1 && seq2_region->getLH(i) < 0.1) {
          seq2_more_info = true;
        }
      }
    }

    // update pos
    pos = end_pos + 1;
  }

  // return result
  if (seq1_more_info) {
    if (seq2_more_info) {
      return 0;
    } else {
      return 1;
    }
  } else if (seq2_more_info) {
    return -1;
  } else {
    return 1;
  }

  // return 0;
}

auto cmaple::SeqRegions::areDiffFrom(
    const std::unique_ptr<SeqRegions>& regions2,
    PositionType seq_length,
    StateType num_states,
    const Params& params) const -> bool {

  assert(seq_length > 0);
  assert(num_states > 0);
  assert(size() > 0);
        
  if (this->size() != regions2->size()) {
    return true;
  }

  // init variables
  PositionType pos = 0;
  const SeqRegions& seq1_regions = *this;
  const SeqRegions& seq2_regions = *regions2;
  size_t iseq1 = 0;
  size_t iseq2 = 0;

  // compare each pair of regions
  while (pos < seq_length) {
    PositionType end_pos = 0;

    // get the next shared segment in the two sequences
    cmaple::SeqRegions::getNextSharedSegment(pos, seq1_regions, seq2_regions,
                                             iseq1, iseq2, end_pos);
    const auto* const seq1_region = &seq1_regions[iseq1];
    const auto* const seq2_region = &seq2_regions[iseq2];

    // if any pair of regions is different in type -> the two sequences are
    // different
    if (seq1_region->type != seq2_region->type) {
      return true;
    }

    // type R/ACGT
    if (seq1_region->type < num_states || seq1_region->type == TYPE_R) {
      // compare plength_from_root and plength_observation
      if (fabs(seq1_region->plength_observation2root -
               seq2_region->plength_observation2root) > params.threshold_prob ||
          fabs(seq1_region->plength_observation2node -
               seq2_region->plength_observation2node) > params.threshold_prob) {
        return true;
      }
    }

    // type O
    if (seq1_region->type == TYPE_O) {
      // compare plength_observation
      if (fabs(seq1_region->plength_observation2node -
               seq2_region->plength_observation2node) > params.threshold_prob) {
        return true;
      }

      // compare likelihood of each state
      for (StateType i = 0; i < num_states; ++i) {
        RealNumType diff = fabs(seq1_region->getLH(i) - seq2_region->getLH(i));

        if (diff > 0) {
          if ((seq1_region->getLH(i) == 0) || (seq2_region->getLH(i) == 0)) {
            return true;
          }

          if (diff > params.thresh_diff_update ||
              (diff > params.threshold_prob &&
               ((diff >
                 params.thresh_diff_fold_update * seq1_region->getLH(i)) ||
                (diff >
                 params.thresh_diff_fold_update * seq2_region->getLH(i))))) {
            return true;
          }
        }
      }
    }

    // update pos
    pos = end_pos + 1;
  }

  return false;
}

auto cmaple::SeqRegions::countSharedSegments(const SeqRegions& seq2_regions,
                                             const size_t seq_length) const
    -> size_t {
  const SeqRegions& seq1_regions = *this;
  size_t count{};
  PositionType pos{};
  size_t iseq1 = 0;
  size_t iseq2 = 0;

  while (pos < static_cast<PositionType>(seq_length)) {
    PositionType end_pos{};

    // get the next shared segment in the two sequences
    cmaple::SeqRegions::getNextSharedSegment(pos, seq1_regions, seq2_regions,
                                             iseq1, iseq2, end_pos);
    ++count;
    pos = end_pos + 1;
  }
  return ++count;
}

void cmaple::merge_N_RACGT(const SeqRegion& reg_racgt,
                           const RealNumType lower_plength,
                           const PositionType end_pos,
                           const RealNumType threshold_prob,
                           SeqRegions& merged_regions) {
  RealNumType plength_observation2node = -1;
  RealNumType plength_observation2root = 0;

  if (reg_racgt.plength_observation2node >= 0) {
    RealNumType new_plength = reg_racgt.plength_observation2node;
    if (lower_plength > 0) {
      new_plength += lower_plength;
    }

    plength_observation2node = new_plength;
  } else {
    if (lower_plength > 0) {
      plength_observation2node = lower_plength;
    } else {
      plength_observation2root = -1;
    }
  }

  // add a new region and try to merge consecutive R regions together
  cmaple::SeqRegions::addNonConsecutiveRRegion(
      merged_regions, reg_racgt.type, plength_observation2node,
      plength_observation2root, end_pos, threshold_prob);
}

void cmaple::merge_RACGT_N(const SeqRegion& reg_n,
                           const RealNumType upper_plength,
                           const PositionType end_pos,
                           const RealNumType threshold_prob,
                           SeqRegions& merged_regions) {
  assert(reg_n.type != TYPE_N && reg_n.type != TYPE_O);
    
  // todo rewrite
  RealNumType plength_observation2node = -1;
  RealNumType plength_observation2root = -1;

  if (reg_n.plength_observation2root >= 0) {
    RealNumType plength_from_root = reg_n.plength_observation2root;
    if (upper_plength > 0) {
      plength_from_root += upper_plength;
    }

    plength_observation2node = reg_n.plength_observation2node;
    plength_observation2root = plength_from_root;
  } else if (reg_n.plength_observation2node >= 0) {
    RealNumType plength_observation = reg_n.plength_observation2node;
    if (upper_plength > 0) {
      plength_observation += upper_plength;
    }

    plength_observation2node = plength_observation;
  } else if (upper_plength > 0) {
    plength_observation2node = upper_plength;
  }

  // add a new region and try to merge consecutive R regions together
  cmaple::SeqRegions::addNonConsecutiveRRegion(
      merged_regions, reg_n.type, plength_observation2node,
      plength_observation2root, end_pos, threshold_prob);
}

auto cmaple::merge_Zero_Distance(const SeqRegion& seq1_region,
                                 const SeqRegion& seq2_region,
                                 const RealNumType total_blength_1,
                                 const RealNumType total_blength_2,
                                 const PositionType end_pos,
                                 const RealNumType threshold_prob,
                                 const StateType num_states,
                                 std::unique_ptr<SeqRegions>& merged_regions)
    -> bool {

  assert(num_states > 0);
  
  // due to 0 distance, the entry will be of same type as entry2
  if ((seq2_region.type < num_states || seq2_region.type == TYPE_R) &&
      total_blength_2 <= 0) {
    if ((seq1_region.type < num_states || seq1_region.type == TYPE_R) &&
        total_blength_1 <= 0) {
      // outError("Sorry! something went wrong. DEBUG: ((seq2_region->type <
      // num_states || seq2_region->type == TYPE_R) && total_blength_2 == 0) &&
      // ((seq1_region->type < num_states || seq1_region->type == TYPE_R) &&
      // total_blength_1 == 0)");
      merged_regions = nullptr;
      return true;
    }

    // add a new region and try to merge consecutive R regions together
    cmaple::SeqRegions::addNonConsecutiveRRegion(
        *merged_regions, seq2_region.type, -1, -1, end_pos, threshold_prob);
    return true;
  }
  // due to 0 distance, the entry will be of same type as entry1
  else if ((seq1_region.type < num_states || seq1_region.type == TYPE_R) &&
           total_blength_1 <= 0) {
    // add a new region and try to merge consecutive R regions together
    cmaple::SeqRegions::addNonConsecutiveRRegion(
        *merged_regions, seq1_region.type, -1, -1, end_pos, threshold_prob);
    return true;
  }

  return false;
}

void cmaple::merge_N_O_TwoLowers(const SeqRegion& seq2_region,
                                 const PositionType end_pos,
                                 const RealNumType plength2,
                                 SeqRegions& merged_regions) {
  assert(seq2_region.type == TYPE_O);
  
  // add merged region into merged_regions
  SeqRegion new_region = SeqRegion::clone(seq2_region);
  new_region.position = end_pos;
  if (seq2_region.plength_observation2node >= 0) {
    if (plength2 > 0) {
      new_region.plength_observation2node += plength2;
    }
  } else {
    if (plength2 > 0) {
      new_region.plength_observation2node = plength2;
    }
  }
  merged_regions.push_back(std::move(new_region));
}

void cmaple::merge_N_RACGT_TwoLowers(const SeqRegion& seq2_region,
                                     const PositionType end_pos,
                                     const RealNumType plength2,
                                     const RealNumType threshold_prob,
                                     SeqRegions& merged_regions) {
  assert(seq2_region.type != TYPE_N && seq2_region.type != TYPE_O);
  
  RealNumType plength_observation2node = -1;

  if (seq2_region.plength_observation2node >= 0) {
    RealNumType new_plength = seq2_region.plength_observation2node;
    if (plength2 > 0) {
      new_plength += plength2;
    }

    plength_observation2node = new_plength;
  } else if (plength2 > 0) {
    plength_observation2node = plength2;
  }

  // add a new region and try to merge consecutive R regions together
  cmaple::SeqRegions::addNonConsecutiveRRegion(merged_regions, seq2_region.type,
                                               plength_observation2node, -1,
                                               end_pos, threshold_prob);
}

void cmaple::merge_identicalRACGT_TwoLowers(
    const SeqRegion& seq1_region,
    const PositionType end_pos,
    RealNumType total_blength_1,
    RealNumType total_blength_2,
    const PositionType pos,
    const RealNumType threshold_prob,
    const ModelBase* model,
    const RealNumType* const cumulative_rate,
    RealNumType& log_lh,
    SeqRegions& merged_regions,
    const bool return_log_lh) {
  assert(seq1_region.type != TYPE_N && seq1_region.type != TYPE_O);
  assert(model);
  assert(cumulative_rate);
    
  // add a new region and try to merge consecutive R regions together
  cmaple::SeqRegions::addNonConsecutiveRRegion(merged_regions, seq1_region.type,
                                               -1, -1, end_pos, threshold_prob);

  if (return_log_lh) {
    // convert total_blength_1 and total_blength_2 to zero if they are -1
    if (total_blength_1 < 0) {
      total_blength_1 = 0;
    }
    if (total_blength_2 < 0) {
      total_blength_2 = 0;
    }

    if (seq1_region.type == TYPE_R) {
      auto prev_log_lh = log_lh;
      log_lh += (total_blength_1 + total_blength_2) *
                (cumulative_rate[end_pos + 1] - cumulative_rate[pos]);
    } else {
      log_lh += model->getDiagonalMutationMatrixEntry(seq1_region.type, end_pos) *
                (total_blength_1 + total_blength_2);
    }
  }
}

void cmaple::calSiteLhs_identicalRACGT(std::vector<RealNumType>&site_lh_contributions,
                               const SeqRegion& seq1_region,
                               const PositionType end_pos,
                               RealNumType total_blength_1,
                               RealNumType total_blength_2,
                               const PositionType pos,
                               const RealNumType threshold_prob,
                               const ModelBase* model,
                               const RealNumType* const cumulative_rate,
                               RealNumType& log_lh,
                               SeqRegions& merged_regions) {
  assert(seq1_region.type != TYPE_N && seq1_region.type != TYPE_O);
  assert(model);
  assert(cumulative_rate);
    
  // add a new region and try to merge consecutive R regions together
  cmaple::SeqRegions::addNonConsecutiveRRegion(merged_regions, seq1_region.type,
                                               -1, -1, end_pos, threshold_prob);

  // compute the (site) lh contributions
  // convert total_blength_1 and total_blength_2 to zero if they are -1
  if (total_blength_1 < 0) {
    total_blength_1 = 0;
  }
  if (total_blength_2 < 0) {
    total_blength_2 = 0;
  }

  RealNumType total_blength = total_blength_1 + total_blength_2;
  if (seq1_region.type == TYPE_R) {
    log_lh +=
        total_blength * (cumulative_rate[end_pos + 1] - cumulative_rate[pos]);

    // compute site lh contributions
    for (PositionType i = pos; i < end_pos + 1; ++i) {
      site_lh_contributions[static_cast<std::vector<RealNumType>
                            ::size_type>(i)] += total_blength
        * (cumulative_rate[i + 1] - cumulative_rate[i]);
    }
  } else {
    log_lh += model->getDiagonalMutationMatrixEntry(seq1_region.type, pos) * total_blength;

    // compute site lh contributions
    site_lh_contributions[static_cast<std::vector<RealNumType>::size_type>(pos)] +=
        model->getDiagonalMutationMatrixEntry(seq1_region.type, pos) * total_blength;
  }
}

auto cmaple::SeqRegions::operator==(const SeqRegions& seqregions_1) const
    -> bool {
  if (size() != seqregions_1.size()) {
    return false;
  }

  for (std::vector<cmaple::SeqRegion>::size_type i = 0; i < size(); ++i) {
    if (!(at(i) == seqregions_1[i])) {
      return false;
    }
  }

  return true;
}

void cmaple::SeqRegions::addNonConsecutiveRRegion(
    SeqRegions& regions,
    const cmaple::StateType new_region_type,
    const cmaple::RealNumType plength_observation2node,
    const cmaple::RealNumType plength_observation2root,
    const cmaple::PositionType end_pos,
    const cmaple::RealNumType threshold_prob) {
  // cannot merge consecutive R regions if no region exists in regions
  if (!regions.empty()) {
    // try to merge consecutive R regions
    SeqRegion& last_region = regions.back();
    if (new_region_type == cmaple::TYPE_R &&
        last_region.type == cmaple::TYPE_R &&
        fabs(last_region.plength_observation2node - plength_observation2node) <
            threshold_prob &&
        fabs(last_region.plength_observation2root - plength_observation2root) <
            threshold_prob) {
      last_region.position = end_pos;
      last_region.plength_observation2node = plength_observation2node;
      last_region.plength_observation2root = plength_observation2root;
      return;
    }
  }

  // if we cannot merge new region into existing R region => just add a new one
  regions.emplace_back(new_region_type, end_pos, plength_observation2node,
                       plength_observation2root);
}

auto cmaple::SeqRegions::simplifyO(cmaple::RealNumType* const partial_lh,
                                   cmaple::StateType ref_state,
                                   cmaple::StateType num_states,
                                   cmaple::RealNumType threshold)
    -> cmaple::StateType {
  // dummy variables
  assert(partial_lh);
  assert(num_states > 0);
  cmaple::RealNumType max_prob = 0;
  cmaple::StateType max_index = 0;
  cmaple::StateType high_prob_count = 0;

  // Check all states one by one
  for (cmaple::StateType i = 0; i < num_states; ++i) {
    // record the state with the highest likelihood
    if (partial_lh[i] > max_prob) {
      max_prob = partial_lh[i];
      max_index = i;
    }

    // count the number of states that have the likelihood greater than a
    // threshold
    if (partial_lh[i] > threshold) {
      ++high_prob_count;
    }
  }

  // if the partial lh concentrates at a specific state -> return new state
  if (high_prob_count == 1) {
    // new state matches with the reference state
    if (max_index == ref_state) {
      return cmaple::TYPE_R;
      // return a nucleotide
    } else {
      return max_index;
    }
  }
  // otherwise, cannot simplify
  else {
    return cmaple::TYPE_O;
  }
}

void cmaple::SeqRegions::addSimplifiedO(
    const cmaple::PositionType end_pos,
    SeqRegion::LHType& new_lh,
    const Alignment* aln,
    const cmaple::RealNumType threshold_prob,
    SeqRegions& merged_regions) {
  assert(aln);
    
  cmaple::StateType new_state = SeqRegions::simplifyO(
      new_lh.data(), aln->ref_seq[static_cast<std::vector<cmaple::StateType>
                                    ::size_type>(end_pos)], aln->num_states, threshold_prob);

  if (new_state == cmaple::TYPE_O) {
    merged_regions.emplace_back(cmaple::TYPE_O, end_pos, 0, 0,
                                std::move(new_lh));
  } else {
    // add a new region and try to merge consecutive R regions together
    SeqRegions::addNonConsecutiveRRegion(merged_regions, new_state, -1, -1,
                                         end_pos, threshold_prob);
  }
}
