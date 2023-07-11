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

// explicit instantiation of templates
template double cmaple::SeqRegions::computeAbsoluteLhAtRoot<2>(const ModelBase*);
template double cmaple::SeqRegions::computeAbsoluteLhAtRoot<4>(const ModelBase*);
template double cmaple::SeqRegions::computeAbsoluteLhAtRoot<20>(const ModelBase*);

template double cmaple::SeqRegions::computeSiteLhAtRoot<2>(std::vector<RealNumType>&, const ModelBase*);
template double cmaple::SeqRegions::computeSiteLhAtRoot<4>(std::vector<RealNumType>&, const ModelBase*);
template double cmaple::SeqRegions::computeSiteLhAtRoot<20>(std::vector<RealNumType>&, const ModelBase*);

template double cmaple::SeqRegions::mergeTwoLowers<2>(std::unique_ptr<SeqRegions>&, const RealNumType, const SeqRegions&, const RealNumType, const AlignmentBase*, const ModelBase*, const RealNumType, const bool) const;
template double cmaple::SeqRegions::mergeTwoLowers<4>(std::unique_ptr<SeqRegions>&, const RealNumType, const SeqRegions&, const RealNumType, const AlignmentBase*, const ModelBase*, const RealNumType, const bool) const;
template double cmaple::SeqRegions::mergeTwoLowers<20>(std::unique_ptr<SeqRegions>&, const RealNumType, const SeqRegions&, const RealNumType, const AlignmentBase*, const ModelBase*, const RealNumType, const bool) const;

template void cmaple::SeqRegions::mergeUpperLower<2>(std::unique_ptr<SeqRegions>&, RealNumType, const SeqRegions&, RealNumType, const AlignmentBase*, const ModelBase*, const RealNumType) const;
template void cmaple::SeqRegions::mergeUpperLower<4>(std::unique_ptr<SeqRegions>&, RealNumType, const SeqRegions&, RealNumType, const AlignmentBase*, const ModelBase*, const RealNumType) const;
template void cmaple::SeqRegions::mergeUpperLower<20>(std::unique_ptr<SeqRegions>&, RealNumType, const SeqRegions&, RealNumType, const AlignmentBase*, const ModelBase*, const RealNumType) const;

template void cmaple::SeqRegions::computeTotalLhAtRoot<2>(std::unique_ptr<SeqRegions>&, const ModelBase*, RealNumType) const;
template void cmaple::SeqRegions::computeTotalLhAtRoot<4>(std::unique_ptr<SeqRegions>&, const ModelBase*, RealNumType) const;
template void cmaple::SeqRegions::computeTotalLhAtRoot<20>(std::unique_ptr<SeqRegions>&, const ModelBase*, RealNumType) const;

template double cmaple::SeqRegions::calculateSiteLhContributions<2>(std::vector<RealNumType>&, std::unique_ptr<SeqRegions>&, const RealNumType, const SeqRegions&, const RealNumType, const AlignmentBase*, const ModelBase*, const RealNumType) const;
template double cmaple::SeqRegions::calculateSiteLhContributions<4>(std::vector<RealNumType>&, std::unique_ptr<SeqRegions>&, const RealNumType, const SeqRegions&, const RealNumType, const AlignmentBase*, const ModelBase*, const RealNumType) const;
template double cmaple::SeqRegions::calculateSiteLhContributions<20>(std::vector<RealNumType>&, std::unique_ptr<SeqRegions>&, const RealNumType, const SeqRegions&, const RealNumType, const AlignmentBase*, const ModelBase*, const RealNumType) const;


template <const StateType num_states>
RealNumType updateLHwithModel(const ModelBase* model,
  const SeqRegion::LHType& prior,
  SeqRegion::LHType& posterior, const RealNumType total_blength)
{
    RealNumType sum_lh = 0;
    RealNumType* mutation_mat_row = model->mutation_mat;
    for (StateType i = 0; i < num_states; ++i, mutation_mat_row += num_states)
    {
        RealNumType tot = 0;
        if (total_blength > 0) // TODO: avoid
        {
            tot += dotProduct<num_states>(&(prior)[0], mutation_mat_row);

            tot *= total_blength;
        }

        tot += prior[i];
        posterior[i] = tot * model->root_freqs[i];
        sum_lh += posterior[i];
    }
    return sum_lh;
}

template <const StateType num_states>
RealNumType updateLHwithMat(const RealNumType* mat_row,
  const SeqRegion::LHType& prior,
  SeqRegion::LHType& posterior, const RealNumType total_blength)
{
    RealNumType sum_lh = 0;
    for (StateType i = 0; i < num_states; ++i, mat_row += num_states)
    {
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
RealNumType updateMultLHwithMat(const RealNumType* mat_row,
  const SeqRegion::LHType& prior,
  SeqRegion::LHType& posterior, const RealNumType total_blength)
{
    RealNumType sum_lh = 0;
    for (StateType i = 0; i < num_states; ++i, mat_row += num_states)
    {
        RealNumType tot = 0;
        if (total_blength > 0) // TODO: avoid
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

cmaple::SeqRegions::SeqRegions(const std::unique_ptr<SeqRegions>& n_regions)
{
  if (!n_regions)
  {
    outError("oops");
  }
    // clone regions one by one
    reserve(n_regions->size());
    for (const auto& region : *n_regions)
      push_back(SeqRegion::clone(region));
} 

int cmaple::SeqRegions::compareWithSample(const SeqRegions& sequence2, PositionType seq_length, StateType num_states) const
{
    ASSERT(seq_length > 0);
    
    // init dummy variables
    bool seq1_more_info = false;
    bool seq2_more_info = false;
    PositionType pos = 0;
    const SeqRegions& seq1_regions = *this;
    const SeqRegions& seq2_regions = sequence2;
    size_t iseq1 = 0;
    size_t iseq2 = 0;

    while (pos < seq_length && (!seq1_more_info || !seq2_more_info))
    {
        PositionType end_pos;
        
        // get the next shared segment in the two sequences
        cmaple::SeqRegions::getNextSharedSegment(pos, seq1_regions, seq2_regions, iseq1, iseq2, end_pos);
        const auto* const seq1_region = &seq1_regions[iseq1];
        const auto* const seq2_region = &seq2_regions[iseq2];
        // The two regions have different types from each other
        if (seq1_region->type != seq2_region->type)
        {
            if (seq1_region->type == TYPE_N)
                seq2_more_info = true;
            else
                if (seq2_region->type == TYPE_N)
                    seq1_more_info = true;
                else if (seq1_region->type == TYPE_O)
                        seq2_more_info = true;
                    else
                        if (seq2_region->type == TYPE_O)
                            seq1_more_info = true;
                        else
                        {
                            seq1_more_info = true;
                            seq2_more_info = true;
                        }
        }
        // Both regions are type O
        else if (seq1_region->type == TYPE_O)
        {
            for (StateType i = 0; i < num_states; ++i)
            {
                if (seq2_region->getLH(i) > 0.1 && seq1_region->getLH(i) < 0.1)
                    seq1_more_info = true;
                else if (seq1_region->getLH(i) > 0.1 && seq2_region->getLH(i) < 0.1)
                    seq2_more_info = true;
            }
        }

        // update pos
        pos = end_pos + 1;
    }

    // return result
    if (seq1_more_info)
        if (seq2_more_info)
            return 0;
        else
            return 1;
    else
        if (seq2_more_info)
            return -1;
        else
            return 1;
    
    return 0;
}

bool cmaple::SeqRegions::areDiffFrom(const SeqRegions& regions2, PositionType seq_length, StateType num_states, const Params& params) const
{
    if (regions2.empty())
        return true;
    
    // init variables
    PositionType pos = 0;
    const SeqRegions& seq1_regions = *this;
    const SeqRegions& seq2_regions = regions2;
    size_t iseq1 = 0;
    size_t iseq2 = 0;
    
    // compare each pair of regions
    while (pos < seq_length)
    {
        PositionType end_pos;
        
        // get the next shared segment in the two sequences
        cmaple::SeqRegions::getNextSharedSegment(pos, seq1_regions, seq2_regions, iseq1, iseq2, end_pos);
        const auto* const seq1_region = &seq1_regions[iseq1];
        const auto* const seq2_region = &seq2_regions[iseq2];

        // if any pair of regions is different in type -> the two sequences are different
        if (seq1_region->type != seq2_region->type)
            return true;
        
        // type R/ACGT
        if (seq1_region->type < num_states || seq1_region->type == TYPE_R)
        {
            // compare plength_from_root and plength_observation
            if (fabs(seq1_region->plength_observation2root - seq2_region->plength_observation2root) > params.threshold_prob
                ||fabs(seq1_region->plength_observation2node - seq2_region->plength_observation2node) > params.threshold_prob)
                return true;
        }
        
        // type O
        if (seq1_region->type == TYPE_O)
        {
            // compare plength_observation
            if (fabs(seq1_region->plength_observation2node - seq2_region->plength_observation2node) > params.threshold_prob)
                return true;
            
            // compare likelihood of each state
            for (StateType i = 0; i < num_states; ++i)
            {
                RealNumType diff = fabs(seq1_region->getLH(i) - seq2_region->getLH(i));
                
                if (diff > 0)
                {
                    if ((seq1_region->getLH(i) == 0) || (seq2_region->getLH(i) == 0))
                        return true;
                    
                    if (diff > params.thresh_diff_update
                        || (diff > params.threshold_prob
                            && ((diff > params.thresh_diff_fold_update * seq1_region->getLH(i))
                                || (diff > params.thresh_diff_fold_update * seq2_region->getLH(i)))))
                        return true;
                }
            }
        }
        
        // update pos
        pos = end_pos + 1;
    }
    
    return false;
}

size_t cmaple::SeqRegions::countSharedSegments(const SeqRegions& seq2_regions, const size_t seq_length) const
{
    const SeqRegions& seq1_regions = *this;
    size_t count{};
    PositionType pos{};
    size_t iseq1 = 0;
    size_t iseq2 = 0;
    
    while (pos < seq_length)
    {
        PositionType end_pos{};
        
        // get the next shared segment in the two sequences
        cmaple::SeqRegions::getNextSharedSegment(pos, seq1_regions, seq2_regions, iseq1, iseq2, end_pos);
        ++count;
        pos = end_pos + 1;
    }
    return ++count;
}

template <const StateType num_states>
void cmaple::merge_N_O(const RealNumType lower_plength, const SeqRegion& reg_o, const ModelBase* model,
               const PositionType end_pos, SeqRegions& merged_target)
{
  RealNumType total_blength = lower_plength;
  if (reg_o.plength_observation2node >= 0)
    total_blength = reg_o.plength_observation2node + (lower_plength > 0 ? lower_plength : 0);

  auto new_lh = std::make_unique<SeqRegion::LHType>(); // = new RealNumType[num_states];
  RealNumType sum_lh = updateLHwithModel<num_states>(model, *reg_o.likelihood, (*new_lh), total_blength);
  // normalize the new partial likelihood
  normalize_arr(new_lh->data(), num_states, sum_lh);

  // add merged region into merged_regions
  merged_target.emplace_back(TYPE_O, end_pos, 0, 0, std::move(new_lh));
}

void cmaple::merge_N_RACGT(const SeqRegion& reg_racgt, const RealNumType lower_plength, const PositionType end_pos,
                   const RealNumType threshold_prob, SeqRegions& merged_regions)
{
  RealNumType plength_observation2node = -1;
  RealNumType plength_observation2root = 0;

  if (reg_racgt.plength_observation2node >= 0)
  {
    RealNumType new_plength = reg_racgt.plength_observation2node;
    if (lower_plength > 0)
      new_plength += lower_plength;

    plength_observation2node = new_plength;
  }
  else
  {
    if (lower_plength > 0)
      plength_observation2node = lower_plength;
    else
      plength_observation2root = -1;
  }

  // add a new region and try to merge consecutive R regions together
  cmaple::SeqRegions::addNonConsecutiveRRegion(merged_regions, reg_racgt.type, plength_observation2node, plength_observation2root, end_pos, threshold_prob);
}

template <const StateType num_states>
void cmaple::merge_O_N(const SeqRegion& reg_o, const RealNumType upper_plength, const PositionType end_pos, const ModelBase* model, SeqRegions& merged_regions)
{
  /*
   NhanLT: total_blength may be wrong if upper_plength == -1 and reg_o.plength_observation2node > 0
  RealNumType total_blength = upper_plength;
  if (reg_o.plength_observation2node > 0)
  {
    total_blength += reg_o.plength_observation2node;
  }*/
    
    RealNumType total_blength = -1;
    
    if (reg_o.plength_observation2node >= 0)
    {
        total_blength = reg_o.plength_observation2node;
        if (upper_plength > 0)
            total_blength += upper_plength;
    }
    else if (upper_plength > 0)
        total_blength = upper_plength;


  if (total_blength > 0)
  {
    auto new_lh = std::make_unique<SeqRegion::LHType>(); // = new RealNumType[num_states];
    RealNumType sum_lh = updateLHwithMat<num_states>(model->transposed_mut_mat, *(reg_o.likelihood), *new_lh, total_blength);

    // normalize the new partial likelihood
    normalize_arr(new_lh->data(), num_states, sum_lh);

    // add merged region into merged_regions
    merged_regions.emplace_back(TYPE_O, end_pos, 0, 0, std::move(new_lh));
  }
  else
  {
    // add merged region into merged_regions
    merged_regions.emplace_back(TYPE_O, end_pos, 0, 0, *(reg_o.likelihood));
  }
}

void cmaple::merge_RACGT_N(const SeqRegion& reg_n, const RealNumType upper_plength, const PositionType end_pos,
  const RealNumType threshold_prob, SeqRegions& merged_regions)
{
  // todo rewrite
  RealNumType plength_observation2node = -1;
  RealNumType plength_observation2root = -1;

  if (reg_n.plength_observation2root >= 0)
  {
    RealNumType plength_from_root = reg_n.plength_observation2root;
    if (upper_plength > 0)
      plength_from_root += upper_plength;

    plength_observation2node = reg_n.plength_observation2node;
    plength_observation2root = plength_from_root;
  }
  else if (reg_n.plength_observation2node >= 0)
  {
    RealNumType plength_observation = reg_n.plength_observation2node;
    if (upper_plength > 0)
      plength_observation += upper_plength;

    plength_observation2node = plength_observation;
  }
  else if (upper_plength > 0)
    plength_observation2node = upper_plength;

  // add a new region and try to merge consecutive R regions together
  cmaple::SeqRegions::addNonConsecutiveRRegion(merged_regions, reg_n.type, plength_observation2node, plength_observation2root, end_pos, threshold_prob);
}

bool cmaple::merge_Zero_Distance(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType total_blength_1, const RealNumType total_blength_2, const PositionType end_pos, const RealNumType threshold_prob, const StateType num_states, std::unique_ptr<SeqRegions>& merged_regions)
{
    // due to 0 distance, the entry will be of same type as entry2
    if ((seq2_region.type < num_states || seq2_region.type == TYPE_R) && total_blength_2 <= 0)
    {
        if ((seq1_region.type < num_states || seq1_region.type == TYPE_R) && total_blength_1 <= 0)
        {
            //outError("Sorry! something went wrong. DEBUG: ((seq2_region->type < num_states || seq2_region->type == TYPE_R) && total_blength_2 == 0) && ((seq1_region->type < num_states || seq1_region->type == TYPE_R) && total_blength_1 == 0)");
            merged_regions = nullptr;
            return true;
        }
        
        // add a new region and try to merge consecutive R regions together
        cmaple::SeqRegions::addNonConsecutiveRRegion(*merged_regions, seq2_region.type, -1, -1, end_pos, threshold_prob);
        return true;
    }
    // due to 0 distance, the entry will be of same type as entry1
    else if ((seq1_region.type < num_states || seq1_region.type == TYPE_R) && total_blength_1 <= 0)
    {
        // add a new region and try to merge consecutive R regions together
        cmaple::SeqRegions::addNonConsecutiveRRegion(*merged_regions, seq1_region.type, -1, -1, end_pos, threshold_prob);
        return true;
    }
    
    return false;
}

template <const StateType num_states>
void cmaple::merge_O_ORACGT(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType total_blength_1, const RealNumType total_blength_2, const PositionType end_pos, const RealNumType threshold_prob, const ModelBase* model, const AlignmentBase* aln, SeqRegions& merged_regions)
{
    auto new_lh = std::make_unique<SeqRegion::LHType>(); // = new RealNumType[num_states];
    auto& new_lh_value = *new_lh;

    // if total_blength_1 > 0 => compute new partial likelihood
    if (total_blength_1 > 0)
      updateLHwithMat<num_states>(model->transposed_mut_mat, *(seq1_region.likelihood), *new_lh, total_blength_1);
    // otherwise, clone the partial likelihood from seq1
    else
      *new_lh = *seq1_region.likelihood;
    
    RealNumType sum_new_lh = 0;

    // seq1 = seq2 = O
    if (seq2_region.type == TYPE_O)
    {
        RealNumType* mutation_mat_row = model->mutation_mat;
        sum_new_lh = updateMultLHwithMat<num_states>(model->mutation_mat, *(seq2_region.likelihood), *new_lh, total_blength_2);
    }
    // seq1 = "O" and seq2 = R or ACGT
    else
    {
        StateType seq2_state = seq2_region.type;
        if (seq2_state == TYPE_R)
            seq2_state = aln->ref_seq[end_pos];
        
        if (total_blength_2 > 0)
        {
            RealNumType* transposed_mut_mat_row = model->transposed_mut_mat + model->row_index[seq2_state];
            sum_new_lh += updateVecWithState<num_states>(new_lh_value.data(), seq2_state, transposed_mut_mat_row, total_blength_2);
        }
        else
            sum_new_lh += resetLhVecExceptState<num_states>(new_lh_value.data(), seq2_state, new_lh_value[seq2_state]);
    }
    
    // normalize the new partial lh
    if (sum_new_lh == 0)
        outError("Sum of new partital lh is zero.");
    
    // normalize the new partial likelihood
    normalize_arr(new_lh->data(), num_states, sum_new_lh);
    cmaple::SeqRegions::addSimplifiedO(end_pos, new_lh_value, aln, threshold_prob, merged_regions);
}

template <const StateType num_states>
void cmaple::merge_RACGT_O(const SeqRegion& seq2_region, const RealNumType total_blength_2, const PositionType end_pos, SeqRegion::LHType& new_lh, const RealNumType threshold_prob, const ModelBase* model, const AlignmentBase* aln, SeqRegions& merged_regions)
{
    RealNumType sum_new_lh = updateMultLHwithMat<num_states>(model->mutation_mat, *(seq2_region.likelihood), new_lh, total_blength_2);
    
    // normalize the new partial likelihood
    normalize_arr(new_lh.data(), num_states, sum_new_lh);
    cmaple::SeqRegions::addSimplifiedO(end_pos, new_lh, aln, threshold_prob, merged_regions);
}

template <const StateType num_states>
void cmaple::merge_RACGT_RACGT(const SeqRegion& seq2_region, const RealNumType total_blength_2, const PositionType end_pos, SeqRegion::LHType& new_lh, const ModelBase* model, const AlignmentBase* aln, SeqRegions& merged_regions)
{
    RealNumType sum_new_lh = 0;
    StateType seq2_state = seq2_region.type;
    
    if (seq2_state == TYPE_R)
        seq2_state = aln->ref_seq[end_pos];
    
    // TODO: this seems a weird operation on `new_lh_value` (since it was just created anew and is all 00000)!
    // please check this makes sense!
    // CHECKED: NHANLT: new_lh_value or new_lh is already initialized or updated in the section "init or update new_lh/new_lh_value" above
    if (total_blength_2 > 0)
    {
        RealNumType* transposed_mut_mat_row = model->transposed_mut_mat + model->row_index[seq2_state];
        sum_new_lh = updateVecWithState<num_states>(new_lh.data(), seq2_state, transposed_mut_mat_row, total_blength_2);
    }
    else
        sum_new_lh = resetLhVecExceptState<num_states>(new_lh.data(), seq2_state, new_lh[seq2_state]);
    
    // normalize the new partial likelihood
    normalize_arr(new_lh.data(), num_states, sum_new_lh);

    // add new region into the merged regions
    merged_regions.emplace_back(TYPE_O, end_pos, 0, 0, std::move(new_lh));
}

template <const StateType num_states>
void cmaple::merge_RACGT_ORACGT(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType total_blength_1, const RealNumType total_blength_2, const RealNumType upper_plength, const PositionType end_pos, const RealNumType threshold_prob, const ModelBase* model, const AlignmentBase* aln, SeqRegions& merged_regions)
{
    StateType seq1_state = seq1_region.type;
    if (seq1_state == TYPE_R)
        seq1_state = aln->ref_seq[end_pos];
    
    auto new_lh = std::make_unique<SeqRegion::LHType>(); // = new RealNumType[num_states];
    auto& new_lh_value = *new_lh;
    
    // init or update new_lh/new_lh_value
    if (seq1_region.plength_observation2root >= 0)
    {
        RealNumType length_to_root = seq1_region.plength_observation2root;
        if (upper_plength > 0)
            length_to_root += upper_plength;
        SeqRegion::LHType root_vec;
        memcpy(root_vec.data(), model->root_freqs, sizeof(RealNumType)* num_states);
        
        RealNumType* transposed_mut_mat_row = model->transposed_mut_mat + model->row_index[seq1_state];
        updateVecWithState<num_states>(root_vec.data(), seq1_state, transposed_mut_mat_row, seq1_region.plength_observation2node);
        
        updateLHwithMat<num_states>(model->transposed_mut_mat, root_vec, *new_lh, length_to_root);
    }
    else
    {
        if (total_blength_1 > 0)
        {
            RealNumType* mutation_mat_row = model->mutation_mat + model->row_index[seq1_state];
            setVecWithState<num_states>(new_lh_value.data(), seq1_state, mutation_mat_row, total_blength_1);
        }
        else
        {
          for (auto& v : new_lh_value) v = 0;
          new_lh_value[seq1_state] = 1;
        }
    }
      
    // seq1 = R/ACGT and seq2 = "O"
    if (seq2_region.type == TYPE_O)
    {
        merge_RACGT_O<num_states>(seq2_region, total_blength_2, end_pos, *new_lh, threshold_prob, model, aln, merged_regions);
    }
    // seq1 = R/ACGT and different from seq2 = R/ACGT
    else
    {
        merge_RACGT_RACGT<num_states>(seq2_region, total_blength_2, end_pos, *new_lh, model, aln, merged_regions);
    }
}

using DoubleState = uint16_t;
static constexpr DoubleState NN = (DoubleState(TYPE_N) << 8) | TYPE_N;
static constexpr DoubleState NO = (DoubleState(TYPE_N) << 8) | TYPE_O;
static constexpr DoubleState RR = (DoubleState(TYPE_R) << 8) | TYPE_R;
static constexpr DoubleState RO = (DoubleState(TYPE_R) << 8) | TYPE_O;
static constexpr DoubleState OO = (DoubleState(TYPE_O) << 8) | TYPE_O;
static constexpr DoubleState ON = (DoubleState(TYPE_O) << 8) | TYPE_N;
template <const StateType num_states>
void cmaple::SeqRegions::mergeUpperLower(std::unique_ptr<SeqRegions>& merged_regions,
                                 RealNumType upper_plength, 
                                 const SeqRegions& lower_regions, 
                                 RealNumType lower_plength, const
                                 AlignmentBase* aln,
                                 const ModelBase* model,
                                 const RealNumType threshold_prob) const
{
    // init variables
    PositionType pos = 0;
    const SeqRegions& seq1_regions = *this;
    const SeqRegions& seq2_regions = lower_regions;
    size_t iseq1 = 0;
    size_t iseq2 = 0;
    const PositionType seq_length = aln->ref_seq.size();
    
    // init merged_regions
    if (merged_regions)
        merged_regions->clear();
    else
        merged_regions = std::make_unique<SeqRegions>();
    
    // avoid realloc of vector data (minimize memory footprint)
    merged_regions->reserve(countSharedSegments(seq2_regions, seq_length)); // avoid realloc of vector data
    const size_t max_elements = merged_regions->capacity(); // remember capacity (may be more than we 'reserved')

    while (pos < seq_length)
    {
        PositionType end_pos;
        
        // get the next shared segment in the two sequences
        cmaple::SeqRegions::getNextSharedSegment(pos, seq1_regions, seq2_regions, iseq1, iseq2, end_pos);
        const auto* const seq1_region = &seq1_regions[iseq1];
        const auto* const seq2_region = &seq2_regions[iseq2];
        const DoubleState s1s2 = (DoubleState(seq1_region->type) << 8) | seq2_region->type;

        // seq1_entry = 'N'
        // seq1_entry = 'N' and seq2_entry = 'N'
        if (s1s2 == NN)
            merged_regions->emplace_back(TYPE_N, end_pos);
        // seq1_entry = 'N' and seq2_entry = O/R/ACGT
        // seq1_entry = 'N' and seq2_entry = O
        else if (s1s2 == NO)
        {
            merge_N_O<num_states>(lower_plength, *seq2_region, model, end_pos, *merged_regions);
        }
        // seq1_entry = 'N' and seq2_entry = R/ACGT
        else if (seq1_region->type == TYPE_N)
        {
            merge_N_RACGT(*seq2_region, lower_plength, end_pos, threshold_prob, *merged_regions);
        }
        // seq2_entry = 'N'
        // seq1_entry = 'O' and seq2_entry = N
        else if (s1s2 == ON)
        {
            merge_O_N<num_states>(*seq1_region, upper_plength, end_pos, model, *merged_regions);
        }
        // seq2_entry = 'N' and seq1_entry = R/ACGT
        else if (seq2_region->type == TYPE_N)
        {
            merge_RACGT_N(*seq1_region, upper_plength, end_pos, threshold_prob, *merged_regions);
        }
        // seq1_entry = seq2_entry = R/ACGT
        // todo: improve condition: a == b & a == RACGT
        else if (seq1_region->type == seq2_region->type && (seq1_region->type < num_states || seq1_region->type == TYPE_R))
        {
            // add a new region and try to merge consecutive R regions together
            addNonConsecutiveRRegion(*merged_regions, seq1_region->type, -1, -1, end_pos, threshold_prob);
        }
        // cases where the new genome list entry will likely be of type "O"
        else
        {
            RealNumType total_blength_1 = upper_plength;
            if (seq1_region->plength_observation2node >= 0)
            {
                total_blength_1 = seq1_region->plength_observation2node;
                if (upper_plength > 0)
                    total_blength_1 += upper_plength;
                
                if (seq1_region->type != TYPE_O && seq1_region->plength_observation2root >= 0)
                    total_blength_1 += seq1_region->plength_observation2root;
            }
            
            RealNumType total_blength_2 = lower_plength;
            if (seq2_region->plength_observation2node >= 0)
            {
                total_blength_2 = seq2_region->plength_observation2node;
                if (lower_plength > 0)
                    total_blength_2 += lower_plength;
            }
            
            // special cases when either total_blength_1 or total_blength_2 is zero
            if (merge_Zero_Distance(*seq1_region, *seq2_region, total_blength_1, total_blength_2, end_pos, threshold_prob, num_states, merged_regions))
            {
                if (!merged_regions) return;
            }
            // seq1_entry = O and seq2_entry = O/R/ACGT
            else if (seq1_region->type == TYPE_O)
            {
                merge_O_ORACGT<num_states>(*seq1_region, *seq2_region, total_blength_1, total_blength_2, end_pos, threshold_prob, model, aln, *merged_regions);
            }
            // seq1_entry = R/ACGT and seq2_entry = O/R/ACGT
            else
            {
                merge_RACGT_ORACGT<num_states>(*seq1_region, *seq2_region, total_blength_1, total_blength_2, upper_plength, end_pos, threshold_prob, model, aln, *merged_regions);
            }
            
            // NHANLT: LOGS FOR DEBUGGING
            /*if (Params::getInstance().debug && merged_regions->at(merged_regions->size()-1).type == TYPE_O)
            {
                SeqRegion::LHType lh =  *merged_regions->at(merged_regions->size()-1).likelihood;
                std::cout << "mergeUpLow " << pos << " " <<  std::setprecision(20) << lh[0] << " " << lh[1] << " " << lh[2] << " " << lh[3] << " " << std::endl;
            }*/
        }

        // update pos
        pos = end_pos + 1;
    }
    
    assert(merged_regions->capacity() == max_elements); // ensure we did the correct reserve, otherwise it was a pessimization
    
    // randomly choose some test-cases for testing
    /*if (Params::getInstance().output_testing && rand() < 2000000)
    {
        std::string output_file(Params::getInstance().output_testing);
        std::ofstream out = std::ofstream(output_file + ".txt", std::ios_base::app);
        out << "// --- New test --- //" << "\t@" << std::setprecision(50) << upper_plength << "\t@" << std::setprecision(50) << lower_plength << std::endl;
        this->writeConstructionCodes("input1", out, num_states);
        lower_regions.writeConstructionCodes("input2", out, num_states);
        merged_regions->writeConstructionCodes("output", out, num_states);
        out.close();
    }*/
}

void cmaple::merge_N_O_TwoLowers(const SeqRegion& seq2_region, const PositionType end_pos, const RealNumType plength2, SeqRegions& merged_regions)
{
    // add merged region into merged_regions
    SeqRegion new_region = SeqRegion::clone(seq2_region);
    new_region.position = end_pos;
    if (seq2_region.plength_observation2node >= 0)
    {
        if (plength2 > 0)
            new_region.plength_observation2node += plength2;
    }
    else
    {
        if (plength2 > 0)
            new_region.plength_observation2node = plength2;
    }
    merged_regions.push_back(std::move(new_region));
}

void cmaple::merge_N_RACGT_TwoLowers(const SeqRegion& seq2_region, const PositionType end_pos, const RealNumType plength2, const RealNumType threshold_prob, SeqRegions& merged_regions)
{
    RealNumType plength_observation2node = -1;
    
    if (seq2_region.plength_observation2node >= 0)
    {
        RealNumType new_plength = seq2_region.plength_observation2node;
        if (plength2 > 0)
            new_plength += plength2;
        
        plength_observation2node = new_plength;
    }
    else if (plength2 > 0)
            plength_observation2node = plength2;
    
    // add a new region and try to merge consecutive R regions together
    cmaple::SeqRegions::addNonConsecutiveRRegion(merged_regions, seq2_region.type, plength_observation2node, -1, end_pos, threshold_prob);
}

void cmaple::merge_identicalRACGT_TwoLowers(const SeqRegion& seq1_region, const PositionType end_pos, RealNumType total_blength_1, RealNumType total_blength_2, const PositionType pos, const RealNumType threshold_prob, const ModelBase* model, RealNumType &log_lh, SeqRegions& merged_regions, const bool return_log_lh)
{
    const RealNumType* const &cumulative_rate = model->cumulative_rate;
    
    // add a new region and try to merge consecutive R regions together
    cmaple::SeqRegions::addNonConsecutiveRRegion(merged_regions, seq1_region.type, -1, -1, end_pos, threshold_prob);
    
    if (return_log_lh)
    {
        // convert total_blength_1 and total_blength_2 to zero if they are -1
        if (total_blength_1 < 0) total_blength_1 = 0;
        if (total_blength_2 < 0) total_blength_2 = 0;
        
        if (seq1_region.type == TYPE_R)
            log_lh += (total_blength_1 + total_blength_2) * (cumulative_rate[end_pos + 1] - cumulative_rate[pos]);
        else
            log_lh += model->diagonal_mut_mat[seq1_region.type] * (total_blength_1 + total_blength_2);
    }
}

template <const StateType num_states>
bool cmaple::merge_O_O_TwoLowers(const SeqRegion& seq2_region, RealNumType total_blength_2, const PositionType end_pos, const AlignmentBase* aln, const ModelBase* model, const RealNumType threshold_prob, RealNumType &log_lh, SeqRegion::LHType& new_lh, std::unique_ptr<SeqRegions>& merged_regions, const bool return_log_lh)
{
    RealNumType sum_lh = updateMultLHwithMat<num_states>(model->mutation_mat, *seq2_region.likelihood, new_lh, total_blength_2);
    
    if (sum_lh == 0)
    {
        merged_regions = nullptr;
        return false;
    }
        
    // normalize the new partial likelihood
    normalize_arr(new_lh.data(), num_states, sum_lh);
    cmaple::SeqRegions::addSimplifiedO(end_pos, new_lh, aln, threshold_prob, *merged_regions);
    
    if (return_log_lh)
        log_lh += log(sum_lh);
    
    // no error
    return true;
}

template <const StateType num_states>
bool cmaple::merge_O_RACGT_TwoLowers(const SeqRegion& seq2_region, RealNumType total_blength_2, const PositionType end_pos, const AlignmentBase* aln, const ModelBase* model, const RealNumType threshold_prob, RealNumType &log_lh, SeqRegion::LHType& new_lh, RealNumType& sum_lh, std::unique_ptr<SeqRegions>& merged_regions, const bool return_log_lh)
{
    StateType seq2_state = seq2_region.type;
    if (seq2_state == TYPE_R)
        seq2_state = aln->ref_seq[end_pos];
    
    if (total_blength_2 > 0)
    {
        RealNumType* transposed_mut_mat_row = model->transposed_mut_mat + model->row_index[seq2_state];
        sum_lh = updateVecWithState<num_states>(new_lh.data(), seq2_state, transposed_mut_mat_row, total_blength_2);
        
        // normalize new partial lh
        // normalize the new partial likelihood
        normalize_arr(new_lh.data(), num_states, sum_lh);
        cmaple::SeqRegions::addSimplifiedO(end_pos, new_lh, aln, threshold_prob, *merged_regions);
        
        if (return_log_lh)
            log_lh += log(sum_lh);
    }
    else
    {
        if (new_lh[seq2_state] == 0)
        {
            merged_regions = nullptr;
            return false;
        }
        
        // add a new region and try to merge consecutive R regions together
        cmaple::SeqRegions::addNonConsecutiveRRegion(*merged_regions, seq2_region.type, -1, -1, end_pos, threshold_prob);
    
        if (return_log_lh)
            log_lh += log(new_lh[seq2_state]);
    }
    
    // no error
    return true;
}

template <const StateType num_states>
bool cmaple::merge_O_ORACGT_TwoLowers(const SeqRegion& seq1_region, const SeqRegion& seq2_region, RealNumType total_blength_1, RealNumType total_blength_2, const PositionType end_pos, const AlignmentBase* aln, const ModelBase* model, const RealNumType threshold_prob, RealNumType &log_lh, std::unique_ptr<SeqRegions>& merged_regions, const bool return_log_lh)
{
    auto new_lh = std::make_unique<SeqRegion::LHType>(); // = new RealNumType[num_states];
    RealNumType sum_lh = 0;
    
    if (total_blength_1 > 0)
    {
      sum_lh = updateLHwithMat<num_states>(model->mutation_mat, *(seq1_region.likelihood), *new_lh, total_blength_1);
    }
    // otherwise, clone the partial likelihood from seq1
    else
        *new_lh = *seq1_region.likelihood;

    // seq1_entry = O and seq2_entry = O
    if (seq2_region.type == TYPE_O)
    {
        return merge_O_O_TwoLowers<num_states>(seq2_region, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, *new_lh, merged_regions, return_log_lh);
    }
    // seq1_entry = O and seq2_entry = R/ACGT
    else
    {
        return merge_O_RACGT_TwoLowers<num_states>(seq2_region, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, *new_lh, sum_lh, merged_regions, return_log_lh);
    }
    
    // no error
    return true;
}

template <const StateType num_states>
bool cmaple::merge_RACGT_O_TwoLowers(const SeqRegion& seq2_region, RealNumType total_blength_2, const PositionType end_pos, const AlignmentBase* aln, const ModelBase* model, const RealNumType threshold_prob, SeqRegion::LHType& new_lh, RealNumType &log_lh, std::unique_ptr<SeqRegions>& merged_regions, const bool return_log_lh)
{
    RealNumType sum_lh = updateMultLHwithMat<num_states>(model->mutation_mat, *(seq2_region.likelihood), new_lh, total_blength_2);
                        
    if (sum_lh == 0)
    {
        merged_regions = nullptr;
        return false;
    }
        
    // normalize the new partial likelihood
    normalize_arr(new_lh.data(), num_states, sum_lh);
    cmaple::SeqRegions::addSimplifiedO(end_pos, new_lh, aln, threshold_prob, *merged_regions);
    
    if (return_log_lh)
        log_lh += log(sum_lh);
    
    // no error
    return true;
}

template <const StateType num_states>
bool cmaple::merge_RACGT_RACGT_TwoLowers(const SeqRegion& seq2_region, RealNumType total_blength_2, const PositionType end_pos, const AlignmentBase* aln, const ModelBase* model, const RealNumType threshold_prob, SeqRegion::LHType& new_lh, RealNumType& sum_lh, RealNumType &log_lh, std::unique_ptr<SeqRegions>& merged_regions, const bool return_log_lh)
{
    StateType seq2_state = seq2_region.type;
    if (seq2_state == TYPE_R)
        seq2_state = aln->ref_seq[end_pos];
    
    if (total_blength_2 > 0)
    {
        RealNumType* transposed_mut_mat_row = model->transposed_mut_mat + model->row_index[seq2_state];
        sum_lh += updateVecWithState<num_states>(new_lh.data(), seq2_state, transposed_mut_mat_row, total_blength_2);
        
        // normalize the new partial likelihood
        normalize_arr(new_lh.data(), num_states, sum_lh);
        cmaple::SeqRegions::addSimplifiedO(end_pos, new_lh, aln, threshold_prob, *merged_regions);
        
        if (return_log_lh)
            log_lh += log(sum_lh);
    }
    else
    {
        // add a new region and try to merge consecutive R regions together
        cmaple::SeqRegions::addNonConsecutiveRRegion(*merged_regions, seq2_region.type, -1, -1, end_pos, threshold_prob);
    
        if (return_log_lh)
            log_lh += log(new_lh[seq2_state]);
    }
    
    // no error
    return true;
}

template <const StateType num_states>
bool cmaple::merge_RACGT_ORACGT_TwoLowers(const SeqRegion& seq1_region, const SeqRegion& seq2_region, RealNumType total_blength_1, RealNumType total_blength_2, const PositionType end_pos, const AlignmentBase* aln, const ModelBase* model, const RealNumType threshold_prob, RealNumType &log_lh, std::unique_ptr<SeqRegions>& merged_regions, const bool return_log_lh)
{
    StateType seq1_state = seq1_region.type;
    if (seq1_state == TYPE_R)
        seq1_state = aln->ref_seq[end_pos];
    
    auto new_lh = std::make_unique<SeqRegion::LHType>(); // = new RealNumType[num_states];
    //auto& new_lh_value = *new_lh;
    RealNumType sum_lh = 0;
    
    if (total_blength_1 > 0)
    {
        RealNumType* transposed_mut_mat_row = model->transposed_mut_mat + model->row_index[seq1_state];
        setVecWithState<num_states>(new_lh->data(), seq1_state, transposed_mut_mat_row, total_blength_1);
    }
    else
        resetLhVecExceptState<num_states>(new_lh->data(), seq1_state, 1);

    // seq1_entry = R/ACGT and seq2_entry = O
    if (seq2_region.type == TYPE_O)
        return merge_RACGT_O_TwoLowers<num_states>(seq2_region, total_blength_2, end_pos, aln, model, threshold_prob, *new_lh, log_lh, merged_regions, return_log_lh);
    
    // otherwise, seq1_entry = R/ACGT and seq2_entry = R/ACGT
    return merge_RACGT_RACGT_TwoLowers<num_states>(seq2_region, total_blength_2, end_pos, aln, model, threshold_prob, *new_lh, sum_lh, log_lh, merged_regions, return_log_lh);
}

template <const StateType num_states>
bool cmaple::merge_notN_notN_TwoLowers(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType plength1, const RealNumType plength2, const PositionType end_pos, const PositionType pos, const AlignmentBase* aln, const ModelBase* model, const RealNumType threshold_prob, RealNumType &log_lh, std::unique_ptr<SeqRegions>& merged_regions, const bool return_log_lh)
{
    RealNumType total_blength_1 = plength1;
    if (seq1_region.plength_observation2node >= 0)
    {
        total_blength_1 = seq1_region.plength_observation2node;
        if (plength1 > 0)
            total_blength_1 += plength1;
    }
    
    RealNumType total_blength_2 = plength2;
    if (seq2_region.plength_observation2node >= 0)
    {
        total_blength_2 = seq2_region.plength_observation2node;
        if (plength2 > 0)
            total_blength_2 += plength2;
    }
    
    // seq1_entry and seq2_entry are identical seq1_entry = R/ACGT
    if (seq1_region.type == seq2_region.type && (seq1_region.type == TYPE_R || seq1_region.type < num_states))
    {
        merge_identicalRACGT_TwoLowers(seq1_region, end_pos, total_blength_1, total_blength_2, pos, threshold_prob, model, log_lh, *merged_regions, return_log_lh);
    }
    // #0 distance between different nucleotides: merge is not possible
    else if (total_blength_1 <= 0 && total_blength_2 <= 0 && (seq1_region.type == TYPE_R || seq1_region.type < num_states) && (seq2_region.type == TYPE_R || seq2_region.type < num_states))
    {
        merged_regions = nullptr;
        return false;
    }
    // seq1_entry = O
    else if (seq1_region.type == TYPE_O)
    {
        return merge_O_ORACGT_TwoLowers<num_states>(seq1_region, seq2_region, total_blength_1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, merged_regions, return_log_lh);
    }
    // seq1_entry = R/ACGT
    else
    {
        return merge_RACGT_ORACGT_TwoLowers<num_states>(seq1_region, seq2_region, total_blength_1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, merged_regions, return_log_lh);
    }
    
    // no error
    return true;
}

template <const StateType num_states>
RealNumType cmaple::SeqRegions::mergeTwoLowers(std::unique_ptr<SeqRegions>& merged_regions, const RealNumType plength1, const SeqRegions& regions2, const RealNumType plength2, const AlignmentBase* aln, const ModelBase* model, const RealNumType threshold_prob, const bool return_log_lh) const
{
    // init variables
    RealNumType log_lh = 0;
    PositionType pos = 0;
    const SeqRegions& seq1_regions = *this;
    const SeqRegions& seq2_regions = regions2;
    size_t iseq1 = 0;
    size_t iseq2 = 0;
    const PositionType seq_length = aln->ref_seq.size();
    
    // init merged_regions
    if (merged_regions)
        merged_regions->clear();
    else
        merged_regions = std::make_unique<SeqRegions>();

    // avoid realloc of vector data (minimize memory footprint)
    merged_regions->reserve(countSharedSegments(seq2_regions, seq_length)); // avoid realloc of vector data
    const size_t max_elements = merged_regions->capacity(); // remember capacity (may be more than we 'reserved')

    while (pos < seq_length)
    {
        PositionType end_pos;

        // get the next shared segment in the two sequences
        cmaple::SeqRegions::getNextSharedSegment(pos, seq1_regions, seq2_regions, iseq1, iseq2, end_pos);
        const auto* const seq1_region = &seq1_regions[iseq1];
        const auto* const seq2_region = &seq2_regions[iseq2];
        const DoubleState s1s2 = (DoubleState(seq1_region->type) << 8) | seq2_region->type;

        // seq1_entry = 'N'
        // seq1_entry = 'N' and seq2_entry = 'N'
        if (s1s2 == NN)
            merged_regions->emplace_back(TYPE_N, end_pos);
        // seq1_entry = 'N' and seq2_entry = O/R/ACGT
        // seq1_entry = 'N' and seq2_entry = O
        else if (s1s2 == NO)
        {
            merge_N_O_TwoLowers(*seq2_region, end_pos, plength2, *merged_regions);
        }
        // seq1_entry = 'N' and seq2_entry = R/ACGT
        else if (seq1_region->type == TYPE_N)
        {
            merge_N_RACGT_TwoLowers(*seq2_region, end_pos, plength2, threshold_prob, *merged_regions);
        }
        // seq2_entry = 'N'
        // seq1_entry = 'O' and seq2_entry = N
        else if (s1s2 == ON)
        {
            // NOTE: merge_N_O_TwoLowers can be used to merge O_N
            merge_N_O_TwoLowers(*seq1_region, end_pos, plength1, *merged_regions);
        }
        // seq1_entry = R/ACGT and seq2_entry = 'N'
        else if (seq2_region->type == TYPE_N)
        {
            // NOTE: merge_N_RACGT can be used to merge RACGT_N
            merge_N_RACGT_TwoLowers(*seq1_region, end_pos, plength1, threshold_prob, *merged_regions);
        }
        // neither seq1_entry nor seq2_entry = N
        else
        {
            if (!merge_notN_notN_TwoLowers<num_states>(*seq1_region, *seq2_region, plength1, plength2, end_pos, pos, aln, model, threshold_prob, log_lh, merged_regions, return_log_lh)) return MIN_NEGATIVE;
        }
        
        // NHANLT: LOGS FOR DEBUGGING
        /*if (Params::getInstance().debug && merged_regions->at(merged_regions->size()-1).type == TYPE_O)
        {
            SeqRegion::LHType lh =  *merged_regions->at(merged_regions->size()-1).likelihood;
            std::cout << "merge2Low " << pos << " " << std::setprecision(20) << lh[0] << " " << lh[1] << " " << lh[2] << " " << lh[3] << " " << std::endl;
        }*/
        
        // update pos
        pos = end_pos + 1;
    }
    
    assert(merged_regions->capacity() == max_elements); // ensure we did the correct reserve, otherwise it was a pessimization

    return log_lh;
}

void calSiteLhs_identicalRACGT(std::vector<RealNumType>& site_lh_contributions, const SeqRegion& seq1_region, const PositionType end_pos, RealNumType total_blength_1, RealNumType total_blength_2, const PositionType pos, const RealNumType threshold_prob, const ModelBase* model, RealNumType &log_lh, SeqRegions& merged_regions)
{
    const RealNumType* const &cumulative_rate = model->cumulative_rate;
    
    // add a new region and try to merge consecutive R regions together
    cmaple::SeqRegions::addNonConsecutiveRRegion(merged_regions, seq1_region.type, -1, -1, end_pos, threshold_prob);
    
    // compute the (site) lh contributions
    // convert total_blength_1 and total_blength_2 to zero if they are -1
    if (total_blength_1 < 0) total_blength_1 = 0;
    if (total_blength_2 < 0) total_blength_2 = 0;
        
    RealNumType total_blength = total_blength_1 + total_blength_2;
    if (seq1_region.type == TYPE_R)
    {
        log_lh += total_blength * (cumulative_rate[end_pos + 1] - cumulative_rate[pos]);
        
        // compute site lh contributions
        for (PositionType i = pos; i < end_pos + 1; ++i)
            site_lh_contributions[i] += total_blength * (cumulative_rate[i + 1] - cumulative_rate[i]);
    }
    else
    {
        log_lh += model->diagonal_mut_mat[seq1_region.type] * total_blength;
        
        // compute site lh contributions
        site_lh_contributions[pos] += model->diagonal_mut_mat[seq1_region.type] * total_blength;
    }
}

template <const StateType num_states>
bool calSiteLhs_O_O(std::vector<RealNumType>& site_lh_contributions, const SeqRegion& seq2_region, RealNumType total_blength_2, const PositionType end_pos, const AlignmentBase* aln, const ModelBase* model, const RealNumType threshold_prob, RealNumType &log_lh, SeqRegion::LHType& new_lh, std::unique_ptr<SeqRegions>& merged_regions)
{
    RealNumType sum_lh = updateMultLHwithMat<num_states>(model->mutation_mat, *seq2_region.likelihood, new_lh, total_blength_2);
    
    if (sum_lh == 0)
    {
        merged_regions = nullptr;
        return false;
    }
        
    // normalize the new partial likelihood
    normalize_arr(new_lh.data(), num_states, sum_lh);
    cmaple::SeqRegions::addSimplifiedO(end_pos, new_lh, aln, threshold_prob, *merged_regions);
    
    // compute (site) lh contributions
    RealNumType lh_contribution = log(sum_lh);
    log_lh += lh_contribution;
    site_lh_contributions[end_pos] += lh_contribution;
    
    // no error
    return true;
}

template <const StateType num_states>
bool calSiteLhs_O_RACGT(std::vector<RealNumType>& site_lh_contributions, const SeqRegion& seq2_region, RealNumType total_blength_2, const PositionType end_pos, const AlignmentBase* aln, const ModelBase* model, const RealNumType threshold_prob, RealNumType &log_lh, SeqRegion::LHType& new_lh, RealNumType& sum_lh, std::unique_ptr<SeqRegions>& merged_regions)
{
    StateType seq2_state = seq2_region.type;
    if (seq2_state == TYPE_R)
        seq2_state = aln->ref_seq[end_pos];
    
    if (total_blength_2 > 0)
    {
        RealNumType* transposed_mut_mat_row = model->transposed_mut_mat + model->row_index[seq2_state];
        sum_lh = updateVecWithState<num_states>(new_lh.data(), seq2_state, transposed_mut_mat_row, total_blength_2);
        
        // normalize new partial lh
        // normalize the new partial likelihood
        normalize_arr(new_lh.data(), num_states, sum_lh);
        cmaple::SeqRegions::addSimplifiedO(end_pos, new_lh, aln, threshold_prob, *merged_regions);
        
        // compute (site) lh contributions
        RealNumType lh_contribution = log(sum_lh);
        log_lh += lh_contribution;
        site_lh_contributions[end_pos] += lh_contribution;
    }
    else
    {
        if (new_lh[seq2_state] == 0)
        {
            merged_regions = nullptr;
            return false;
        }
        
        // add a new region and try to merge consecutive R regions together
        cmaple::SeqRegions::addNonConsecutiveRRegion(*merged_regions, seq2_region.type, -1, -1, end_pos, threshold_prob);
        
        // compute (site) lh contributions
        RealNumType lh_contribution = log(new_lh[seq2_state]);
        log_lh += lh_contribution;
        site_lh_contributions[end_pos] += lh_contribution;
    }
    
    // no error
    return true;
}

template <const StateType num_states>
bool calSiteLhs_O_ORACGT(std::vector<RealNumType>& site_lh_contributions, const SeqRegion& seq1_region, const SeqRegion& seq2_region, RealNumType total_blength_1, RealNumType total_blength_2, const PositionType end_pos, const AlignmentBase* aln, const ModelBase* model, const RealNumType threshold_prob, RealNumType &log_lh, std::unique_ptr<SeqRegions>& merged_regions)
{
    auto new_lh = std::make_unique<SeqRegion::LHType>(); // = new RealNumType[num_states];
    RealNumType sum_lh = 0;
    
    if (total_blength_1 > 0)
    {
      sum_lh = updateLHwithMat<num_states>(model->mutation_mat, *(seq1_region.likelihood), *new_lh, total_blength_1);
    }
    // otherwise, clone the partial likelihood from seq1
    else
        *new_lh = *seq1_region.likelihood;

    // seq1_entry = O and seq2_entry = O
    if (seq2_region.type == TYPE_O)
    {
        return calSiteLhs_O_O<num_states>(site_lh_contributions, seq2_region, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, *new_lh, merged_regions);
    }
    // seq1_entry = O and seq2_entry = R/ACGT
    else
    {
        return calSiteLhs_O_RACGT<num_states>(site_lh_contributions, seq2_region, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, *new_lh, sum_lh, merged_regions);
    }
    
    // no error
    return true;
}

template <const StateType num_states>
bool calSiteLhs_RACGT_O(std::vector<RealNumType>& site_lh_contributions, const SeqRegion& seq2_region, RealNumType total_blength_2, const PositionType end_pos, const AlignmentBase* aln, const ModelBase* model, const RealNumType threshold_prob, SeqRegion::LHType& new_lh, RealNumType &log_lh, std::unique_ptr<SeqRegions>& merged_regions)
{
    RealNumType sum_lh = updateMultLHwithMat<num_states>(model->mutation_mat, *(seq2_region.likelihood), new_lh, total_blength_2);
                        
    if (sum_lh == 0)
    {
        merged_regions = nullptr;
        return false;
    }
        
    // normalize the new partial likelihood
    normalize_arr(new_lh.data(), num_states, sum_lh);
    cmaple::SeqRegions::addSimplifiedO(end_pos, new_lh, aln, threshold_prob, *merged_regions);
    
    // compute (site) lh contributions
    RealNumType lh_contribution = log(sum_lh);
    log_lh += lh_contribution;
    site_lh_contributions[end_pos] += lh_contribution;
    
    // no error
    return true;
}

template <const StateType num_states>
bool calSiteLhs_RACGT_RACGT(std::vector<RealNumType>& site_lh_contributions, const SeqRegion& seq2_region, RealNumType total_blength_2, const PositionType end_pos, const AlignmentBase* aln, const ModelBase* model, const RealNumType threshold_prob, SeqRegion::LHType& new_lh, RealNumType& sum_lh, RealNumType &log_lh, std::unique_ptr<SeqRegions>& merged_regions)
{
    StateType seq2_state = seq2_region.type;
    if (seq2_state == TYPE_R)
        seq2_state = aln->ref_seq[end_pos];
    
    if (total_blength_2 > 0)
    {
        RealNumType* transposed_mut_mat_row = model->transposed_mut_mat + model->row_index[seq2_state];
        sum_lh += updateVecWithState<num_states>(new_lh.data(), seq2_state, transposed_mut_mat_row, total_blength_2);
        
        // normalize the new partial likelihood
        normalize_arr(new_lh.data(), num_states, sum_lh);
        cmaple::SeqRegions::addSimplifiedO(end_pos, new_lh, aln, threshold_prob, *merged_regions);
        
        // compute (site) lh contributions
        RealNumType lh_contribution = log(sum_lh);
        log_lh += lh_contribution;
        site_lh_contributions[end_pos] += lh_contribution;
    }
    else
    {
        // add a new region and try to merge consecutive R regions together
        cmaple::SeqRegions::addNonConsecutiveRRegion(*merged_regions, seq2_region.type, -1, -1, end_pos, threshold_prob);
    
        // compute (site) lh contributions
        RealNumType lh_contribution = log(new_lh[seq2_state]);
        log_lh += lh_contribution;
        site_lh_contributions[end_pos] += lh_contribution;
    }
    
    // no error
    return true;
}

template <const StateType num_states>
bool calSiteLhs_RACGT_ORACGT(std::vector<RealNumType>& site_lh_contributions, const SeqRegion& seq1_region, const SeqRegion& seq2_region, RealNumType total_blength_1, RealNumType total_blength_2, const PositionType end_pos, const AlignmentBase* aln, const ModelBase* model, const RealNumType threshold_prob, RealNumType &log_lh, std::unique_ptr<SeqRegions>& merged_regions)
{
    StateType seq1_state = seq1_region.type;
    if (seq1_state == TYPE_R)
        seq1_state = aln->ref_seq[end_pos];
    
    auto new_lh = std::make_unique<SeqRegion::LHType>(); // = new RealNumType[num_states];
    //auto& new_lh_value = *new_lh;
    RealNumType sum_lh = 0;
    
    if (total_blength_1 > 0)
    {
        RealNumType* transposed_mut_mat_row = model->transposed_mut_mat + model->row_index[seq1_state];
        setVecWithState<num_states>(new_lh->data(), seq1_state, transposed_mut_mat_row, total_blength_1);
    }
    else
        resetLhVecExceptState<num_states>(new_lh->data(), seq1_state, 1);

    // seq1_entry = R/ACGT and seq2_entry = O
    if (seq2_region.type == TYPE_O)
        return calSiteLhs_RACGT_O<num_states>(site_lh_contributions, seq2_region, total_blength_2, end_pos, aln, model, threshold_prob, *new_lh, log_lh, merged_regions);
    
    // otherwise, seq1_entry = R/ACGT and seq2_entry = R/ACGT
    return calSiteLhs_RACGT_RACGT<num_states>(site_lh_contributions, seq2_region, total_blength_2, end_pos, aln, model, threshold_prob, *new_lh, sum_lh, log_lh, merged_regions);
}

template <const StateType num_states>
bool calSiteLhs_notN_notN(std::vector<RealNumType>& site_lh_contributions, const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType plength1, const RealNumType plength2, const PositionType end_pos, const PositionType pos, const AlignmentBase* aln, const ModelBase* model, const RealNumType threshold_prob, RealNumType &log_lh, std::unique_ptr<SeqRegions>& merged_regions)
{
    RealNumType total_blength_1 = plength1;
    if (seq1_region.plength_observation2node >= 0)
    {
        total_blength_1 = seq1_region.plength_observation2node;
        if (plength1 > 0)
            total_blength_1 += plength1;
    }
    
    RealNumType total_blength_2 = plength2;
    if (seq2_region.plength_observation2node >= 0)
    {
        total_blength_2 = seq2_region.plength_observation2node;
        if (plength2 > 0)
            total_blength_2 += plength2;
    }
    
    // seq1_entry and seq2_entry are identical seq1_entry = R/ACGT
    if (seq1_region.type == seq2_region.type && (seq1_region.type == TYPE_R || seq1_region.type < num_states))
    {
        calSiteLhs_identicalRACGT(site_lh_contributions, seq1_region, end_pos, total_blength_1, total_blength_2, pos, threshold_prob, model, log_lh, *merged_regions);
    }
    // #0 distance between different nucleotides: merge is not possible
    else if (total_blength_1 == 0 && total_blength_2 == 0 && (seq1_region.type == TYPE_R || seq1_region.type < num_states) && (seq2_region.type == TYPE_R || seq2_region.type < num_states))
    {
        merged_regions = nullptr;
        return false;
    }
    // seq1_entry = O
    else if (seq1_region.type == TYPE_O)
    {
        return calSiteLhs_O_ORACGT<num_states>(site_lh_contributions, seq1_region, seq2_region, total_blength_1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, merged_regions);
    }
    // seq1_entry = R/ACGT
    else
    {
        return calSiteLhs_RACGT_ORACGT<num_states>(site_lh_contributions, seq1_region, seq2_region, total_blength_1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, merged_regions);
    }
    
    // no error
    return true;
}

template <const StateType num_states>
RealNumType cmaple::SeqRegions::calculateSiteLhContributions(std::vector<RealNumType>& site_lh_contributions, std::unique_ptr<SeqRegions>& merged_regions, const RealNumType plength1, const SeqRegions& regions2, const RealNumType plength2, const AlignmentBase* aln, const ModelBase* model, const RealNumType threshold_prob) const
{
    // init variables
    RealNumType log_lh = 0;
    PositionType pos = 0;
    const SeqRegions& seq1_regions = *this;
    const SeqRegions& seq2_regions = regions2;
    size_t iseq1 = 0;
    size_t iseq2 = 0;
    const PositionType seq_length = aln->ref_seq.size();
    ASSERT(site_lh_contributions.size() == seq_length);
    
    // init merged_regions
    if (merged_regions)
        merged_regions->clear();
    else
        merged_regions = std::make_unique<SeqRegions>();

    // avoid realloc of vector data (minimize memory footprint)
    merged_regions->reserve(countSharedSegments(seq2_regions, seq_length)); // avoid realloc of vector data
    const size_t max_elements = merged_regions->capacity(); // remember capacity (may be more than we 'reserved')

    while (pos < seq_length)
    {
        PositionType end_pos;

        // get the next shared segment in the two sequences
        cmaple::SeqRegions::getNextSharedSegment(pos, seq1_regions, seq2_regions, iseq1, iseq2, end_pos);
        const auto* const seq1_region = &seq1_regions[iseq1];
        const auto* const seq2_region = &seq2_regions[iseq2];
        const DoubleState s1s2 = (DoubleState(seq1_region->type) << 8) | seq2_region->type;

        // seq1_entry = 'N'
        // seq1_entry = 'N' and seq2_entry = 'N'
        if (s1s2 == NN)
            merged_regions->emplace_back(TYPE_N, end_pos);
        // seq1_entry = 'N' and seq2_entry = O/R/ACGT
        // seq1_entry = 'N' and seq2_entry = O
        else if (s1s2 == NO)
        {
            merge_N_O_TwoLowers(*seq2_region, end_pos, plength2, *merged_regions);
        }
        // seq1_entry = 'N' and seq2_entry = R/ACGT
        else if (seq1_region->type == TYPE_N)
        {
            merge_N_RACGT_TwoLowers(*seq2_region, end_pos, plength2, threshold_prob, *merged_regions);
        }
        // seq2_entry = 'N'
        // seq1_entry = 'O' and seq2_entry = N
        else if (s1s2 == ON)
        {
            // NOTE: merge_N_O_TwoLowers can be used to merge O_N
            merge_N_O_TwoLowers(*seq1_region, end_pos, plength1, *merged_regions);
        }
        // seq1_entry = R/ACGT and seq2_entry = 'N'
        else if (seq2_region->type == TYPE_N)
        {
            // NOTE: merge_N_RACGT can be used to merge RACGT_N
            merge_N_RACGT_TwoLowers(*seq1_region, end_pos, plength1, threshold_prob, *merged_regions);
        }
        // neither seq1_entry nor seq2_entry = N
        else
        {
            if (!calSiteLhs_notN_notN<num_states>(site_lh_contributions, *seq1_region, *seq2_region, plength1, plength2, end_pos, pos, aln, model, threshold_prob, log_lh, merged_regions))
            {
                if (cmaple::verbose_mode >= cmaple::VB_DEBUG)
                    outWarning("calculateSiteLhContributions() returns MIN_NEGATIVE!");
                return MIN_NEGATIVE;
            }
        }
        
        // NHANLT: LOGS FOR DEBUGGING
        /*if (Params::getInstance().debug && merged_regions->at(merged_regions->size()-1).type == TYPE_O)
        {
            SeqRegion::LHType lh =  *merged_regions->at(merged_regions->size()-1).likelihood;
            std::cout << "merge2Low " << pos << " " << std::setprecision(20) << lh[0] << " " << lh[1] << " " << lh[2] << " " << lh[3] << " " << std::endl;
        }*/
        
        // update pos
        pos = end_pos + 1;
    }
    
    assert(merged_regions->capacity() == max_elements); // ensure we did the correct reserve, otherwise it was a pessimization

    return log_lh;
}

template <const StateType num_states>
RealNumType cmaple::SeqRegions::computeAbsoluteLhAtRoot(const ModelBase* model)
{
    // dummy variables
    RealNumType log_lh = 0;
    RealNumType log_factor = 1;
    PositionType start_pos = 0;
    const SeqRegions& regions = *this;
    const std::vector< std::vector<PositionType> > &cumulative_base = model->cumulative_base;
    
    // browse regions one by one to compute the likelihood of each region
    for (const SeqRegion& region : regions)
    {
        // type R
        if (region.type == TYPE_R)
        {
            for (StateType i = 0; i < num_states; ++i)
                log_lh += model->root_log_freqs[i] * (cumulative_base[region.position + 1][i] - cumulative_base[start_pos][i]);
        }
        // type ACGT
        else if (region.type < num_states)
            log_lh += model->root_log_freqs[region.type];
        // type O
        else if (region.type == TYPE_O)
        {
            RealNumType tot = 0;
            tot += dotProduct<num_states>(&(*region.likelihood)[0], model->root_freqs);
            log_factor *= tot;
        }
        
        // maintain start_pos
        start_pos = region.position + 1;
        
        // NhanLT: avoid underflow on log_factor
        if (log_factor <= MIN_CARRY_OVER)
        {
            //if (log_factor < MIN_POSITIVE)
            //    return MIN_NEGATIVE;
            
            //lh_cost += log(total_factor);
            //total_factor = 1.0;
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
RealNumType cmaple::SeqRegions::computeSiteLhAtRoot(std::vector<RealNumType>& site_lh_contributions, const ModelBase* model)
{
    // dummy variables
    RealNumType log_lh = 0;
    RealNumType log_factor = 1;
    PositionType start_pos = 0;
    const SeqRegions& regions = *this;
    const std::vector< std::vector<PositionType> > &cumulative_base = model->cumulative_base;
    
    // browse regions one by one to compute the likelihood of each region
    for (const SeqRegion& region : regions)
    {
        // type R
        if (region.type == TYPE_R)
        {
            for (StateType i = 0; i < num_states; ++i)
                log_lh += model->root_log_freqs[i] * (cumulative_base[region.position + 1][i] - cumulative_base[start_pos][i]);
            
            // calculate site lhs
            for (PositionType pos = start_pos; pos < region.position + 1; ++pos)
            {
                for (StateType i = 0; i < num_states; ++i)
                    site_lh_contributions[pos] += model->root_log_freqs[i] * (cumulative_base[pos + 1][i] - cumulative_base[pos][i]);
            }
        }
        // type ACGT
        else if (region.type < num_states)
        {
            RealNumType lh_contribution = model->root_log_freqs[region.type];
            log_lh += lh_contribution;
            
            // calculate site lhs
            site_lh_contributions[start_pos] += lh_contribution;
        }
        // type O
        else if (region.type == TYPE_O)
        {
            RealNumType tot = 0;
            tot += dotProduct<num_states>(&(*region.likelihood)[0], model->root_freqs);
            log_factor *= tot;
            
            // calculate site lhs
            site_lh_contributions[start_pos] += log(tot);
        }
        
        // maintain start_pos
        start_pos = region.position + 1;
        
        // NhanLT: avoid underflow on log_factor
        if (log_factor <= MIN_CARRY_OVER)
        {
            //if (log_factor < MIN_POSITIVE)
            //    return MIN_NEGATIVE;
            
            //lh_cost += log(total_factor);
            //total_factor = 1.0;
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
void cmaple::SeqRegions::computeTotalLhAtRoot(std::unique_ptr<SeqRegions>& total_lh, const ModelBase* model, RealNumType blength) const
{
    if (total_lh)
        total_lh->clear();
    else
        total_lh = std::make_unique<SeqRegions>();
    
    total_lh->reserve(size()); // avoid realloc of vector data
    for (const SeqRegion& elem : (*this))
    {
        const auto* const region = &elem;
        // type N
        if (region->type == TYPE_N)
        {
            total_lh->emplace_back(region->type, region->position, region->plength_observation2node, region->plength_observation2root);
        }
        else
        {
            // type O
            if (region->type == TYPE_O)
            {
                // compute total blength
                RealNumType total_blength = blength;
                if (region->plength_observation2node >= 0)
                {
                    total_blength = region->plength_observation2node;
                    if (blength > 0)
                        total_blength += blength;
                }
                
                // init new likelihood
                auto new_lh = std::make_unique<SeqRegion::LHType>(); // = new RealNumType[num_states];
                auto& new_lh_value = *new_lh;
                RealNumType sum_lh = updateLHwithModel<num_states>(model, *region->likelihood, (*new_lh), total_blength);
                // normalize the new partial likelihood
                normalize_arr(new_lh->data(), num_states, sum_lh);
                
                // add new region to the total_lh_regions
                total_lh->emplace_back(region->type, region->position, region->plength_observation2node, region->plength_observation2root, std::move(new_lh));
            }
            // other types: R or A/C/G/T
            else
            {
                // add new region to the total_lh_regions
                SeqRegion& new_region = total_lh->emplace_back(region->type, region->position, region->plength_observation2node, region->plength_observation2root);
                
                if (new_region.plength_observation2node >= 0)
                {
                    if (blength > 0)
                        new_region.plength_observation2node += blength;

                    new_region.plength_observation2root = 0;
                }
                else if (blength > 0)
                {
                    new_region.plength_observation2node = blength;
                    new_region.plength_observation2root = 0;
                }
            }
        }
    }
}

bool cmaple::SeqRegions::operator==(const SeqRegions& seqregions_1) const
{
    if (size() != seqregions_1.size())
        return false;
    
    for (PositionType i = 0; i < size(); ++i)
        if (!(at(i) == seqregions_1[i]))
            return false;
    
    return true;
}

void cmaple::SeqRegions::writeConstructionCodes(const std::string regions_name, std::ofstream& out, const StateType num_states) const
{
    const SeqRegions& regions = *this;
    
    // browse regions one by one to export the construction codes
    for (const SeqRegion& region : regions)
        region.writeConstructionCodes(regions_name, out, num_states);
}
