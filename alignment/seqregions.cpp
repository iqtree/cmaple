//
//  regions.cpp
//  alignment
//
//  Created by NhanLT on 31/3/2022.
//

#include "seqregions.h"
#include <cassert>


RealNumType updateLHwithModel(int num_states, const Model& model,
  const SeqRegion::LHType& prior,
  SeqRegion::LHType& posterior, const RealNumType total_blength)
{
    RealNumType sum_lh = 0;
    RealNumType* mutation_mat_row = model.mutation_mat;
    for (StateType i = 0; i < num_states; ++i, mutation_mat_row += num_states)
    {
        RealNumType tot = 0;
        if (total_blength > 0) // TODO: avoid
        {
            for (StateType j = 0; j < num_states; ++j)
                tot += mutation_mat_row[j] * prior[j];

            tot *= total_blength;
        }

        tot += prior[i];
        posterior[i] = tot * model.root_freqs[i];
        sum_lh += posterior[i];
    }
    return sum_lh;
}


RealNumType updateLHwithMat(int num_states, const RealNumType* mat_row,
  const SeqRegion::LHType& prior,
  SeqRegion::LHType& posterior, const RealNumType total_blength)
{
    RealNumType sum_lh = 0;
    for (StateType i = 0; i < num_states; ++i, mat_row += num_states)
    {
      RealNumType tot = 0;
      for (StateType j = 0; j < num_states; ++j)
          tot += mat_row[j] * prior[j];

      tot *= total_blength;
      tot += prior[i];
      posterior[i] = tot;
      sum_lh += tot;
    }
    return sum_lh;
}


RealNumType updateMultLHwithMat(int num_states, const RealNumType* mat_row,
  const SeqRegion::LHType& prior,
  SeqRegion::LHType& posterior, const RealNumType total_blength)
{
    RealNumType sum_lh = 0;
    for (StateType i = 0; i < num_states; ++i, mat_row += num_states)
    {
        RealNumType tot = 0;
        if (total_blength > 0) // TODO: avoid
        {
            for (StateType j = 0; j < num_states; ++j)
                tot += mat_row[j] * prior[j];

            tot *= total_blength;
        }
        tot += prior[i];
        posterior[i] *= tot;
        sum_lh += posterior[i];
    }
    return sum_lh;
}

SeqRegions::SeqRegions(SeqRegions* n_regions)
{
  if (!n_regions)
  {
    outError("oops");
  }
    // clone regions one by one
    reserve(n_regions->size());
    for (const auto& region : *n_regions)
      emplace_back(SeqRegion::clone(region));
} 

SeqRegions::~SeqRegions()
{
}

void SeqRegions::deleteRegions()
{
  clear();
}

int SeqRegions::compareWithSample(const SeqRegions& sequence2, PositionType seq_length, StateType num_states) const
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
    PositionType end_pos;

    while (pos < seq_length && (!seq1_more_info || !seq2_more_info))
    {
        // get the next shared segment in the two sequences
        SeqRegions::getNextSharedSegment(pos, seq1_regions, seq2_regions, iseq1, iseq2, end_pos);
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

bool SeqRegions::areDiffFrom(const SeqRegions& regions2, PositionType seq_length, StateType num_states, const Params* params) const
{
    if (regions2.empty())
        return true;
    
    // init variables
    PositionType pos = 0;
    const SeqRegions& seq1_regions = *this;
    const SeqRegions& seq2_regions = regions2;
    size_t iseq1 = 0;
    size_t iseq2 = 0;
    PositionType end_pos;
    
    // compare each pair of regions
    while (pos < seq_length)
    {
        // get the next shared segment in the two sequences
        SeqRegions::getNextSharedSegment(pos, seq1_regions, seq2_regions, iseq1, iseq2, end_pos);
        const auto* const seq1_region = &seq1_regions[iseq1];
        const auto* const seq2_region = &seq2_regions[iseq2];

        // if any pair of regions is different in type -> the two sequences are different
        if (seq1_region->type != seq2_region->type)
            return true;
        
        // type R/ACGT
        if (seq1_region->type < num_states || seq1_region->type == TYPE_R)
        {
            // compare plength_from_root and plength_observation
            if (fabs(seq1_region->plength_observation2root - seq2_region->plength_observation2root) > params->threshold_prob
                ||fabs(seq1_region->plength_observation2node - seq2_region->plength_observation2node) > params->threshold_prob)
                return true;
        }
        
        // type O
        if (seq1_region->type == TYPE_O)
        {
            // compare plength_observation
            if (fabs(seq1_region->plength_observation2node - seq2_region->plength_observation2node) > params->threshold_prob)
                return true;
            
            // compare likelihood of each state
            for (StateType i = 0; i < num_states; ++i)
            {
                RealNumType diff = fabs(seq1_region->getLH(i) - seq2_region->getLH(i));
                
                if (diff > 0)
                {
                    if ((seq1_region->getLH(i) == 0) || (seq2_region->getLH(i) == 0))
                        return true;
                    
                    if (diff > params->thresh_diff_update
                        || (diff > params->threshold_prob
                            && ((diff > params->thresh_diff_fold_update * seq1_region->getLH(i))
                                || (diff > params->thresh_diff_fold_update * seq2_region->getLH(i)))))
                        return true;
                }
            }
        }
        
        // update pos
        pos = end_pos + 1;
    }
    
    return false;
}

size_t SeqRegions::countSharedSegments(const SeqRegions& seq2_regions, const size_t seq_length) const
{
    const SeqRegions& seq1_regions = *this;
    size_t count{};
    PositionType pos{}, end_pos{};
    size_t iseq1 = 0;
    size_t iseq2 = 0;
    while (pos < seq_length)
    {
        // get the next shared segment in the two sequences
        SeqRegions::getNextSharedSegment(pos, seq1_regions, seq2_regions, iseq1, iseq2, end_pos);
        ++count;
        pos = end_pos + 1;
    }
    return ++count;
}

void SeqRegions::mergeUpperLower(SeqRegions* &merged_regions, 
                                 RealNumType upper_plength, 
                                 const SeqRegions& lower_regions, 
                                 RealNumType lower_plength, const
                                 Alignment& aln, 
                                 const Model& model,
                                 RealNumType threshold_prob) const
{
    // init variables
    PositionType pos = 0;
    const SeqRegions& seq1_regions = *this;
    const SeqRegions& seq2_regions = lower_regions;
    size_t iseq1 = 0;
    size_t iseq2 = 0;
    PositionType end_pos;
    StateType num_states = aln.num_states;
    PositionType seq_length = aln.ref_seq.size();
    
    // init merged_regions
    if (merged_regions)
        merged_regions->deleteRegions();
    else
        merged_regions = new SeqRegions();
    
    // avoid realloc of vector data (minimize memory footprint)
    merged_regions->reserve(countSharedSegments(seq2_regions, seq_length)); // avoid realloc of vector data
    const size_t max_elements = merged_regions->capacity(); // remember capacity (may be more than we 'reserved')

    while (pos < seq_length)
    {
        // get the next shared segment in the two sequences
        SeqRegions::getNextSharedSegment(pos, seq1_regions, seq2_regions, iseq1, iseq2, end_pos);
        const auto* const seq1_region = &seq1_regions[iseq1];
        const auto* const seq2_region = &seq2_regions[iseq2];

        // seq1_entry = 'N'
        if (seq1_region->type == TYPE_N)
        {
            // seq1_entry = 'N' and seq2_entry = 'N'
            if (seq2_region->type == TYPE_N)
                merged_regions->emplace_back(TYPE_N, end_pos);
            // seq1_entry = 'N' and seq2_entry = O/R/ACGT
            else
            {
                // seq1_entry = 'N' and seq2_entry = O
                if (seq2_region->type == TYPE_O)
                {
                    RealNumType total_blength = lower_plength;
                    if (seq2_region->plength_observation2node >= 0)
                        total_blength = seq2_region->plength_observation2node + (lower_plength > 0 ? lower_plength: 0);
                        
                    auto new_lh = std::make_unique<SeqRegion::LHType>(); // = new RealNumType[num_states];
                    RealNumType sum_lh = updateLHwithModel(num_states, model, *seq2_region->likelihood, (*new_lh), total_blength);
                    // normalize the new partial likelihood
                    normalize_arr(new_lh->data(), num_states, sum_lh);

                    // add merged region into merged_regions
                    merged_regions->emplace_back(TYPE_O, end_pos, 0, 0, std::move(new_lh));
                }
                // seq1_entry = 'N' and seq2_entry = R/ACGT
                else
                {
                    RealNumType plength_observation2node = -1;
                    RealNumType plength_observation2root = 0;
                    
                    if (seq2_region->plength_observation2node >= 0)
                    {
                        RealNumType new_plength = seq2_region->plength_observation2node ;
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
                    addNonConsecutiveRRegion(merged_regions, seq2_region->type, plength_observation2node, plength_observation2root, end_pos, threshold_prob);
                }
            }
        }
        // seq2_entry = 'N'
        else if (seq2_region->type == TYPE_N)
        {
            // seq1_entry = 'O' and seq2_entry = N
            if (seq1_region->type == TYPE_O)
            {
                RealNumType total_blength = -1;
                
                if (seq1_region->plength_observation2node >= 0)
                {
                    total_blength = seq1_region->plength_observation2node;
                    if (upper_plength > 0)
                        total_blength += upper_plength;
                }
                else if (upper_plength > 0)
                    total_blength = upper_plength;
                
                if (total_blength > 0)
                {
                    auto new_lh = std::make_unique<SeqRegion::LHType>(); // = new RealNumType[num_states];
                    RealNumType sum_lh = updateLHwithMat(num_states, model.transposed_mut_mat, *(seq1_region->likelihood), *new_lh, total_blength);

                    // normalize the new partial likelihood
                    normalize_arr(new_lh->data(), num_states, sum_lh);
                    
                    // add merged region into merged_regions
                    merged_regions->emplace_back(TYPE_O, end_pos, 0, 0, std::move(new_lh));
                }
                else
                {
                    // add merged region into merged_regions
                    merged_regions->emplace_back(TYPE_O, end_pos, 0, 0, *(seq1_region->likelihood));
                }
            }
            // seq2_entry = 'N' and seq1_entry = R/ACGT
            else
            {
                RealNumType plength_observation2node = -1;
                RealNumType plength_observation2root = -1;
                
                if (seq1_region->plength_observation2root >= 0)
                {
                    RealNumType plength_from_root = seq1_region->plength_observation2root;
                    if (upper_plength > 0)
                        plength_from_root += upper_plength;
                    
                    plength_observation2node = seq1_region->plength_observation2node;
                    plength_observation2root = plength_from_root;
                }
                else if (seq1_region->plength_observation2node >= 0)
                {
                    RealNumType plength_observation = seq1_region->plength_observation2node;
                    if (upper_plength > 0)
                        plength_observation += upper_plength;
                    
                    plength_observation2node = plength_observation;
                }
                else if (upper_plength > 0)
                        plength_observation2node = upper_plength;
                
                // add a new region and try to merge consecutive R regions together
                addNonConsecutiveRRegion(merged_regions, seq1_region->type, plength_observation2node, plength_observation2root, end_pos, threshold_prob);
            }
        }
        // seq1_entry = seq2_entry = R/ACGT
        else if (seq1_region->type == seq2_region->type && (seq1_region->type < num_states || seq1_region->type == TYPE_R))
        {
            // add a new region and try to merge consecutive R regions together
            addNonConsecutiveRRegion(merged_regions, seq1_region->type, -1, -1, end_pos, threshold_prob);
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
            
            // due to 0 distance, the entry will be of same type as entry2
            if ((seq2_region->type < num_states || seq2_region->type == TYPE_R) && total_blength_2 <= 0)
            {
                if ((seq1_region->type < num_states || seq1_region->type == TYPE_R) && total_blength_1 <= 0)
                {
                    //outError("Sorry! something went wrong. DEBUG: ((seq2_region->type < num_states || seq2_region->type == TYPE_R) && total_blength_2 == 0) && ((seq1_region->type < num_states || seq1_region->type == TYPE_R) && total_blength_1 == 0)");
                    delete merged_regions;
                    merged_regions = NULL;
                    return;
                }
                
                // add a new region and try to merge consecutive R regions together
                addNonConsecutiveRRegion(merged_regions, seq2_region->type, -1, -1, end_pos, threshold_prob);
            }
            // due to 0 distance, the entry will be of same type as entry1
            else if ((seq1_region->type < num_states || seq1_region->type == TYPE_R) && total_blength_1 <= 0)
            {
                // add a new region and try to merge consecutive R regions together
                addNonConsecutiveRRegion(merged_regions, seq1_region->type, -1, -1, end_pos, threshold_prob);
            }
            // seq1_entry = O
            else if (seq1_region->type == TYPE_O)
            {
                auto new_lh = std::make_unique<SeqRegion::LHType>(); // = new RealNumType[num_states];
                auto& new_lh_value = *new_lh;

                // if total_blength_1 > 0 => compute new partial likelihood
                if (total_blength_1 > 0)
                {
                  updateLHwithMat(num_states, model.transposed_mut_mat, *(seq1_region->likelihood), *new_lh, total_blength_1);
                }
                // otherwise, clone the partial likelihood from seq1
                else
                  *new_lh = *seq1_region->likelihood;
                
                RealNumType sum_new_lh = 0;

                // seq2 = O
                if (seq2_region->type == TYPE_O)
                {
                    RealNumType* mutation_mat_row = model.mutation_mat;
                    sum_new_lh = updateMultLHwithMat(num_states, model.mutation_mat, *(seq2_region->likelihood), *new_lh, total_blength_2);
                }
                // seq1 = "O" and seq2 = ACGT
                else
                {
                    StateType seq2_state = seq2_region->type;
                    if (seq2_state == TYPE_R)
                        seq2_state = aln.ref_seq[end_pos];
                    

                    if (total_blength_2 > 0)
                    {
                        RealNumType* transposed_mut_mat_row = model.transposed_mut_mat + model.row_index[seq2_state];
                        for (StateType i = 0; i < num_states; ++i)
                        {
                            if (i == seq2_state)
                              new_lh_value[i] *= (1.0 + transposed_mut_mat_row[i] * total_blength_2);
                            else
                              new_lh_value[i] *= (transposed_mut_mat_row[i] * total_blength_2);
                                
                            sum_new_lh += new_lh_value[i];
                        }
                    }
                    else
                    {
                        for (StateType i = 0; i < num_states; ++i)
                        {
                            if (i != seq2_state)
                              new_lh_value[i] = 0;
                            
                            sum_new_lh += new_lh_value[i];
                        }
                    }
                }
                
                // normalize the new partial lh
                if (sum_new_lh == 0)
                    outError("Sum of new partital lh is zero.");
                
                // normalize the new partial likelihood
                normalize_arr(new_lh->data(), num_states, sum_new_lh);

                StateType new_state = simplifyO(new_lh_value.data(), aln.ref_seq[end_pos], num_states, threshold_prob);

                if (new_state == TYPE_O)
                    merged_regions->emplace_back(TYPE_O, end_pos, 0, 0, std::move(new_lh));
                else
                {
                    // add a new region and try to merge consecutive R regions together
                    addNonConsecutiveRRegion(merged_regions, new_state, -1, -1, end_pos, threshold_prob);
                }
            }
            // seq1_entry = R/ACGT
            else
            {
                StateType seq1_state = seq1_region->type;
                if (seq1_state == TYPE_R)
                    seq1_state = aln.ref_seq[end_pos];
                
                auto new_lh = std::make_unique<SeqRegion::LHType>(); // = new RealNumType[num_states];
                auto& new_lh_value = *new_lh;
                RealNumType sum_new_lh = 0;
                
                if (seq1_region->plength_observation2root >= 0)
                {
                    RealNumType length_to_root = seq1_region->plength_observation2root;
                    if (upper_plength > 0)
                        length_to_root += upper_plength;
                    SeqRegion::LHType root_vec;
                    memcpy(root_vec.data(), model.root_freqs, sizeof(RealNumType)* num_states);
                    
                    RealNumType* transposed_mut_mat_row = model.transposed_mut_mat + model.row_index[seq1_state];
                    for (StateType i = 0; i < num_states; ++i)
                    {
                        if (i == seq1_state)
                            root_vec[i] *= (1.0 + transposed_mut_mat_row[i] * seq1_region->plength_observation2node);
                        else
                            root_vec[i] *= transposed_mut_mat_row[i] * seq1_region->plength_observation2node;
                    }
                    
                    updateLHwithMat(num_states, model.transposed_mut_mat, root_vec, *new_lh, length_to_root);
                }
                else
                {
                    if (total_blength_1 > 0)
                    {
                        RealNumType* mutation_mat_row = model.mutation_mat + model.row_index[seq1_state];
                        for (StateType i = 0; i < num_states; ++i)
                        {
                            if (i == seq1_state)
                              new_lh_value[i] = 1.0 + mutation_mat_row[i] * total_blength_1;
                            else
                              new_lh_value[i] = mutation_mat_row[i] * total_blength_1;
                        }
                    }
                    else
                    {
                      for (auto& v : new_lh_value) v = 0;
                      new_lh_value[seq1_state] = 1;
                    }
                }
                  
                // seq2 = "O" and seq1 = ACGT
                if (seq2_region->type == TYPE_O)
                {
                  sum_new_lh = updateMultLHwithMat(num_states, model.mutation_mat, *(seq2_region->likelihood), *new_lh, total_blength_2);
                  // normalize the new partial likelihood
                  normalize_arr(new_lh->data(), num_states, sum_new_lh);
                  StateType new_state = simplifyO(new_lh_value.data(), aln.ref_seq[end_pos], num_states, threshold_prob);
                    
                  if (new_state == TYPE_O)
                      merged_regions->emplace_back(TYPE_O, end_pos, 0, 0, std::move(new_lh));
                  else
                  {
                      // add a new region and try to merge consecutive R regions together
                      addNonConsecutiveRRegion(merged_regions, new_state, -1, -1, end_pos, threshold_prob);
                  }
                }
                // seq1 = ACGT and different from seq2 = R/ACGT
                else
                {
                    StateType seq2_state = seq2_region->type;
                    
                    if (seq2_state == TYPE_R)
                        seq2_state = aln.ref_seq[end_pos];
                    
                    if (total_blength_2 > 0)
                    {
                        RealNumType* transposed_mut_mat_row = model.transposed_mut_mat + model.row_index[seq2_state];
                        for (StateType i = 0; i < num_states; ++i)
                        {
                            if (i == seq2_state)
                                new_lh_value[i] *= 1.0 + transposed_mut_mat_row[i] * total_blength_2;
                            else
                              new_lh_value[i] *= transposed_mut_mat_row[i] * total_blength_2;
                            
                            sum_new_lh += new_lh_value[i];
                        }
                    }
                    else
                    {
                        for (StateType i = 0; i < num_states; ++i)
                        {
                            if (i != seq2_state)
                              new_lh_value[i] = 0;
                        }
                        sum_new_lh = new_lh_value[seq2_state];
                    }
                    
                    // normalize the new partial likelihood
                    normalize_arr(new_lh->data(), num_states, sum_new_lh);

                    // add new region into the merged regions
                    merged_regions->emplace_back(TYPE_O, end_pos, 0, 0, std::move(new_lh));
                }
            }
        }

        // update pos
        pos = end_pos + 1;
    }
    
    assert(merged_regions->capacity() == max_elements); // ensure we did the correct reserve, otherwise it was a pessimization
}

StateType SeqRegions::simplifyO(RealNumType* const partial_lh, StateType ref_state, StateType num_states, RealNumType threshold_prob) const
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
        if (partial_lh[i] > threshold_prob)
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

RealNumType SeqRegions::mergeTwoLowers(SeqRegions* &merged_regions, RealNumType plength1, const SeqRegions* const regions2, RealNumType plength2, const Alignment& aln, const
                                       Model& model, RealNumType threshold_prob, RealNumType* cumulative_rate, bool return_log_lh)
{
    // init variables
    RealNumType log_lh = 0;
    PositionType pos = 0;
    StateType num_states = aln.num_states;
    const SeqRegions& seq1_regions = *this;
    const SeqRegions& seq2_regions = *regions2;
    size_t iseq1 = 0;
    size_t iseq2 = 0;
    PositionType end_pos;
    PositionType seq_length = aln.ref_seq.size();
    
    // init merged_regions
    if (merged_regions)
        merged_regions->deleteRegions();
    else
        merged_regions = new SeqRegions();

    // avoid realloc of vector data (minimize memory footprint)
    merged_regions->reserve(countSharedSegments(seq2_regions, seq_length)); // avoid realloc of vector data
    const size_t max_elements = merged_regions->capacity(); // remember capacity (may be more than we 'reserved')

    while (pos < seq_length)
    {
        // get the next shared segment in the two sequences
        // get the next shared segment in the two sequences
        SeqRegions::getNextSharedSegment(pos, seq1_regions, seq2_regions, iseq1, iseq2, end_pos); 
        const auto* const seq1_region = &seq1_regions[iseq1];
        const auto* const seq2_region = &seq2_regions[iseq2];

        // seq1_entry = 'N'
        if (seq1_region->type == TYPE_N)
        {
            // seq1_entry = 'N' and seq2_entry = 'N'
            if (seq2_region->type == TYPE_N)
                merged_regions->emplace_back(TYPE_N, end_pos);
            // seq1_entry = 'N' and seq2_entry = O/R/ACGT
            else
            {
                // seq1_entry = 'N' and seq2_entry = O
                if (seq2_region->type == TYPE_O)
                {
                    // add merged region into merged_regions
                    SeqRegion new_region = SeqRegion::clone(*seq2_region);
                    new_region.position = end_pos;
                    if (seq2_region->plength_observation2node >= 0)
                    {
                        if (plength2 > 0)
                            new_region.plength_observation2node += plength2;
                    }
                    else
                    {
                        if (plength2 > 0)
                            new_region.plength_observation2node = plength2;
                    }
                    merged_regions->emplace_back(std::move(new_region));
                }
                // seq1_entry = 'N' and seq2_entry = R/ACGT
                else
                {
                    RealNumType plength_observation2node = -1;
                    
                    if (seq2_region->plength_observation2node >= 0)
                    {
                        RealNumType new_plength = seq2_region->plength_observation2node;
                        if (plength2 > 0)
                            new_plength += plength2;
                        
                        plength_observation2node = new_plength;
                    }
                    else if (plength2 > 0)
                            plength_observation2node = plength2;
                    
                    // add a new region and try to merge consecutive R regions together
                    addNonConsecutiveRRegion(merged_regions, seq2_region->type, plength_observation2node, -1, end_pos, threshold_prob);
                }
            }
        }
        // seq2_entry = 'N'
        else if (seq2_region->type == TYPE_N)
        {
            // seq1_entry = 'O' and seq2_entry = N
            if (seq1_region->type == TYPE_O)
            {
                // add merged region into merged_regions
                auto new_region = SeqRegion::clone(*seq1_region);
                new_region.position = end_pos;
                if (seq1_region->plength_observation2node >= 0)
                {
                    if (plength1 > 0)
                        new_region.plength_observation2node += plength1;
                }
                else
                {
                    if (plength1 > 0)
                        new_region.plength_observation2node = plength1;
                }
                merged_regions->emplace_back(std::move(new_region));
            }
            // seq1_entry = 'N' and seq2_entry = R/ACGT
            else
            {
                RealNumType plength_observation2node = -1;
                
                if (seq1_region->plength_observation2node >= 0)
                {
                    RealNumType new_plength = seq1_region->plength_observation2node;
                    if (plength1 > 0)
                        new_plength += plength1;
                    
                    plength_observation2node = new_plength;
                }
                else if (plength1 > 0)
                        plength_observation2node = plength1;
                
                // add a new region and try to merge consecutive R regions together
                addNonConsecutiveRRegion(merged_regions, seq1_region->type, plength_observation2node, -1, end_pos, threshold_prob);
            }
        }
        // neither seq1_entry nor seq2_entry = N
        else
        {
            RealNumType total_blength_1 = plength1;
            if (seq1_region->plength_observation2node >= 0)
            {
                total_blength_1 = seq1_region->plength_observation2node;
                if (plength1 > 0)
                    total_blength_1 += plength1;
            }
            
            RealNumType total_blength_2 = plength2;
            if (seq2_region->plength_observation2node >= 0)
            {
                total_blength_2 = seq2_region->plength_observation2node;
                if (plength2 > 0)
                    total_blength_2 += plength2;
            }
            
            // seq1_entry and seq2_entry are identical seq1_entry = R/ACGT
            if (seq1_region->type == seq2_region->type && (seq1_region->type == TYPE_R || seq1_region->type < num_states))
            {
                // add a new region and try to merge consecutive R regions together
                addNonConsecutiveRRegion(merged_regions, seq1_region->type, -1, -1, end_pos, threshold_prob);
                
                if (return_log_lh)
                {
                    // convert total_blength_1 and total_blength_2 to zero if they are -1
                    if (total_blength_1 < 0) total_blength_1 = 0;
                    if (total_blength_2 < 0) total_blength_2 = 0;
                    
                    if (seq1_region->type == TYPE_R)
                        log_lh += (total_blength_1 + total_blength_2) * (cumulative_rate[end_pos + 1] - cumulative_rate[pos]);
                    else
                        log_lh += model.diagonal_mut_mat[seq1_region->type] * (total_blength_1 + total_blength_2);
                }
            }
            // #0 distance between different nucleotides: merge is not possible
            else if (total_blength_1 == 0 && total_blength_2 == 0 && (seq1_region->type == TYPE_R || seq1_region->type < num_states) && (seq2_region->type == TYPE_R || seq2_region->type < num_states))
            {
                delete merged_regions;
                merged_regions = NULL;
                return MIN_NEGATIVE;
            }
            // seq1_entry = O
            else if (seq1_region->type == TYPE_O)
            {
                auto new_lh = std::make_unique<SeqRegion::LHType>(); // = new RealNumType[num_states];
                auto& new_lh_value = *new_lh;
                RealNumType sum_lh = 0;
                if (total_blength_1 > 0)
                {
                  sum_lh = updateLHwithMat(num_states, model.mutation_mat, *(seq1_region->likelihood), *new_lh, total_blength_1);
                }
                // otherwise, clone the partial likelihood from seq1
                else
                    *new_lh = *seq1_region->likelihood;

                // seq1_entry = O and seq2_entry = O
                if (seq2_region->type == TYPE_O)
                {
                    sum_lh = updateMultLHwithMat(num_states, model.mutation_mat, *(seq2_region->likelihood), *new_lh, total_blength_2);
                    
                    if (sum_lh == 0)
                    {
                        delete merged_regions;
                        merged_regions = NULL;
                        return MIN_NEGATIVE;
                    }
                        
                    // normalize the new partial likelihood
                    normalize_arr(new_lh->data(), num_states, sum_lh);

                    StateType new_state = simplifyO(new_lh_value.data(), aln.ref_seq[end_pos], num_states, threshold_prob);

                    if (new_state == TYPE_O)
                        merged_regions->emplace_back(new_state, end_pos, 0, 0, std::move(new_lh));
                    else
                    {
                        // add a new region and try to merge consecutive R regions together
                        addNonConsecutiveRRegion(merged_regions, new_state, -1, -1, end_pos, threshold_prob);
                    }
                    
                    if (return_log_lh)
                        log_lh += log(sum_lh);
                }
                // seq1_entry = O and seq2_entry = ACGT
                else
                {
                    StateType seq2_state = seq2_region->type;
                    if (seq2_state == TYPE_R)
                        seq2_state = aln.ref_seq[end_pos];
                    
                    if (total_blength_2 > 0)
                    {
                        RealNumType* transposed_mut_mat_row = model.transposed_mut_mat + model.row_index[seq2_state];
                        for (StateType i = 0; i < num_states; ++i)
                        {
                            if (seq2_state == i)
                                new_lh_value[i] *= (1 + transposed_mut_mat_row[i] * total_blength_2);
                            else
                              new_lh_value[i] *= (transposed_mut_mat_row[i] * total_blength_2);
                            
                            sum_lh += new_lh_value[i];
                        }
                        
                        // normalize new partial lh
                        // normalize the new partial likelihood
                        normalize_arr(new_lh->data(), num_states, sum_lh);

                        StateType new_state = simplifyO(new_lh_value.data(), aln.ref_seq[end_pos],  num_states, threshold_prob);

                        if (new_state == TYPE_O)
                            merged_regions->emplace_back(new_state, end_pos, 0, 0, std::move(new_lh));
                        else
                        {
                            // add a new region and try to merge consecutive R regions together
                            addNonConsecutiveRRegion(merged_regions, new_state, -1, -1, end_pos, threshold_prob);
                        }
                        
                        if (return_log_lh)
                            log_lh += log(sum_lh);
                    }
                    else
                    {
                        if (new_lh_value[seq2_state] == 0)
                        {
                            delete merged_regions;
                            merged_regions = NULL;
                            return MIN_NEGATIVE;
                        }
                        
                        // add a new region and try to merge consecutive R regions together
                        addNonConsecutiveRRegion(merged_regions, seq2_region->type, -1, -1, end_pos, threshold_prob);
                    
                        if (return_log_lh)
                            log_lh += log(new_lh_value[seq2_state]);
                    }
                }
            }
            // seq1_entry = R/ACGT
            else
            {
                StateType seq1_state = seq1_region->type;
                if (seq1_state == TYPE_R)
                    seq1_state = aln.ref_seq[end_pos];
                
                auto new_lh = std::make_unique<SeqRegion::LHType>(); // = new RealNumType[num_states];
                auto& new_lh_value = *new_lh;
                RealNumType sum_lh = 0;
                
                if (total_blength_1 > 0)
                {
                    RealNumType* transposed_mut_mat_row = model.transposed_mut_mat + model.row_index[seq1_state];
                    
                    for (StateType i = 0; i < num_states; ++i)
                    {
                        if (seq1_state == i)
                          new_lh_value[i] = 1 + transposed_mut_mat_row[i] * total_blength_1;
                        else
                          new_lh_value[i] = transposed_mut_mat_row[i] * total_blength_1;
                    }
                }
                else
                {
                    for (StateType i = 0; i < num_states; ++i)
                      new_lh_value[i] = 0;
                    new_lh_value[seq1_state] = 1;
                }

                // seq1_entry = ACGT and seq2_entry = O
                if (seq2_region->type == TYPE_O)
                {
                    sum_lh = updateMultLHwithMat(num_states, model.mutation_mat, *(seq2_region->likelihood), *new_lh, total_blength_2);
                                        
                    if (sum_lh == 0)
                    {
                        delete merged_regions;
                        merged_regions = NULL;
                        return MIN_NEGATIVE;
                    }
                        
                    // normalize the new partial likelihood
                    normalize_arr(new_lh->data(), num_states, sum_lh);
                    
                    StateType new_state = simplifyO(new_lh_value.data(), aln.ref_seq[end_pos], num_states, threshold_prob);

                    if (new_state == TYPE_O)
                        merged_regions->emplace_back(new_state, end_pos, 0, 0, std::move(new_lh));
                    else
                    {
                        // add a new region and try to merge consecutive R regions together
                        addNonConsecutiveRRegion(merged_regions, new_state, -1, -1, end_pos, threshold_prob);
                    }
                    
                    if (return_log_lh)
                        log_lh += log(sum_lh);
                }
                // seq1_entry = ACGT and seq2_entry = R/ACGT
                else
                {
                    StateType seq2_state = seq2_region->type;
                    if (seq2_state == TYPE_R)
                        seq2_state = aln.ref_seq[end_pos];
                    
                    if (total_blength_2 > 0)
                    {
                        RealNumType* transposed_mut_mat_row = model.transposed_mut_mat + model.row_index[seq2_state];
                        for (StateType i = 0; i < num_states; ++i)
                        {
                            if (seq2_state == i)
                              new_lh_value[i] *= (1 + transposed_mut_mat_row[i] * total_blength_2);
                            else
                              new_lh_value[i] *= (transposed_mut_mat_row[i] * total_blength_2);
                            
                            sum_lh += new_lh_value[i];
                        }
                        
                        // normalize the new partial likelihood
                        normalize_arr(new_lh->data(), num_states, sum_lh);
                        StateType new_state = simplifyO(new_lh_value.data(), aln.ref_seq[end_pos], num_states, threshold_prob);

                        if (new_state == TYPE_O)
                            merged_regions->emplace_back(new_state, end_pos, 0, 0, std::move(new_lh));
                        else
                        {
                            // add a new region and try to merge consecutive R regions together
                            addNonConsecutiveRRegion(merged_regions, new_state, -1, -1, end_pos, threshold_prob);
                        }
                        
                        if (return_log_lh)
                            log_lh += log(sum_lh);
                    }
                    else
                    {
                        // add a new region and try to merge consecutive R regions together
                        addNonConsecutiveRRegion(merged_regions, seq2_region->type, -1, -1, end_pos, threshold_prob);
                    
                        if (return_log_lh)
                            log_lh += log(new_lh_value[seq2_state]);
                    }
                }
            }
        }
        // update pos
        pos = end_pos + 1;
    }
    
    assert(merged_regions->capacity() == max_elements); // ensure we did the correct reserve, otherwise it was a pessimization

    return log_lh;
}

RealNumType SeqRegions::computeAbsoluteLhAtRoot(const Alignment& aln, const Model& model, vector< vector<PositionType> > &cumulative_base)
{
    // dummy variables
    RealNumType log_lh = 0;
    RealNumType log_factor = 1;
    const StateType num_states = aln.num_states;
    PositionType start_pos = 0;
    const SeqRegions& regions = *this;
    
    // browse regions one by one to compute the likelihood of each region
    for (const SeqRegion& region : regions)
    {
        // type R
        if (region.type == TYPE_R)
        {
            for (StateType i = 0; i < num_states; ++i)
                log_lh += model.root_log_freqs[i] * (cumulative_base[region.position + 1][i] - cumulative_base[start_pos][i]);
        }
        // type ACGT
        else if (region.type < num_states)
            log_lh += model.root_log_freqs[region.type];
        // type O
        else if (region.type == TYPE_O)
        {
            RealNumType tot = 0;
            for (StateType i = 0; i < num_states; ++i)
                tot += model.root_freqs[i] * region.getLH(i);
            log_factor *= tot;
        }
        
        // maintain start_pos
        start_pos = region.position + 1;
    }

    // update log_lh
    log_lh += log(log_factor);
    
    // return the absolute likelihood
    return log_lh;
}

SeqRegions* SeqRegions::computeTotalLhAtRoot(StateType num_states, const Model& model, RealNumType blength)
{
    SeqRegions* total_lh = new SeqRegions();
    total_lh->reserve(this->size()); // avoid realloc of vector data
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
                RealNumType sum_lh = updateLHwithModel(num_states, model, *region->likelihood, (*new_lh), total_blength);
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
    
    return total_lh;
}
