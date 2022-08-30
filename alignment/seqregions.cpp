//
//  regions.cpp
//  alignment
//
//  Created by NhanLT on 31/3/2022.
//

#include "seqregions.h"

SeqRegions::SeqRegions()
{
    // do nothing
}

SeqRegions::~SeqRegions()
{
    deleteRegions();
}

void SeqRegions::deleteRegions()
{
    for (iterator it = begin(); it != end(); ++it)
        delete (*it);
    clear();
}

void SeqRegions::copyRegions(SeqRegions* n_regions, StateType num_states)
{
    // delete the current regions
    deleteRegions();
    
    // clone regions one by one
    resize(n_regions->size());
    for (PositionType i = 0; i < (PositionType) n_regions->size(); ++i)
    {
        SeqRegion* region = n_regions->at(i);
        at(i) = new SeqRegion(region, num_states, true);
    }
}

void SeqRegions::mergeRegionR(StateType num_states, RealNumType threshold)
{
    // dummy variables
    PositionType start_R_index = -1;
    
    // browse regions one by one to detect identical sequences of R regions
    for (PositionType i = 0; i < (PositionType) size(); ++i)
    {
        // get the current region
        SeqRegion* current_region = at(i);
        
        // if start_R_index is not set -> check if the current region is an R region
        if (start_R_index == -1)
        {
            if (current_region->type == TYPE_R)
                start_R_index = i;
        }
        // otherwise, check if the current region is an R region and identical to the starting R region
        else
        {
            SeqRegion* start_R_region = at(start_R_index);
            
            // if the current region is NOT identical to the starting region -> merging sequence of identical R regions starting from start_R_index to the region before the current region
            if (!(current_region->type == TYPE_R && fabs(start_R_region->plength_observation2node - current_region->plength_observation2node) < threshold && fabs(start_R_region->plength_observation2root - current_region->plength_observation2root) < threshold))
             {
                 // remove identical R regions
                 for (PositionType j = 0; j < i - start_R_index - 1; ++j)
                     delete at(start_R_index + j);
                 erase(begin() + start_R_index, begin() + i - 1);
                 
                 // recompute i
                 i -= i - start_R_index - 1;
                 
                 // reset start_R_index
                 start_R_index = -1;
             }
        }
    }
    
    // merge the last R regions (if any)
    if  (start_R_index != -1)
    {
        PositionType num_regions = size();
        
        // remove identical R regions
        for (PositionType j = 0; j < num_regions - 2 - start_R_index; ++j)
            delete at(start_R_index + j);
        erase(begin() + start_R_index, end() - 1);
    }
}

int SeqRegions::compareWithSample(SeqRegions* sequence2, PositionType seq_length, StateType num_states)
{
    ASSERT(seq_length > 0);
    
    // init dummy variables
    bool seq1_more_info = false;
    bool seq2_more_info = false;
    PositionType pos = 0;
    SeqRegion **seq1_region_pointer = &front();
    SeqRegion **seq2_region_pointer = &sequence2->front();
    SeqRegion *seq1_region = (*seq1_region_pointer);
    SeqRegion *seq2_region = (*seq2_region_pointer);
    PositionType end_pos;

    while (pos < seq_length && (!seq1_more_info || !seq2_more_info))
    {
        // get the next shared segment in the two sequences
        SeqRegions::getNextSharedSegment(pos, seq1_region, seq2_region, seq1_region_pointer, seq2_region_pointer, end_pos);
        
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
                if (seq2_region->likelihood[i] > 0.1 && seq1_region->likelihood[i] < 0.1)
                    seq1_more_info = true;
                else if (seq1_region->likelihood[i] > 0.1 && seq2_region->likelihood[i] < 0.1)
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

bool SeqRegions::areDiffFrom(SeqRegions* regions2, PositionType seq_length, StateType num_states, Params* params)
{
    if (!regions2 || regions2->size() == 0)
        return true;
    
    // init variables
    PositionType pos = 0;
    SeqRegion **seq1_region_pointer = &front();
    SeqRegion **seq2_region_pointer = &regions2->front();
    SeqRegion *seq1_region = (*seq1_region_pointer);
    SeqRegion *seq2_region = (*seq2_region_pointer);
    PositionType end_pos;
    
    // compare each pair of regions
    while (pos < seq_length)
    {
        // get the next shared segment in the two sequences
        SeqRegions::getNextSharedSegment(pos, seq1_region, seq2_region, seq1_region_pointer, seq2_region_pointer, end_pos);
        
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
                RealNumType diff = fabs(seq1_region->likelihood[i] - seq2_region->likelihood[i]);
                
                if (diff > 0)
                {
                    if ((seq1_region->likelihood[i] == 0) || (seq2_region->likelihood[i] == 0))
                        return true;
                    
                    if (diff > params->thresh_diff_update
                        || (diff > params->threshold_prob
                            && ((diff / seq1_region->likelihood[i] > params->thresh_diff_update)
                                || (diff / seq2_region->likelihood[i] > params->thresh_diff_update))))
                        return true;
                }
            }
        }
        
        // update pos
        pos = end_pos + 1;
    }
    
    return false;
}

void SeqRegions::mergeUpperLower(SeqRegions* &merged_regions, RealNumType upper_plength, SeqRegions* lower_regions, RealNumType lower_plength, Alignment* aln, Model* model, RealNumType threshold_prob)
{
    // init variables
    PositionType pos = 0;
    SeqRegion **seq1_region_pointer = &front();
    SeqRegion **seq2_region_pointer = &lower_regions->front();
    SeqRegion *seq1_region = (*seq1_region_pointer);
    SeqRegion *seq2_region = (*seq2_region_pointer);
    PositionType end_pos;
    StateType num_states = aln->num_states;
    PositionType seq_length = aln->ref_seq.size();
    
    // init merged_regions
    if (merged_regions)
        merged_regions->deleteRegions();
    else
        merged_regions = new SeqRegions();
                
    while (pos <  seq_length)
    {
        // get the next shared segment in the two sequences
        SeqRegions::getNextSharedSegment(pos, seq1_region, seq2_region, seq1_region_pointer, seq2_region_pointer, end_pos);
        
        // seq1_entry = 'N'
        if (seq1_region->type == TYPE_N)
        {
            // seq1_entry = 'N' and seq2_entry = 'N'
            if (seq2_region->type == TYPE_N)
                merged_regions->push_back(new SeqRegion(seq1_region->type, end_pos));
            // seq1_entry = 'N' and seq2_entry = O/R/ACGT
            else
            {
                // seq1_entry = 'N' and seq2_entry = O
                if (seq2_region->type == TYPE_O)
                {
                    RealNumType total_blength = lower_plength;
                    if (seq2_region->plength_observation2node >= 0)
                        total_blength = seq2_region->plength_observation2node + (lower_plength > 0 ? lower_plength: 0);
                        
                    RealNumType* new_lh = new RealNumType[num_states];
                    RealNumType sum_lh = 0;
                    
                    for (StateType i = 0; i < num_states; ++i)
                    {
                        RealNumType tot = 0;
                        
                        if (total_blength > 0)
                        {
                            for (StateType j = 0; j < num_states; ++j)
                                tot += model->mutation_mat[model->row_index[i] + j] * seq2_region->likelihood[j];
                            
                            tot *= total_blength;
                        }
                        
                        tot += seq2_region->likelihood[i];
                        new_lh[i] = tot * model->root_freqs[i];
                        sum_lh += new_lh[i];
                    }

                    // normalize the new partial likelihood
                    RealNumType inverse_sum_lh = 1.0 / sum_lh;
                    for (StateType i = 0; i < num_states; ++i)
                        new_lh[i] *= inverse_sum_lh;
                    
                    // add merged region into merged_regions
                    merged_regions->push_back(new SeqRegion(TYPE_O, end_pos, 0, 0, new_lh));
                }
                // seq1_entry = 'N' and seq2_entry = R/ACGT
                else
                {
                    if (seq2_region->plength_observation2node >= 0)
                    {
                        RealNumType new_plength = seq2_region->plength_observation2node ;
                        if (lower_plength > 0)
                            new_plength += lower_plength;
                        merged_regions->push_back(new SeqRegion(seq2_region->type, end_pos, new_plength, 0));
                    }
                    else
                    {
                        if (lower_plength > 0)
                            merged_regions->push_back(new SeqRegion(seq2_region->type, end_pos, lower_plength, 0));
                        else
                            merged_regions->push_back(new SeqRegion(seq2_region->type, end_pos));
                    }
                    
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
                    RealNumType* new_lh = new RealNumType[num_states];
                    RealNumType sum_lh = 0;
                    
                    for (StateType i = 0; i < num_states; ++i)
                    {
                        RealNumType tot = 0;
                        
                        for (StateType j = 0; j < num_states; ++j)
                            tot += model->transposed_mut_mat[model->row_index[i] + j] * seq1_region->likelihood[j];
                        
                        tot *= total_blength;
                        tot += seq1_region->likelihood[i];
                        new_lh[i] = tot;
                        sum_lh += new_lh[i];
                    }

                    // normalize the new partial likelihood
                    RealNumType inverse_sum_lh = 1.0 / sum_lh;
                    for (StateType i = 0; i < num_states; ++i)
                        new_lh[i] *= inverse_sum_lh;
                    
                    // add merged region into merged_regions
                    merged_regions->push_back(new SeqRegion(TYPE_O, end_pos, 0, 0, new_lh));
                }
                else
                {
                    RealNumType* new_lh = new RealNumType[num_states];
                    memcpy(new_lh, seq1_region->likelihood, sizeof(RealNumType) * num_states);
                    // add merged region into merged_regions
                    merged_regions->push_back(new SeqRegion(seq1_region->type, end_pos, 0, 0, new_lh));
                }
            }
            // seq2_entry = 'N' and seq1_entry = R/ACGT
            else
            {
                if (seq1_region->plength_observation2root >= 0)
                {
                    RealNumType plength_from_root = seq1_region->plength_observation2root;
                    if (upper_plength > 0)
                        plength_from_root += upper_plength;
                    merged_regions->push_back(new SeqRegion(seq1_region->type, end_pos, seq1_region->plength_observation2node, plength_from_root));
                }
                else if (seq1_region->plength_observation2node >= 0)
                {
                    RealNumType plength_observation = seq1_region->plength_observation2node;
                    if (upper_plength > 0)
                        plength_observation += upper_plength;
                    merged_regions->push_back(new SeqRegion(seq1_region->type, end_pos, plength_observation));
                }
                else
                {
                    if (upper_plength > 0)
                        merged_regions->push_back(new SeqRegion(seq1_region->type, end_pos, upper_plength));
                    else
                        merged_regions->push_back(new SeqRegion(seq1_region->type, end_pos));
                }
            }
        }
        // seq1_entry = seq2_entry = R/ACGT
        else if (seq1_region->type == seq2_region->type && (seq1_region->type < num_states || seq1_region->type == TYPE_R))
            merged_regions->push_back(new SeqRegion(seq1_region->type, end_pos));
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
                
                merged_regions->push_back(new SeqRegion(seq2_region->type, end_pos));
            }
            // due to 0 distance, the entry will be of same type as entry1
            else if ((seq1_region->type < num_states || seq1_region->type == TYPE_R) && total_blength_1 <= 0)
            {
                merged_regions->push_back(new SeqRegion(seq1_region->type, end_pos));
            }
            // seq1_entry = O
            else if (seq1_region->type == TYPE_O)
            {
                RealNumType* new_lh = new RealNumType[num_states];
                
                // if total_blength_1 > 0 => compute new partial likelihood
                if (total_blength_1 > 0)
                {
                    for (StateType i = 0; i < num_states; ++i)
                    {
                        RealNumType tot = 0;
                        
                        for (StateType j = 0; j < num_states; ++j)
                            tot += model->transposed_mut_mat[model->row_index[i] + j] * seq1_region->likelihood[j];
                        
                        tot *= total_blength_1;
                        tot += seq1_region->likelihood[i];
                        new_lh[i] = tot;
                    }
                }
                // otherwise, clone the partial likelihood from seq1
                else
                    memcpy(new_lh, seq1_region->likelihood, sizeof(RealNumType) * num_states);
                
                RealNumType sum_new_lh = 0;

                // seq2 = O
                if (seq2_region->type == TYPE_O)
                {
                    for (StateType i = 0; i < num_states; ++i)
                    {
                        RealNumType tot = 0;
                        
                        if (total_blength_2 > 0)
                        {
                            for (StateType j = 0; j < num_states; ++j)
                                tot += model->mutation_mat[model->row_index[i] + j] * seq2_region->likelihood[j];
        
                            tot *= total_blength_2;
                        }
                        tot += seq2_region->likelihood[i];
                        new_lh[i] *= tot;
                        sum_new_lh += new_lh[i];
                    }
                }
                // seq1 = "O" and seq2 = ACGT
                else
                {
                    StateType seq2_state = seq2_region->type;
                    if (seq2_state == TYPE_R)
                        seq2_state = aln->ref_seq[end_pos];
                    
                    if (total_blength_2 > 0)
                    {
                        StateType mutation_index = model->row_index[seq2_state];
                        for (StateType i = 0; i < num_states; ++i)
                        {
                            if (i == seq2_state)
                                new_lh[i] *= (1.0 + model->transposed_mut_mat[mutation_index + i] * total_blength_2);
                            else
                                new_lh[i] *= (model->transposed_mut_mat[mutation_index + i] * total_blength_2);
                                
                            sum_new_lh += new_lh[i];
                        }
                    }
                    else
                    {
                        for (StateType i = 0; i < num_states; ++i)
                        {
                            if (i != seq2_state)
                                new_lh[i] = 0;
                            
                            sum_new_lh += new_lh[i];
                        }
                    }
                }
                
                // normalize the new partial lh
                if (sum_new_lh == 0)
                    outError("Sum of new partital lh is zero.");
                
                RealNumType inverse_sum_lh = 1.0 / sum_new_lh;
                for (StateType i = 0; i < num_states; ++i)
                    new_lh[i] *= inverse_sum_lh;
                
                StateType new_state = simplifyO(new_lh, aln->ref_seq[end_pos], num_states, threshold_prob);

                if (new_state == TYPE_O)
                    merged_regions->push_back(new SeqRegion(TYPE_O, end_pos, 0, 0, new_lh));
                else
                {
                    delete[] new_lh;
                    merged_regions->push_back(new SeqRegion(new_state, end_pos));
                }
            }
            // seq1_entry = R/ACGT
            else
            {
                StateType seq1_state = seq1_region->type;
                if (seq1_state == TYPE_R)
                    seq1_state = aln->ref_seq[end_pos];
                
                RealNumType* new_lh = new RealNumType[num_states];
                RealNumType sum_new_lh = 0;
                
                if (seq1_region->plength_observation2root >= 0)
                {
                    RealNumType length_to_root = seq1_region->plength_observation2root;
                    if (upper_plength > 0)
                        length_to_root += upper_plength;
                    RealNumType* root_vec = new RealNumType[num_states];
                    memcpy(root_vec, model->root_freqs, sizeof(RealNumType) * num_states);
                    
                    StateType mutation_index = model->row_index[seq1_state];
                    for (StateType i = 0; i < num_states; ++i)
                    {
                        if (i == seq1_state)
                            root_vec[i] *= (1.0 + model->transposed_mut_mat[mutation_index + i] * seq1_region->plength_observation2node);
                        else
                            root_vec[i] *= model->transposed_mut_mat[mutation_index + i] * seq1_region->plength_observation2node;
                    }
                    
                    for (StateType i = 0; i < num_states; ++i)
                    {
                        RealNumType tot = 0;
                        
                        for (StateType j = 0; j < num_states; ++j)
                            tot += model->transposed_mut_mat[model->row_index[i] + j] * root_vec[j];
                        
                        tot *= length_to_root;
                        tot += root_vec[i];
                        new_lh[i] = tot;
                    }
                }
                else
                {
                    if (total_blength_1 > 0)
                    {
                        StateType mutation_index = model->row_index[seq1_state];
                        for (StateType i = 0; i < num_states; ++i)
                        {
                            if (i == seq1_state)
                                new_lh[i] = 1.0 + model->mutation_mat[mutation_index + i] * total_blength_1;
                            else
                                new_lh[i] = model->mutation_mat[mutation_index + i] * total_blength_1;
                        }
                    }
                    else
                    {
                        for (StateType i = 0; i < num_states; ++i)
                            new_lh[i] = 0;
                        new_lh[seq1_state] = 1;
                    }
                }
                  
                // seq2 = "O" and seq1 = ACGT
                if (seq2_region->type == TYPE_O)
                {
                    for (StateType i = 0; i < num_states; ++i)
                    {
                        RealNumType tot = 0;
                        
                        if (total_blength_2 > 0)
                        {
                            for (StateType j = 0; j < num_states; ++j)
                                tot += model->mutation_mat[model->row_index[i] + j] * seq2_region->likelihood[j];
                            
                            tot *= total_blength_2;
                        }
                        
                        tot += seq2_region->likelihood[i];
                        new_lh[i] *= tot;
                        sum_new_lh += new_lh[i];
                    }
                    
                    // normalize the new partial lh
                    RealNumType inverse_sum_lh = 1.0 / sum_new_lh;
                    for (StateType i = 0; i < num_states; ++i)
                        new_lh[i] *= inverse_sum_lh;
                    
                    StateType new_state = simplifyO(new_lh, aln->ref_seq[end_pos], num_states, threshold_prob);
                    
                    if (new_state == TYPE_O)
                        merged_regions->push_back(new SeqRegion(TYPE_O, end_pos, 0, 0, new_lh));
                    else
                        merged_regions->push_back(new SeqRegion(new_state, end_pos));
                }
                // seq1 = ACGT and different from seq2 = R/ACGT
                else
                {
                    StateType seq2_state = seq2_region->type;
                    
                    if (seq2_state == TYPE_R)
                        seq2_state = aln->ref_seq[end_pos];
                    
                    if (total_blength_2 > 0)
                    {
                        StateType mutation_index = model->row_index[seq2_state];
                        for (StateType i = 0; i < num_states; ++i)
                        {
                            if (i == seq2_state)
                                new_lh[i] *= 1.0 + model->transposed_mut_mat[mutation_index + i] * total_blength_2;
                            else
                                new_lh[i] *= model->transposed_mut_mat[mutation_index + i] * total_blength_2;
                            
                            sum_new_lh += new_lh[i];
                        }
                    }
                    else
                    {
                        for (StateType i = 0; i < num_states; ++i)
                        {
                            if (i != seq2_state)
                                new_lh[i] = 0;
                                
                            sum_new_lh += new_lh[i];
                        }
                    }
                    
                    // normalize the new partial lh
                    RealNumType inverse_sum_lh = 1.0 / sum_new_lh;
                    for (StateType i = 0; i < num_states; ++i)
                        new_lh[i] *= inverse_sum_lh;

                    // add new region into the merged regions
                    merged_regions->push_back(new SeqRegion(TYPE_O, end_pos, 0, 0, new_lh));
                }
            }
        }

        // update pos
        pos = end_pos + 1;
    }

    // try to merge consecutive and similar 'R' regions together
    merged_regions->mergeRegionR(num_states, threshold_prob);
}

StateType SeqRegions::simplifyO(RealNumType* &partial_lh, StateType ref_state, StateType num_states, RealNumType threshold_prob)
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

RealNumType SeqRegions::mergeTwoLowers(SeqRegions* &merged_regions, RealNumType plength1, SeqRegions* regions2, RealNumType plength2, Alignment* aln, Model* model, RealNumType threshold_prob, RealNumType* cumulative_rate, bool return_log_lh)
{
    // init variables
    RealNumType log_lh = 0;
    PositionType pos = 0;
    StateType num_states = aln->num_states;
    SeqRegion **seq1_region_pointer = &front();
    SeqRegion **seq2_region_pointer = &regions2->front();
    SeqRegion *seq1_region = (*seq1_region_pointer);
    SeqRegion *seq2_region = (*seq2_region_pointer);
    PositionType end_pos;
    PositionType seq_length = aln->ref_seq.size();
    
    // init merged_regions
    if (merged_regions)
        merged_regions->deleteRegions();
    else
        merged_regions = new SeqRegions();
                
    while (pos < seq_length)
    {
        // get the next shared segment in the two sequences
        SeqRegions::getNextSharedSegment(pos, seq1_region, seq2_region, seq1_region_pointer, seq2_region_pointer, end_pos);
        
        // seq1_entry = 'N'
        if (seq1_region->type == TYPE_N)
        {
            // seq1_entry = 'N' and seq2_entry = 'N'
            if (seq2_region->type == TYPE_N)
                merged_regions->push_back(new SeqRegion(seq1_region->type, end_pos));
            // seq1_entry = 'N' and seq2_entry = O/R/ACGT
            else
            {
                // seq1_entry = 'N' and seq2_entry = O
                if (seq2_region->type == TYPE_O)
                {
                    // add merged region into merged_regions
                    SeqRegion* new_region = new SeqRegion(seq2_region, num_states);
                    new_region->position = end_pos;
                    if (seq2_region->plength_observation2node >= 0)
                    {
                        if (plength2 > 0)
                            new_region->plength_observation2node += plength2;
                    }
                    else
                    {
                        if (plength2 > 0)
                            new_region->plength_observation2node = plength2;
                    }
                    merged_regions->push_back(new_region);
                }
                // seq1_entry = 'N' and seq2_entry = R/ACGT
                else
                {
                    if (seq2_region->plength_observation2node >= 0)
                    {
                        RealNumType new_plength = seq2_region->plength_observation2node;
                        if (plength2 > 0)
                            new_plength += plength2;
                        merged_regions->push_back(new SeqRegion(seq2_region->type, end_pos, new_plength));
                    }
                    else
                    {
                        if (plength2 > 0)
                            merged_regions->push_back(new SeqRegion(seq2_region->type, end_pos, plength2));
                        else
                            merged_regions->push_back(new SeqRegion(seq2_region->type, end_pos));
                    }
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
                SeqRegion* new_region = new SeqRegion(seq1_region, num_states);
                new_region->position = end_pos;
                if (seq1_region->plength_observation2node >= 0)
                {
                    if (plength1 > 0)
                        new_region->plength_observation2node += plength1;
                }
                else
                {
                    if (plength1 > 0)
                        new_region->plength_observation2node = plength1;
                }
                merged_regions->push_back(new_region);
            }
            // seq1_entry = 'N' and seq2_entry = R/ACGT
            else
            {
                if (seq1_region->plength_observation2node >= 0)
                {
                    RealNumType new_plength = seq1_region->plength_observation2node;
                    if (plength1 > 0)
                        new_plength += plength1;
                    merged_regions->push_back(new SeqRegion(seq1_region->type, end_pos, new_plength));
                }
                else
                {
                    if (plength1 > 0)
                        merged_regions->push_back(new SeqRegion(seq1_region->type, end_pos, plength1));
                    else
                        merged_regions->push_back(new SeqRegion(seq1_region->type, end_pos));
                        
                }
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
                merged_regions->push_back(new SeqRegion(seq1_region->type, end_pos));
                
                if (return_log_lh)
                {
                    // convert total_blength_1 and total_blength_2 to zero if they are -1
                    if (total_blength_1 < 0) total_blength_1 = 0;
                    if (total_blength_2 < 0) total_blength_2 = 0;
                    
                    if (seq1_region->type == TYPE_R)
                        log_lh += (total_blength_1 + total_blength_2) * (cumulative_rate[end_pos + 1] - cumulative_rate[pos]);
                    else
                        log_lh += model->diagonal_mut_mat[seq1_region->type] * (total_blength_1 + total_blength_2);
                }
            }
            // #0 distance between different nucleotides: merge is not possible
            else if (total_blength_1 == 0 && total_blength_2 == 0 && (seq1_region->type == TYPE_R || seq1_region->type < num_states) && (seq2_region->type == TYPE_R || seq2_region->type < num_states))
            {
                delete merged_regions;
                merged_regions = NULL;
                return -DBL_MAX;
            }
            // seq1_entry = O
            else if (seq1_region->type == TYPE_O)
            {
                RealNumType* new_lh = new RealNumType[num_states];
                RealNumType sum_lh = 0;
                
                if (total_blength_1 > 0)
                {
                    for (StateType i = 0; i < num_states; ++i)
                    {
                        RealNumType tot = 0;
                        
                        for (StateType j = 0; j < num_states; ++j)
                            tot += model->mutation_mat[model->row_index[i] + j] * seq1_region->likelihood[j];
                        
                        tot *= total_blength_1;
                        tot += seq1_region->likelihood[i];
                        new_lh[i] = tot;
                    }
                }
                // otherwise, clone the partial likelihood from seq1
                else
                    memcpy(new_lh, seq1_region->likelihood, sizeof(RealNumType) * num_states);

                // seq1_entry = O and seq2_entry = O
                if (seq2_region->type == TYPE_O)
                {
                    for (StateType i = 0; i < num_states; ++i)
                    {
                        RealNumType tot = 0;
                        
                        if (total_blength_2 > 0)
                        {
                            for (StateType j = 0; j < num_states; ++j)
                                tot += model->mutation_mat[model->row_index[i] + j] * seq2_region->likelihood[j];
                            
                            tot *= total_blength_2;
                        }
                        tot += seq2_region->likelihood[i];
                        new_lh[i] *= tot;
                        sum_lh += new_lh[i];
                    }
                    
                    if (sum_lh == 0)
                    {
                        delete merged_regions;
                        merged_regions = NULL;
                        return -DBL_MAX;
                    }
                        
                    // normalize new partial lh
                    RealNumType inverse_sum_lh = 1.0 / sum_lh;
                    for (StateType i = 0; i < num_states; ++i)
                        new_lh[i] *= inverse_sum_lh;
                    
                    StateType new_state = simplifyO(new_lh, aln->ref_seq[end_pos], num_states, threshold_prob);

                    if (new_state == TYPE_O)
                        merged_regions->push_back(new SeqRegion(new_state, end_pos, 0, 0, new_lh));
                    else
                        merged_regions->push_back(new SeqRegion(new_state, end_pos));
                    
                    if (return_log_lh)
                        log_lh += log(sum_lh);
                }
                // seq1_entry = O and seq2_entry = ACGT
                else
                {
                    StateType seq2_state = seq2_region->type;
                    if (seq2_state == TYPE_R)
                        seq2_state = aln->ref_seq[end_pos];
                    
                    if (total_blength_2 > 0)
                    {
                        StateType mutation_index = model->row_index[seq2_state];
                        for (StateType i = 0; i < num_states; ++i)
                        {
                            if (seq2_state == i)
                                new_lh[i] *= (1 + model->transposed_mut_mat[mutation_index + i] * total_blength_2);
                            else
                                new_lh[i] *= (model->transposed_mut_mat[mutation_index + i] * total_blength_2);
                            
                            sum_lh += new_lh[i];
                        }
                        
                        // normalize new partial lh
                        RealNumType inverse_sum_lh = 1.0 / sum_lh;
                        for (StateType i = 0; i < num_states; ++i)
                            new_lh[i] *= inverse_sum_lh;
                        
                        StateType new_state = simplifyO(new_lh, aln->ref_seq[end_pos],  num_states, threshold_prob);

                        if (new_state == TYPE_O)
                            merged_regions->push_back(new SeqRegion(new_state, end_pos, 0, 0, new_lh));
                        else
                            merged_regions->push_back(new SeqRegion(new_state, end_pos));
                        
                        if (return_log_lh)
                            log_lh += log(sum_lh);
                    }
                    else
                    {
                        if (new_lh[seq2_state] == 0)
                        {
                            {
                                delete merged_regions;
                                merged_regions = NULL;
                                return -DBL_MAX;
                            }
                        }
                        
                        merged_regions->push_back(new SeqRegion(seq2_region->type, end_pos));
                    
                        if (return_log_lh)
                            log_lh += log(new_lh[seq2_state]);
                    }
                }
            }
            // seq1_entry = R/ACGT
            else
            {
                StateType seq1_state = seq1_region->type;
                if (seq1_state == TYPE_R)
                    seq1_state = aln->ref_seq[end_pos];
                
                RealNumType* new_lh = new RealNumType[num_states];
                RealNumType sum_lh = 0;
                
                if (total_blength_1 > 0)
                {
                    StateType mutation_index = model->row_index[seq1_state];
                    for (StateType i = 0; i < num_states; ++i)
                    {
                        if (seq1_state == i)
                            new_lh[i] = 1 + model->transposed_mut_mat[mutation_index + i] * total_blength_1;
                        else
                            new_lh[i] = model->transposed_mut_mat[mutation_index + i] * total_blength_1;
                    }
                }
                else
                {
                    for (StateType i = 0; i < num_states; ++i)
                        new_lh[i] = 0;
                    new_lh[seq1_state] = 1;
                }

                // seq1_entry = ACGT and seq2_entry = O
                if (seq2_region->type == TYPE_O)
                {
                    for (StateType i = 0; i < num_states; ++i)
                    {
                        RealNumType tot = 0;
                        
                        if (total_blength_2 > 0)
                        {
                            for (StateType j = 0; j < num_states; ++j)
                                tot += model->mutation_mat[model->row_index[i] + j] * seq2_region->likelihood[j];
                            
                            tot *= total_blength_2;
                        }
                        tot += seq2_region->likelihood[i];
                        new_lh[i] *= tot;
                        sum_lh += new_lh[i];
                    }
                    
                    if (sum_lh == 0)
                    {
                        delete merged_regions;
                        merged_regions = NULL;
                        return -DBL_MAX;
                    }
                        
                    // normalize new partial lh
                    RealNumType inverse_sum_lh = 1.0 / sum_lh;
                    for (StateType i = 0; i < num_states; ++i)
                        new_lh[i] *= inverse_sum_lh;
                    
                    StateType new_state = simplifyO(new_lh, aln->ref_seq[end_pos], num_states, threshold_prob);

                    if (new_state == TYPE_O)
                        merged_regions->push_back(new SeqRegion(new_state, end_pos, 0, 0, new_lh));
                    else
                        merged_regions->push_back(new SeqRegion(new_state, end_pos));
                    
                    if (return_log_lh)
                        log_lh += log(sum_lh);
                }
                // seq1_entry = ACGT and seq2_entry = R/ACGT
                else
                {
                    StateType seq2_state = seq2_region->type;
                    if (seq2_state == TYPE_R)
                        seq2_state = aln->ref_seq[end_pos];
                    
                    if (total_blength_2 > 0)
                    {
                        StateType mutation_index = model->row_index[seq2_state];
                        for (StateType i = 0; i < num_states; ++i)
                        {
                            if (seq2_state == i)
                                new_lh[i] *= (1 + model->transposed_mut_mat[mutation_index + i] * total_blength_2);
                            else
                                new_lh[i] *= (model->transposed_mut_mat[mutation_index + i] * total_blength_2);
                            
                            sum_lh += new_lh[i];
                        }
                        
                        // normalize new partial lh
                        RealNumType inverse_sum_lh = 1.0 / sum_lh;
                        for (StateType i = 0; i < num_states; ++i)
                            new_lh[i] *= inverse_sum_lh;
                        
                        StateType new_state = simplifyO(new_lh, aln->ref_seq[end_pos], num_states, threshold_prob);

                        if (new_state == TYPE_O)
                            merged_regions->push_back(new SeqRegion(new_state, end_pos, 0, 0, new_lh));
                        else
                            merged_regions->push_back(new SeqRegion(new_state, end_pos));
                        
                        if (return_log_lh)
                            log_lh += log(sum_lh);
                    }
                    else
                    {
                        merged_regions->push_back(new SeqRegion(seq2_region->type, end_pos));
                    
                        if (return_log_lh)
                            log_lh += log(new_lh[seq2_state]);
                    }
                }
            }
        }
        // update pos
        pos = end_pos + 1;
    }

    // try to merge consecutive and similar 'R' regions together
    merged_regions->mergeRegionR(num_states, threshold_prob);
    
    return log_lh;
}

RealNumType SeqRegions::computeAbsoluteLhAtRoot(Alignment* aln, Model* model, vector< vector<PositionType> > &cumulative_base)
{
    // dummy variables
    RealNumType log_lh = 0;
    RealNumType log_factor = 1;
    StateType num_states = aln->num_states;
    PositionType start_pos = 0;
    
    // browse regions one by one to compute the likelihood of each region
    for (PositionType region_index = 0; region_index < (PositionType) size(); ++region_index)
    {
        SeqRegion* region = at(region_index);
        
        // type R
        if (region->type == TYPE_R)
        {
            for (StateType i = 0; i < num_states; ++i)
                log_lh += model->root_log_freqs[i] * (cumulative_base[region->position + 1][i] - cumulative_base[start_pos][i]);
        }
        // type ACGT
        else if (region->type < num_states)
            log_lh += model->root_log_freqs[region->type];
        // type O
        else if (region->type == TYPE_O)
        {
            RealNumType tot = 0;
            for (StateType i = 0; i < num_states; ++i)
                tot += model->root_freqs[i] * region->likelihood[i];
            log_factor *= tot;
        }
        
        // maintain start_pos
        start_pos = region->position + 1;
    }

    // update log_lh
    log_lh += log(log_factor);
    
    // return the absolute likelihood
    return log_lh;
}

SeqRegions* SeqRegions::computeTotalLhAtRoot(StateType num_states, Model* model, RealNumType blength)
{
    SeqRegions* total_lh = new SeqRegions();
    
    for (SeqRegion* region : (*this))
    {
        // type N
        if (region->type == TYPE_N)
        {
            SeqRegion* new_region = new SeqRegion(region, num_states, false);
            total_lh->push_back(new_region);
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
                RealNumType* new_likelihood = new RealNumType[num_states];
                RealNumType sum_likelihood = 0;
                
                for (StateType i = 0; i < num_states; ++i)
                {
                    RealNumType tot = 0.0;
                    
                    if (total_blength > 0)
                    {
                        for (StateType j = 0; j < num_states; ++j)
                            tot += model->mutation_mat[model->row_index[i] + j] * region->likelihood[j];
                        tot *= total_blength;
                    }
                    
                    tot += region->likelihood[i];
                    new_likelihood[i] = tot * model->root_freqs[i];
                    sum_likelihood += new_likelihood[i];
                }
                
                // normalize likelihood
                sum_likelihood = 1.0 / sum_likelihood;
                for (StateType i = 0; i < num_states; ++i)
                    new_likelihood[i] *= sum_likelihood;
                
                // add new region to the total_lh_regions
                SeqRegion* new_region = new SeqRegion(region, num_states, false);
                new_region->likelihood = new_likelihood;
                total_lh->push_back(new_region);
            }
            // other types: R or A/C/G/T
            else
            {
                // add new region to the total_lh_regions
                SeqRegion* new_region = new SeqRegion(region, num_states, false);
                
                if (new_region->plength_observation2node >= 0)
                {
                    if (blength > 0)
                        new_region->plength_observation2node += blength;

                    new_region->plength_observation2root = 0;
                }
                else if (blength > 0)
                {
                    new_region->plength_observation2node = blength;
                    new_region->plength_observation2root = 0;
                }
                
                total_lh->push_back(new_region);
            }
        }
    }
    
    return total_lh;
}
