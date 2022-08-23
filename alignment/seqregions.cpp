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
                 for (PositionType j = 1; j < i - start_R_index; ++j)
                     delete at(start_R_index + j);
                 erase(begin() + start_R_index + 1, begin() + i);
                 
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
        for (PositionType j = 1; j < num_regions - 1 - start_R_index; ++j)
            delete at(start_R_index + j);
        erase(begin() + start_R_index + 1, end());
    }
}

void SeqRegions::move2NextRegion(SeqRegions* sequence, PositionType region_index, SeqRegion* &region, PositionType &current_pos, PositionType &end_pos, PositionType seq_length)
{
    // get the current region
    region = sequence->at(region_index);
    
    // get the current position and end position
    current_pos = region->position;
    end_pos = (region_index < (PositionType) sequence->size() - 1 ? sequence->at(region_index + 1)->position : seq_length) - 1;
}

void SeqRegions::getNextSharedSegment(PositionType current_pos, PositionType seq_length, SeqRegions* sequence1, SeqRegions* sequence2, PositionType &seq1_index, PositionType &seq2_index, SeqRegion* &seq1_region, SeqRegion* &seq2_region, PositionType &seq1_end_pos, PositionType &seq2_end_pos, PositionType &length)
{
    PositionType seq1_pos, seq2_pos, end_pos;
    
    // move to the next region in sequence 1
    if (current_pos > seq1_end_pos)
    {
        ++seq1_index;
        move2NextRegion(sequence1, seq1_index, seq1_region, seq1_pos, seq1_end_pos, seq_length);
    }
    
    // move to the next region in sequence 2
    if (current_pos > seq2_end_pos)
    {
        ++seq2_index;
        move2NextRegion(sequence2, seq2_index, seq2_region, seq2_pos, seq2_end_pos, seq_length);
    }
    
    // compute the end_pos for the shared segment
    end_pos = seq1_end_pos < seq2_end_pos ? seq1_end_pos : seq2_end_pos;
    length = end_pos + 1 - current_pos;
}

int SeqRegions::compareWithSample(SeqRegions* sequence2, PositionType seq_length, StateType num_states)
{
    ASSERT(seq_length > 0);
    
    // init dummy variables
    PositionType seq1_index = -1;
    PositionType seq2_index = -1;
    bool seq1_more_info = false;
    bool seq2_more_info = false;
    PositionType pos = 0;
    SeqRegion *seq1_region, *seq2_region;
    PositionType seq1_end = -1, seq2_end = -1;
    PositionType length;

    while (pos < seq_length && (!seq1_more_info || !seq2_more_info))
    {
        // get the next shared segment in the two sequences
        getNextSharedSegment(pos, seq_length, this, sequence2, seq1_index, seq2_index, seq1_region, seq2_region, seq1_end, seq2_end, length);
        
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
        pos += length;
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
    PositionType seq1_index = -1;
    PositionType seq2_index = -1;
    PositionType pos = 0;
    SeqRegion *seq1_region, *seq2_region;
    PositionType seq1_end = -1, seq2_end = -1;
    PositionType length;
    const RealNumType THRESH_PROB = params->threshold_prob;
    const RealNumType THRESH_DIFF_UPDATE = params->thresh_diff_update;
    
    // compare each pair of regions
    while (pos < seq_length)
    {
        // get the next shared segment in the two sequences
        SeqRegions::getNextSharedSegment(pos, seq_length, this, regions2, seq1_index, seq2_index, seq1_region, seq2_region, seq1_end, seq2_end, length);
        
        // if any pair of regions is different in type -> the two sequences are different
        if (seq1_region->type != seq2_region->type)
            return true;
        
        // type R/ACGT
        if (seq1_region->type < num_states || seq1_region->type == TYPE_R)
        {
            // compare plength_from_root and plength_observation
            if (fabs(seq1_region->plength_observation2root - seq2_region->plength_observation2root) > THRESH_PROB
                ||fabs(seq1_region->plength_observation2node - seq2_region->plength_observation2node) > THRESH_PROB)
                return true;
        }
        
        // type O
        if (seq1_region->type == TYPE_O)
        {
            // compare plength_observation
            if (fabs(seq1_region->plength_observation2node - seq2_region->plength_observation2node) > THRESH_PROB)
                return true;
            
            // compare likelihood of each state
            for (StateType i = 0; i < num_states; ++i)
            {
                RealNumType diff = fabs(seq1_region->likelihood[i] - seq2_region->likelihood[i]);
                
                if (diff > 0)
                {
                    if ((seq1_region->likelihood[i] == 0) || (seq2_region->likelihood[i] == 0))
                        return true;
                    
                    if (diff > THRESH_DIFF_UPDATE
                        || (diff > THRESH_PROB
                            && ((diff / seq1_region->likelihood[i] > THRESH_DIFF_UPDATE)
                                || (diff / seq2_region->likelihood[i] > THRESH_DIFF_UPDATE))))
                        return true;
                }
            }
        }
        
        // update pos
        pos += length;
    }
    
    return false;
}

void SeqRegions::mergeUpperLower(SeqRegions* &merged_regions, RealNumType upper_plength, SeqRegions* lower_regions, RealNumType lower_plength, Alignment* aln, Model* model, RealNumType threshold_prob)
{
    // init variables
    PositionType seq1_index = -1;
    PositionType seq2_index = -1;
    PositionType pos = 0;
    SeqRegion *seq1_region, *seq2_region;
    PositionType seq1_end = -1, seq2_end = -1;
    PositionType length;
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
        SeqRegions::getNextSharedSegment(pos, seq_length, this, lower_regions, seq1_index, seq2_index, seq1_region, seq2_region, seq1_end, seq2_end, length);
        
        // seq1_entry = 'N'
        if (seq1_region->type == TYPE_N)
        {
            // seq1_entry = 'N' and seq2_entry = 'N'
            if (seq2_region->type == TYPE_N)
                merged_regions->push_back(new SeqRegion(seq1_region->type, pos));
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
                    
                    StateType start_index = 0;
                    for (StateType i = 0; i < num_states; ++i)
                    {
                        RealNumType tot = 0;
                        
                        if (total_blength > 0)
                        {
                            for (StateType j = 0; j < num_states; ++j)
                                tot += model->mutation_mat[start_index + j] * seq2_region->likelihood[j];
                            
                            tot *= total_blength;
                        }
                        
                        tot += seq2_region->likelihood[i];
                        new_lh[i] = tot * model->root_freqs[i];
                        sum_lh += new_lh[i];
                        
                        // update start_index
                        start_index += num_states;
                    }

                    // normalize the new partial likelihood
                    RealNumType inverse_sum_lh = 1.0 / sum_lh;
                    for (StateType i = 0; i < num_states; ++i)
                        new_lh[i] *= inverse_sum_lh;
                    
                    // add merged region into merged_regions
                    merged_regions->push_back(new SeqRegion(TYPE_O, pos, 0, 0, new_lh));
                }
                // seq1_entry = 'N' and seq2_entry = R/ACGT
                else
                {
                    if (seq2_region->plength_observation2node >= 0)
                    {
                        RealNumType new_plength = seq2_region->plength_observation2node ;
                        if (lower_plength > 0)
                            new_plength += lower_plength;
                        merged_regions->push_back(new SeqRegion(seq2_region->type, pos, new_plength, 0));
                    }
                    else
                    {
                        if (lower_plength > 0)
                            merged_regions->push_back(new SeqRegion(seq2_region->type, pos, lower_plength, 0));
                        else
                            merged_regions->push_back(new SeqRegion(seq2_region->type, pos));
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
                    
                    StateType start_index = 0;
                    for (StateType i = 0; i < num_states; ++i)
                    {
                        RealNumType tot = 0;
                        
                        for (StateType j = 0; j < num_states; ++j)
                            tot += model->transposed_mut_mat[start_index + j] * seq1_region->likelihood[j];
                        
                        tot *= total_blength;
                        tot += seq1_region->likelihood[i];
                        new_lh[i] = tot;
                        sum_lh += new_lh[i];
                        
                        // update start_index
                        start_index += num_states;
                    }

                    // normalize the new partial likelihood
                    RealNumType inverse_sum_lh = 1.0 / sum_lh;
                    for (StateType i = 0; i < num_states; ++i)
                        new_lh[i] *= inverse_sum_lh;
                    
                    // add merged region into merged_regions
                    merged_regions->push_back(new SeqRegion(TYPE_O, pos, 0, 0, new_lh));
                }
                else
                {
                    RealNumType* new_lh = new RealNumType[num_states];
                    memcpy(new_lh, seq1_region->likelihood, sizeof(RealNumType) * num_states);
                    // add merged region into merged_regions
                    merged_regions->push_back(new SeqRegion(seq1_region->type, pos, 0, 0, new_lh));
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
                    merged_regions->push_back(new SeqRegion(seq1_region->type, pos, seq1_region->plength_observation2node, plength_from_root));
                }
                else if (seq1_region->plength_observation2node >= 0)
                {
                    RealNumType plength_observation = seq1_region->plength_observation2node;
                    if (upper_plength > 0)
                        plength_observation += upper_plength;
                    merged_regions->push_back(new SeqRegion(seq1_region->type, pos, plength_observation));
                }
                else
                {
                    if (upper_plength > 0)
                        merged_regions->push_back(new SeqRegion(seq1_region->type, pos, upper_plength));
                    else
                        merged_regions->push_back(new SeqRegion(seq1_region->type, pos));
                }
            }
        }
        // seq1_entry = seq2_entry = R/ACGT
        else if (seq1_region->type == seq2_region->type && (seq1_region->type < num_states || seq1_region->type == TYPE_R))
            merged_regions->push_back(new SeqRegion(seq1_region->type, pos));
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
                
                merged_regions->push_back(new SeqRegion(seq2_region->type, pos));
            }
            // due to 0 distance, the entry will be of same type as entry1
            else if ((seq1_region->type < num_states || seq1_region->type == TYPE_R) && total_blength_1 <= 0)
            {
                merged_regions->push_back(new SeqRegion(seq1_region->type, pos));
            }
            // seq1_entry = O
            else if (seq1_region->type == TYPE_O)
            {
                RealNumType* new_lh = new RealNumType[num_states];
                
                // if total_blength_1 > 0 => compute new partial likelihood
                if (total_blength_1 > 0)
                {
                    StateType start_index = 0;
                    for (StateType i = 0; i < num_states; ++i)
                    {
                        RealNumType tot = 0;
                        
                        for (StateType j = 0; j < num_states; ++j)
                            tot += model->transposed_mut_mat[start_index + j] * seq1_region->likelihood[j];
                        
                        tot *= total_blength_1;
                        tot += seq1_region->likelihood[i];
                        new_lh[i] = tot;
                        
                        // update start_index
                        start_index += num_states;
                    }
                }
                // otherwise, clone the partial likelihood from seq1
                else
                    memcpy(new_lh, seq1_region->likelihood, sizeof(RealNumType) * num_states);
                
                RealNumType sum_new_lh = 0;

                // seq2 = O
                if (seq2_region->type == TYPE_O)
                {
                    StateType start_index = 0;
                    for (StateType i = 0; i < num_states; ++i)
                    {
                        RealNumType tot = 0;
                        
                        if (total_blength_2 > 0)
                        {
                            for (StateType j = 0; j < num_states; ++j)
                                tot += model->mutation_mat[start_index + j] * seq2_region->likelihood[j];
        
                            tot *= total_blength_2;
                        }
                        tot += seq2_region->likelihood[i];
                        new_lh[i] *= tot;
                        sum_new_lh += new_lh[i];
                        
                        // update start_index
                        start_index += num_states;
                    }
                }
                // seq1 = "O" and seq2 = ACGT
                else
                {
                    StateType seq2_state = seq2_region->type;
                    if (seq2_state == TYPE_R)
                        seq2_state = aln->ref_seq[pos];
                    
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
                
                StateType new_state = simplifyO(new_lh, aln->ref_seq[pos], num_states, threshold_prob);

                if (new_state == TYPE_O)
                    merged_regions->push_back(new SeqRegion(TYPE_O, pos, 0, 0, new_lh));
                else
                {
                    delete[] new_lh;
                    merged_regions->push_back(new SeqRegion(new_state, pos));
                }
            }
            // seq1_entry = R/ACGT
            else
            {
                StateType seq1_state = seq1_region->type;
                if (seq1_state == TYPE_R)
                    seq1_state = aln->ref_seq[pos];
                
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
                    
                    StateType start_index = 0;
                    for (StateType i = 0; i < num_states; ++i)
                    {
                        RealNumType tot = 0;
                        
                        for (StateType j = 0; j < num_states; ++j)
                            tot += model->transposed_mut_mat[start_index + j] * root_vec[j];
                        
                        tot *= length_to_root;
                        tot += root_vec[i];
                        new_lh[i] = tot;
                        
                        // update start_index
                        start_index += num_states;
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
                        {
                            if (i == seq1_state)
                                new_lh[i] = 1;
                            else
                                new_lh[i] = 0;
                        }
                    }
                }
                  
                // seq2 = "O" and seq1 = ACGT
                if (seq2_region->type == TYPE_O)
                {
                    StateType start_index = 0;
                    for (StateType i = 0; i < num_states; ++i)
                    {
                        RealNumType tot = 0;
                        
                        if (total_blength_2 > 0)
                        {
                            for (StateType j = 0; j < num_states; ++j)
                                tot += model->mutation_mat[start_index + j] * seq2_region->likelihood[j];
                            
                            tot *= total_blength_2;
                        }
                        
                        tot += seq2_region->likelihood[i];
                        new_lh[i] *= tot;
                        sum_new_lh += new_lh[i];
                        
                        // update start_index
                        start_index += num_states;
                    }
                    
                    // normalize the new partial lh
                    RealNumType inverse_sum_lh = 1.0 / sum_new_lh;
                    for (StateType i = 0; i < num_states; ++i)
                        new_lh[i] *= inverse_sum_lh;
                    
                    StateType new_state = simplifyO(new_lh, aln->ref_seq[pos], num_states, threshold_prob);
                    
                    if (new_state == TYPE_O)
                        merged_regions->push_back(new SeqRegion(TYPE_O, pos, 0, 0, new_lh));
                    else
                        merged_regions->push_back(new SeqRegion(new_state, pos));
                }
                // seq1 = ACGT and different from seq2 = R/ACGT
                else
                {
                    StateType seq2_state = seq2_region->type;
                    
                    if (seq2_state == TYPE_R)
                        seq2_state = aln->ref_seq[pos];
                    
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
                    merged_regions->push_back(new SeqRegion(TYPE_O, pos, 0, 0, new_lh));
                }
            }
        }

        // update pos
        pos += length;
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
    PositionType seq1_index = -1;
    PositionType seq2_index = -1;
    PositionType pos = 0;
    StateType num_states = aln->num_states;
    SeqRegion *seq1_region, *seq2_region;
    PositionType seq1_end = -1, seq2_end = -1;
    PositionType length;
    PositionType seq_length = aln->ref_seq.size();
    
    // init merged_regions
    if (merged_regions)
        merged_regions->deleteRegions();
    else
        merged_regions = new SeqRegions();
                
    while (pos < seq_length)
    {
        // get the next shared segment in the two sequences
        SeqRegions::getNextSharedSegment(pos, seq_length, this, regions2, seq1_index, seq2_index, seq1_region, seq2_region, seq1_end, seq2_end, length);
        
        // seq1_entry = 'N'
        if (seq1_region->type == TYPE_N)
        {
            // seq1_entry = 'N' and seq2_entry = 'N'
            if (seq2_region->type == TYPE_N)
                merged_regions->push_back(new SeqRegion(seq1_region->type, pos));
            // seq1_entry = 'N' and seq2_entry = O/R/ACGT
            else
            {
                // seq1_entry = 'N' and seq2_entry = O
                if (seq2_region->type == TYPE_O)
                {
                    // add merged region into merged_regions
                    SeqRegion* new_region = new SeqRegion(seq2_region, num_states);
                    new_region->position = pos;
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
                        merged_regions->push_back(new SeqRegion(seq2_region->type, pos, new_plength));
                    }
                    else
                    {
                        if (plength2 > 0)
                            merged_regions->push_back(new SeqRegion(seq2_region->type, pos, plength2));
                        else
                            merged_regions->push_back(new SeqRegion(seq2_region->type, pos));
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
                new_region->position = pos;
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
                    merged_regions->push_back(new SeqRegion(seq1_region->type, pos, new_plength));
                }
                else
                {
                    if (plength1 > 0)
                        merged_regions->push_back(new SeqRegion(seq1_region->type, pos, plength1));
                    else
                        merged_regions->push_back(new SeqRegion(seq1_region->type, pos));
                        
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
                merged_regions->push_back(new SeqRegion(seq1_region->type, pos));
                
                if (return_log_lh)
                {
                    // convert total_blength_1 and total_blength_2 to zero if they are -1
                    if (total_blength_1 < 0) total_blength_1 = 0;
                    if (total_blength_2 < 0) total_blength_2 = 0;
                    
                    if (seq1_region->type == TYPE_R)
                        log_lh += (total_blength_1 + total_blength_2) * (cumulative_rate[pos + length] - cumulative_rate[pos]);
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
                    StateType start_index = 0;
                    for (StateType i = 0; i < num_states; ++i)
                    {
                        RealNumType tot = 0;
                        
                        for (StateType j = 0; j < num_states; ++j)
                            tot += model->mutation_mat[start_index + j] * seq1_region->likelihood[j];
                        
                        tot *= total_blength_1;
                        tot += seq1_region->likelihood[i];
                        new_lh[i] = tot;
                        
                        // update start_index
                        start_index += num_states;
                    }
                }
                // otherwise, clone the partial likelihood from seq1
                else
                    memcpy(new_lh, seq1_region->likelihood, sizeof(RealNumType) * num_states);

                // seq1_entry = O and seq2_entry = O
                if (seq2_region->type == TYPE_O)
                {
                    StateType start_index = 0;
                    for (StateType i = 0; i < num_states; ++i)
                    {
                        RealNumType tot = 0;
                        
                        if (total_blength_2 > 0)
                        {
                            for (StateType j = 0; j < num_states; ++j)
                                tot += model->mutation_mat[start_index + j] * seq2_region->likelihood[j];
                            
                            tot *= total_blength_2;
                        }
                        tot += seq2_region->likelihood[i];
                        new_lh[i] *= tot;
                        sum_lh += new_lh[i];
                        
                        // update start_index
                        start_index += num_states;
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
                    
                    StateType new_state = simplifyO(new_lh, aln->ref_seq[pos], num_states, threshold_prob);

                    if (new_state == TYPE_O)
                        merged_regions->push_back(new SeqRegion(new_state, pos, 0, 0, new_lh));
                    else
                        merged_regions->push_back(new SeqRegion(new_state, pos));
                    
                    if (return_log_lh)
                        log_lh += log(sum_lh);
                }
                // seq1_entry = O and seq2_entry = ACGT
                else
                {
                    StateType seq2_state = seq2_region->type;
                    if (seq2_state == TYPE_R)
                        seq2_state = aln->ref_seq[pos];
                    
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
                        
                        StateType new_state = simplifyO(new_lh, aln->ref_seq[pos],  num_states, threshold_prob);

                        if (new_state == TYPE_O)
                            merged_regions->push_back(new SeqRegion(new_state, pos, 0, 0, new_lh));
                        else
                            merged_regions->push_back(new SeqRegion(new_state, pos));
                        
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
                        
                        merged_regions->push_back(new SeqRegion(seq2_region->type, pos));
                    
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
                    seq1_state = aln->ref_seq[pos];
                
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
                    {
                        if (seq1_state == i)
                            new_lh[i] = 1;
                        else
                            new_lh[i] = 0;
                    }
                }

                // seq1_entry = ACGT and seq2_entry = O
                if (seq2_region->type == TYPE_O)
                {
                    StateType start_index = 0;
                    for (StateType i = 0; i < num_states; ++i)
                    {
                        RealNumType tot = 0;
                        
                        if (total_blength_2 > 0)
                        {
                            for (StateType j = 0; j < num_states; ++j)
                                tot += model->mutation_mat[start_index + j] * seq2_region->likelihood[j];
                            
                            tot *= total_blength_2;
                        }
                        tot += seq2_region->likelihood[i];
                        new_lh[i] *= tot;
                        sum_lh += new_lh[i];
                        
                        // update start_index
                        start_index += num_states;
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
                    
                    StateType new_state = simplifyO(new_lh, aln->ref_seq[pos], num_states, threshold_prob);

                    if (new_state == TYPE_O)
                        merged_regions->push_back(new SeqRegion(new_state, pos, 0, 0, new_lh));
                    else
                        merged_regions->push_back(new SeqRegion(new_state, pos));
                    
                    if (return_log_lh)
                        log_lh += log(sum_lh);
                }
                // seq1_entry = ACGT and seq2_entry = R/ACGT
                else
                {
                    StateType seq2_state = seq2_region->type;
                    if (seq2_state == TYPE_R)
                        seq2_state = aln->ref_seq[pos];
                    
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
                        
                        StateType new_state = simplifyO(new_lh, aln->ref_seq[pos], num_states, threshold_prob);

                        if (new_state == TYPE_O)
                            merged_regions->push_back(new SeqRegion(new_state, pos, 0, 0, new_lh));
                        else
                            merged_regions->push_back(new SeqRegion(new_state, pos));
                        
                        if (return_log_lh)
                            log_lh += log(sum_lh);
                    }
                    else
                    {
                        merged_regions->push_back(new SeqRegion(seq2_region->type, pos));
                    
                        if (return_log_lh)
                            log_lh += log(new_lh[seq2_state]);
                    }
                }
            }
        }
        // update pos
        pos += length;
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
    PositionType seq_length = aln->ref_seq.size();
    
    // browse regions one by one to compute the likelihood of each region
    for (PositionType region_index = 0; region_index < (PositionType) size(); ++region_index)
    {
        SeqRegion* region = at(region_index);
        
        // type R
        if (region->type == TYPE_R)
        {
            PositionType start_pos = region->position;
            PositionType end_pos = (region_index == ((PositionType) size()) - 1 ? seq_length - 1 : at(region_index + 1)->position - 1);
            
            for (StateType i = 0; i < num_states; ++i)
                log_lh += model->root_log_freqs[i] * (cumulative_base[end_pos + 1][i] - cumulative_base[start_pos][i]);
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
                
                StateType start_index = 0;
                for (StateType i = 0; i < num_states; ++i)
                {
                    RealNumType tot = 0.0;
                    
                    if (total_blength > 0)
                    {
                        for (StateType j = 0; j < num_states; ++j)
                            tot += model->mutation_mat[start_index + j] * region->likelihood[j];
                        tot *= total_blength;
                    }
                    
                    tot += region->likelihood[i];
                    new_likelihood[i] = tot * model->root_freqs[i];
                    sum_likelihood += new_likelihood[i];
                    
                    // update start_index
                    start_index += num_states;
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

// this implementation derives from appendProbNode
RealNumType SeqRegions::calculateSubTreePlacementCost(Alignment* aln, Model* model, RealNumType* cumulative_rate, SeqRegions* child_regions, RealNumType blength)
{
    // init dummy variables
    RealNumType lh_cost = 0;
    PositionType seq1_index = -1;
    PositionType seq2_index = -1;
    PositionType pos = 0;
    RealNumType total_factor = 1;
    StateType num_states = aln->num_states;
    SeqRegion *seq1_region, *seq2_region;
    PositionType seq1_end = -1, seq2_end = -1;
    PositionType length;
    RealNumType minimum_carry_over = DBL_MIN * 1e50;
    RealNumType total_blength = blength;
    PositionType seq_length = aln->ref_seq.size();
    
    while (pos < seq_length)
    {
        // get the next shared segment in the two sequences
        SeqRegions::getNextSharedSegment(pos, seq_length, this, child_regions, seq1_index, seq2_index, seq1_region, seq2_region, seq1_end, seq2_end, length);
        
        // 1. e1.type = N || e2.type = N
        if (seq2_region->type == TYPE_N || seq1_region->type == TYPE_N)
        {
            pos += length;
            continue;
        }
        // e1.type != N && e2.type != N
        else
        {
            // total_blength will be here the total length from the root or from the upper node, down to the down node.
            if (seq1_region->plength_observation2root >= 0)
                total_blength = seq1_region->plength_observation2root + (blength >= 0 ? blength : 0);
            else if (seq1_region->plength_observation2node >= 0)
                total_blength = seq1_region->plength_observation2node + (blength >= 0 ? blength : 0);
            else
                total_blength = blength;
                
            if (seq2_region->plength_observation2node >= 0)
                total_blength = (total_blength > 0 ? total_blength : 0) + seq2_region->plength_observation2node;
            
            // 2. e1.type = R
            if (seq1_region->type == TYPE_R)
            {
                // 2.1. e1.type = R and e2.type = R
                if (seq2_region->type == TYPE_R)
                {
                    if (seq1_region->plength_observation2root >= 0)
                        total_blength += seq1_region->plength_observation2node;
                    
                    if (total_blength > 0)
                        lh_cost += total_blength * (cumulative_rate[pos + length] - cumulative_rate[pos]);
                }
                // 2.2. e1.type = R and e2.type = O
                else if (seq2_region->type == TYPE_O)
                {
                    RealNumType tot = 0;
                    StateType seq1_state = aln->ref_seq[pos];
                    
                    if (seq1_region->plength_observation2root >= 0)
                    {
                        StateType mutation_index = model->row_index[seq1_state];
                        StateType start_index = 0;
                        for (StateType i = 0; i < num_states; ++i)
                        {
                            RealNumType tot2;
                            
                            if (seq1_state == i)
                                tot2 = model->root_freqs[i] * (1.0 + model->transposed_mut_mat[mutation_index + i] * seq1_region->plength_observation2node);
                            else
                                tot2 = model->root_freqs[i] * (model->transposed_mut_mat[mutation_index + i] * seq1_region->plength_observation2node);
                                
                            RealNumType tot3 = 0;
                            if (total_blength > 0)
                            {
                                for (StateType j = 0; j < num_states; ++j)
                                    tot3 += model->mutation_mat[start_index + j] * seq2_region->likelihood[j];
                            }
                            
                            tot += tot2 * (seq2_region->likelihood[i] + total_blength * tot3);
                            
                            // update start_index
                            start_index += num_states;
                        }
                        
                        tot *= model->inverse_root_freqs[seq1_state];
                    }
                    else
                    {
                        if (total_blength > 0)
                        {
                            StateType mutation_index = model->row_index[seq1_state];
                            for (StateType j = 0; j < num_states; ++j)
                                tot += model->mutation_mat[mutation_index + j] * seq2_region->likelihood[j];
                            
                            tot *= total_blength;
                        }
                        
                        tot += seq2_region->likelihood[seq1_state];
                    }
                        
                    total_factor *= tot;
                }
                // 2.3. e1.type = R and e2.type = A/C/G/T
                else
                {
                    StateType seq1_state = aln->ref_seq[pos];
                    StateType seq2_state = seq2_region->type;
                    
                    if (seq1_region->plength_observation2root >= 0)
                    {
                        if (total_blength > 0)
                        {
                            RealNumType seq1_state_evolves_seq2_state = model->mutation_mat[model->row_index[seq1_state] + seq2_state] * total_blength * (1.0 + model->diagonal_mut_mat[seq1_state] * seq1_region->plength_observation2node);
                            
                            RealNumType seq2_state_evolves_seq1_state = model->freqi_freqj_qij[model->row_index[seq2_state] + seq1_state] * seq1_region->plength_observation2node * (1.0 + model->diagonal_mut_mat[seq2_state] * total_blength);
                            
                            total_factor *= seq1_state_evolves_seq2_state + seq2_state_evolves_seq1_state;
                        }
                        else
                        {
                            total_factor *= model->freqi_freqj_qij[model->row_index[seq2_state] + seq1_state] * seq1_region->plength_observation2node;
                        }
                    }
                    else
                    {
                        if (total_blength > 0)
                            total_factor *= model->mutation_mat[model->row_index[seq1_state] + seq2_state] * total_blength;
                        else
                            return -DBL_MAX;
                    }
                }
            }
            // 3. e1.type = O
            else if (seq1_region->type == TYPE_O)
            {
                // 3.1. e1.type = O and e2.type = O
                if (seq2_region->type == TYPE_O)
                {
                    RealNumType tot = 0;
                    
                    if (total_blength > 0)
                    {
                        StateType start_index = 0;
                        for (StateType i = 0; i < num_states; ++i)
                        {
                            RealNumType tot2 = 0;
                            
                            for (StateType j = 0; j < num_states; ++j)
                                tot2 += model->mutation_mat[start_index + j] * seq2_region->likelihood[j];
                            
                            tot += seq1_region->likelihood[i] * (seq2_region->likelihood[i] + total_blength * tot2);
                            
                            // update start_index
                            start_index += num_states;
                        }
                    }
                    else
                    {
                        for (StateType i = 0; i < num_states; ++i)
                            tot += seq1_region->likelihood[i] * seq2_region->likelihood[i];
                    }
                    
                    total_factor *= tot;
                }
                // 3.2. e1.type = O and e2.type = R or A/C/G/T
                else
                {
                    StateType seq2_state = seq2_region->type;
                    if (seq2_state == TYPE_R)
                        seq2_state = aln->ref_seq[pos];
                    
                    if (total_blength > 0)
                    {
                        RealNumType tot2 = 0;
                        StateType mutation_index = model->row_index[seq2_state];
                        for (StateType j = 0; j < num_states; ++j)
                            tot2 += seq1_region->likelihood[j] * model->transposed_mut_mat[mutation_index + j];
                        
                        total_factor *= seq1_region->likelihood[seq2_state] + total_blength * tot2;
                    }
                    else
                        total_factor *= seq1_region->likelihood[seq2_state];
                }
            }
            // 4. e1.type = A/C/G/T
            else
            {
                int seq1_state = seq1_region->type;
                
                // 4.1. e1.type =  e2.type
                if (seq1_region->type == seq2_region->type)
                {
                    if (seq1_region->plength_observation2root >= 0)
                        total_blength += seq1_region->plength_observation2node;
                    
                    if (total_blength > 0)
                        lh_cost += model->diagonal_mut_mat[seq1_state] * total_blength;
                }
                // e1.type = A/C/G/T and e2.type = O/A/C/G/T
                else
                {
                    // 4.2. e1.type = A/C/G/T and e2.type = O
                    if (seq2_region->type == TYPE_O)
                    {
                        RealNumType tot = 0.0;
                        
                        if (seq1_region->plength_observation2root >= 0)
                        {
                            StateType mutation_index = model->row_index[seq1_state];
                            StateType start_index = 0;
                            for (StateType i = 0; i < num_states; ++i)
                            {
                                RealNumType tot2;
                                
                                if (seq1_state == i)
                                    tot2 = model->root_freqs[i] * (1.0 + model->transposed_mut_mat[mutation_index + i] * seq1_region->plength_observation2node);
                                else
                                    tot2 = model->root_freqs[i] * (model->transposed_mut_mat[mutation_index + i] * seq1_region->plength_observation2node);
                                    
                                RealNumType tot3 = 0;
                                for (StateType j = 0; j < num_states; ++j)
                                    tot3 += model->mutation_mat[start_index + j] * seq2_region->likelihood[j];
                                tot += tot2 * (seq2_region->likelihood[i] + total_blength * tot3);
                                
                                // update start_index
                                start_index += num_states;
                            }
                            
                            total_factor *= (tot * model->inverse_root_freqs[seq1_state]);
                        }
                        else
                        {
                            StateType mutation_index = model->row_index[seq1_state];
                            for (StateType j = 0; j < num_states; ++j)
                                tot += model->mutation_mat[mutation_index + j] * seq2_region->likelihood[j];
                            
                            tot *= total_blength;
                            tot += seq2_region->likelihood[seq1_state];
                            total_factor *= tot;
                        }
                    }
                    // 4.3. e1.type = A/C/G/T and e2.type = R or A/C/G/T
                    else
                    {
                        StateType seq2_state = seq2_region->type;
                        if (seq2_state == TYPE_R)
                            seq2_state = aln->ref_seq[pos];
                        
                        if (seq1_region->plength_observation2root >= 0)
                        {
                            if (total_blength > 0)
                            {
                                RealNumType seq1_state_evolves_seq2_state = model->mutation_mat[model->row_index[seq1_state] + seq2_state] * total_blength * (1.0 + model->diagonal_mut_mat[seq1_state] * seq1_region->plength_observation2node);
                                
                                RealNumType seq2_state_evolves_seq1_state = model->freqi_freqj_qij[model->row_index[seq2_state] + seq1_state] * seq1_region->plength_observation2node * (1.0 + model->diagonal_mut_mat[seq2_state] * total_blength);
                                
                                total_factor *= seq1_state_evolves_seq2_state + seq2_state_evolves_seq1_state;
                            }
                            else
                            {
                                total_factor *= model->freqi_freqj_qij[model->row_index[seq2_state] + seq1_state] * seq1_region->plength_observation2node;
                            }
                        }
                        else
                        {
                            if (total_blength > 0)
                                total_factor *= model->mutation_mat[model->row_index[seq1_state] + seq2_state] * total_blength;
                            else
                                return -DBL_MAX;
                        }
                    }
                }
            }
        }
         
        // approximately update lh_cost and total_factor
        if (total_factor <= minimum_carry_over)
        {
            if (total_factor < DBL_MIN)
                return -DBL_MAX;
            lh_cost += log(total_factor);
            total_factor = 1.0;
        }
        
        // update pos
        pos += length;
    }
    
    return lh_cost + log(total_factor);
}

#define MIN_CARRY_OVER 1e-250
// this implementation derives from appendProb
RealNumType SeqRegions::calculateSamplePlacementCost(Alignment* aln, Model* model, RealNumType* cumulative_rate, SeqRegions* child_regions, RealNumType blength)
{
    // init dummy variables
    RealNumType lh_cost = 0;
    PositionType seq1_index = -1;
    PositionType seq2_index = -1;
    PositionType pos = 0;
    RealNumType total_factor = 1;
    StateType num_states = aln->num_states;
    SeqRegion *seq1_region, *seq2_region;
    PositionType seq1_end = -1, seq2_end = -1;
    PositionType length;
    //RealNumType minimum_carry_over = DBL_MIN * 1e50;
    if (blength < 0) blength = 0;
    RealNumType total_blength = blength;
    PositionType seq_length = aln->ref_seq.size();
    
    while (pos < seq_length)
    {
        // get the next shared segment in the two sequences
        SeqRegions::getNextSharedSegment(pos, seq_length, this, child_regions, seq1_index, seq2_index, seq1_region, seq2_region, seq1_end, seq2_end, length);
        
        // 1. e1.type = N || e2.type = N
        if (seq2_region->type == TYPE_N || seq1_region->type == TYPE_N)
        {
            pos += length;
            continue;
        }
        // e1.type != N && e2.type != N
        else {
            // A,C,G,T
            // R -> same as the reference
            // N -> gaps
            // O -> vector of 4 probability to observe A, C, G, T
            // 2. e1.type = R
            if (seq1_region->type == TYPE_R)
            {
                // 2.1. e1.type = R and e2.type = R
                if (seq2_region->type == TYPE_R)
                {
                    if (seq1_region->plength_observation2node < 0 && seq1_region->plength_observation2root < 0)
                        lh_cost += blength * (cumulative_rate[pos + length] - cumulative_rate[pos]);
                    else
                    {
                        total_blength = blength + seq1_region->plength_observation2node;
                        if (seq1_region->plength_observation2root < 0)
                            lh_cost += total_blength * (cumulative_rate[pos + length] - cumulative_rate[pos]);
                        else
                            // here contribution from root frequency gets added and subtracted so it's ignored
                            lh_cost += (total_blength + seq1_region->plength_observation2root) * (cumulative_rate[pos + length] - cumulative_rate[pos]);
                    }
                }
                // 2.2. e1.type = R and e2.type = O
                else if (seq2_region->type == TYPE_O)
                {
                    StateType seq1_state = aln->ref_seq[pos];
                    if (seq1_region->plength_observation2root >= 0)
                    {
                        total_blength = seq1_region->plength_observation2root + blength;
                        
                        if (seq2_region->likelihood[seq1_state] > 0.1)
                        {
                            total_blength += seq1_region->plength_observation2node;
                            
                            // here contribution from root frequency can also be also ignored
                            lh_cost += model->diagonal_mut_mat[seq1_state] * total_blength;
                        }
                        else
                        {
                            RealNumType tot = 0;
                            StateType mutation_index = model->row_index[seq1_state];
                            StateType start_index = 0;
                            for (StateType i = 0; i < num_states; ++i)
                            {
                                RealNumType tot2;
                                
                                if (seq1_state == i)
                                    tot2 = model->root_freqs[i] * (1.0 + model->transposed_mut_mat[mutation_index + i] * seq1_region->plength_observation2node);
                                else
                                    tot2 = model->root_freqs[i] * (model->transposed_mut_mat[mutation_index + i] * seq1_region->plength_observation2node);
                                    
                                RealNumType tot3 = 0;
                                if (seq2_region->likelihood[i] > 0.1)
                                    tot3 = 1;
                                
                                for (StateType j = 0; j < num_states; ++j)
                                    if (seq2_region->likelihood[j] > 0.1)
                                        tot3 += model->mutation_mat[start_index + j]; // TODO
                                tot3 *= total_blength;
                                
                                tot += tot2 * tot3;
                                
                                // update start_index
                                start_index += num_states;
                            }
                            
                            total_factor *= tot * model->inverse_root_freqs[seq1_state];
                        }
                    }
                    else
                    {
                        if (seq2_region->likelihood[seq1_state] > 0.1)
                        {
                            if (seq1_region->plength_observation2node >= 0)
                                lh_cost += model->diagonal_mut_mat[seq1_state] * (blength + seq1_region->plength_observation2node);
                            else
                                lh_cost += model->diagonal_mut_mat[seq1_state] * blength;
                        }
                        else
                        {
                            RealNumType tot = 0;
                            StateType mutation_index = model->row_index[seq1_state];
                            for (StateType i = 0; i < num_states; ++i)
                                if (seq2_region->likelihood[i] > 0.1)
                                    tot += model->mutation_mat[mutation_index + i];
                            
                            if (seq1_region->plength_observation2node >= 0)
                                total_factor *= tot * (blength + seq1_region->plength_observation2node);
                            else
                                total_factor *= tot * blength;
                        }
                    }
                }
                // 2.3. e1.type = R and e2.type = A/C/G/T
                else
                {
                    StateType seq1_state = aln->ref_seq[pos];
                    StateType seq2_state = seq2_region->type;
                    
                    if (seq1_region->plength_observation2root >= 0)
                    {
                        RealNumType seq1_state_evolves_seq2_state = model->mutation_mat[model->row_index[seq1_state] + seq2_state] * blength * (1.0 + model->diagonal_mut_mat[seq1_state] * seq1_region->plength_observation2node);
                        
                        RealNumType seq2_state_evolves_seq1_state = model->freqi_freqj_qij[model->row_index[seq2_state] + seq1_state] * seq1_region->plength_observation2node * (1.0 + model->diagonal_mut_mat[seq2_state] * (blength + seq1_region->plength_observation2root));
                                                                                                                                                                                                    
                        total_factor *= seq1_state_evolves_seq2_state + seq2_state_evolves_seq1_state;
                    }
                    else
                    {
                        total_factor *= model->mutation_mat[model->row_index[seq1_state] + seq2_state] * (blength + (seq1_region->plength_observation2node < 0 ? 0 : seq1_region->plength_observation2node));
                    }
                }
            }
            // 3. e1.type = O
            else if (seq1_region->type == TYPE_O)
            {
                RealNumType blength13 = blength;
                if (seq1_region->plength_observation2node >= 0)
                {
                    blength13 = seq1_region->plength_observation2node;
                    if (blength > 0)
                        blength13 += blength;
                }
                    
                // 3.1. e1.type = O and e2.type = O
                if (seq2_region->type == TYPE_O)
                {
                    RealNumType tot = 0;
                    
                    StateType start_index = 0;
                    for (StateType i = 0; i < num_states; ++i)
                    {
                        RealNumType tot2 = 0;
                        
                        for (StateType j = 0; j < num_states; ++j)
                            if (seq2_region->likelihood[j] > 0.1)
                                tot2 += model->mutation_mat[start_index + j];
                        
                        tot2 *= blength13;
                        
                        if (seq2_region->likelihood[i] > 0.1)
                            tot2 += 1;
                        
                        tot += tot2 * seq1_region->likelihood[i];
                        
                        // update start_index
                        start_index += num_states;
                    }
                        
                    total_factor *= tot;
                }
                // 3.2. e1.type = O and e2.type = R or A/C/G/T
                else
                {
                    StateType seq2_state = seq2_region->type;
                    if (seq2_state == TYPE_R)
                        seq2_state = aln->ref_seq[pos];
                    
                    RealNumType tot2 = 0;
                    StateType mutation_index = model->row_index[seq2_state];
                    for (StateType j = 0; j < num_states; ++j)
                        tot2 += model->transposed_mut_mat[mutation_index + j] * seq1_region->likelihood[j];
                    
                    total_factor *= seq1_region->likelihood[seq2_state] + blength13 * tot2;
                }
            }
            // 4. e1.type = A/C/G/T
            else
            {
                int seq1_state = seq1_region->type;
                
                // 4.1. e1.type =  e2.type
                if (seq1_region->type == seq2_region->type)
                {
                    RealNumType total_blength = blength;
                    total_blength += (seq1_region->plength_observation2node < 0 ? 0 : seq1_region->plength_observation2node);
                    total_blength += (seq1_region->plength_observation2root < 0 ? 0 : seq1_region->plength_observation2root);

                    lh_cost += model->diagonal_mut_mat[seq1_state] * total_blength;
                }
                // e1.type = A/C/G/T and e2.type = O/A/C/G/T
                else
                {
                    // 4.2. e1.type = A/C/G/T and e2.type = O
                    if (seq2_region->type == TYPE_O)
                    {
                        RealNumType tot = 0.0;
                        
                        if (seq1_region->plength_observation2root >= 0)
                        {
                            RealNumType blength15 = blength + seq1_region->plength_observation2root;
                            
                            if (seq2_region->likelihood[seq1_state] > 0.1)
                                lh_cost += model->diagonal_mut_mat[seq1_state] * (blength15 + seq1_region->plength_observation2node);
                            else
                            {
                                StateType mutation_index = model->row_index[seq1_state];
                                StateType start_index = 0;
                                for (StateType i = 0; i < num_states; ++i)
                                {
                                    RealNumType tot2;
                                    if (seq1_state == i)
                                        tot2 = model->root_freqs[i] * (1.0 + model->transposed_mut_mat[mutation_index + i] * seq1_region->plength_observation2node);
                                    else
                                        tot2 = model->root_freqs[i] * (model->transposed_mut_mat[mutation_index + i] * seq1_region->plength_observation2node);
                                        
                                    RealNumType tot3 = 0;
                                    for (StateType j = 0; j < num_states; ++j)
                                        if (seq2_region->likelihood[j] > 0.1)
                                            tot3 += model->mutation_mat[start_index + j];
                                    
                                    if (seq2_region->likelihood[i] > 0.1)
                                        tot += tot2 * (1.0 + blength15 * tot3);
                                    else
                                        tot += tot2 * blength15 * tot3;
                                    
                                    // update start_index
                                    start_index += num_states;
                                }
                                
                                total_factor *= (tot * model->inverse_root_freqs[seq1_state]);
                            }
                        }
                        else
                        {
                            RealNumType tmp_blength = blength + (seq1_region->plength_observation2node < 0 ? 0 : seq1_region->plength_observation2node);
                            if (seq2_region->likelihood[seq1_state] > 0.1)
                                lh_cost += model->diagonal_mut_mat[seq1_state] * tmp_blength;
                            else
                            {
                                StateType start_index = model->row_index[seq1_state];
                                for (StateType j = 0; j < num_states; ++j)
                                    if (seq2_region->likelihood[j] > 0.1)
                                        tot += model->mutation_mat[start_index + j];
                                
                                total_factor *= tot * tmp_blength;
                            }
                        }
                    }
                    // 4.3. e1.type = A/C/G/T and e2.type = R or A/C/G/T
                    else
                    {
                        StateType seq2_state = seq2_region->type;
                        if (seq2_state == TYPE_R)
                            seq2_state = aln->ref_seq[pos];
                        
                        if (seq1_region->plength_observation2root >= 0)
                        {
                            // here we ignore contribution of non-parsimonious mutational histories
                            RealNumType seq1_state_evoloves_seq2_state = model->mutation_mat[model->row_index[seq1_state] + seq2_state] * (blength + seq1_region->plength_observation2root) * (1.0 + model->diagonal_mut_mat[seq1_state] * seq1_region->plength_observation2node);
                            
                            RealNumType seq2_state_evolves_seq2_state = model->freqi_freqj_qij[model->row_index[seq2_state] + seq1_state] * seq1_region->plength_observation2node * (1.0 + model->diagonal_mut_mat[seq2_state] * (blength + seq1_region->plength_observation2root));
                            
                            total_factor *= seq1_state_evoloves_seq2_state + seq2_state_evolves_seq2_state;
                        }
                        else
                        {
                            RealNumType tmp_blength = blength + (seq1_region->plength_observation2node < 0 ? 0 : seq1_region->plength_observation2node);
                            
                            total_factor *= model->mutation_mat[model->row_index[seq1_state] + seq2_state] * tmp_blength;
                        }
                    }
                }
            }
        }
         
        // approximately update lh_cost and total_factor
        if (total_factor <= MIN_CARRY_OVER)
        {
            if (total_factor < DBL_MIN)
                return -DBL_MAX;
            lh_cost += log(total_factor);
            total_factor = 1.0;
        }
        
        // update pos
        pos += length;
    }
    
    return lh_cost + log(total_factor);
}
