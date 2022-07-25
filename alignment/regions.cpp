//
//  regions.cpp
//  alignment
//
//  Created by NhanLT on 31/3/2022.
//

#include "regions.h"

Regions::Regions()
{
    // do nothing
}

Regions::~Regions()
{
    deleteRegions();
}

void Regions::deleteRegions()
{
    for (iterator it = begin(); it != end(); it++)
        delete (*it);
    clear();
}

Region* Regions::getRegion(PositionType id)
{
    return at(id);
}

void Regions::copyRegions(Regions* n_regions, StateType num_states)
{
    // delete the current regions
    deleteRegions();
    
    // clone regions one by one
    resize(n_regions->size());
    for (PositionType i = 0; i < n_regions->size(); i++)
    {
        Region* region = n_regions->getRegion(i);
        at(i) = new Region(region, num_states, true);
    }
}

void Regions::mergeRegionR(StateType num_states, double threshold)
{
    // clone the current regions
    Regions* tmp_regions = new Regions();
    tmp_regions->copyRegions(this, num_states);
    
    // delete the current regions
    deleteRegions();
    
    // browse regions one by one and try to merge consecutive regions of type R
    Region* region = tmp_regions->getRegion(0);
    
    for (PositionType i = 0; i < tmp_regions->size() - 1; i++)
    {
        Region* next_region = tmp_regions->getRegion(i + 1);
        
        // merge two consecutive regions if they are both 'R' and "similar" to each other (regarding plength_observation, and plength_from_root)
        if (!(region->type == TYPE_R && next_region->type == TYPE_R && fabs(region->plength_observation - next_region->plength_observation) < threshold && fabs(region->plength_from_root - next_region->plength_from_root) < threshold))
        {
            // add the current region into the current regions
            push_back(new Region(region, num_states, true));
            
            // move to the next entry
            region = tmp_regions->getRegion(i + 1);
        }
    }
    
    // add the last region into new_regions
    push_back(new Region(region, num_states, true));
    
    // delete tmp_regions
    delete tmp_regions;
}

void Regions::move2NextRegion(Regions* sequence, PositionType region_index, Region* &region, PositionType &current_pos, PositionType &end_pos, PositionType seq_length)
{
    ASSERT(region_index < sequence->size());
    
    // get the current region
    region = sequence->getRegion(region_index);
    
    // get the current position and end position
    current_pos = region->position;
    PositionType length = (region_index < sequence->size() - 1 ? sequence->getRegion(region_index + 1)->position : seq_length) - current_pos;
    end_pos = current_pos + length - 1;
}

void Regions::getNextSharedSegment(PositionType current_pos, PositionType seq_length, Regions* sequence1, Regions* sequence2, PositionType &seq1_index, PositionType &seq2_index, Region* &seq1_region, Region* &seq2_region, PositionType &seq1_end_pos, PositionType &seq2_end_pos, PositionType &length)
{
    PositionType seq1_pos, seq2_pos, end_pos;
    
    // move to the next region in sequence 1
    if (current_pos > seq1_end_pos)
    {
        seq1_index++;
        move2NextRegion(sequence1, seq1_index, seq1_region, seq1_pos, seq1_end_pos, seq_length);
    }
    
    // move to the next region in sequence 2
    if (current_pos > seq2_end_pos)
    {
        seq2_index++;
        move2NextRegion(sequence2, seq2_index, seq2_region, seq2_pos, seq2_end_pos, seq_length);
    }
    
    // compute the end_pos for the shared segment
    end_pos = seq1_end_pos < seq2_end_pos ? seq1_end_pos : seq2_end_pos;
    length = end_pos + 1 - current_pos;
}

int Regions::compareWithSample(Regions* sequence2, PositionType seq_length, StateType num_states)
{
    ASSERT(seq_length > 0);
    
    // init dummy variables
    PositionType seq1_index = -1;
    PositionType seq2_index = -1;
    bool seq1_more_info = false;
    bool seq2_more_info = false;
    PositionType pos = 0;
    Region *seq1_region, *seq2_region;
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
            for (StateType i = 0; i < num_states; i++)
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

bool Regions::areDiffFrom(Regions* regions2, PositionType seq_length, StateType num_states, Params* params)
{
    if (!regions2 || regions2->size() == 0)
        return true;
    
    // init variables
    PositionType seq1_index = -1;
    PositionType seq2_index = -1;
    PositionType pos = 0;
    Region *seq1_region, *seq2_region;
    PositionType seq1_end = -1, seq2_end = -1;
    PositionType length;
    
    // compare each pair of regions
    while (pos < seq_length)
    {
        // get the next shared segment in the two sequences
        Regions::getNextSharedSegment(pos, seq_length, this, regions2, seq1_index, seq2_index, seq1_region, seq2_region, seq1_end, seq2_end, length);
        
        // if any pair of regions is different in type -> the two sequences are different
        if (seq1_region->type != seq2_region->type)
            return true;
        
        // type R/ACGT
        if (seq1_region->type < num_states || seq1_region->type == TYPE_R)
        {
            // compare plength_from_root and plength_observation
            if (fabs(seq1_region->plength_from_root - seq2_region->plength_from_root) > params->threshold_prob
                ||fabs(seq1_region->plength_observation - seq2_region->plength_observation) > params->threshold_prob)
                return true;
        }
        
        // type O
        if (seq1_region->type == TYPE_O)
        {
            // compare plength_observation
            if (fabs(seq1_region->plength_observation - seq2_region->plength_observation) > params->threshold_prob)
                return true;
            
            // compare likelihood of each state
            for (StateType i = 0; i < num_states; i++)
            {
                double diff = fabs(seq1_region->likelihood[i] - seq2_region->likelihood[i]);
                
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
        pos += length;
    }
    
    return false;
}

void Regions::mergeUpperLower(Regions* &merged_regions, double upper_plength, Regions* lower_regions, double lower_plength, Alignment* aln, Model* model, double threshold_prob)
{
    // init variables
    PositionType seq1_index = -1;
    PositionType seq2_index = -1;
    PositionType pos = 0;
    Region *seq1_region, *seq2_region;
    PositionType seq1_end = -1, seq2_end = -1;
    PositionType length;
    StateType num_states = aln->num_states;
    PositionType seq_length = aln->ref_seq.size();
    
    // init merged_regions
    if (merged_regions) delete merged_regions;
    merged_regions = new Regions();
                
    while (pos <  seq_length)
    {
        // get the next shared segment in the two sequences
        Regions::getNextSharedSegment(pos, seq_length, this, lower_regions, seq1_index, seq2_index, seq1_region, seq2_region, seq1_end, seq2_end, length);
        
        // seq1_entry = 'N'
        if (seq1_region->type == TYPE_N)
        {
            // seq1_entry = 'N' and seq2_entry = 'N'
            if (seq2_region->type == TYPE_N)
                merged_regions->push_back(new Region(seq1_region->type, pos));
            // seq1_entry = 'N' and seq2_entry = O/R/ACGT
            else
            {
                // seq1_entry = 'N' and seq2_entry = O
                if (seq2_region->type == TYPE_O)
                {
                    double total_blength = lower_plength;
                    if (seq2_region->plength_observation >= 0)
                        total_blength = seq2_region->plength_observation + (lower_plength > 0 ? lower_plength: 0);
                        
                    double* new_lh = new double[num_states];
                    double sum_lh = 0;
                    
                    for (StateType i = 0; i < num_states; i++)
                    {
                        double tot = 0;
                        
                        if (total_blength > 0)
                        {
                            for (StateType j = 0; j < num_states; j++)
                                tot += model->mutation_mat[i * num_states + j] * seq2_region->likelihood[j];
                            
                            tot *= total_blength;
                        }
                        
                        tot += seq2_region->likelihood[i];
                        new_lh[i] = tot * model->root_freqs[i];
                        sum_lh += new_lh[i];
                    }

                    // normalize the new partial likelihood
                    for (StateType i = 0; i < num_states; i++)
                        new_lh[i] /= sum_lh;
                    
                    // add merged region into merged_regions
                    merged_regions->push_back(new Region(TYPE_O, pos, 0, 0, new_lh));
                }
                // seq1_entry = 'N' and seq2_entry = R/ACGT
                else
                {
                    if (seq2_region->plength_observation >= 0)
                    {
                        double new_plength = seq2_region->plength_observation ;
                        if (lower_plength > 0)
                            new_plength += lower_plength;
                        merged_regions->push_back(new Region(seq2_region->type, pos, new_plength, 0));
                    }
                    else
                    {
                        if (lower_plength > 0)
                            merged_regions->push_back(new Region(seq2_region->type, pos, lower_plength, 0));
                        else
                            merged_regions->push_back(new Region(seq2_region->type, pos));
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
                double total_blength = -1;
                
                if (seq1_region->plength_observation >= 0)
                {
                    total_blength = seq1_region->plength_observation;
                    if (upper_plength > 0)
                        total_blength += upper_plength;
                }
                else if (upper_plength > 0)
                    total_blength = upper_plength;
                
                if (total_blength > 0)
                {
                    double* new_lh = new double[num_states];
                    double sum_lh = 0;
                    
                    for (StateType i = 0; i < num_states; i++)
                    {
                        double tot = 0;
                        
                        for (StateType j = 0; j < num_states; j++)
                            tot += model->mutation_mat[j * num_states + i] * seq1_region->likelihood[j];
                        
                        tot *= total_blength;
                        tot += seq1_region->likelihood[i];
                        new_lh[i] = tot;
                        sum_lh += new_lh[i];
                    }

                    // normalize the new partial likelihood
                    for (StateType i = 0; i < num_states; i++)
                        new_lh[i] /= sum_lh;
                    
                    // add merged region into merged_regions
                    merged_regions->push_back(new Region(TYPE_O, pos, 0, 0, new_lh));
                }
                else
                {
                    double* new_lh = new double[num_states];
                    memcpy(new_lh, seq1_region->likelihood, sizeof(double) * num_states);
                    // add merged region into merged_regions
                    merged_regions->push_back(new Region(seq1_region->type, pos, 0, 0, new_lh));
                }
            }
            // seq2_entry = 'N' and seq1_entry = R/ACGT
            else
            {
                if (seq1_region->plength_from_root >= 0)
                {
                    double plength_from_root = seq1_region->plength_from_root;
                    if (upper_plength > 0)
                        plength_from_root += upper_plength;
                    merged_regions->push_back(new Region(seq1_region->type, pos, seq1_region->plength_observation, plength_from_root));
                }
                else if (seq1_region->plength_observation >= 0)
                {
                    double plength_observation = seq1_region->plength_observation;
                    if (upper_plength > 0)
                        plength_observation += upper_plength;
                    merged_regions->push_back(new Region(seq1_region->type, pos, plength_observation));
                }
                else
                {
                    if (upper_plength > 0)
                        merged_regions->push_back(new Region(seq1_region->type, pos, upper_plength));
                    else
                        merged_regions->push_back(new Region(seq1_region->type, pos));
                }
            }
        }
        // seq1_entry = seq2_entry = R/ACGT
        else if (seq1_region->type == seq2_region->type && (seq1_region->type < num_states || seq1_region->type == TYPE_R))
            merged_regions->push_back(new Region(seq1_region->type, pos));
        // cases where the new genome list entry will likely be of type "O"
        else
        {
            double total_blength_1 = upper_plength;
            if (seq1_region->plength_observation >= 0)
            {
                total_blength_1 = seq1_region->plength_observation;
                if (upper_plength > 0)
                    total_blength_1 += upper_plength;
                
                if (seq1_region->type != TYPE_O && seq1_region->plength_from_root >= 0)
                    total_blength_1 += seq1_region->plength_from_root;
            }
            
            double total_blength_2 = lower_plength;
            if (seq2_region->plength_observation >= 0)
            {
                total_blength_2 = seq2_region->plength_observation;
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
                
                merged_regions->push_back(new Region(seq2_region->type, pos));
            }
            // due to 0 distance, the entry will be of same type as entry1
            else if ((seq1_region->type < num_states || seq1_region->type == TYPE_R) && total_blength_1 <= 0)
            {
                merged_regions->push_back(new Region(seq1_region->type, pos));
            }
            // seq1_entry = O
            else if (seq1_region->type == TYPE_O)
            {
                double* new_lh = new double[num_states];
                
                // if total_blength_1 > 0 => compute new partial likelihood
                if (total_blength_1 > 0)
                {
                    for (StateType i = 0; i < num_states; i++)
                    {
                        double tot = 0;
                        
                        for (StateType j = 0; j < num_states; j++)
                            tot += model->mutation_mat[j * num_states + i] * seq1_region->likelihood[j];
                        
                        tot *= total_blength_1;
                        tot += seq1_region->likelihood[i];
                        new_lh[i] = tot;
                    }
                }
                // otherwise, clone the partial likelihood from seq1
                else
                    memcpy(new_lh, seq1_region->likelihood, sizeof(double) * num_states);
                
                double sum_new_lh = 0;

                // seq2 = O
                if (seq2_region->type == TYPE_O)
                {
                    for (StateType i = 0; i < num_states; i++)
                    {
                        double tot = 0;
                        
                        if (total_blength_2 > 0)
                        {
                            for (StateType j = 0; j < num_states; j++)
                                tot += model->mutation_mat[i * num_states + j] * seq2_region->likelihood[j];
        
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
                        seq2_state = aln->ref_seq[pos];
                    
                    if (total_blength_2 > 0)
                    {
                        for (StateType i = 0; i < num_states; i++)
                        {
                            StateType mutation_index = i * num_states + seq2_state;
                            if (i == seq2_state)
                                new_lh[i] *= (1.0 + model->mutation_mat[mutation_index] * total_blength_2);
                            else
                                new_lh[i] *= (model->mutation_mat[mutation_index] * total_blength_2);
                                
                            sum_new_lh += new_lh[i];
                        }
                    }
                    else
                    {
                        for (StateType i = 0; i < num_states; i++)
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
                for (StateType i = 0; i < num_states; i++)
                    new_lh[i] /= sum_new_lh;
                
                StateType new_state = simplifyO(new_lh, aln->ref_seq[pos], num_states, threshold_prob);

                if (new_state == TYPE_O)
                    merged_regions->push_back(new Region(TYPE_O, pos, 0, 0, new_lh));
                else
                {
                    delete[] new_lh;
                    merged_regions->push_back(new Region(new_state, pos));
                }
            }
            // seq1_entry = R/ACGT
            else
            {
                StateType seq1_state = seq1_region->type;
                if (seq1_state == TYPE_R)
                    seq1_state = aln->ref_seq[pos];
                
                double* new_lh = new double[num_states];
                double sum_new_lh = 0;
                
                if (seq1_region->plength_from_root >= 0)
                {
                    double length_to_root = seq1_region->plength_from_root;
                    if (upper_plength > 0)
                        length_to_root += upper_plength;
                    double* root_vec = new double[num_states];
                    memcpy(root_vec, model->root_freqs, sizeof(double) * num_states);
                    
                    for (StateType i = 0; i < num_states; i++)
                    {
                        StateType mutation_index = i * num_states + seq1_state;
                        
                        if (i == seq1_state)
                            root_vec[i] *= (1.0 + model->mutation_mat[mutation_index] * seq1_region->plength_observation);
                        else
                            root_vec[i] *= model->mutation_mat[mutation_index] * seq1_region->plength_observation;
                    }
                        
                    for (StateType i = 0; i < num_states; i++)
                    {
                        double tot = 0;
                        
                        for (StateType j = 0; j < num_states; j++)
                            tot += model->mutation_mat[j * num_states + i] * root_vec[j];
                        
                        tot *= length_to_root;
                        tot += root_vec[i];
                        new_lh[i] = tot;
                    }
                }
                else
                {
                    if (total_blength_1 > 0)
                    {
                        for (StateType i = 0; i < num_states; i++)
                        {
                            StateType mutation_index = seq1_state * num_states + i;
                            
                            if (i == seq1_state)
                                new_lh[i] = 1.0 + model->mutation_mat[mutation_index] * total_blength_1;
                            else
                                new_lh[i] = model->mutation_mat[mutation_index] * total_blength_1;
                        }
                    }
                    else
                    {
                        for (StateType i = 0; i < num_states; i++)
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
                    for (StateType i = 0; i < num_states; i++)
                    {
                        double tot = 0;
                        
                        if (total_blength_2 > 0)
                        {
                            for (StateType j = 0; j < num_states; j++)
                                tot += model->mutation_mat[i * num_states + j] * seq2_region->likelihood[j];
                            
                            tot *= total_blength_2;
                        }
                        
                        tot += seq2_region->likelihood[i];
                        new_lh[i] *= tot;
                        sum_new_lh += new_lh[i];
                    }
                    
                    // normalize the new partial lh
                    for (StateType i = 0; i < num_states; i++)
                        new_lh[i] /= sum_new_lh;
                    
                    StateType new_state = simplifyO(new_lh, aln->ref_seq[pos], num_states, threshold_prob);
                    
                    if (new_state == TYPE_O)
                        merged_regions->push_back(new Region(TYPE_O, pos, 0, 0, new_lh));
                    else
                        merged_regions->push_back(new Region(new_state, pos));
                }
                // seq1 = ACGT and different from seq2 = R/ACGT
                else
                {
                    StateType seq2_state = seq2_region->type;
                    
                    if (seq2_state == TYPE_R)
                        seq2_state = aln->ref_seq[pos];
                    
                    if (total_blength_2 > 0)
                    {
                        for (StateType i = 0; i < num_states; i++)
                        {
                            StateType mutation_index = i * num_states + seq2_state;
                            
                            if (i == seq2_state)
                                new_lh[i] *= 1.0 + model->mutation_mat[mutation_index] * total_blength_2;
                            else
                                new_lh[i] *= model->mutation_mat[mutation_index] * total_blength_2;
                            
                            sum_new_lh += new_lh[i];
                        }
                    }
                    else
                    {
                        for (StateType i = 0; i < num_states; i++)
                        {
                            if (i != seq2_state)
                                new_lh[i] = 0;
                                
                            sum_new_lh += new_lh[i];
                        }
                    }
                    
                    // normalize the new partial lh
                    for (StateType i = 0; i < num_states; i++)
                        new_lh[i] /= sum_new_lh;

                    // add new region into the merged regions
                    merged_regions->push_back(new Region(TYPE_O, pos, 0, 0, new_lh));
                }
            }
        }

        // update pos
        pos += length;
    }

    // try to merge consecutive and similar 'R' regions together
    merged_regions->mergeRegionR(num_states, threshold_prob);
}

StateType Regions::simplifyO(double* &partial_lh, StateType ref_state, StateType num_states, double threshold_prob)
{
    // dummy variables
    ASSERT(partial_lh);
    double max_prob = 0;
    double max_index = 0;
    StateType high_prob_count = 0;
    
    // Check all states one by one
    for (StateType i = 0; i < num_states; i++)
    {
        // record the state with the highest likelihood
        if (partial_lh[i] > max_prob)
        {
            max_prob = partial_lh[i];
            max_index = i;
        }
        
        // count the number of states that have the likelihood greater than a threshold
        if (partial_lh[i] > threshold_prob)
            high_prob_count++;
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

double Regions::mergeTwoLowers(Regions* &merged_regions, double plength1, Regions* regions2, double plength2, Alignment* aln, Model* model, double threshold_prob, double* cumulative_rate, bool return_log_lh)
{
    // init variables
    double log_lh = 0;
    PositionType seq1_index = -1;
    PositionType seq2_index = -1;
    PositionType pos = 0;
    StateType num_states = aln->num_states;
    Region *seq1_region, *seq2_region;
    PositionType seq1_end = -1, seq2_end = -1;
    PositionType length;
    PositionType seq_length = aln->ref_seq.size();
    
    // init merged_regions
    if (merged_regions) delete merged_regions;
    merged_regions = new Regions();
                
    while (pos < seq_length)
    {
        // get the next shared segment in the two sequences
        Regions::getNextSharedSegment(pos, seq_length, this, regions2, seq1_index, seq2_index, seq1_region, seq2_region, seq1_end, seq2_end, length);
        
        // seq1_entry = 'N'
        if (seq1_region->type == TYPE_N)
        {
            // seq1_entry = 'N' and seq2_entry = 'N'
            if (seq2_region->type == TYPE_N)
                merged_regions->push_back(new Region(seq1_region->type, pos));
            // seq1_entry = 'N' and seq2_entry = O/R/ACGT
            else
            {
                // seq1_entry = 'N' and seq2_entry = O
                if (seq2_region->type == TYPE_O)
                {
                    // add merged region into merged_regions
                    Region* new_region = new Region(seq2_region, num_states);
                    new_region->position = pos;
                    if (seq2_region->plength_observation >= 0)
                    {
                        if (plength2 > 0)
                            new_region->plength_observation += plength2;
                    }
                    else
                    {
                        if (plength2 > 0)
                            new_region->plength_observation = plength2;
                    }
                    merged_regions->push_back(new_region);
                }
                // seq1_entry = 'N' and seq2_entry = R/ACGT
                else
                {
                    if (seq2_region->plength_observation >= 0)
                    {
                        double new_plength = seq2_region->plength_observation;
                        if (plength2 > 0)
                            new_plength += plength2;
                        merged_regions->push_back(new Region(seq2_region->type, pos, new_plength));
                    }
                    else
                    {
                        if (plength2 > 0)
                            merged_regions->push_back(new Region(seq2_region->type, pos, plength2));
                        else
                            merged_regions->push_back(new Region(seq2_region->type, pos));
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
                Region* new_region = new Region(seq1_region, num_states);
                new_region->position = pos;
                if (seq1_region->plength_observation >= 0)
                {
                    if (plength1 > 0)
                        new_region->plength_observation += plength1;
                }
                else
                {
                    if (plength1 > 0)
                        new_region->plength_observation = plength1;
                }
                merged_regions->push_back(new_region);
            }
            // seq1_entry = 'N' and seq2_entry = R/ACGT
            else
            {
                if (seq1_region->plength_observation >= 0)
                {
                    double new_plength = seq1_region->plength_observation;
                    if (plength1 > 0)
                        new_plength += plength1;
                    merged_regions->push_back(new Region(seq1_region->type, pos, new_plength));
                }
                else
                {
                    if (plength1 > 0)
                        merged_regions->push_back(new Region(seq1_region->type, pos, plength1));
                    else
                        merged_regions->push_back(new Region(seq1_region->type, pos));
                        
                }
            }
        }
        // neither seq1_entry nor seq2_entry = N
        else
        {
            double total_blength_1 = plength1;
            if (seq1_region->plength_observation >= 0)
            {
                total_blength_1 = seq1_region->plength_observation;
                if (plength1 > 0)
                    total_blength_1 += plength1;
            }
            
            double total_blength_2 = plength2;
            if (seq2_region->plength_observation >= 0)
            {
                total_blength_2 = seq2_region->plength_observation;
                if (plength2 > 0)
                    total_blength_2 += plength2;
            }
            
            // seq1_entry and seq2_entry are identical seq1_entry = R/ACGT
            if (seq1_region->type == seq2_region->type && (seq1_region->type == TYPE_R || seq1_region->type < num_states))
            {
                merged_regions->push_back(new Region(seq1_region->type, pos));
                
                if (return_log_lh)
                {
                    // convert total_blength_1 and total_blength_2 to zero if they are -1
                    if (total_blength_1 < 0) total_blength_1 = 0;
                    if (total_blength_2 < 0) total_blength_2 = 0;
                    
                    if (seq1_region->type == TYPE_R)
                        log_lh += (total_blength_1 + total_blength_2) * (cumulative_rate[pos + length - 1] - (pos == 0 ? 0 : cumulative_rate[pos - 1]));
                    else
                        log_lh += model->mutation_mat[seq1_region->type * (num_states + 1)] * (total_blength_1 + total_blength_2);
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
                double* new_lh = new double[num_states];
                double sum_lh = 0;
                
                if (total_blength_1 > 0)
                {
                    for (StateType i = 0; i < num_states; i++)
                    {
                        double tot = 0;
                        
                        for (StateType j = 0; j < num_states; j++)
                            tot += model->mutation_mat[i * num_states + j] * seq1_region->likelihood[j];
                        
                        tot *= total_blength_1;
                        tot += seq1_region->likelihood[i];
                        new_lh[i] = tot;
                    }
                }
                // otherwise, clone the partial likelihood from seq1
                else
                    memcpy(new_lh, seq1_region->likelihood, sizeof(double) * num_states);

                // seq1_entry = O and seq2_entry = O
                if (seq2_region->type == TYPE_O)
                {
                    for (StateType i = 0; i < num_states; i++)
                    {
                        double tot = 0;
                        
                        if (total_blength_2 > 0)
                        {
                            for (StateType j = 0; j < num_states; j++)
                                tot += model->mutation_mat[i * num_states + j] * seq2_region->likelihood[j];
                            
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
                    for (StateType i = 0; i < num_states; i++)
                        new_lh[i] /= sum_lh;
                    
                    StateType new_state = simplifyO(new_lh, aln->ref_seq[pos], num_states, threshold_prob);

                    if (new_state == TYPE_O)
                        merged_regions->push_back(new Region(new_state, pos, 0, 0, new_lh));
                    else
                        merged_regions->push_back(new Region(new_state, pos));
                    
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
                        for (StateType i = 0; i < num_states; i++)
                        {
                            StateType mutation_index = i * num_states + seq2_state;
                            
                            if (seq2_state == i)
                                new_lh[i] *= (1 + model->mutation_mat[mutation_index] * total_blength_2);
                            else
                                new_lh[i] *= (model->mutation_mat[mutation_index] * total_blength_2);
                            
                            sum_lh += new_lh[i];
                        }
                        
                        // normalize new partial lh
                        for (StateType i = 0; i < num_states; i++)
                            new_lh[i] /= sum_lh;
                        
                        StateType new_state = simplifyO(new_lh, aln->ref_seq[pos],  num_states, threshold_prob);

                        if (new_state == TYPE_O)
                            merged_regions->push_back(new Region(new_state, pos, 0, 0, new_lh));
                        else
                            merged_regions->push_back(new Region(new_state, pos));
                        
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
                        
                        merged_regions->push_back(new Region(seq2_region->type, pos));
                    
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
                
                double* new_lh = new double[num_states];
                double sum_lh = 0;
                
                if (total_blength_1 > 0)
                {
                    for (StateType i = 0; i < num_states; i++)
                    {
                        StateType mutation_index = i * num_states + seq1_state;
                        
                        if (seq1_state == i)
                            new_lh[i] = 1 + model->mutation_mat[mutation_index] * total_blength_1;
                        else
                            new_lh[i] = model->mutation_mat[mutation_index] * total_blength_1;
                    }
                }
                else
                {
                    for (StateType i = 0; i < num_states; i++)
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
                    for (StateType i = 0; i < num_states; i++)
                    {
                        double tot = 0;
                        
                        if (total_blength_2 > 0)
                        {
                            for (StateType j = 0; j < num_states; j++)
                                tot += model->mutation_mat[i * num_states + j] * seq2_region->likelihood[j];
                            
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
                    for (StateType i = 0; i < num_states; i++)
                        new_lh[i] /= sum_lh;
                    
                    StateType new_state = simplifyO(new_lh, aln->ref_seq[pos], num_states, threshold_prob);

                    if (new_state == TYPE_O)
                        merged_regions->push_back(new Region(new_state, pos, 0, 0, new_lh));
                    else
                        merged_regions->push_back(new Region(new_state, pos));
                    
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
                        for (StateType i = 0; i < num_states; i++)
                        {
                            StateType mutation_index = i * num_states + seq2_state;
                            
                            if (seq2_state == i)
                                new_lh[i] *= (1 + model->mutation_mat[mutation_index] * total_blength_2);
                            else
                                new_lh[i] *= (model->mutation_mat[mutation_index] * total_blength_2);
                            
                            sum_lh += new_lh[i];
                        }
                        
                        // normalize new partial lh
                        for (StateType i = 0; i < num_states; i++)
                            new_lh[i] /= sum_lh;
                        
                        StateType new_state = simplifyO(new_lh, aln->ref_seq[pos], num_states, threshold_prob);

                        if (new_state == TYPE_O)
                            merged_regions->push_back(new Region(new_state, pos, 0, 0, new_lh));
                        else
                            merged_regions->push_back(new Region(new_state, pos));
                        
                        if (return_log_lh)
                            log_lh += log(sum_lh);
                    }
                    else
                    {
                        merged_regions->push_back(new Region(seq2_region->type, pos));
                    
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

double Regions::computeAbsoluteLhAtRoot(Alignment* aln, Model* model, vector<vector<PositionType>> &cumulative_base)
{
    // dummy variables
    double log_lh = 0;
    double log_factor = 1;
    StateType num_states = aln->num_states;
    PositionType seq_length = aln->ref_seq.size();
    
    // browse regions one by one to compute the likelihood of each region
    for (PositionType region_index = 0; region_index < size(); region_index++)
    {
        Region* region = getRegion(region_index);
        
        // type R
        if (region->type == TYPE_R)
        {
            PositionType start_pos = region->position;
            PositionType end_pos = (region_index == size() - 1 ? seq_length - 1 : getRegion(region_index + 1)->position - 1);
            
            for (StateType i = 0; i < num_states; i++)
                log_lh += model->root_log_freqs[i] * (cumulative_base[end_pos][i] - (start_pos == 0 ? 0 : cumulative_base[start_pos - 1][i]));
        }
        // type ACGT
        else if (region->type < num_states)
            log_lh += model->root_log_freqs[region->type];
        // type O
        else if (region->type == TYPE_O)
        {
            double tot = 0;
            for (StateType i = 0; i < num_states; i++)
                tot += model->root_freqs[i] * region->likelihood[i];
                log_factor *= tot;
        }
    }

    // update log_lh
    log_lh += log(log_factor);
    
    // return the absolute likelihood
    return log_lh;
}

Regions* Regions::computeTotalLhAtRoot(StateType num_states, Model* model, double blength)
{
    Regions* total_lh = new Regions();
    
    for (Region* region : (*this))
    {
        // type N
        if (region->type == TYPE_N)
        {
            Region* new_region = new Region(region, num_states, false);
            total_lh->push_back(new_region);
        }
        else
        {
            // type O
            if (region->type == TYPE_O)
            {
                // compute total blength
                double total_blength = blength;
                if (region->plength_observation >= 0)
                {
                    total_blength = region->plength_observation;
                    if (blength > 0)
                        total_blength += blength;
                }
                
                // init new likelihood
                double* new_likelihood = new double[num_states];
                double sum_likelihood = 0;
                
                for (StateType i = 0; i < num_states; i++)
                {
                    double tot = 0.0;
                    
                    if (total_blength > 0)
                    {
                        for (StateType j = 0; j < num_states; j++)
                            tot += model->mutation_mat[i * num_states + j] * region->likelihood[j];
                        tot *= total_blength;
                    }
                    
                    tot += region->likelihood[i];
                    new_likelihood[i] = tot * model->root_freqs[i];
                    sum_likelihood += new_likelihood[i];
                }
                
                // normalize likelihood
                sum_likelihood = 1 / sum_likelihood;
                for (StateType i = 0; i < num_states; i++)
                    new_likelihood[i] *= sum_likelihood;
                
                // add new region to the total_lh_regions
                Region* new_region = new Region(region, num_states, false);
                new_region->likelihood = new_likelihood;
                total_lh->push_back(new_region);
            }
            // other types: R or A/C/G/T
            else
            {
                // add new region to the total_lh_regions
                Region* new_region = new Region(region, num_states, false);
                
                if (new_region->plength_observation >= 0)
                {
                    if (blength > 0)
                        new_region->plength_observation += blength;

                    new_region->plength_from_root = 0;
                }
                else if (blength > 0)
                {
                    new_region->plength_observation = blength;
                    new_region->plength_from_root = 0;
                }
                
                total_lh->push_back(new_region);
            }
        }
    }
    
    return total_lh;
}

// this implementation derives from appendProbNode
/*double Regions::calculatePlacementCost(Regions* parent_regions, Regions* child_regions, double blength)
{
    // init dummy variables
    double lh_cost = 0;
    PositionType seq1_index = -1;
    PositionType seq2_index = -1;
    PositionType pos = 0;
    double total_factor = 1;
    StateType num_states = tree->aln->num_states;
    Region *seq1_region, *seq2_region;
    PositionType seq1_end = -1, seq2_end = -1;
    PositionType length;
    double minimum_carry_over = DBL_MIN * 1e50;
    double total_blength = blength;
    PositionType seq_length = tree->aln->ref_seq.size();
    
    while (pos < seq_length)
    {
        // get the next shared segment in the two sequences
        Regions::getNextSharedSegment(pos, seq_length, parent_regions, child_regions, seq1_index, seq2_index, seq1_region, seq2_region, seq1_end, seq2_end, length);
        
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
            if (seq1_region->plength_from_root >= 0)
                total_blength = seq1_region->plength_from_root + (blength >= 0 ? blength : 0);
            else if (seq1_region->plength_observation >= 0)
                total_blength = seq1_region->plength_observation + (blength >= 0 ? blength : 0);
            else
                total_blength = blength;
                
            if (seq2_region->plength_observation >= 0)
                total_blength = (total_blength > 0 ? total_blength : 0) + seq2_region->plength_observation;
            
            // 2. e1.type = R
            if (seq1_region->type == TYPE_R)
            {
                // 2.1. e1.type = R and e2.type = R
                if (seq2_region->type == TYPE_R)
                {
                    if (seq1_region->plength_from_root >= 0)
                        total_blength += seq1_region->plength_observation;
                    
                    if (total_blength > 0)
                        lh_cost += total_blength * (cumulative_rate[pos + length - 1] - (pos == 0 ? 0 : cumulative_rate[pos - 1]));
                }
                // 2.2. e1.type = R and e2.type = O
                else if (seq2_region->type == TYPE_O)
                {
                    double tot = 0;
                    StateType seq1_state = tree->aln->ref_seq[pos];
                    
                    if (seq1_region->plength_from_root >= 0)
                    {
                        for (StateType i = 0; i < num_states; i++)
                        {
                            double tot2;
                            StateType mutation_index = i * num_states + seq1_state;
                            
                            if (seq1_state == i)
                                tot2 = model->root_freqs[i] * (1.0 + model->mutation_mat[mutation_index] * seq1_region->plength_observation);
                            else
                                tot2 = model->root_freqs[i] * (model->mutation_mat[mutation_index] * seq1_region->plength_observation);
                                
                            double tot3 = 0;
                            if (total_blength > 0)
                            {
                                for (StateType j = 0; j < num_states; j++)
                                    tot3 += model->mutation_mat[i * num_states + j] * seq2_region->likelihood[j];
                            }
                            
                            tot += tot2 * (seq2_region->likelihood[i] + total_blength * tot3);
                        }
                        
                        tot /= model->root_freqs[seq1_state];
                    }
                    else
                    {
                        if (total_blength > 0)
                        {
                            for (StateType j = 0; j < num_states; j++)
                                tot += model->mutation_mat[seq1_state * num_states + j] * seq2_region->likelihood[j];
                            
                            tot *= total_blength;
                        }
                        
                        tot += seq2_region->likelihood[seq1_state];
                    }
                        
                    total_factor *= tot;
                }
                // 2.3. e1.type = R and e2.type = A/C/G/T
                else
                {
                    StateType seq1_state = tree->aln->ref_seq[pos];
                    StateType seq2_state = seq2_region->type;
                    
                    if (seq1_region->plength_from_root >= 0)
                    {
                        if (total_blength > 0)
                            total_factor *= ((model->root_freqs[seq1_state] * model->mutation_mat[seq1_state * num_states + seq2_state] * total_blength * (1.0 + model->mutation_mat[seq1_state * (num_states + 1)] * seq1_region->plength_observation) + model->root_freqs[seq2_state] * model->mutation_mat[seq2_state * num_states + seq1_state] * seq1_region->plength_observation * (1.0 + model->mutation_mat[seq2_state * (num_states + 1)] * total_blength)) / model->root_freqs[seq1_state]);
                        else
                            total_factor *= ((model->root_freqs[seq2_state] * model->mutation_mat[seq2_state * num_states + seq1_state] * seq1_region->plength_observation) / model->root_freqs[seq1_state]);
                    }
                    else
                    {
                        if (total_blength > 0)
                            total_factor *= model->mutation_mat[seq1_state * num_states + seq2_state] * total_blength;
                        else
                            return -DBL_MAX;
                    }
                }
            }
            // 3. e1.type = O
            else if (seq1_region->type == TYPE_O)
            {
                double blength13 = blength;
                if (seq1_region->plength_observation >= 0)
                {
                    blength13 = seq1_region->plength_observation;
                    if (blength > 0)
                        blength13 += blength;
                }
                    
                // 3.1. e1.type = O and e2.type = O
                if (seq2_region->type == TYPE_O)
                {
                    double tot = 0;
                    
                    for (StateType i = 0; i < num_states; i++)
                    {
                        double tot2 = 0;
                        
                        for (StateType j = 0; j < num_states; j++)
                            if (seq2_region->likelihood[j] > 0.1)
                                tot2 += model->mutation_mat[i * num_states + j];
                        
                        tot2 *= blength13;
                        
                        if (seq2_region->likelihood[i] > 0.1)
                            tot2 += 1;
                        
                        tot += tot2 * seq1_region->likelihood[i];
                    }
                        
                    total_factor *= tot;
                }
                // 3.2. e1.type = O and e2.type = R or A/C/G/T
                else
                {
                    StateType seq2_state = seq2_region->type;
                    if (seq2_state == TYPE_R)
                        seq2_state = tree->aln->ref_seq[pos];
                    
                    double tot2 = 0;
                    if (total_blength > 0)
                    {
                        for (StateType j = 0; j < num_states; j++)
                            tot2 += model->mutation_mat[j * num_states + seq2_state] * seq1_region->likelihood[j];
                    }
                    
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
                    if (seq1_region->plength_from_root >= 0)
                        total_blength += seq1_region->plength_observation;
                    
                    if (total_blength > 0)
                        lh_cost += model->mutation_mat[seq1_state * (num_states + 1)] * total_blength;
                }
                // e1.type = A/C/G/T and e2.type = O/A/C/G/T
                else
                {
                    // 4.2. e1.type = A/C/G/T and e2.type = O
                    if (seq2_region->type == TYPE_O)
                    {
                        double tot = 0.0;
                        
                        if (seq1_region->plength_from_root >= 0)
                        {
                            for (StateType i = 0; i < num_states; i++)
                            {
                                double tot2;
                                StateType mutation_index = i * num_states + seq1_state;
                                if (seq1_state == i)
                                    tot2 = model->root_freqs[i] * (1.0 + model->mutation_mat[mutation_index] * seq1_region->plength_observation);
                                else
                                    tot2 = model->root_freqs[i] * (model->mutation_mat[mutation_index] * seq1_region->plength_observation);
                                    
                                double tot3 = 0;
                                for (StateType j = 0; j < num_states; j++)
                                    tot3 += model->mutation_mat[i * num_states + j] * seq2_region->likelihood[j];
                                tot += tot2 * (seq2_region->likelihood[i] + total_blength * tot3);
                            }
                            
                            total_factor *= (tot / model->root_freqs[seq1_state]);
                        }
                        else
                        {
                            for (StateType j = 0; j < num_states; j++)
                                tot += model->mutation_mat[seq1_state * num_states + j] * seq2_region->likelihood[j];
                            
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
                            seq2_state = tree->aln->ref_seq[pos];
                        
                        if (seq1_region->plength_from_root >= 0)
                        {
                            if (total_blength > 0)
                                total_factor *= ((model->root_freqs[seq1_state] * model->mutation_mat[seq1_state * num_states + seq2_state] * total_blength * (1.0 + model->mutation_mat[seq1_state * (num_states + 1)] * seq1_region->plength_observation) + model->root_freqs[seq2_state] * model->mutation_mat[seq2_state * num_states + seq1_state] * seq1_region->plength_observation * (1.0 + model->mutation_mat[seq2_state * (num_states + 1)] * total_blength)) / model->root_freqs[seq1_state]);
                            else
                                total_factor *= ((model->root_freqs[seq2_state] * model->mutation_mat[seq2_state * num_states + seq1_state] * seq1_region->plength_observation) / model->root_freqs[seq1_state]);
                        }
                        else
                        {
                            if (total_blength > 0)
                                total_factor *= model->mutation_mat[seq1_state * num_states + seq2_state] * total_blength;
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
}*/

// this implementation derives from appendProb
double Regions::calculatePlacementCost(Alignment* aln, Model* model, double* cumulative_rate, Regions* child_regions, double blength)
{
    // init dummy variables
    double lh_cost = 0;
    PositionType seq1_index = -1;
    PositionType seq2_index = -1;
    PositionType pos = 0;
    double total_factor = 1;
    StateType num_states = aln->num_states;
    Region *seq1_region, *seq2_region;
    PositionType seq1_end = -1, seq2_end = -1;
    PositionType length;
    double minimum_carry_over = DBL_MIN * 1e50;
    if (blength < 0) blength = 0;
    double total_blength = blength;
    PositionType seq_length = aln->ref_seq.size();
    
    while (pos < seq_length)
    {
        // get the next shared segment in the two sequences
        Regions::getNextSharedSegment(pos, seq_length, this, child_regions, seq1_index, seq2_index, seq1_region, seq2_region, seq1_end, seq2_end, length);
        
        // 1. e1.type = N || e2.type = N
        if (seq2_region->type == TYPE_N || seq1_region->type == TYPE_N)
        {
            pos += length;
            continue;
        }
        // e1.type != N && e2.type != N
        else
        {
            // 2. e1.type = R
            if (seq1_region->type == TYPE_R)
            {
                // 2.1. e1.type = R and e2.type = R
                if (seq2_region->type == TYPE_R)
                {
                    if (seq1_region->plength_observation < 0 && seq1_region->plength_from_root < 0)
                        lh_cost += blength * (cumulative_rate[pos + length - 1] - (pos == 0 ? 0 : cumulative_rate[pos - 1]));
                    else
                    {
                        total_blength = blength + seq1_region->plength_observation;
                        if (seq1_region->plength_from_root < 0)
                            lh_cost += total_blength * (cumulative_rate[pos + length - 1] - (pos == 0 ? 0 : cumulative_rate[pos - 1]));
                        else
                            // here contribution from root frequency gets added and subtracted so it's ignored
                            lh_cost += (total_blength + seq1_region->plength_from_root) * (cumulative_rate[pos + length - 1] - (pos == 0 ? 0 : cumulative_rate[pos - 1]));
                    }
                }
                // 2.2. e1.type = R and e2.type = O
                else if (seq2_region->type == TYPE_O)
                {
                    StateType seq1_state = aln->ref_seq[pos];
                    if (seq1_region->plength_from_root >= 0)
                    {
                        total_blength = seq1_region->plength_from_root + (blength > 0 ? blength : 0);
                        
                        if (seq2_region->likelihood[seq1_state] > 0.1)
                        {
                            total_blength += seq1_region->plength_observation;
                            
                            // here contribution from root frequency can also be also ignored
                            lh_cost += model->mutation_mat[seq1_state * (num_states + 1)] * total_blength;
                        }
                        else
                        {
                            double tot = 0;
                            for (StateType i = 0; i < num_states; i++)
                            {
                                double tot2;
                                StateType mutation_index = i * num_states + seq1_state;
                                
                                if (seq1_state == i)
                                    tot2 = model->root_freqs[i] * (1.0 + model->mutation_mat[mutation_index] * seq1_region->plength_observation);
                                else
                                    tot2 = model->root_freqs[i] * (model->mutation_mat[mutation_index] * seq1_region->plength_observation);
                                    
                                double tot3 = 0;
                                for (StateType j = 0; j < num_states; j++)
                                    if (seq2_region->likelihood[j] > 0.1)
                                        tot3 += model->mutation_mat[i * num_states + j];
                                tot3 *= total_blength;
                                
                                if (seq2_region->likelihood[i] > 0.1)
                                    tot3 += 1;
                                
                                tot += tot2 * tot3;
                            }
                            
                            total_factor *= tot /model->root_freqs[seq1_state];
                        }
                    }
                    else
                    {
                        if (seq2_region->likelihood[seq1_state] > 0.1)
                        {
                            if (seq1_region->plength_observation >= 0)
                                lh_cost += model->mutation_mat[seq1_state * (num_states + 1)] * (blength + seq1_region->plength_observation);
                            else
                                lh_cost += model->mutation_mat[seq1_state * (num_states + 1)] * blength;
                        }
                        else
                        {
                            double tot = 0;
                            for (StateType i = 0; i < num_states; i++)
                                if (seq2_region->likelihood[i] > 0.1)
                                    tot += model->mutation_mat[seq1_state * num_states + i];
                            
                            if (seq1_region->plength_observation >= 0)
                                total_factor *= tot * (blength + seq1_region->plength_observation);
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
                    
                    if (seq1_region->plength_from_root >= 0)
                    {
                        total_factor *= ((model->root_freqs[seq1_state] * model->mutation_mat[seq1_state * num_states + seq2_state] * blength * (1.0 + model->mutation_mat[seq1_state * (num_states + 1)] * seq1_region->plength_observation) + model->root_freqs[seq2_state] * model->mutation_mat[seq2_state * num_states + seq1_state] * seq1_region->plength_observation * (1.0 + model->mutation_mat[seq2_state * (num_states + 1)] * (blength + seq1_region->plength_from_root))) / model->root_freqs[seq1_state]);
                    }
                    else
                        total_factor *= model->mutation_mat[seq1_state * num_states + seq2_state] * (blength + (seq1_region->plength_observation < 0 ? 0 : seq1_region->plength_observation));
                }
            }
            // 3. e1.type = O
            else if (seq1_region->type == TYPE_O)
            {
                double blength13 = blength;
                if (seq1_region->plength_observation >= 0)
                {
                    blength13 = seq1_region->plength_observation;
                    if (blength > 0)
                        blength13 += blength;
                }
                    
                // 3.1. e1.type = O and e2.type = O
                if (seq2_region->type == TYPE_O)
                {
                    double tot = 0;
                    
                    for (StateType i = 0; i < num_states; i++)
                    {
                        double tot2 = 0;
                        
                        for (StateType j = 0; j < num_states; j++)
                            if (seq2_region->likelihood[j] > 0.1)
                                tot2 += model->mutation_mat[i * num_states + j];
                        
                        tot2 *= blength13;
                        
                        if (seq2_region->likelihood[i] > 0.1)
                            tot2 += 1;
                        
                        tot += tot2 * seq1_region->likelihood[i];
                    }
                        
                    total_factor *= tot;
                }
                // 3.2. e1.type = O and e2.type = R or A/C/G/T
                else
                {
                    StateType seq2_state = seq2_region->type;
                    if (seq2_state == TYPE_R)
                        seq2_state = aln->ref_seq[pos];
                    
                    double tot2 = 0;
                    for (StateType j = 0; j < num_states; j++)
                        tot2 += model->mutation_mat[j * num_states + seq2_state] * seq1_region->likelihood[j];
                    
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
                    double total_blength = blength;
                    total_blength += (seq1_region->plength_observation < 0 ? 0 : seq1_region->plength_observation);
                    total_blength += (seq1_region->plength_from_root < 0 ? 0 : seq1_region->plength_from_root);

                    lh_cost += model->mutation_mat[seq1_state * (num_states + 1)] * total_blength;
                }
                // e1.type = A/C/G/T and e2.type = O/A/C/G/T
                else
                {
                    // 4.2. e1.type = A/C/G/T and e2.type = O
                    if (seq2_region->type == TYPE_O)
                    {
                        double tot = 0.0;
                        
                        if (seq1_region->plength_from_root >= 0)
                        {
                            double blength15 = blength + seq1_region->plength_from_root;
                            
                            if (seq2_region->likelihood[seq1_state] > 0.1)
                                lh_cost += model->mutation_mat[seq1_state * (num_states + 1)] * (blength15 + seq1_region->plength_observation);
                            else
                            {
                                for (StateType i = 0; i < num_states; i++)
                                {
                                    double tot2;
                                    StateType mutation_index = i * num_states + seq1_state;
                                    if (seq1_state == i)
                                        tot2 = model->root_freqs[i] * (1.0 + model->mutation_mat[mutation_index] * seq1_region->plength_observation);
                                    else
                                        tot2 = model->root_freqs[i] * (model->mutation_mat[mutation_index] * seq1_region->plength_observation);
                                        
                                    double tot3 = 0;
                                    for (StateType j = 0; j < num_states; j++)
                                        if (seq2_region->likelihood[j] > 0.1)
                                            tot3 += model->mutation_mat[i * num_states + j];
                                    
                                    if (seq2_region->likelihood[i] > 0.1)
                                        tot += tot2 * (1.0 + blength15 * tot3);
                                    else
                                        tot += tot2 * blength15 * tot3;
                                }
                                
                                total_factor *= (tot / model->root_freqs[seq1_state]);
                            }
                        }
                        else
                        {
                            double tmp_blength = blength + (seq1_region->plength_observation < 0 ? 0 : seq1_region->plength_observation);
                            if (seq2_region->likelihood[seq1_state] > 0.1)
                                lh_cost += model->mutation_mat[seq1_state * (num_states + 1)] * tmp_blength;
                            else
                            {
                                for (StateType j = 0; j < num_states; j++)
                                    if (seq2_region->likelihood[j] > 0.1)
                                        tot += model->mutation_mat[seq1_state * num_states + j];
                                
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
                        
                        if (seq1_region->plength_from_root >= 0)
                        {
                            // here we ignore contribution of non-parsimonious mutational histories
                            total_factor *= ((model->root_freqs[seq1_state] * model->mutation_mat[seq1_state * num_states + seq2_state] * (blength + seq1_region->plength_from_root) * (1.0 + model->mutation_mat[seq1_state * (num_states + 1)] * seq1_region->plength_observation) + model->root_freqs[seq2_state] * model->mutation_mat[seq2_state * num_states + seq1_state] * seq1_region->plength_observation * (1.0 + model->mutation_mat[seq2_state * (num_states + 1)] * (blength + seq1_region->plength_from_root))) / model->root_freqs[seq1_state]);
                        }
                        else
                        {
                            double tmp_blength = blength + (seq1_region->plength_observation < 0 ? 0 : seq1_region->plength_observation);
                            
                            total_factor *= model->mutation_mat[seq1_state * num_states + seq2_state] * tmp_blength;
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
