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
