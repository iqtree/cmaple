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
}

/**
*  get a Region
*/
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
