//
//  regions.h
//  alignment
//
//  Created by Nhan Ly-Trong on 24/01/2022.
//

#include "region.h"

#ifndef REGIONS_H
#define REGIONS_H

class Regions: public vector<Region*> {
private:
    /**
    *  Delete all region items one by one
    */
    void deleteRegions();
    
public:
    /**
    *  Regions constructor
    */
    Regions();
    
    /**
    *  Regions deconstructor
    */
    ~Regions();
    
    /**
    *  copy regions
    */
    void copyRegions(Regions* n_regions, StateType num_states);
    
    /**
    *  get a Region
    */
    Region* getRegion(PositionType id);
};
#endif
