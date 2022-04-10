/*
 *  cmaple.h
 *
 *  Created on: Mar 13, 2021
 *      Author: Nhan Ly-Trong
 */


#include "alignment/alignment.h"

#ifndef CMAPLE_H
#define CMAPLE_H

class CMaple {
private:
   
    
public:
    Params *params;
    Alignment *aln;
    
    /**
    *  CMaple constructor
    */
    CMaple();
    
    /**
    *  CMaple constructor
    */
    CMaple(Params *params);
    
    /**
    *  CMaple deconstructor
    */
    ~CMaple();
    
    /**
            Extract Diff file from and alignment file
            @param NULL
     */
    void extractDiffFile();
};

// Run CMaple
void runCMaple(Params &params);
#endif
