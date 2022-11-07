/*
 *  cmaple.h
 *
 *  Created on: Mar 13, 2021
 *      Author: Nhan Ly-Trong
 */

#include "tree/tree.h"
#include "utils/timeutil.h"
#include "alignment/seqregions.h"

#ifndef CMAPLE_H
#define CMAPLE_H

/** CMaple class contains all main methods*/
class CMaple {
private:
    
    /**
        Build an Initial Tree
     */
    void buildInitialTree();
    
    /**
        Optimize the current tree
     */
    void optimizeTree();
    
    /**
        Optimize the tree topology
     */
    void optimizeTreeTopology(bool short_range_search = false);
    
    /**
        Optimize the branch lengths of the current tree
     */
    void optimizeBranchLengthsOfTree();
    
    /**
        export output files
     */
    void exportOutput(std::string filename);
    
public:
    /** Model-related parameters */
    RealNumType *cumulative_rate;
    std::vector< std::vector<PositionType> > cumulative_base;
    
    /** The phylogenetic tree */
    Tree* tree;
    
    /** tree-related params */
    RealNumType default_blength;
    RealNumType min_blength, max_blength, min_blength_mid, min_blength_sensitivity;
    
    /**
    *  CMaple constructor
    */
    CMaple();
    
    /**
    *  CMaple constructor
     @param params user-specified parameters
    */
    CMaple(Params *params);
    
    /**
    *  CMaple deconstructor
    */
    ~CMaple();
    
    /**
        Load input data
     */
    void loadInput();
    
    /**
        Prepare for the inference
     */
    void preInference();
    
    /**
        Do the inference
     */
    void doInference();
    
    /**
        Complete the inference
     */
    void postInference();
    
    /**
        Temporarily method for testing
     */
    void tmpTestingMethod();
};

/** Method to run CMaple */
void runCMaple(Params &params);
#endif
