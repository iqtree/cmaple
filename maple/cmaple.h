/*
 *  cmaple.h
 *
 *  Created on: Mar 13, 2021
 *      Author: Nhan Ly-Trong
 */

#include "tree/tree.h"
#include "utils/timeutil.h"
#include "alignment/seqregions.h"

#pragma once

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
    void exportOutput(const std::string &filename);
    
    /**
        calculate branch supports
     */
    void calculateBranchSupports();
    
public:
    
    /** The phylogenetic tree */
    Tree tree;
    
    /**
    *  CMaple constructor
    */
    CMaple() = default;
    
    /**
    *  CMaple constructor
     @param params user-specified parameters
    */
    CMaple(Params params):tree(std::move(params)){};
    
    /**
    *  CMaple destructor
    */
    ~CMaple() = default;
    
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
};

/** Method to run CMaple */
void runCMaple(Params &params);
