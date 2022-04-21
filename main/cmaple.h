/*
 *  cmaple.h
 *
 *  Created on: Mar 13, 2021
 *      Author: Nhan Ly-Trong
 */

#include "tree/tree.h"
#include "utils/timeutil.h"

#ifndef CMAPLE_H
#define CMAPLE_H

class CMaple {
private:
    
    /**
        Extract Diff file from and alignment file
        @param NULL
     */
    void extractDiffFile();
    
    /**
        compute cumulative rate of the ref genome
     */
    void computeCumulativeRate();
    
    /**
        Check to redo the inference
     */
    void checkRedoInference();
    
    /**
        Compute the thresholds for approximations
     */
    void computeThresholds();
    
    /**
        Build an Initial Tree
     */
    void buildInitialTree();
    
public:
    // model-related params
    double *cumulative_rate;
    double* pseu_mutation_count;
    
    // The phylogenetic tree
    Tree* tree;
    // tree-related params
    double default_blength;
    
    // thresholds for approximations
    double threshold_prob2, threshold_prob4;
    
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
        Temporarily method for testing
     */
    void tmpTestingMethod();
};

// Run CMaple
void runCMaple(Params &params);
#endif
