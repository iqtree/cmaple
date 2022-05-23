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
    
    /**
        Add the root node to the tree
     */
    void addRoot();
    
    /**
        Update the mutation matrix periodically from the empirical count of mutations
     */
    void updateMutationMatFromEmpiricalCount();
    
    /**
           seek a position for a sample placement
     */
    void seekPlacement(Node* start_node, string seq_name, vector<Region*> regions, Node* &selected_node, double &best_lh_diff , bool &is_mid_branch, double &best_up_lh_diff, double &best_down_lh_diff, Node* &best_child, bool &adjust_blen);
    
    /**
           examine a position for placing a new sample
     */
    void examinePosition(Node* node, double parent_lh, int failure_count, string seq_name, vector<Region*> regions, Node* &selected_node, double &best_lh_diff , bool &is_mid_branch, double &best_up_lh_diff, double &best_down_lh_diff, Node* &best_child, bool &adjust_blen);
    
    /**
        Traverse tree from a node (node = root by default)
     */
    void traverseTree(Node* node = NULL);
    
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
