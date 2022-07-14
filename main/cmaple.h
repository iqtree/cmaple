/*
 *  cmaple.h
 *
 *  Created on: Mar 13, 2021
 *      Author: Nhan Ly-Trong
 */

#include "tree/tree.h"
#include "utils/timeutil.h"
#include <queue>
#include "alignment/regions.h"

#ifndef CMAPLE_H
#define CMAPLE_H

class CMaple {
private:
    
    /**
        Extract Diff file from and alignment file
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
        Update the mutation matrix periodically from the empirical count of mutations
     */
    void updateMutationMatEmpirical();
    
    /**
        update pseudocounts from new sample to improve the estimate of the substitution rates
        @param node_regions the genome list at the node where the appending happens;
        @param sample_regions the genome list for the new sample.
     */
    void updatePesudoCount(Regions* node_regions, Regions* sample_regions);
    
    /**
        iteratively update partial_lh starting from the nodes in node_queue
        @param node_queue queue of nodes;
     */
    void updatePartialLh(queue<Node*> &node_queue);
    
    /**
            calculate the placement cost
            @param parent_regions, child_regions: vector of regions of the parent and the new sample
     */
    double calculatePlacementCost(Regions* parent_regions, Regions* child_regions, double blength);
    
    /**
           seek a position for a sample placement starting at the start_node
     */
    void seekPlacement(Node* start_node, string seq_name, Regions* regions, Node* &selected_node, double &best_lh_diff , bool &is_mid_branch, double &best_up_lh_diff, double &best_down_lh_diff, Node* &best_child);
    
    /**
           place a new sample on the tree
     */
    void placeNewSample(Node* selected_node, Regions* sample, string seq_name, double best_lh_diff , bool is_mid_branch, double best_up_lh_diff, double best_down_lh_diff, Node* best_child);
    
    /**
        merge two likelihood vectors, one from above and one from below
     */
    void mergeUpperLower(Regions* &merged_regions, Regions* upper_regions, double upper_plength, Regions* lower_regions, double lower_plength);
    
    /**
        merge two lower likelihood vectors
     */
    double mergeTwoLowers(Regions* &merged_regions, Regions* regions1, double plength1, Regions* regions2, double plength2, bool return_log_lh = false);
    
    /**
        compute total lh/upper left_right for root node
        @param lower_regions the lower lh vector
        @param blength the branch length; (-1 by default).
     */
    Regions* computeTotalLhAtRoot(Regions* lower_regions, double blength = -1);
    
    /**
    *  compute the total likelihood vector for a node.
    */
    Regions* computeTotalLhAtNode(Node* node, bool update = true);
    
    /**
        compute the likelihood by merging the lower lh with root frequencies
     */
    double computeAbsoluteLhAtRoot(Regions* region);
    
    /**
        get partial_lh of a node
     */
    Regions* getPartialLhAtNode(Node* node);
    
    /**
        convert an entry 'O' into a normal nucleotide if its probability dominated others
     */
    StateType simplifyO(double* &partial_lh, StateType ref_state);
    
    /**
        increase the length of a 0-length branch to resolve the inconsistency when updating regions in updatePartialLh()
     */
    void updateZeroBlength(queue<Node*> &node_queue, Node* node);
    
public:
    // model-related params
    double *cumulative_rate;
    vector<vector<PositionType>> cumulative_base;
    double* pseu_mutation_count;
    
    // The phylogenetic tree
    Tree* tree;
    // tree-related params
    double default_blength;
    double min_blength, max_blength;
    
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
