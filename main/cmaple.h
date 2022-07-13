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
            calculate the placement cost
            @param parent_regions, child_regions: vector of regions of the parent and the new sample
     */
    double calculatePlacementCost(Regions* parent_regions, Regions* child_regions, double blength);
    
    /**
           seek a position for a sample placement
     */
    void seekPlacement(Node* start_node, string seq_name, Regions* regions, Node* &selected_node, double &best_lh_diff , bool &is_mid_branch, double &best_up_lh_diff, double &best_down_lh_diff, Node* &best_child, bool &adjust_blen);
    
    /**
           place a new sample on the tree
     */
    void placeNewSample(Node* selected_node, Regions* sample, string seq_name, double best_lh_diff , bool is_mid_branch, double best_up_lh_diff, double best_down_lh_diff, Node* best_child, bool adjust_blen);
    
    /**
        compare two sequences regarding the amount of information
        @param sequence1, sequence2
        @return 0: if the two sequences are incomparable; 1: if sequence1 is more or equally informative than/to sequence2; -1; if sequence1 is less informative than sequence2
     */
    int compareSequences(Regions* sequence1, Regions* sequence2);
    
    /**
        move to the next region in the vector of regions
        @param sequence: a vector of regions; region_index: the index of the next region
        @return region: a region; current_pos: the current position; end_pos: the end position
     */
    void move2NextRegion(Regions* sequence, PositionType region_index, Region* &region, PositionType &current_pos, PositionType &end_pos);
    
    /**
        get the shared segment between the next regions of two sequences
        @param current_pos: current site posisition; sequence1, sequence2: vectors of regions; seq1_index, seq2_index: the indexes of the current regions; seq1_end_pos, seq2_end_pos: the end positions of the current regions
        @return seq1_region, seq2_region: the regions contains the shared segment; length: length of the shared segment
     */
    void getNextSharedSegment(PositionType current_pos, Regions* sequence1, Regions* sequence2, PositionType &seq1_index, PositionType &seq2_index, Region* &seq1_region, Region* &seq2_region, PositionType &seq1_end_pos, PositionType &seq2_end_pos, PositionType &length);
    
    /**
        merge two likelihood vectors, one from above and one from below
     */
    void mergeLhUpDown(Regions* &merged_regions, Regions* upper_regions, double upper_plength, Regions* lower_regions, double lower_plength);
    
    /**
        merge two lower likelihood vectors
     */
    double mergeLhTwoLower(Regions* &merged_regions, Regions* regions1, double plength1, Regions* regions2, double plength2, bool return_log_lh = false);
    
    /**
    *  get lower likelihood sequence at tip (converting a vector of Mutations into a vector of Regions)
    */
    void getLowerLhSeqAtTip(Regions* &output_regions, vector<Mutation*> mutations, double blength);
    
    /**
        compute total lh/upper left_right for root node
     */
    Regions* computeTotalLhForRoot(Regions* lower_regions, double blength = 0);
    
    /**
    *  calculate the total likelihood regions for a node.
    */
    Regions* computeTotalLhAtNode(Node* node, bool update = true);
    
    /**
        convert an entry 'O' into a normal nucleotide if its probability dominated others
     */
    StateType simplify(double* &partial_lh, StateType ref_state);
    
    /**
        merge consecutive R regions
     */
    void mergeRegionR(Regions* &regions);
    
    /**
        compute the likelihood by merging the lower lh with root frequencies
     */
    double computeLhAtRoot(Regions* region);
    
    /**
        get partial_lh of a node
     */
    Regions* getPartialLhAtNode(Node* node);
    
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
        Check if two regions represent the same partial likelihoods or not -> be used to stop traversing the tree further for updating partial likelihoods
     */
    bool areDiffRegions(Regions* regions1, Regions* regions2);
    
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
