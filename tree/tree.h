#include "node.h"
#include "alignment/alignment.h"
#include "model/model.h"

#ifndef TREE_H
#define TREE_H

class Tree {
private:
    /**
           traverse the intial tree from root to re-calculate all lower likelihoods regarding the latest/final estimated model parameters
     */
    void refreshAllLowerLhs(RealNumType *cumulative_rate, RealNumType default_blength, RealNumType max_blength, RealNumType min_blength);
    
    /**
           traverse the intial tree from root to re-calculate all non-lower likelihoods regarding the latest/final estimated model parameters
     */
    void refreshAllNonLowerLhs(RealNumType *cumulative_rate, RealNumType default_blength, RealNumType max_blength, RealNumType min_blength);
    
    /**
        try to improve a subtree rooted at node with SPR moves
        @return total improvement
     */
    RealNumType improveSubTree(Node* node, RealNumType *cumulative_rate, RealNumType default_blength, RealNumType max_blength, RealNumType min_blength);
    
public:
    // parameters
    Params *params;
    
    // alignment
    Alignment *aln;
    
    // evolutionary model
    Model* model;
    
    // root node of the tree
    Node* root;
    
    /**
        constructor
     */
    Tree();
    
    /**
    *  constructor
    */
    Tree(Params *params, Node* n_root = NULL);
    
    /**
        deconstructor
     */
    ~Tree();
    
    /**
    *  add a new node to become a sibling ò the current node
    */
    void addNode(Node* new_node, Node* current_node, int &new_id, bool insert_to_right = true);
    
    /**
        export tree string
     */
    string exportTreeString(Node* node = NULL);
    
    /**
        iteratively update partial_lh starting from the nodes in node_stack
        @param node_stack stack of nodes;
     */
    void updatePartialLh(stack<Node*> &node_stack, RealNumType* cumulative_rate, RealNumType default_blength, RealNumType min_blength, RealNumType max_blength);
    
    /**
        seek a position for a sample placement starting at the start_node
     */
    void seekSamplePlacement(Node* start_node, string seq_name, SeqRegions* sample_regions, Node* &selected_node, RealNumType &best_lh_diff , bool &is_mid_branch, RealNumType &best_up_lh_diff, RealNumType &best_down_lh_diff, Node* &best_child, RealNumType* cumulative_rate, RealNumType default_blength, RealNumType min_blength_mid);
    
    /**
        seek a position for placing a subtree/sample starting at the start_node
     */
    void seekPlacement(Node* &best_node, RealNumType &best_lh_diff, bool &is_mid_branch, RealNumType &best_up_lh_diff, RealNumType &best_down_lh_diff, Node* &best_child, Node* start_node, RealNumType &removed_blength, RealNumType* cumulative_rate, RealNumType default_blength, RealNumType min_blength_mid, bool search_subtree_placement = true, SeqRegions* sample_regions = NULL);

    /**
        place a new sample on the tree
     */
    void placeNewSample(Node* selected_node, SeqRegions* sample, string seq_name, RealNumType best_lh_diff , bool is_mid_branch, RealNumType best_up_lh_diff, RealNumType best_down_lh_diff, Node* best_child, RealNumType* cumulative_rate, vector< vector<PositionType> > &cumulative_base, RealNumType default_blength, RealNumType max_blength, RealNumType min_blength);
    
    /**
        traverse the intial tree from root to re-calculate all likelihoods regarding the latest/final estimated model parameters
     */
    void refreshAllLhs(RealNumType *cumulative_rate, RealNumType default_blength, RealNumType max_blength, RealNumType min_blength);
    
    /**
        traverse the tree to set the dirty flag (which is used to prevent traversing the same part of the tree multiple times) of all nodes to TRUE
     */
    void setAllNodeDirty();
    
    /**
        try to improve the entire tree with SPR moves
        @return total improvement
     */
    RealNumType improveEntireTree(RealNumType *cumulative_rate, RealNumType default_blength, RealNumType max_blength, RealNumType min_blength);
};

#endif
