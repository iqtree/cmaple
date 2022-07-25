#include "node.h"
#include "alignment/alignment.h"
#include "model/model.h"

#ifndef TREE_H
#define TREE_H

class Tree {
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
    void updatePartialLh(stack<Node*> &node_stack, double* cumulative_rate, double default_blength, double min_blength, double max_blength);
    
    /**
           seek a position for a sample placement starting at the start_node
     */
    void seekPlacement(Node* start_node, string seq_name, Regions* sample_regions, Node* &selected_node, double &best_lh_diff , bool &is_mid_branch, double &best_up_lh_diff, double &best_down_lh_diff, Node* &best_child, double* cumulative_rate, double default_blength, double min_blength_mid);

    /**
           place a new sample on the tree
     */
    void placeNewSample(Node* selected_node, Regions* sample, string seq_name, double best_lh_diff , bool is_mid_branch, double best_up_lh_diff, double best_down_lh_diff, Node* best_child, double* cumulative_rate, vector<vector<PositionType>> &cumulative_base, double default_blength, double max_blength, double min_blength);
};

#endif
