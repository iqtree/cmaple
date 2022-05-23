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
        get overall likelihood sequence at root node that take into account the root frequencies
        @param input_regions the input sequence (vector of regions); blength the branch length; threshold_prob4 the threshold for approximation
    */
    vector<Region*> getOverallLhSeqAtRoot(vector<Region*> input_regions, double blength, double threshold_prob4);
    
    /**
    *  get lower likelihood sequence at tip (converting a vector of Mutations into a vector of Regions)
    */
    vector<Region*> getLowerLhSeqAtTip(vector<Mutation*> mutations, double blength);
    
    /**
    *  add a new node to become a sibling ò the current node
    */
    void addNode(Node* new_node, Node* current_node, int &new_id, bool insert_to_right = true);

};

#endif
