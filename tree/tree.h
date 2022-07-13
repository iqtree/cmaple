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

};

#endif
