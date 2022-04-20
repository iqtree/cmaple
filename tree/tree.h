#include "node.h"

#ifndef TREE_H
#define TREE_H

class Tree {
public:
    // root node of the tree
    Node* root;
    
    /**
        constructor
     */
    Tree();
    
    /**
        constructor
        @param node the root node
     */
    Tree(Node* node);
    
    /**
        deconstructor
     */
    ~Tree();

};

#endif
