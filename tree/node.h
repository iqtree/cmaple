#include "alignment/region.h"
#include "alignment/sequence.h"

#ifndef NODE_H
#define NODE_H

class Node {
public:
    // node's id
    int id;
    
    // sequence name
    string seq_name;
    
    // sequence
    vector<Region*> sequence;
    
    // length of branch connecting to parent
    double length;
    
    // parent node
    Node* parent;
    
    // vector of children: 0 - left; 1 - right; >1 : others
    vector<Node*> children;

    // flexible attributes
    map<string,string> attributes;
    
    /**
        constructor
     */
    Node();
    
    /**
        constructor
        @param n_seq_name sequence name of this node
        @param n_sequence the sequence (vector of mutations)
     */
    Node(string n_seq_name, Sequence sequence);

};

#endif
