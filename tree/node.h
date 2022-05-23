#include "alignment/region.h"
#include "alignment/sequence.h"
#include "alignment/alignment.h"

#ifndef NODE_H
#define NODE_H

class Node {
public:
    // node's id
    int id;
    
    // sequence name
    string seq_name;
    
    // overall likelihood
    vector<Region*> overall_lh;
    
    // overall likelihood mid branch
    vector<Region*> overall_lh_mid_branch;
    
    // likelihood from other nodes
    vector<Region*> lh_from_others;
    
    // lower likelihood sequence
    vector<Region*> lower_lh_seq;
    
    // length of branch connecting to parent
    double length;
    
    // next node in the circle of neighbors. For tips, circle = NULL
    Node* relative;
    
    // next node in the phylo tree
    Node* neighbor;
    
    // vector of less informative sequences
    vector<string> less_info_seqs;

    // flexible attributes
    map<string,string> attributes;
    
    /**
        constructor
     */
    Node();
    
    /**
        constructor
        @param n_seq_name the sequence name
     */
    Node(int n_id = -1, string n_seq_name = "");
    
    /**
        deconstructor
     */
    ~Node();
    
    /**
        TRUE if this node is a leaf
     */
    bool isLeave();
};

#endif

#define FOR_RELATIVE(node, relative_node) \
for(relative_node = node->relative; relative_node && relative_node != node; relative_node = relative_node->relative)
