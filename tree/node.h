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
    
    // lower likelihood sequence
    vector<Region*> lower_lh_seq;
    
    // overall likelihood sequence
    vector<Region*> overall_lh_seq;
    
    // length of branch connecting to parent
    double length;
    
    // next node in the circle of neighbors. For tips, circle = NULL
    Node* circle;
    
    // next node in the phylo tree
    Node* next;
    
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
        @param sequence the sequence
        @param aln the alignment
     */
    Node(Sequence* sequence, Alignment* aln);
    
    /**
        deconstructor
     */
    ~Node();
};

#endif
