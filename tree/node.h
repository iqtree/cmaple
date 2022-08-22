#include "alignment/regions.h"
#include "alignment/sequence.h"
#include "alignment/alignment.h"
#include "model/model.h"

#ifndef NODE_H
#define NODE_H

/** A node of the tree following the cyclic tree structure */
class Node {
public:
    /** node's id */
    int id;
    
    // sequence name
    string seq_name;
    
    /** likelihood of the subtree rooted at this node.
     It is computed from next->partial_lh, next->next->partial_lh and so on */
    Regions* partial_lh = NULL;
    
    /**
     combined likelihood vector at this node used to quickly calculate the placement likelihood.
     In essence, it is
     this->partial_lh combined with this->neighbor->partial_lh  extended by this branch length
    TODO: rename this
     */
    Regions* total_lh = NULL;
    
    /**
     combined likelihood vector at the midpoint of this node and its neighbor.
     It's used to quickly calculate the placement likelihood at the midpoint.
     In essence, it is
     this->partial_lh extended by length/2 combined with this->neighbor->partial_lh extended by length/2
     */
    Regions* mid_branch_lh = NULL;
    
    // length of branch connecting to parent
    RealNumType length;
    
    // next node in the circle of neighbors. For tips, circle = NULL
    Node* next;
    
    // next node in the phylo tree
    Node* neighbor;
    
    // vector of less informative sequences
    vector<string> less_info_seqs;
    
    // flag tp prevent traversing the same part of the tree multiple times.
    // TODO: change to computed or something similar
    bool dirty;
    
    /**
        TRUE if this node is the top node in a phylogenetic node structure
     */
    bool is_top;

    // flexible string attributes
    map<string,string> str_attributes;
    
    // flexible RealNumType attributes
    map<string,RealNumType> real_number_attributes;
    
    /**
        constructor
        @param n_seq_name the sequence name
     */
    Node(int n_id, string n_seq_name = "");
    
    /**
        constructor
        @param n_seq_name the sequence name
     */
    Node(string n_seq_name);
    
    /**
        constructor
        @param is_top_node TRUE if this is the first node in the next circle visited from root
     */
    Node(bool is_top_node = false);
    
    /**
        deconstructor
     */
    ~Node();
    
    /**
        TRUE if this node is a leaf
     */
    bool isLeave();
    
    /**
        get the top node
     */
    Node* getTopNode();
    
    /**
        export string: name + branch length
     */
    string exportString();
    
    /**
        get partial_lh of a node
     */
    Regions* getPartialLhAtNode(Alignment* aln, Model* model, RealNumType threshold_prob, RealNumType* cumulative_rate);
    
    /**
        increase the length of a 0-length branch (connecting this node to its parent) to resolve the inconsistency when updating regions in updatePartialLh()
     */
    void updateZeroBlength(stack<Node*> &node_stack, Alignment* aln, Model* model, RealNumType threshold_prob, RealNumType* cumulative_rate, RealNumType default_blength, RealNumType min_blength, RealNumType max_blength);
    
    /**
    *  compute the total likelihood vector for a node.
    */
    Regions* computeTotalLhAtNode(Alignment* aln, Model* model, RealNumType threshold_prob, RealNumType* cumulative_rate, bool is_root, bool update = true);
};

#endif

#define FOR_NEXT(node, next_node) \
for(next_node = node->next; next_node && next_node != node; next_node = next_node->next)

#define FOR_NEIGHBOR(node, neighbor_node) \
if (node->next) \
for(neighbor_node = node->next->neighbor; neighbor_node && neighbor_node != node->neighbor; neighbor_node = neighbor_node->neighbor->next->neighbor)
