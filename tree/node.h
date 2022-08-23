#include "alignment/seqregions.h"
#include "alignment/sequence.h"
#include "alignment/alignment.h"
#include "model/model.h"

#ifndef NODE_H
#define NODE_H

/** A node of the tree following the cyclic tree structure */
class Node {
public:
    /**
        Node's id
    */
    int id;
    
    /**
        Sequence name
     */
    string seq_name;
    
    /**
        Likelihood of the subtree rooted at this node.
        It is computed from next->partial_lh, next->next->partial_lh and so on
     */
    SeqRegions* partial_lh = NULL;
    
    /**
        Combined likelihood vector at this node used to quickly calculate the placement likelihood.
        In essence, it is
        this->partial_lh combined with this->neighbor->partial_lh  extended by this branch length
     */
    SeqRegions* total_lh = NULL;
    
    /**
        Combined likelihood vector at the midpoint of this node and its neighbor.
        It's used to quickly calculate the placement likelihood at the midpoint.
        In essence, it is
        this->partial_lh extended by length/2 combined with this->neighbor->partial_lh extended by length/2
     */
    SeqRegions* mid_branch_lh = NULL;
    
    /**
        Length of branch connecting to parent
     */
    RealNumType length;
    
    /**
        Next node in the circle of neighbors. For tips, circle = NULL
     */
    Node* next;
    
    /**
        Neighbor node in the phylo tree
     */
    Node* neighbor;
    
    /**
        Vector of sequences that are less informative than the sequence of this node
     */
    vector<string> less_info_seqs;
    
    /**
        Flag to prevent traversing the same part of the tree multiple times.
     */
    bool outdated;
    
    /**
        TRUE if this node is the top node in a phylogenetic node structure
     */
    bool is_top;

    /**
        Flexible string attributes
     
     */
    map<string,string> str_attributes;
    
    /**
        Constructor
        @param n_seq_name the sequence name
     */
    Node(int n_id, string n_seq_name = "");
    
    /**
        Constructor
        @param n_seq_name the sequence name
     */
    Node(string n_seq_name);
    
    /**
        Constructor
        @param is_top_node TRUE if this is the first node in the next circle visited from root
     */
    Node(bool is_top_node = false);
    
    /**
        Deconstructor
     */
    ~Node();
    
    /**
        TRUE if this node is a leaf
     */
    bool isLeave();
    
    /**
        Get the top node of a phylo-node
     */
    Node* getTopNode();
    
    /**
        Export string: name + branch length
     */
    string exportString();
    
    /**
        Get/(or compute) partial_lh of a node
     */
    SeqRegions* getPartialLhAtNode(Alignment* aln, Model* model, RealNumType threshold_prob, RealNumType* cumulative_rate);
    
    /**
        Increase the length of a 0-length branch (connecting this node to its parent) to resolve the inconsistency when updating regions in updatePartialLh()
     */
    void updateZeroBlength(stack<Node*> &node_stack, Alignment* aln, Model* model, RealNumType threshold_prob, RealNumType* cumulative_rate, RealNumType default_blength, RealNumType min_blength, RealNumType max_blength);
    
    /**
        Compute the total likelihood vector for a node.
    */
    SeqRegions* computeTotalLhAtNode(Alignment* aln, Model* model, RealNumType threshold_prob, RealNumType* cumulative_rate, bool is_root, bool update = true);
};

/** An extension of node storing more dummy data used for browsing all nodes in a stack  */
class ExtendedNode {
public:
    /**
        Pointer to a node
     */
    Node* node;
    
    /**
        Count of the number of failures when traversing until the current node
     */
    short int failure_count;
    
    /**
        Cache the likelihood difference computed at the parent node
     */
    RealNumType likelihood_diff;
    
    /**
        Constructor
     */
    ExtendedNode();
    
    /**
        Constructor
     */
    ExtendedNode(Node* n_node, short int n_failure_count, RealNumType n_lh_diff);
    
    /**
        Deconstructor
     */
    ~ExtendedNode();
};

#endif

#define FOR_NEXT(node, next_node) \
for(next_node = node->next; next_node && next_node != node; next_node = next_node->next)

#define FOR_NEIGHBOR(node, neighbor_node) \
if (node->next) \
for(neighbor_node = node->next->neighbor; neighbor_node && neighbor_node != node->neighbor; neighbor_node = neighbor_node->neighbor->next->neighbor)
