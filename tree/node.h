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
    TRUE if the likelihoods were updated due to some SPR moves
    */
    bool outdated;
    
    /**
    TRUE if we applied SPR move on the upper branch of this node
    */
    bool SPR_applied;
    
    /**
        TRUE if this node is the top node in a phylogenetic node structure
     */
    bool is_top;

    
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
    std::vector<std::string> less_info_seqs;
    
    /**
        Sequence name
     */
    std::string seq_name;
    /**
        Flexible string attributes
     
     */
    //std::map<std::string, std::string> str_attributes;
    
    /**
        Constructor
        @param n_seq_name the sequence name
     */
    Node(int n_id, std::string n_seq_name = "");
    
    /**
        Constructor
        @param n_seq_name the sequence name
     */
    Node(std::string n_seq_name);
    
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
        Get the other next node (not the top)
     */
    Node* getOtherNextNode();
    
    /**
        Export string: name + branch length
     */
    std::string exportString(bool binary = false);
    
    /**
        Get/(or compute) partial_lh of a node
     */
    SeqRegions* getPartialLhAtNode(const Alignment& aln, const Model& model, RealNumType threshold_prob);
    
    /**
        Compute the total likelihood vector for a node.
    */
    SeqRegions* computeTotalLhAtNode(const Alignment& aln, const Model& model, RealNumType threshold_prob, bool is_root, bool update = true, RealNumType blength = -1);
};

/** An extension of node storing more dummy data used for browsing all nodes in a stack  */
class TraversingNode {
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
    TraversingNode();
    
    /**
        Constructor
     */
    TraversingNode(Node* n_node, short int n_failure_count, RealNumType n_lh_diff);
    
    /**
        Deconstructor
     */
    ~TraversingNode();
};

/** An extension of node storing more dummy data used for updating nodes in the tree  */
class UpdatingNode: public TraversingNode {
public:
    /**
        an updated regions from the direction where we come from (taking into account the removal of the given subtree),
     */
    SeqRegions* incoming_regions;
    
    /**
        a branch length separating the node from this updated regions (useful for the fact that the removal of the subtree changes the branch length at the removal node)
     */
    RealNumType branch_length;
    
    /**
        a flag says if the updated regions passed needs still updating, or if it has become identical to the pre-existing genome list in the tree (which usually happens after a while)
     */
    bool need_updating;
    
    /**
        TRUE to delete incoming_regions with the deconstructor
     */
    bool delete_regions;
    
    /**
        Constructor
     */
    UpdatingNode();
    
    /**
        Constructor
     */
    UpdatingNode(Node* n_node, SeqRegions* n_incoming_regions, RealNumType n_branch_length, bool n_need_updating, RealNumType n_lh_diff, short int n_failure_count, bool n_delete_regions);
    
    /**
        Deconstructor
     */
    ~UpdatingNode();
};

#endif

#define FOR_NEXT(node, next_node) \
for(next_node = node->next; next_node && next_node != node; next_node = next_node->next)

#define FOR_NEIGHBOR(node, neighbor_node) \
if (node->next) \
for(neighbor_node = node->next->neighbor; neighbor_node && neighbor_node != node->neighbor; neighbor_node = neighbor_node->neighbor->next->neighbor)
