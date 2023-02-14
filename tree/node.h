#include "alignment/seqregions.h"
#include "alignment/sequence.h"
#include "alignment/alignment.h"
#include "model/model.h"

#ifndef NODE_H
#define NODE_H

// ########### BEGIN OF NEW DATA STRUCTURES FOR PHYLOGENETIC NODES ###########
/** Index for a mininode of a phylonode */
enum MiniIndex : uint32_t
{
    TOP,
    LEFT,
    RIGHT
};

/** Holds a space efficient index into a vector<PhyloNode> and subindex for the MiniNode inside the Phylonode */
struct Index
{
    Index(uint32_t vec_index, MiniIndex mi)
    {
        setVectorIndex(vec_index);
        setMiniIndex(mi);
    }

    MiniIndex getMiniIndex() const
    {
        return MiniIndex(value_ >> 30);
    }
    
    void setMiniIndex(MiniIndex mi)
    {
        value_ = (value_ & keep_value) | (uint32_t(mi) << 30);
    }

    uint32_t getVectorIndex() const
    {
        // kill the first two bits (they belong to the mini index)
        return value_ & keep_value;
    }
    
    void setVectorIndex(uint32_t vecindex)
    {
        value_ = (value_ & keep_mini) | vecindex;
    }
    
    /** constructor */
    Index() = default;
    
    /** deconstructor */
    ~Index() = default;
    
private:
    // stores two ints: 2bit for mini index, and 30 bit for vector index
    // i.e.  MMVVVVV...VVVV
    uint32_t value_;
    // keep_value is: 0011111...111
    static constexpr uint32_t keep_value = uint32_t(~0) >> 2;
    // keep_mini is: 1100000...000
    static constexpr uint32_t keep_mini = 3 << 30;
};

/** An internal node of the tree containing 3 minonodes*/
struct InternalNode
{
    // There partial_lh(s) for the three mininodes, which represent the lower; upper left; upper right likelihoods
    std::array<std::unique_ptr<SeqRegions>,3> partial_lh;
    
    // We need to keep track of the index of the neighbor miniNode
    std::array<Index,3> neighbor_index;
    
    /** constructor */
    InternalNode() = default;
    
    /// Move CTor
    InternalNode(InternalNode&& internal) = default;
    
    /// Move Assignment
    InternalNode& operator=(InternalNode&& internal) = default;
    
    /** deconstructor */
    ~InternalNode()
    {
        std::cout << "~InternalNode()" << std::endl;
    };
}; // 40 (includes 4 byte padding at the end)

/** A leaf node of the tree*/
struct LeafNode
{
    // vector is wasteful here.. (24 bytes)
    std::vector<std::string> less_info_seqs; // leafs only (contains seq_names)
    // better use 2x4 bytes and store sequence names in a separate (global) vector
    //uint32_t index_start;
    //uint32_t index_end;   // but this requires some sorting... so lets not go there for now
    // a quick improvement would be to implement our own vector
    // which only takes 16 bytes (8 bytes for a pointer to heap memory + 4byte for size + 4 byte for capacity)

    /// however we can quickly optimize the seq_name to just a char* and store sequence names as a really long contatenated global string
    // std::string seq_name;  // too wasteful
    //char* seq_name;
    // and we can do even better when just storing an index into a combined string or a vector<string> (which would need to be provided externally)
    uint32_t seq_name_index;
    
    //MiniNode leaf; // this leads to padding... so better list the members individually
    // .. and in the right order
    Index neighbor_index;
    std::unique_ptr<SeqRegions> partial_lh;
    
    /** constructor */
    LeafNode() = default;
    
    /// Move CTor
    LeafNode(LeafNode&& leaf) = default;
    
    /// Move Assignment
    LeafNode& operator=(LeafNode&& leaf) = default;
    
    /** deconstructor */
    ~LeafNode()
    {
        std::cout << "~LeafNode()" << std::endl;
    };
}; // size: 40 bytes. no padding! :) -- we could bring this down to 32 bytes if need be

struct OtherLh
{
    // lh[0] ~ total_lh
    // lh[1] ~ mid_branch_lh;
    SeqRegions lh[2];
    
    /** constructor */
    OtherLh() = default;
};

/** An intermediate data structure to store either InternalNode or LeafNode */
union MyVariant
{
    InternalNode internal;
    LeafNode leaf;
    
    /** constructor */
    MyVariant():internal() {};
    
    /** constructor */
    MyVariant(LeafNode&& leaf_):leaf(std::move(leaf_)) {};
    
    /** constructor */
    MyVariant(InternalNode&& internal_):internal(std::move(internal_)) {};
    
    /** deconstructor */
    ~MyVariant() {};
    
    /** switch from  InternalNode to LeafNode */
    void switch2LeafNode();
    
    /** switch from  LeafNode to InternalNode */
    void switch2InternalNode();
};

/** An node in a phylogenetic tree, which could be either an internal or a leaf */
class PhyloNode
{
    std::unique_ptr<OtherLh> other_lh_ = std::make_unique<OtherLh>(OtherLh());
    
    const bool is_internal_; // is our MyVariant an internal node or a leaf?
    
    /// common members for Internal and Leaf nodes
    bool outdated_;
    float length_; // using float allows it to fit into the 6 bytes padding after the two bools.
    // .. using a double would make PyhloNode 8 bytes larger

    MyVariant data_;    
  public:
    /** constructor */
    PhyloNode()
    {
      // could also be part of a separate class test, but we need to make sure that this is enforced and noone accidentally
      // changes it
      std::static_assert(sizeof(PhyloNode) == 64);  // make sure it fits on a cacheline
    };
    
    /** constructor */
    PhyloNode(LeafNode&& leaf): is_internal_{false}, data_(std::move(leaf))
    {};
    
    /** constructor */
    PhyloNode(InternalNode&& internal): is_internal_{true}, data(std::move(internal))
    {};
    
    /** destructor */
    ~PhyloNode()
    {
        std::cout << "~PhyloNode()" << std::endl;
        
        // delete MyVariant (data)
        if (is_internal_)
            data.internal.~InternalNode();
        else
            data.leaf.~LeafNode();
    };
    
    /**
        Get total_lh
     */
    SeqRegions& getTotalLh();
    
    /**
        Set total_lh
     */
    void setTotalLh(SeqRegions&& total_lh);
    
    /**
        Get mid_branch_lh
     */
    SeqRegions& getMidBranchLh();
    
    /**
        Set mid_branch_lh
     */
    void setMidBranchLh(SeqRegions&& mid_branch_lh);
    
    /**
        Get partial_lh
     */
    SeqRegions& getPartialLh(const Index index);
    
    /**
        Set partial_lh
     */
    void setPartialLh(const Index index, SeqRegions&& partial_lh_);
    
    /**
        Get the index of the neighbor node
     */
    Index getNeighborIndex(const Index index) const;
    
    /**
        Set the index of the neighbor node
     */
    void setNeighborIndex(const Index index, const Index neighbor_index_);
    
    /**
        Get the list of less-informative-sequences
     */
    std::vector<std::string>& getLessInfoSeqs();
    
    /**
        Add less-informative-sequence
     */
    void addLessInfoSeqs(std::string&& seq_name);
    
    /**
        Get the index of the sequence name
     */
    uint32_t getSeqNameIndex() const;
    
    /**
        Set the index of the sequence name
     */
    void setSeqNameIndex(uint32_t seq_name_index_);
    
    /**
        Update MyVariant to LeafNode
     */
    void updateMyVariant(LeafNode&& leaf);
    
    /**
        Update MyVariant to InternalNode
     */
    void updateMyVariant(InternalNode&& internal);

};

// just for testing
/** operator<< */
std::ostream& operator<<(std::ostream& os, const Index& index);

// ########### END OF NEW DATA STRUCTURES FOR PHYLOGENETIC NODES ###########

/** A node of the tree following the cyclic tree structure */
class Node {
public:
    /**
    Flag to prevent traversing the same part of the tree multiple times.
    */
    bool outdated;

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
    double length;
    
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
