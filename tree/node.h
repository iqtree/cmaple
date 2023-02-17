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
    RIGHT,
    UNDEFINED // to handle case return node = null
};

/** Holds a space efficient index into a vector<PhyloNode> and subindex for the MiniNode inside the Phylonode */
struct Index
{
    /**
        constructor
        it's required by the constructor of LeafNode
     */
    Index() = default;
    
    /**
        constructor
     */
    Index(uint32_t vec_index, MiniIndex mini_index)
    {
        setVectorIndex(vec_index);
        setMiniIndex(mini_index);
    }

    /**
        return the mini-index
     */
    MiniIndex getMiniIndex() const
    {
        return mini_index_;
    }
    
    /**
        return the mini-index of the other next node (i.e., getFlipMiniIndex(LEFT) = RIGHT; getFlipMiniIndex(RIGHT) = LEFT)
     */
    MiniIndex getFlipMiniIndex() const
    {
        return MiniIndex(3 - mini_index_);
    }
    
    /**
        set the mini-index
     */
    void setMiniIndex(MiniIndex mini_index)
    {
        mini_index_ = mini_index;
    }

    /**
        return the index of a phylonode
     */
    uint32_t getVectorIndex() const
    {
        return vector_index_;
    }
    
    /**
        set the index of a phylonode
     */
    void setVectorIndex(uint32_t vec_index)
    {
        vector_index_ = vec_index;
    }
    
private:
    uint32_t vector_index_ : 30;
    MiniIndex mini_index_ : 2;
};

/** An internal node of the tree containing 3 minonodes*/
struct InternalNode
{
    // There partial_lh(s) for the three mininodes, which represent the lower; upper left; upper right likelihoods
    std::array<std::unique_ptr<SeqRegions>,3> partial_lh3_;
    
    // We need to keep track of the index of the neighbor miniNode
    std::array<Index,3> neighbor_index3_;
    
}; // 40 (includes 4 byte padding at the end)

/** A leaf node of the tree*/
struct LeafNode
{
    // vector is wasteful here.. (24 bytes)
    std::vector<uint32_t> less_info_seqs_; // leafs only (contains seq_names)
    // better use 2x4 bytes and store sequence names in a separate (global) vector
    //uint32_t index_start;
    //uint32_t index_end;   // but this requires some sorting... so lets not go there for now
    // a quick improvement would be to implement our own vector
    // which only takes 16 bytes (8 bytes for a pointer to heap memory + 4byte for size + 4 byte for capacity)

    /// however we can quickly optimize the seq_name to just a char* and store sequence names as a really long contatenated global string
    // and we can do even better when just storing an index into a combined string or a vector<string> (which would need to be provided externally)
    uint32_t seq_name_index_;
    
    // index to connect to its neighbor node
    Index neighbor_index_;
    
    // The partial_lh (lower likelihood)
    std::unique_ptr<SeqRegions> partial_lh_;
  
    /** constructor */
    LeafNode() {};
    
    /** constructor */
    LeafNode(uint32_t new_seq_name_index):seq_name_index_(new_seq_name_index) {};
}; // size: 40 bytes. no padding! :) -- we could bring this down to 32 bytes if need be

/** An node in a phylogenetic tree, which could be either an internal or a leaf */
class PhyloNode
{
private:
    // store the total_lh and mid_branch_lh
    struct OtherLh
    {
        std::unique_ptr<SeqRegions> total_lh;
        std::unique_ptr<SeqRegions> mid_branch_lh;
    };
    
    /** An intermediate data structure to store either InternalNode or LeafNode */
    union MyVariant
    {
        InternalNode internal_;
        LeafNode leaf_;
        
        /**
            constructor
            the compiler complains if I don't explicitly declare this function, it's required by the constructor/destructor of PhyloNode
         */
        MyVariant():internal_(std::move(InternalNode())) {};
        
        /** constructor */
        MyVariant(LeafNode&& leaf):leaf_(std::move(leaf)) {};
        
        /** constructor */
        MyVariant(InternalNode&& internal):internal_(std::move(internal)) {};
        
        /** move constructor */
        MyVariant(MyVariant&& myvariant, const bool is_internal)
        {
            if (is_internal)
                internal_ = std::move(myvariant.internal_);
            else
                leaf_ = std::move(myvariant.leaf_);
        };
        
        /**
            destructor
            the compiler complains if I don't explicitly declare this function, it's required by the constructor/destructor of PhyloNode
         */
        ~MyVariant() {};
    };
    
    // total_lh and mid_branch_lh
    std::unique_ptr<OtherLh> other_lh_ = std::make_unique<OtherLh>(OtherLh());
    
    // NOTES: we can save 1 byte by using Bit Fields to store is_internal_ and outdated_ together
    // is our MyVariant an internal node or a leaf?
    const bool is_internal_;
    
    /// flag to avoid traversing a clade multiple times during topology optimization
    bool outdated_;
    
    // branch length
    double length_ = 0; // using float allows it to fit into the 6 bytes padding after the two bools.
    // .. using a double would make PyhloNode 8 bytes larger
    
    // NOTES: we still have 2-byte gap here

    // store our node (leaf or internal)
    MyVariant data_;
    
  public:
    /** constructor */
    PhyloNode(): is_internal_{true}
    {
        // could also be part of a separate class test, but we need to make sure that this is enforced and noone accidentally
        // changes it
        static_assert(sizeof(PhyloNode) <= 64);  // make sure it fits on a cacheline
    };
    
    /** constructor */
    PhyloNode(LeafNode&& leaf): is_internal_{false}, data_(std::move(leaf)) {};
    
    /** constructor */
    PhyloNode(InternalNode&& internal): is_internal_{true}, data_(std::move(internal)) {};
    
    /** move constructor */
    PhyloNode(PhyloNode&& node): is_internal_{node.isInternal()}, data_(std::move(node.getNode()), is_internal_),
    other_lh_(std::move(node.getOtherLh())), outdated_(node.isOutdated()), length_(node.getUpperLength()) {};
    
    /** destructor */
    ~PhyloNode()
    {
        if (is_internal_)
            data_.internal_.~InternalNode();
        else
            data_.leaf_.~LeafNode();
    }
    
    /**
        Get other_lh_
     */
    std::unique_ptr<OtherLh>& getOtherLh();
    
    /**
        Get total_lh
     */
    std::unique_ptr<SeqRegions>& getTotalLh();
    
    /**
        Set total_lh
     */
    void setTotalLh(std::unique_ptr<SeqRegions>&& total_lh);
    
    /**
        Get mid_branch_lh
     */
    std::unique_ptr<SeqRegions>& getMidBranchLh();
    
    /**
        Set mid_branch_lh
     */
    void setMidBranchLh(std::unique_ptr<SeqRegions>&& mid_branch_lh);
    
    /**
        TRUE if it's an internal node
     */
    bool isInternal() const;
    
    /**
        Get outdated_
     */
    bool isOutdated() const;
    
    /**
        Set outdated_
     */
    void setOutdated(bool new_outdated);
    
    /**
        Get length of the upper branch
     */
    RealNumType getUpperLength() const;
    
    /**
        Set length of the upper branch
     */
    void setUpperLength(const RealNumType new_length);
    
    /**
        Get length of the corresponding branch connecting this (mini-)node; it could be an upper/lower branch
     */
    RealNumType getCorrespondingLength(const MiniIndex mini_index, std::vector<PhyloNode>& nodes) const;
    
    /**
        Set length of the corresponding branch connecting this (mini-)node; it could be an upper/lower branch
     */
    void setCorrespondingLength(const MiniIndex mini_index, std::vector<PhyloNode>& nodes, const RealNumType new_length);
    
    /**
        TRUE if it's a leaf or the top of the three mini-nodes of an internal node
     */
    bool isTop(const MiniIndex mini_index) const;
    
    /**
        Update LeafNode
     */
    void setNode(LeafNode&& leaf);
    
    /**
        Update InternalNode
     */
    void setNode(InternalNode&& internal);
    
    /**
        Get MyVariant data_
     */
    MyVariant& getNode();
    
    /**
        Get partial_lh
     */
    std::unique_ptr<SeqRegions>& getPartialLh(const MiniIndex mini_index);
    
    /**
        Set partial_lh
     */
    void setPartialLh(const MiniIndex mini_index, std::unique_ptr<SeqRegions>&& partial_lh);
    
    /**
        Get the index of the neighbor node
     */
    Index getNeighborIndex(const MiniIndex mini_index) const;
    
    /**
        Set the index of the neighbor node
     */
    void setNeighborIndex(const MiniIndex mini_index, const Index neighbor_index_);
    
    /**
        Get the list of less-informative-sequences
     */
    const std::vector<uint32_t>& getLessInfoSeqs() const;
    
    /**
        Add less-informative-sequence
     */
    void addLessInfoSeqs(uint32_t seq_name_index);
    
    /**
        Get the index of the sequence name
     */
    uint32_t getSeqNameIndex() const;
    
    /**
        Set the index of the sequence name
     */
    void setSeqNameIndex(const uint32_t seq_name_index_);
    
    /**
        Update the total likelihood vector for a node.
    */
    void updateTotalLhAtNode(PhyloNode& neighbor, const Alignment& aln, const Model& model, const RealNumType threshold_prob, const bool is_root, const RealNumType blength = -1);
    
    /**
        Compute the total likelihood vector for a node.
    */
    std::unique_ptr<SeqRegions> computeTotalLhAtNode(PhyloNode& neighbor, const Alignment& aln, const Model& model, const RealNumType threshold_prob, const bool is_root, const RealNumType blength = -1);
    
    /**
        Get a vector of the indexes of neighbors
    */
    std::vector<Index> getNeighborIndexes(MiniIndex mini_index) const;
    
    /**
        Export string: name + branch length
     */
    const std::string exportString(const bool binary, const Alignment& aln) const;
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
    double length = 0;
    
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
        destructor
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
    
    
};

/** An extension of node storing more dummy data used for browsing all nodes in a stack  */
class TraversingNode {
public:
    /**
        Index to a node
     */
    Index node_index_;
    
    /**
        Count of the number of failures when traversing until the current node
     */
    short int failure_count_;
    
    /**
        Cache the likelihood difference computed at the parent node
     */
    RealNumType likelihood_diff_;
    
    /**
        Constructor
     */
    TraversingNode();
    
    /**
        Constructor
     */
    TraversingNode(const Index node_index, const short int n_failure_count, const RealNumType n_lh_diff):node_index_(node_index), failure_count_(n_failure_count), likelihood_diff_(n_lh_diff) {};
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
        TRUE to delete incoming_regions with the destructor
     */
    bool delete_regions;
    
    /**
        Constructor
     */
    UpdatingNode();
    
    /**
        Constructor
     */
    UpdatingNode(Index node_index, SeqRegions* n_incoming_regions, RealNumType n_branch_length, bool n_need_updating, RealNumType n_lh_diff, short int n_failure_count, bool n_delete_regions);
    
    /**
        destructor
     */
    ~UpdatingNode();
};

#endif

#define FOR_NEXT(node, next_node) \
for(next_node = node->next; next_node && next_node != node; next_node = next_node->next)

#define FOR_NEIGHBOR(node, neighbor_node) \
if (node->next) \
for(neighbor_node = node->next->neighbor; neighbor_node && neighbor_node != node->neighbor; neighbor_node = neighbor_node->neighbor->next->neighbor)
