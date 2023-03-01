#include "internal.h"
#include "leaf.h"

#pragma once

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
        MyVariant():internal_(InternalNode()) {};
        
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
    /// TRUE if the likelihoods were updated due to some SPR moves
    bool outdated_;
    
    // TRUE if we applied SPR move on the upper branch of this node
    bool spr_applied_;
    
    // branch length
    float length_ = 0; // using float allows it to fit into the 6 bytes padding after the two bools.
    // .. using a double would make PyhloNode 8 bytes larger
    
    // NOTES: we still have 2-byte gap here

    // store our node (leaf or internal)
    MyVariant data_;
    
  public:
    /** constructor */
    PhyloNode(): is_internal_{true} {};
    
    /** constructor */
    PhyloNode(LeafNode&& leaf): is_internal_{false}, data_(std::move(leaf)) {};
    
    /** constructor */
    // unused
    // PhyloNode(InternalNode&& internal): is_internal_{true}, data_(std::move(internal)) {};
    
    /** move constructor */
    PhyloNode(PhyloNode&& node): is_internal_{node.isInternal()}, data_(std::move(node.getNode()), is_internal_),
    other_lh_(std::move(node.getOtherLh())), outdated_(node.isOutdated()), spr_applied_(node.isSPRApplied()), length_(node.getUpperLength()) {};
    
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
        Get spr_applied_
     */
    const bool isSPRApplied() const;
    
    /**
        Set spr_applied_
     */
    void setSPRApplied(bool spr_applied);
    
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
    void setNeighborIndex(const MiniIndex mini_index, const Index neighbor_index);
    
    /**
        Get the list of less-informative-sequences
     */
    const std::vector<NumSeqsType>& getLessInfoSeqs() const;
    
    /**
        Add less-informative-sequence
     */
    void addLessInfoSeqs(NumSeqsType seq_name_index);
    
    /**
        Get the index of the sequence name
     */
    NumSeqsType getSeqNameIndex() const;
    
    /**
        Set the index of the sequence name
     */
    void setSeqNameIndex(const NumSeqsType seq_name_index);
    
    /**
        Compute the total likelihood vector for a node.
    */
    void computeTotalLhAtNode(std::unique_ptr<SeqRegions>& total_lh, PhyloNode& neighbor, const Alignment& aln, const Model& model, const RealNumType threshold_prob, const bool is_root, const RealNumType blength = -1);
    
    /**
        Get the index of the node likelihood
     */
    const NumSeqsType getNodelhIndex() const;
    
    /**
        Set the index of the node likelihood
     */
    void setNodeLhIndex(const NumSeqsType node_lh_index);
    
    /**
        Get a vector of the indexes of neighbors
    */
    // std::vector<Index> getNeighborIndexes(MiniIndex mini_index) const;
    
    /**
        Export string: name + branch length
     */
    const std::string exportString(const bool binary, const Alignment& aln) const;
};

/** An intermediate data structure to store data for calculating aLRT-SH  */
struct NodeLh
{
    /**
        Get aLRT_
     */
    const RealNumType getaLRT() const;
    
    /**
        Set aLRT_
     */
    void setaLRT(const RealNumType aLRT);
    
    /**
        Get lh_contribution_
     */
    const RealNumType getLhContribution() const;
    
    /**
        Set lh_contribution_
     */
    void setLhContribution(const RealNumType lh_contribution);
    
    /**
        Get site_lh_contributions_
     */
    std::vector<RealNumType>& getSiteLhContributions();
    
    /*
        Constructor
     */
    NodeLh(RealNumType lh_contribution):lh_contribution_(lh_contribution) {};
    
private:
    /*
        The aLRT value of the upper branch of this node
     */
    RealNumType aLRT_;
    
    /*
        The likelihood contribution of this node to the total likelihood
     */
    RealNumType lh_contribution_;
    
    /*
        The site-likelihood contribution of this node to the total likelihood
     */
    std::vector<RealNumType> site_lh_contributions_;
};

// just for testing
/** operator<< */
std::ostream& operator<<(std::ostream& os, const Index& index);

