#include "internal.h"
#include "leaf.h"

#pragma once

namespace cmaple {
/** An node in a phylogenetic tree, which could be either an internal or a leaf
 */
class PhyloNode {
 private:
  // store the total_lh and mid_branch_lh
  struct OtherLh {
    std::unique_ptr<SeqRegions> total_lh;
    std::unique_ptr<SeqRegions> mid_branch_lh;
  };

  /** An intermediate data structure to store either InternalNode or LeafNode */
  union MyVariant {
    InternalNode internal_;
    LeafNode leaf_;

    MyVariant() = delete;  ///< not needed (would incur overhead)

    /** constructor */
    MyVariant(LeafNode&& leaf) : leaf_(std::move(leaf)){}

    /** constructor */
    MyVariant(InternalNode&& internal) : internal_(std::move(internal)){}

    /** move-like constructor */
    MyVariant(MyVariant&& myvariant, const bool is_internal) {
      // note: none of our members are constructed, so using them (even as move
      // target) is a bad idea
      //       we need to construct them first
      if (is_internal)
        new (&internal_) InternalNode(std::move(myvariant.internal_));
      else
        new (&leaf_) LeafNode(std::move(myvariant.leaf_));
    }

    /**
     destructor
     the compiler complains if I don't explicitly declare this function, it's
     required by the constructor/destructor of PhyloNode
     */
    ~MyVariant(){}
  };

  // total_lh and mid_branch_lh
  // std::unique_ptr<OtherLh> other_lh_ = cmaple::make_unique<OtherLh>(OtherLh());
  std::unique_ptr<OtherLh> other_lh_ = std::unique_ptr<OtherLh>(new OtherLh());

  // NOTES: we can save 1 byte by using Bit Fields to store is_internal_ and
  // outdated_ together is our MyVariant an internal node or a leaf?
  const bool is_internal_;

  /// flag to avoid traversing a clade multiple times during topology
  /// optimization TRUE if the likelihoods were updated due to some SPR moves
  bool outdated_;

  // Count the number of SPR moves applied at the upper branch of this node
  uint8_t spr_count_ = 0;

  // branch length
  float length_ = 0;  // Using float allows it to fit into the 5 bytes padding
                      // after the two bools + uint8.
                      // Using a double would make PyhloNode 8 bytes larger

  // store our node (leaf or internal)
  MyVariant data_;

 public:
  /** constructor */
  PhyloNode() = delete;

  /** constructor */
  PhyloNode(LeafNode&& leaf)
      : is_internal_{false},
        outdated_{true},
        spr_count_{0},
        data_(std::move(leaf)){}

  /** constructor */
  PhyloNode(InternalNode&& internal) noexcept
      : is_internal_{true},
        outdated_{true},
        spr_count_{0},
        data_(std::move(internal)){}

  /** move constructor */
  PhyloNode(PhyloNode&& node) noexcept
      : is_internal_{node.is_internal_},
        data_(std::move(node.data_), node.is_internal_),
        other_lh_(std::move(node.other_lh_)),
        outdated_(node.outdated_),
        spr_count_(node.spr_count_),
        length_(node.length_){}

  /** destructor */
  ~PhyloNode() {
    // could also be part of a separate class test, but we need to make sure
    // that this is enforced and noone accidentally changes it
    // make sure it fits on a cacheline
    static_assert(sizeof(PhyloNode) <= 64, "PhyloNode does not fit on a cacheline");
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
   Get spr_count_
   */
  const uint8_t getSPRCount() const;

  /**
   Set spr_count_
   */
  void setSPRCount(uint8_t spr_count);

  /**
   Get length of the upper branch
   */
  cmaple::RealNumType getUpperLength() const;

  /**
   Set length of the upper branch
   */
  void setUpperLength(const cmaple::RealNumType new_length);

  /**
   Get length of the corresponding branch connecting this (mini-)node; it could
   be an upper/lower branch
   */
  cmaple::RealNumType getCorrespondingLength(
      const cmaple::MiniIndex mini_index,
      std::vector<PhyloNode>& nodes) const;

  /**
   Set length of the corresponding branch connecting this (mini-)node; it could
   be an upper/lower branch
   */
  void setCorrespondingLength(const cmaple::MiniIndex mini_index,
                              std::vector<PhyloNode>& nodes,
                              const cmaple::RealNumType new_length);

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
  std::unique_ptr<SeqRegions>& getPartialLh(const cmaple::MiniIndex mini_index);

  /**
   Set partial_lh
   */
  void setPartialLh(const cmaple::MiniIndex mini_index,
                    std::unique_ptr<SeqRegions>&& partial_lh);

  /**
   Get the index of the neighbor node
   */
  cmaple::Index getNeighborIndex(const cmaple::MiniIndex mini_index) const;

  /**
   Set the index of the neighbor node
   */
  void setNeighborIndex(const cmaple::MiniIndex mini_index,
                        const cmaple::Index neighbor_index);

  /**
   Get the list of less-informative-sequences
   */
  std::vector<cmaple::NumSeqsType>& getLessInfoSeqs();

  /**
   Add less-informative-sequence
   */
  void addLessInfoSeqs(cmaple::NumSeqsType seq_name_index);

  /**
   Get the index of the sequence name
   */
  cmaple::NumSeqsType getSeqNameIndex() const;

  /**
   Set the index of the sequence name
   */
  void setSeqNameIndex(const cmaple::NumSeqsType seq_name_index);

  /**
   Compute the total likelihood vector for a node.
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void computeTotalLhAtNode(std::unique_ptr<SeqRegions>& total_lh,
                            std::unique_ptr<SeqRegions>& this_node_mutations,
                            PhyloNode& neighbor,
                            const Alignment* aln,
                            const ModelBase* model,
                            const cmaple::RealNumType threshold_prob,
                            const bool is_root,
                            const cmaple::RealNumType blength = -1);

  /**
   Get the index of the node likelihood
   */
  const cmaple::NumSeqsType getNodelhIndex() const;

  /**
   Set the index of the node likelihood
   */
  void setNodeLhIndex(const cmaple::NumSeqsType node_lh_index);
    
    /**
     Integrate mutations (local ref) to all lh vectors
     @throw std::logic\_error if unexpected values/behaviors found during the
     operations
     */
    template <const cmaple::StateType num_states>
    auto integrateMutAllRegions(std::unique_ptr<SeqRegions>& mutations,
        const Alignment* aln, const bool inverse = false) -> void;

  /**
   Get a vector of the indexes of neighbors
   */
  // std::vector<Index> getNeighborIndexes(cmaple::MiniIndex mini_index) const;

  /**
   Export string: name + branch length
   */
  const std::string exportString(const bool binary,
                                 const std::vector<std::string>& seq_names,
                                 const bool print_internal_id,
                                 const bool show_branch_supports,
                                 const bool print_sprta_less_info_seq,
                                 const std::string& annotation_str);
};

/** An intermediate data structure to store data for calculating aLRT-SH  */
struct NodeLh {
  /**
   Set neighbor_2_lh_diff_
   */
  void setLhDiff2(const cmaple::RealNumType lh_diff);

  /**
   Get neighbor_2_lh_diff_
   */
  const cmaple::RealNumType getLhDiff2() const;

  /**
   Set neighbor_3_lh_diff_
   */
  void setLhDiff3(const cmaple::RealNumType lh_diff);

  /**
   Get neighbor_3_lh_diff_
   */
  const cmaple::RealNumType getLhDiff3() const;

  /**
   Get half of aLRT (~aLRT / 2)
   */
  const cmaple::RealNumType getHalf_aLRT() const;

  /**
   Get lh_contribution_
   */
  const cmaple::RealNumType getLhContribution() const;

  /**
   Set lh_contribution_
   */
  void setLhContribution(const cmaple::RealNumType lh_contribution);

  /**
   Get aLRT_SH
   */
  const cmaple::RealNumType get_aLRT_SH() const;

  /**
   Set aLRT_SH
   */
  void set_aLRT_SH(const cmaple::RealNumType aLRT_SH);

  /*
   Constructor
   */
  NodeLh(cmaple::RealNumType lh_contribution)
      : lh_contribution_(lh_contribution){}

 private:
  /*
   The likelihood different between NNI neighbor 2 and the ML tree (T1)
   */
  cmaple::RealNumType neighbor_2_lh_diff_;

  /*
   The likelihood different between NNI neighbor 3 and the ML tree (T1)
   */
  cmaple::RealNumType neighbor_3_lh_diff_;

  /*
   The aLRT_SH value of the upper branch
   */
  cmaple::RealNumType aLRT_SH_;

  /*
   The likelihood contribution of this node to the total likelihood
   */
  cmaple::RealNumType lh_contribution_;
};

template <const StateType num_states>
void cmaple::PhyloNode::computeTotalLhAtNode(
    std::unique_ptr<SeqRegions>& total_lh,
    std::unique_ptr<SeqRegions>& this_node_mutations,
    PhyloNode& neighbor,
    const Alignment* aln,
    const ModelBase* model,
    const RealNumType threshold_prob,
    const bool is_root,
    const RealNumType blength) {
  assert(aln);
  assert(model);
    
  // if node is root
  if (is_root) {
    getPartialLh(TOP)->computeTotalLhAtRoot<num_states>(total_lh, model,
                                                        blength);
    // if node is normal nodes
  } else {
      std::unique_ptr<SeqRegions>& ori_upper_lr_regions =
        neighbor.getPartialLh(getNeighborIndex(TOP).getMiniIndex());
      // 1. create a new upper_lr_regions that integrate the mutations, if any
      std::unique_ptr<SeqRegions> mut_integrated_upper_lr_regions =
          (this_node_mutations && this_node_mutations->size())
          ? ori_upper_lr_regions
            ->integrateMutations<num_states>(this_node_mutations, aln)
          : nullptr;
      // 2. create the pointer that points to the appropriate upper_lr_regions
      const std::unique_ptr<SeqRegions>* upper_lr_regions_ptr =
          (this_node_mutations && this_node_mutations->size())
          ? &mut_integrated_upper_lr_regions
          : &ori_upper_lr_regions;
      // 3. create a reference from that pointer
      auto& upper_lr_regions = *upper_lr_regions_ptr;
      
    std::unique_ptr<SeqRegions>& lower_regions = getPartialLh(TOP);
    upper_lr_regions->mergeUpperLower<num_states>(total_lh, getUpperLength(),
                                      *lower_regions, blength, aln, model,
                                      threshold_prob);
  }
}

template <const StateType num_states>
auto cmaple::PhyloNode::integrateMutAllRegions(
    std::unique_ptr<SeqRegions>& mutations,
    const Alignment* aln, const bool inverse) -> void
{
    // lower regions
    setPartialLh(TOP, getPartialLh(TOP)->integrateMutations<num_states>(mutations, aln, inverse));
    
    // upper left/right regions
    if (isInternal())
    {
        setPartialLh(LEFT, getPartialLh(LEFT)->integrateMutations<num_states>(mutations, aln, inverse));
        setPartialLh(RIGHT, getPartialLh(RIGHT)->integrateMutations<num_states>(mutations, aln, inverse));
    }
    
    // mid-branch regions
    if (getMidBranchLh())
    {
        if (getUpperLength() > 0)
        {
            setMidBranchLh(getMidBranchLh()->integrateMutations<num_states>(mutations, aln, inverse));
        }
        else
        {
            setMidBranchLh(nullptr);
        }
    }
}

}  // namespace cmaple

// just for testing
/** operator<< */
std::ostream& operator<<(std::ostream& os, const cmaple::Index& index);
