#include "phylonode.h"
using namespace std;
using namespace cmaple;

std::ostream& operator<<(std::ostream& os, const Index& index) {
  os << index.getVectorIndex() << " " << index.getMiniIndex();
  return os;
}

std::unique_ptr<cmaple::PhyloNode::OtherLh>& cmaple::PhyloNode::getOtherLh() {
  return other_lh_;
}

std::unique_ptr<SeqRegions>& cmaple::PhyloNode::getTotalLh() {
  return other_lh_->total_lh;
}

void cmaple::PhyloNode::setTotalLh(std::unique_ptr<SeqRegions>&& total_lh) {
  other_lh_->total_lh = std::move(total_lh);
}

std::unique_ptr<SeqRegions>& cmaple::PhyloNode::getMidBranchLh() {
  return other_lh_->mid_branch_lh;
}

void cmaple::PhyloNode::setMidBranchLh(
    std::unique_ptr<SeqRegions>&& mid_branch_lh) {
  other_lh_->mid_branch_lh = std::move(mid_branch_lh);
}

auto cmaple::PhyloNode::isInternal() const -> bool {
  return is_internal_;
}

auto cmaple::PhyloNode::isOutdated() const -> bool {
  return outdated_;
}

void cmaple::PhyloNode::setOutdated(bool new_outdated) {
  outdated_ = new_outdated;
}

auto cmaple::PhyloNode::getSPRCount() const -> const uint8_t {
  return spr_count_;
}

void cmaple::PhyloNode::setSPRCount(uint8_t spr_count) {
  spr_count_ = spr_count;
}

auto cmaple::PhyloNode::getUpperLength() const -> RealNumType {
  return static_cast<RealNumType>(length_);
}

void cmaple::PhyloNode::setUpperLength(const double new_length) {
  length_ = new_length;
}

auto cmaple::PhyloNode::getCorrespondingLength(
    const MiniIndex mini_index,
    std::vector<PhyloNode>& nodes) const -> RealNumType {
  // if it's a left/right mini node -> return the length of the lower left/right
  // branch
  if (is_internal_ && mini_index != TOP) {
    Index neighbor_index = getNeighborIndex(mini_index);
    return nodes[neighbor_index.getVectorIndex()].getUpperLength();
  }

  // by default, return the length of the upper branch
  return static_cast<RealNumType>(length_);
}

void cmaple::PhyloNode::setCorrespondingLength(const MiniIndex mini_index,
                                               std::vector<PhyloNode>& nodes,
                                               const RealNumType new_length) {
  // if it's a left/right mini node -> update the length of the lower left/right
  // branch
  if (is_internal_ && mini_index != TOP) {
    Index neighbor_index = getNeighborIndex(mini_index);
    return nodes[neighbor_index.getVectorIndex()].setUpperLength(new_length);
  }
  // by default, return the length of the upper branch
  else {
    length_ = new_length;
  }
}

void cmaple::PhyloNode::setNode(LeafNode&& leaf) {
  assert(!is_internal_ && "Wrong data type! data should be a leaf node");

  // update leaf
  data_.leaf_ = std::move(leaf);
}

void cmaple::PhyloNode::setNode(InternalNode&& internal) {
  assert(is_internal_ && "Wrong data type! data should be an internal node");

  // update internal
  data_.internal_ = std::move(internal);
}

auto cmaple::PhyloNode::getNode() -> cmaple::PhyloNode::MyVariant& {
  return data_;
}

std::unique_ptr<SeqRegions>& cmaple::PhyloNode::getPartialLh(
    const MiniIndex mini_index) {
  // if it's an internal node
  if (is_internal_) {
    // return the corresponding partial_lh based on the mini-index
    return data_.internal_.partial_lh3_[mini_index];
  }

  // if it's a leaf
  return data_.leaf_.partial_lh_;
}

void cmaple::PhyloNode::setPartialLh(const MiniIndex mini_index,
                                     std::unique_ptr<SeqRegions>&& partial_lh) {
  // if it's an internal node
  if (is_internal_) {
    // update the corresponding partial_lh based on the mini-index
    data_.internal_.partial_lh3_[mini_index] = std::move(partial_lh);
  }
  // if it's a leaf
  else {
    // update the partial_lh
    data_.leaf_.partial_lh_ = std::move(partial_lh);
  }
}

auto cmaple::PhyloNode::getNeighborIndex(const MiniIndex mini_index) const
    -> Index {
  // if it's an internal node
  if (is_internal_) {
    // return the corresponding neighbor_index based on the mini-index
    return data_.internal_.neighbor_index3_[mini_index];
  }

  // if it's a leaf
  return data_.leaf_.neighbor_index_;
}

void cmaple::PhyloNode::setNeighborIndex(const MiniIndex mini_index,
                                         const Index neighbor_index_) {
  // if it's an internal node
  if (is_internal_) {
    // update the corresponding neighbor_index based on the mini-index
    data_.internal_.neighbor_index3_[mini_index] = neighbor_index_;
  }
  // if it's a leaf
  else {
    data_.leaf_.neighbor_index_ = neighbor_index_;
  }
}

std::vector<NumSeqsType>& cmaple::PhyloNode::getLessInfoSeqs() {
  assert(!is_internal_);
  return data_.leaf_.less_info_seqs_;
}

void cmaple::PhyloNode::addLessInfoSeqs(NumSeqsType seq_name_index) {
  assert(!is_internal_);
  data_.leaf_.less_info_seqs_.push_back(seq_name_index);
}

auto cmaple::PhyloNode::getSeqNameIndex() const -> NumSeqsType {
  assert(!is_internal_);
  return data_.leaf_.seq_name_index_;
}

void cmaple::PhyloNode::setSeqNameIndex(NumSeqsType seq_name_index_) {
  assert(!is_internal_);
  data_.leaf_.seq_name_index_ = seq_name_index_;
}

/*std::vector<Index> cmaple::PhyloNode::getNeighborIndexes(MiniIndex mini_index)
const
{
    std::vector<Index> neighbor_indexes;
    neighbor_indexes.reserve(2);

    // if it's an internal -> return the indexes of its 2 neighbors
    if (isInternal())
    {
        // I don't use a loop here because I want to keep the order of neighbors
same as those in the previous implementation of pointers switch (mini_index)
        {
            case TOP:
            {
                neighbor_indexes.push_back(getNeighborIndex(RIGHT));
                neighbor_indexes.push_back(getNeighborIndex(LEFT));
                break;
            }
            case LEFT:
            {
                neighbor_indexes.push_back(getNeighborIndex(TOP));
                neighbor_indexes.push_back(getNeighborIndex(RIGHT));
                break;
            }
            case RIGHT:
            {
                neighbor_indexes.push_back(getNeighborIndex(LEFT));
                neighbor_indexes.push_back(getNeighborIndex(TOP));
                break;
            }
            default:
                break;
        }
    }
    // if it's a leaf -> return an empty vector

    return neighbor_indexes;
}*/

const std::string cmaple::PhyloNode::exportString(
    const bool binary,
    const std::vector<std::string>& seq_names,
    const bool show_branch_supports) {
  if (!isInternal()) {
    string length_str = getUpperLength() <= 0
                            ? "0"
                            : convertDoubleToString(getUpperLength(), 12);
    // without minor sequences -> simply return node's name and its branch
    // length
    std::vector<NumSeqsType>& less_info_seqs = getLessInfoSeqs();
    const std::vector<NumSeqsType>::size_type num_less_info_seqs = less_info_seqs.size();
    if (num_less_info_seqs == 0) {
      return seq_names[getSeqNameIndex()] + ":" + length_str;
      // with minor sequences -> return minor sequences' names with zero
      // branch lengths
    } else {
      string branch_support = show_branch_supports ? "0" : "";

      string output = "";
      // export less informative sequences in binary tree format
      if (binary) {
        output.resize(num_less_info_seqs, '(');
        output += seq_names[getSeqNameIndex()];
        
        // add the first less-info-seq
        output += ":0," + seq_names[less_info_seqs[0]] + ":0)";
        
        // add the remaining less-info-seqs
        for (std::vector<NumSeqsType>::size_type  i = 1; i < less_info_seqs.size(); ++i) {
          output += branch_support + ":0," + seq_names[less_info_seqs[i]] + ":0)";
        }
        
        output += branch_support + ":" + length_str;
      }
      // export less informative sequences in mutifurcating tree format
      else {
        output = "(" + seq_names[getSeqNameIndex()] + ":0";
        
        for (auto minor_seq_name_index : less_info_seqs) {
          output += "," + seq_names[minor_seq_name_index] + ":0";
        }
        
        output += ")" + branch_support + ":" + length_str;
      }
      return output;
    }
  }

  return "";
}

auto cmaple::PhyloNode::getNodelhIndex() const -> const NumSeqsType {
  assert(isInternal());
  return data_.internal_.node_lh_index_;
}

void cmaple::PhyloNode::setNodeLhIndex(const NumSeqsType node_lh_index) {
  assert(isInternal());
  data_.internal_.node_lh_index_ = node_lh_index;
}

void cmaple::NodeLh::setLhDiff2(const RealNumType lh_diff) {
  neighbor_2_lh_diff_ = lh_diff;
}

auto cmaple::NodeLh::getLhDiff2() const -> const RealNumType {
  return neighbor_2_lh_diff_;
}

void cmaple::NodeLh::setLhDiff3(const RealNumType lh_diff) {
  neighbor_3_lh_diff_ = lh_diff;
}

auto cmaple::NodeLh::getLhDiff3() const -> const RealNumType {
  return neighbor_3_lh_diff_;
}

auto cmaple::NodeLh::getHalf_aLRT() const -> const RealNumType {
  /*  aLRT = 2(LT1 - max(LT2x,LT2y)); where LT2x, LT2y are the log lh of the two
   NNI neighbors of T1
   <=> aLRT = 2(LT1 - max(diff_x + LT1, diff_y + LT1)); where diff_x = LT2x -
   LT1; and similarly to diff_y
   <=> aLRT = 2(LT1 - max(diff_x, diff_y) - LT1)
   <=> aLRT = 2(-max(diff_x, diff_y)) = 2 * (-max_nni_neighbor_lh_diff)
   <=> (half aLRT) = aLRT / 2 = -max_nni_neighbor_lh_diff
   */
  return neighbor_2_lh_diff_ > neighbor_3_lh_diff_ ? -neighbor_2_lh_diff_
                                                   : -neighbor_3_lh_diff_;
}

auto cmaple::NodeLh::getLhContribution() const -> const RealNumType {
  return lh_contribution_;
}

void cmaple::NodeLh::setLhContribution(const RealNumType lh_contribution) {
  lh_contribution_ = lh_contribution;
}

auto cmaple::NodeLh::get_aLRT_SH() const -> const RealNumType {
  return aLRT_SH_;
}

void cmaple::NodeLh::set_aLRT_SH(const RealNumType aLRT_SH) {
  aLRT_SH_ = aLRT_SH;
}
