#include "phylonode.h"
using namespace std;

std::ostream& operator<<(std::ostream& os, const Index& index)
{
    os << index.getVectorIndex() << " " << index.getMiniIndex();
    return os;
}

std::unique_ptr<PhyloNode::OtherLh>& PhyloNode::getOtherLh()
{
    return other_lh_;
}

std::unique_ptr<SeqRegions>& PhyloNode::getTotalLh()
{
    return other_lh_->total_lh;
}

void PhyloNode::setTotalLh(std::unique_ptr<SeqRegions>&& total_lh)
{
    other_lh_->total_lh = std::move(total_lh);
}

std::unique_ptr<SeqRegions>& PhyloNode::getMidBranchLh()
{
    return other_lh_->mid_branch_lh;
}

void PhyloNode::setMidBranchLh(std::unique_ptr<SeqRegions>&& mid_branch_lh)
{
    other_lh_->mid_branch_lh = std::move(mid_branch_lh);
}

bool PhyloNode::isInternal() const
{
    return is_internal_;
}

bool PhyloNode::isOutdated() const
{
    return outdated_;
}

void PhyloNode::setOutdated(bool new_outdated)
{
    outdated_ = new_outdated;
}

const bool PhyloNode::isSPRApplied() const
{
    return spr_applied_;
}

void PhyloNode::setSPRApplied(bool spr_applied)
{
    spr_applied_ = spr_applied;
}

RealNumType PhyloNode::getUpperLength() const
{
    return length_;
}

void PhyloNode::setUpperLength(const double new_length)
{
    length_ = new_length;
}

RealNumType PhyloNode::getCorrespondingLength(const MiniIndex mini_index, std::vector<PhyloNode>& nodes) const
{
    // if it's a left/right mini node -> return the length of the lower left/right branch
    if (is_internal_ && mini_index != TOP)
    {
        Index neighbor_index = getNeighborIndex(mini_index);
        return nodes[neighbor_index.getVectorIndex()].getUpperLength();
    }
    
    // by default, return the length of the upper branch
    return length_;
}

void PhyloNode::setCorrespondingLength(const MiniIndex mini_index, std::vector<PhyloNode>& nodes, const RealNumType new_length)
{
    // if it's a left/right mini node -> update the length of the lower left/right branch
    if (is_internal_ && mini_index != TOP)
    {
        Index neighbor_index = getNeighborIndex(mini_index);
        return nodes[neighbor_index.getVectorIndex()].setUpperLength(new_length);
    }
    // by default, return the length of the upper branch
    else
        length_ = new_length;
}

void PhyloNode::setNode(LeafNode&& leaf)
{
    ASSERT(!is_internal_ && "Wrong data type! data should be a leaf node");
    
    // update leaf
    data_.leaf_ = std::move(leaf);
}

void PhyloNode::setNode(InternalNode&& internal)
{
    ASSERT(is_internal_ && "Wrong data type! data should be an internal node");
    
    // update internal
    data_.internal_ = std::move(internal);
}

PhyloNode::MyVariant& PhyloNode::getNode()
{
    return data_;
}

std::unique_ptr<SeqRegions>& PhyloNode::getPartialLh(const MiniIndex mini_index)
{
    // if it's an internal node
    if (is_internal_)
    {
        // return the corresponding partial_lh based on the mini-index
        return data_.internal_.partial_lh3_[mini_index] ;
    }
    
    // if it's a leaf
    return data_.leaf_.partial_lh_;
}

void PhyloNode::setPartialLh(const MiniIndex mini_index, std::unique_ptr<SeqRegions>&& partial_lh)
{
    // if it's an internal node
    if (is_internal_)
    {
        // update the corresponding partial_lh based on the mini-index
        data_.internal_.partial_lh3_[mini_index] = std::move(partial_lh);
    }
    // if it's a leaf
    else
    {
        // update the partial_lh
        data_.leaf_.partial_lh_ = std::move(partial_lh);
    }
}

Index PhyloNode::getNeighborIndex(const MiniIndex mini_index) const
{
    // if it's an internal node
    if (is_internal_)
    {
        // return the corresponding neighbor_index based on the mini-index
        return data_.internal_.neighbor_index3_[mini_index] ;
    }
    
    // if it's a leaf
    return data_.leaf_.neighbor_index_;
}

void PhyloNode::setNeighborIndex(const MiniIndex mini_index, const Index neighbor_index_)
{
    // if it's an internal node
    if (is_internal_)
    {
        // update the corresponding neighbor_index based on the mini-index
        data_.internal_.neighbor_index3_[mini_index]  = neighbor_index_;
    }
    // if it's a leaf
    else
    {
        data_.leaf_.neighbor_index_ = neighbor_index_;
    }
}

std::vector<NumSeqsType>& PhyloNode::getLessInfoSeqs()
{
    ASSERT(!is_internal_);
    return data_.leaf_.less_info_seqs_;
}

void PhyloNode::addLessInfoSeqs(NumSeqsType seq_name_index)
{
    ASSERT(!is_internal_);
    data_.leaf_.less_info_seqs_.push_back(seq_name_index);
}

NumSeqsType PhyloNode::getSeqNameIndex() const
{
    ASSERT(!is_internal_);
    return data_.leaf_.seq_name_index_;
}

void PhyloNode::setSeqNameIndex(NumSeqsType seq_name_index_)
{
    ASSERT(!is_internal_);
    data_.leaf_.seq_name_index_ = seq_name_index_;
}

/*std::vector<Index> PhyloNode::getNeighborIndexes(MiniIndex mini_index) const
{
    std::vector<Index> neighbor_indexes;
    neighbor_indexes.reserve(2);
    
    // if it's an internal -> return the indexes of its 2 neighbors
    if (isInternal())
    {
        // I don't use a loop here because I want to keep the order of neighbors same as those in the previous implementation of pointers
        switch (mini_index)
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

void PhyloNode::computeTotalLhAtNode(std::unique_ptr<SeqRegions>& total_lh, PhyloNode& neighbor, const Alignment& aln, const std::unique_ptr<Model>& model, const RealNumType threshold_prob, const bool is_root, const RealNumType blength)
{
    // if node is root
    if (is_root)
        getPartialLh(TOP)->computeTotalLhAtRoot(total_lh, aln.num_states, model, blength);
    // if not is normal nodes
    else
    {
        std::unique_ptr<SeqRegions>& lower_regions = getPartialLh(TOP);
        neighbor.getPartialLh(getNeighborIndex(TOP).getMiniIndex())->mergeUpperLower(total_lh, getUpperLength(), *lower_regions, blength, aln, model, threshold_prob);
    }
}

const std::string PhyloNode::exportString(const bool binary, const Alignment& aln, const bool show_branch_supports)
{
    if (!isInternal())
    {
        string length_str = getUpperLength() < 0 ? "0" : convertDoubleToString(getUpperLength(), 20);
        // without minor sequences -> simply return node's name and its branch length
        std::vector<NumSeqsType>& less_info_seqs = getLessInfoSeqs();
        const PositionType num_less_info_seqs = less_info_seqs.size();
        if (num_less_info_seqs == 0)
            return aln.data[getSeqNameIndex()].seq_name + ":" + length_str;
        // with minor sequences -> return minor sequences' names with zero branch lengths
        else
        {
            string branch_support = show_branch_supports ? "0" : "";
            
            string output = "(" + aln.data[getSeqNameIndex()].seq_name + ":0";
            // export less informative sequences in binary tree format
            if (binary)
            {
                string closing_brackets = "";
                for (PositionType i = 0; i < num_less_info_seqs - 1; i++)
                {
                    output += ",(" + aln.data[less_info_seqs[i]].seq_name + ":0";
                    closing_brackets += ")" + branch_support + ":0";
                }
                output += "," + aln.data[less_info_seqs[num_less_info_seqs - 1]].seq_name + ":0" + closing_brackets;
            }
            // export less informative sequences in mutifurcating tree format
            else
            {
                for (auto minor_seq_name_index : less_info_seqs)
                    output += "," + aln.data[minor_seq_name_index].seq_name + ":0";
            }
            output += ")" + branch_support + ":" + length_str;
            return output;
        }
    }
    
    return "";
}

const NumSeqsType PhyloNode::getNodelhIndex() const
{
    ASSERT(isInternal());
    return data_.internal_.node_lh_index_;
}

void PhyloNode::setNodeLhIndex(const NumSeqsType node_lh_index)
{
    ASSERT(isInternal());
    data_.internal_.node_lh_index_ = node_lh_index;
}

void NodeLh::setLhDiff2(const RealNumType lh_diff)
{
    neighbor_2_lh_diff_ = lh_diff;
}

const RealNumType NodeLh::getLhDiff2() const
{
    return neighbor_2_lh_diff_;
}

void NodeLh::setLhDiff3(const RealNumType lh_diff)
{
    neighbor_3_lh_diff_ = lh_diff;
}

const RealNumType NodeLh::getLhDiff3() const
{
    return neighbor_3_lh_diff_;
}

const RealNumType NodeLh::getHalf_aLRT() const
{
    /*  aLRT = 2(LT1 - max(LT2x,LT2y)); where LT2x, LT2y are the log lh of the two NNI neighbors of T1
     <=> aLRT = 2(LT1 - max(diff_x + LT1, diff_y + LT1)); where diff_x = LT2x - LT1; and similarly to diff_y
     <=> aLRT = 2(LT1 - max(diff_x, diff_y) - LT1)
     <=> aLRT = 2(-max(diff_x, diff_y)) = 2 * (-max_nni_neighbor_lh_diff)
     <=> (half aLRT) = aLRT / 2 = -max_nni_neighbor_lh_diff
     */
    return neighbor_2_lh_diff_ > neighbor_3_lh_diff_ ? -neighbor_2_lh_diff_ : -neighbor_3_lh_diff_;
}

const RealNumType NodeLh::getLhContribution() const
{
    return lh_contribution_;
}

void NodeLh::setLhContribution(const RealNumType lh_contribution)
{
    lh_contribution_ = lh_contribution;
}

const RealNumType NodeLh::get_aLRT_SH() const
{
    return aLRT_SH_;
}

void NodeLh::set_aLRT_SH(const RealNumType aLRT_SH)
{
    aLRT_SH_ = aLRT_SH;
}
