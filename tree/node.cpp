#include "node.h"
using namespace std;

// ########### BEGIN OF NEW DATA STRUCTURES FOR PHYLOGENETIC NODES ###########
std::ostream& operator<<(std::ostream& os, const Index& index)
{
    os << index.getVectorIndex() << " " << index.getMiniIndex();
    return os;
}

SeqRegions& PhyloNode::getTotalLh()
{
    return other_lh_->getTotalLh();
}

void PhyloNode::setTotalLh(SeqRegions&& total_lh)
{
    other_lh_->setTotalLh(std::move(total_lh));
}

SeqRegions& PhyloNode::getMidBranchLh()
{
    return other_lh_->getMidBranchLh();
}

void PhyloNode::setMidBranchLh(SeqRegions&& mid_branch_lh)
{
    other_lh_->setMidBranchLh(std::move(mid_branch_lh));
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

float PhyloNode::getLength() const
{
    return length_;
}

void PhyloNode::setLength(float new_length)
{
    length_ = new_length;
}

bool PhyloNode::isTop(const MiniIndex mini_index) const
{
    return (!is_internal_
            || mini_index == TOP);
}

void PhyloNode::setNode(LeafNode&& leaf)
{
    ASSERT(!is_internal_ && "Wrong data type! data should be a leaf node");
    
    // update leaf
    data_.leaf = std::move(leaf);
}

void PhyloNode::setNode(InternalNode&& internal)
{
    ASSERT(is_internal_ && "Wrong data type! data should be an internal node");
    
    // update internal
    data_.internal = std::move(internal);
}

std::unique_ptr<SeqRegions>& PhyloNode::getPartialLh(const MiniIndex mini_index)
{
    // if it's an internal node
    if (is_internal_)
    {
        // return the corresponding partial_lh based on the mini-index
        return data_.internal.partial_lh3_[mini_index] ;
    }
    
    // if it's a leaf
    return data_.leaf.partial_lh_;
}

void PhyloNode::setPartialLh(const MiniIndex mini_index, std::unique_ptr<SeqRegions>& partial_lh)
{
    // if it's an internal node
    if (is_internal_)
    {
        // update the corresponding partial_lh based on the mini-index
        data_.internal.partial_lh3_[mini_index] = std::move(partial_lh);
    }
    // if it's a leaf
    else
    {
        // update the partial_lh
        data_.leaf.partial_lh_ = std::move(partial_lh);
    }
}

Index PhyloNode::getNeighborIndex(const MiniIndex mini_index) const
{
    // if it's an internal node
    if (is_internal_)
    {
        // return the corresponding neighbor_index based on the mini-index
        return data_.internal.neighbor_index3_[mini_index] ;
    }
    
    // if it's a leaf
    return data_.leaf.neighbor_index_;
}

void PhyloNode::setNeighborIndex(const MiniIndex mini_index, const Index neighbor_index_)
{
    // if it's an internal node
    if (is_internal_)
    {
        // update the corresponding neighbor_index based on the mini-index
        data_.internal.neighbor_index3_[mini_index]  = neighbor_index_;
    }
    // if it's a leaf
    else
    {
        data_.leaf.neighbor_index_ = neighbor_index_;
    }
}

const std::vector<std::string>& PhyloNode::getLessInfoSeqs() const
{
    ASSERT(!is_internal_);
    return data_.leaf.less_info_seqs_;
}

void PhyloNode::addLessInfoSeqs(std::string&& seq_name)
{
    ASSERT(!is_internal_);
    data_.leaf.less_info_seqs_.push_back(std::move(seq_name));
}

uint32_t PhyloNode::getSeqNameIndex() const
{
    ASSERT(!is_internal_);
    return data_.leaf.seq_name_index_;
}

void PhyloNode::setSeqNameIndex(uint32_t seq_name_index_)
{
    ASSERT(!is_internal_);
    data_.leaf.seq_name_index_ = seq_name_index_;
}

// ########### END OF NEW DATA STRUCTURES FOR PHYLOGENETIC NODES ###########

Node::Node(bool is_top_node)
{
    seq_name = "";
    length = 0;
    is_top = is_top_node;
    next = NULL;
    neighbor = NULL;
    partial_lh = NULL;
    mid_branch_lh = NULL;
    total_lh = NULL;
    outdated = false;
}

Node::Node(int n_id, string n_seq_name)
{
    seq_name = std::move(n_seq_name);
    length = 0;
    is_top = false;
    next = NULL;
    neighbor = NULL;
    partial_lh = NULL;
    mid_branch_lh = NULL;
    total_lh = NULL;
    outdated = false;
}

Node::Node(string n_seq_name)
{
    seq_name = std::move(n_seq_name);
    length = 0;
    is_top = true;
    next = NULL;
    neighbor = NULL;
    partial_lh = NULL;
    mid_branch_lh = NULL;
    total_lh = NULL;
    outdated = false;
}

Node::~Node()
{
    if (partial_lh) delete partial_lh;
    if (total_lh) delete total_lh;
    if (mid_branch_lh) delete mid_branch_lh;
}

bool Node::isLeave()
{
    return next == NULL;
}

Node* Node::getTopNode()
{
    if (this->is_top)
        return this;
    
    Node* next_node;
    Node* node = this;
    
    FOR_NEXT(node, next_node)
    {
        if (next_node->is_top)
            return next_node;
    }
    
    return NULL;
}

Node* Node::getOtherNextNode()
{
    Node* next_node;
    Node* node = this;
    
    FOR_NEXT(node, next_node)
    {
        if (!next_node->is_top)
            return next_node;
    }
    
    return NULL;
}

string Node::exportString(bool binary)
{
    if (isLeave())
    {
        string length_str = length < 0 ? "0" : convertDoubleToString(length);
        // without minor sequences -> simply return node's name and its branch length
        PositionType num_less_info_seqs = less_info_seqs.size();
        if (num_less_info_seqs == 0)
            return seq_name + ":" + length_str;
        // with minor sequences -> return minor sequences' names with zero branch lengths
        else
        {
            string output = "(" + seq_name + ":0";
            // export less informative sequences in binary tree format
            if (binary)
            {
                string closing_brackets = "";
                for (PositionType i = 0; i < num_less_info_seqs - 1; i++)
                {
                    output += ",(" + less_info_seqs[i] + ":0";
                    closing_brackets += "):0";
                }
                output += "," + less_info_seqs[num_less_info_seqs - 1] + ":0" + closing_brackets;
            }
            // export less informative sequences in mutifurcating tree format
            else
            {
                for (string minor_seq_name : less_info_seqs)
                    output += "," + std::move(minor_seq_name) + ":0";
            }
            output += "):" + length_str;
            return output;
        }
    }
    
    return "";
}

SeqRegions* Node::getPartialLhAtNode(const Alignment& aln, const Model& model, RealNumType threshold_prob)
{
    // if partial_lh has not yet computed (~NULL) -> compute it from next nodes
    if (!partial_lh)
    {        
        // the phylonode is an internal node
        if (next)
        {
            // if node is a top node -> partial_lh is the lower lh regions
            if (is_top)
            {
                // extract the two lower vectors of regions
                Node* next_node_1 = next;
                SeqRegions* regions1 = next_node_1->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
                Node* next_node_2 = next_node_1->next;
                SeqRegions* regions2 = next_node_2->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
                
                // compute partial_lh
                regions1->mergeTwoLowers(partial_lh, next_node_1->length, *regions2, next_node_2->length, aln, model, threshold_prob);
            }
            // otherwise -> partial_lh is the upper left/right regions
            else
            {
                // extract the upper and the lower vectors of regions
                SeqRegions* upper_regions, *lower_regions;
                RealNumType upper_blength, lower_blength;
                
                Node* next_node_1 = next;
                Node* next_node_2 = next_node_1->next;
                
                if (next_node_1->is_top)
                {
                    upper_regions = next_node_1->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
                    upper_blength = next_node_1->length;
                    lower_regions = next_node_2->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
                    lower_blength = next_node_2->length;
                }
                else
                {
                    upper_regions = next_node_2->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
                    upper_blength = next_node_2->length;
                    lower_regions = next_node_1->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
                    lower_blength = next_node_1->length;
                }
                
                // compute partial_lh
                upper_regions->mergeUpperLower(partial_lh, upper_blength, *lower_regions, lower_blength, aln, model, threshold_prob);
            }
        }
        // the phylonode is a tip, partial_lh must be already computed
        else
            outError("Something went wrong! Lower likelihood regions has not been computed at tip!");
    }
    
    // return
    return partial_lh;
}

SeqRegions* Node::computeTotalLhAtNode(const Alignment& aln, const Model& model, RealNumType threshold_prob, bool is_root, bool update, RealNumType blength)
{
    SeqRegions* new_regions = NULL;
    
    // if node is root
    if (is_root)
        new_regions = getPartialLhAtNode(aln, model, threshold_prob)->computeTotalLhAtRoot(aln.num_states, model, blength);
    // if not is normal nodes
    else
    {
        SeqRegions* lower_regions = getPartialLhAtNode(aln, model, threshold_prob);
        neighbor->getPartialLhAtNode(aln, model, threshold_prob)->mergeUpperLower(new_regions, length, *lower_regions, blength, aln, model, threshold_prob);
    }
    
    // update if necessary
    if (update)
    {
        if (total_lh) delete total_lh;
        total_lh = new_regions;
    }
    
    return new_regions;
}

TraversingNode::TraversingNode()
{
    node = NULL;
    failure_count = 0;
    likelihood_diff = 0;
}

TraversingNode::TraversingNode(Node* n_node, short int n_failure_count, RealNumType n_lh_diff)
{
    node = n_node;
    failure_count = n_failure_count;
    likelihood_diff = n_lh_diff;
}

TraversingNode::~TraversingNode()
{
    // do nothing
}

UpdatingNode::UpdatingNode():TraversingNode()
{
    incoming_regions = NULL;
    branch_length = 0;
    need_updating = false;
    delete_regions = false;
}

UpdatingNode::UpdatingNode(Node* n_node, SeqRegions* n_incoming_regions, RealNumType n_branch_length, bool n_need_updating, RealNumType n_lh_diff, short int n_failure_count, bool n_delete_regions):TraversingNode(n_node, n_failure_count, n_lh_diff)
{
    incoming_regions = n_incoming_regions;
    branch_length = n_branch_length;
    need_updating = n_need_updating;
    delete_regions = n_delete_regions;
}

UpdatingNode::~UpdatingNode()
{
    if (delete_regions && incoming_regions) delete incoming_regions;
}
