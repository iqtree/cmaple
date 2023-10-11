#include "updatingnode.h"
using namespace std;
using namespace cmaple;

const std::unique_ptr<SeqRegions>& cmaple::UpdatingNode::getIncomingRegions() const
{
    // return incoming_regions_ if it's not null
    if (incoming_regions_) return incoming_regions_;
    
    // by deafult, return incoming_regions_ref_
    return incoming_regions_ref_;
}

const RealNumType cmaple::UpdatingNode::getBranchLength() const
{
    return branch_length_;
}

const bool cmaple::UpdatingNode::needUpdate() const
{
    return need_updating_;
}

void cmaple::UpdatingNode::setUpdate(const bool need_update)
{
    need_updating_ = need_update;
}
