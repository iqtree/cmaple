#include "updatingnode.h"
using namespace std;
using namespace cmaple;

const std::unique_ptr<SeqRegions>& cmaple::UpdatingNode::getIncomingRegions()
    const {
  // return incoming_regions_ if it's not null
  if (incoming_regions_) {
    return incoming_regions_;
  }

  // by deafult, return incoming_regions_ref_
  return incoming_regions_ref_;
}

auto cmaple::UpdatingNode::getBranchLength() const -> const RealNumType {
  return branch_length_;
}

auto cmaple::UpdatingNode::needUpdate() const -> const bool {
  return need_updating_;
}

void cmaple::UpdatingNode::setUpdate(const bool need_update) {
  need_updating_ = need_update;
}
