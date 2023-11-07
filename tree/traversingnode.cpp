#include "traversingnode.h"
using namespace std;
using namespace cmaple;

auto cmaple::TraversingNode::getIndex() const -> const Index {
  return node_index_;
}

auto cmaple::TraversingNode::getFailureCount() const -> const short int {
  return failure_count_;
}

void cmaple::TraversingNode::setFailureCount(const short int failure_count)
{
    failure_count_ = failure_count;
}

void cmaple::TraversingNode::increaseFailureCount()
{
    ++failure_count_;
}

auto cmaple::TraversingNode::getLhDiff() const -> const RealNumType {
  return likelihood_diff_;
}
