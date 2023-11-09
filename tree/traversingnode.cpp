#include "traversingnode.h"
using namespace std;
using namespace cmaple;

const Index cmaple::TraversingNode::getIndex() const
{
    return node_index_;
}

const short int cmaple::TraversingNode::getFailureCount() const
{
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

const RealNumType cmaple::TraversingNode::getLhDiff() const
{
    return likelihood_diff_;
}
