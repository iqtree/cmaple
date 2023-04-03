#include "traversingnode.h"
using namespace std;

const Index TraversingNode::getIndex() const
{
    return node_index_;
}

const short int TraversingNode::getFailureCount() const
{
    return failure_count_;
}

void TraversingNode::setFailureCount(const short int failure_count)
{
    failure_count_ = failure_count;
}

void TraversingNode::increaseFailureCount()
{
    ++failure_count_;
}

const RealNumType TraversingNode::getLhDiff() const
{
    return likelihood_diff_;
}
