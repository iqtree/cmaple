#include "phylonode.h"

#pragma once

namespace cmaple
{
    /** An extension of node storing more dummy data used for browsing all nodes in a stack  */
    class TraversingNode {
    private:
        /**
         Index to a node
         */
        const cmaple::Index node_index_;
        
        /**
         Count of the number of failures when traversing until the current node
         */
        short int failure_count_;
        
        /**
         Cache the likelihood difference computed at the parent node
         */
        cmaple::RealNumType likelihood_diff_;
        
    public:
        /**
         Constructor
         */
        TraversingNode(const cmaple::Index node_index, const short int n_failure_count, const cmaple::RealNumType n_lh_diff):node_index_(node_index), failure_count_(n_failure_count), likelihood_diff_(n_lh_diff) {};
        
        /**
         Get node's index
         */
        const cmaple::Index getIndex() const;
        
        /**
         Get failure_count_
         */
        const short int getFailureCount() const;
        
        /**
         Set failure_count_
         */
        void setFailureCount(const short int failure_count);
        
        /**
         Increase failure_count_ by 1
         */
        void increaseFailureCount();
        
        /**
         Get RealNumType
         */
        const cmaple::RealNumType getLhDiff() const;
    };
}

