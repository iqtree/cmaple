#include "traversingnode.h"

#pragma once

namespace cmaple {
/** An extension of node storing more dummy data used for assising root position  */
class RootCandidate : public cmaple::TraversingNode {
 private:
    
    /**
     an updated regions from the direction where we come from,
     */
    std::unique_ptr<SeqRegions> incoming_regions_;

    /**
     a branch length separating the node from this updated regions
     */
    const cmaple::RealNumType branch_length_;

  /**
   The likelihood that needs to be deducted when changing the root position
   */
  const cmaple::RealNumType lh_deducted_;

 public:
  /**
   Constructor
   */
    RootCandidate(const cmaple::Index node_index,
               std::unique_ptr<SeqRegions>&& incoming_regions,
               const cmaple::RealNumType branch_length,
               const cmaple::RealNumType lh_deducted,
               const cmaple::RealNumType lh_diff,
               const short int failure_count)
      : TraversingNode(node_index, failure_count, lh_diff),
    incoming_regions_(std::move(incoming_regions)),
    branch_length_(branch_length),
    lh_deducted_(lh_deducted){}

  /**
   Get lh_deducted_
   */
    const cmaple::RealNumType getLhDeducted() const
    {
        return lh_deducted_;
    };
    
    /**
     Get Incoming_regions_
     */
    const std::unique_ptr<SeqRegions>& getIncomingRegions() const
    {
        return incoming_regions_;
    };

    /**
     Get branch_length_
     */
    const cmaple::RealNumType getBranchLength() const
    {
        return branch_length_;
    };
};
}  // namespace cmaple
