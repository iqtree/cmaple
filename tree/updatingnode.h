#include "traversingnode.h"

#pragma once

namespace cmaple {
/** An extension of node storing more dummy data used for updating nodes in the
 * tree  */
class UpdatingNode : public cmaple::TraversingExtNode {
 private:
  /**
   an updated regions from the direction where we come from (taking into account
   the removal of the given subtree),
   */
  std::unique_ptr<SeqRegions> incoming_regions_;

  /**
   a reference to an updated regions from the direction where we come from
   (taking into account the removal of the given subtree),
   */
  std::unique_ptr<SeqRegions>& incoming_regions_ref_;

  /**
   a branch length separating the node from this updated regions (useful for the
   fact that the removal of the subtree changes the branch length at the removal
   node)
   */
  const cmaple::RealNumType branch_length_;

  /**
   a flag says if the updated regions passed needs still updating, or if it has
   become identical to the pre-existing genome list in the tree (which usually
   happens after a while)
   */
  bool need_updating_;

 public:
  /**
   Constructor
   */
  UpdatingNode(const cmaple::Index node_index,
               std::unique_ptr<SeqRegions>&& subtree_regions,
               std::unique_ptr<SeqRegions>&& incoming_regions,
               std::unique_ptr<SeqRegions>& incoming_regions_ref,
               const cmaple::RealNumType branch_length,
               const bool need_updating,
               const cmaple::RealNumType lh_diff,
               const short int failure_count)
      : TraversingExtNode(node_index, failure_count, lh_diff, std::move(subtree_regions)),
        incoming_regions_(std::move(incoming_regions)),
        incoming_regions_ref_(incoming_regions_ref),
        branch_length_(branch_length),
        need_updating_(need_updating){}

  /**
   Get Incoming_regions_[ref_]
   */
  const std::unique_ptr<SeqRegions>& getIncomingRegions() const;

  /**
   Get branch_length_
   */
  const cmaple::RealNumType getBranchLength() const;

  /**
   Get need_updating_
   */
  const bool needUpdate() const;

  /**
   Set need_updating_
   */
  void setUpdate(const bool need_update);
};
}  // namespace cmaple
