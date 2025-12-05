#include "seqregionsmem.h"

using namespace cmaple;

/** SeqRegionsWithCount **/

cmaple::SeqRegionsWithCount::SeqRegionsWithCount(
    std::unique_ptr<SeqRegions>&& lower_regions,
    NumSeqsType count):
    lower_regions_(std::move(lower_regions)),
    count_(count) {}

auto cmaple::SeqRegionsWithCount::getCount() const -> const NumSeqsType&
{
    return count_;
}

auto cmaple::SeqRegionsWithCount::increaseCount() -> void
{
    ++count_;
}

auto cmaple::SeqRegionsWithCount::descreaseCount() -> void
{
    assert(count_ >= 0);
    --count_;
    
    // delete the seqregion if nowhere uses it
    if (!count_)
    {
        lower_regions_ = nullptr;
    }
}

auto cmaple::SeqRegionsWithCount::getSeqRegions()
    -> std::unique_ptr<SeqRegions>&
{
    return lower_regions_;
}

auto cmaple::SeqRegionsWithCount::setSeqRegions(
    std::unique_ptr<SeqRegions>&& lower_regions) -> void
{
    lower_regions_ = std::move(lower_regions);
}
