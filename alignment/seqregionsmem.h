#pragma once

#include "seqregions.h"

namespace cmaple {

/** A class stores sequence regions and the number of times that regions are being used
    This class is used to store the temporary lower regions of a sample or a subtree during their placement search
    When considering a specific placement, this class stores the sample/subtree regions wrt the local ref at that placement
 */

    class SeqRegionsWithCount
    {
        private:
            /**
             The lower regions (with branch-mutations integrated)
             */
            std::unique_ptr<SeqRegions> lower_regions_ = nullptr;

            /**
             The number of times this lower regions is being used
             */
            NumSeqsType count_ = 0;
        
            /**
             The vector of local refs were integrated to build
             the current lower regions from its original regions
             */
            std::vector<cmaple::Index> local_ref_vec_;
        
        public:
        
            /**
             *  constructor
             */
            SeqRegionsWithCount() = default;
            
            /**
             *  constructor
             */
            SeqRegionsWithCount(std::unique_ptr<SeqRegions>&& lower_regions,
                                NumSeqsType count = 0);
        
            /**
             Get count_
             */
            auto getCount() const -> const NumSeqsType&;
            
            /**
             Increase count_ by 1
             */
            auto increaseCount() -> void;
            
            /**
             Decrease count_ by 1
             */
            auto descreaseCount() -> void;
            
            /**
             Get lower_regions_
             */
            auto getSeqRegions() -> std::unique_ptr<SeqRegions>&;
            
            /**
             Set lower_regions_
             */
            auto setSeqRegions(std::unique_ptr<SeqRegions>&& lower_regions) -> void;
        
            /**
             Get local_ref_vec_
             */
            auto getLocalRefVec() -> std::vector<cmaple::Index>&;
        
            /**
             Set local_ref_vec_
             */
            auto setLocalRefVec(
                const std::vector<cmaple::Index>& local_ref_vec) -> void;
        
            /**
             Add a ref index into local_ref_vec_
             */
            auto addLocalRefIndex(const cmaple::Index& local_ref_index) -> void;
    };

/** A dedicated memory class to store a vector of SeqRegionsWithCounts
 */

    class SeqRegionsMem: public std::deque<SeqRegionsWithCount>
    {
        public:
        
            /**
             Get seqregions that integrates mutations
             @param[in] in_seqregions the input seq_regions
             @param[in] aln the alignment
             @param[in] mutations the mutations
             @param[in] inverse TRUE to de-integrate mutations. Otherwise, integrate muations
             return a pointer to the seqregions that integrates the mutations
             */
            template <const cmaple::StateType num_states>
            auto getMutIntegratedSeqRegions(
                SeqRegionsWithCount* const in_seqregions,
                const Alignment* aln,
                std::unique_ptr<SeqRegions>& mutations,
                const bool inverse = false)
                -> SeqRegionsWithCount*;
        
            /**
             Same as the above but with a local ref index added
             to keep track of the list of local refs are being used
             */
            template <const cmaple::StateType num_states>
            auto getMutIntegratedSeqRegions(
                SeqRegionsWithCount* const in_seqregions,
                const Alignment* aln,
                const cmaple::Index& local_ref_index,
                std::unique_ptr<SeqRegions>& mutations,
                const bool inverse = false)
            -> SeqRegionsWithCount*;
    };

template <const cmaple::StateType num_states>
auto cmaple::SeqRegionsMem::getMutIntegratedSeqRegions(
    SeqRegionsWithCount* const in_seqregions,
    const Alignment* aln,
    std::unique_ptr<SeqRegions>& mutations,
    const bool inverse)
    -> SeqRegionsWithCount*
{
    // if no mutations -> reuse the input seqregions,
    // and increase its number of use
    if (!mutations)
    {
        in_seqregions->increaseCount();
        return in_seqregions;
    }
    
    assert(mutations);
    
    // otherwise, create a new seqregions and store it into the memory
    // create a new seqregions
    std::unique_ptr<SeqRegions> mut_integrated_seqregions =
    in_seqregions->getSeqRegions()->integrateMutations<num_states>(
                                    mutations, aln, inverse);
    
    // check if any existing slot is empty
    int avai_slot_id = -1;
    for (size_t i = 0; i < size(); ++i)
    {
        if (!at(i).getCount())
        {
            avai_slot_id = i;
            break;
        }
    }
    // if no slot available, create a new one
    if (avai_slot_id == -1)
    {
        push_back(SeqRegionsWithCount());
        avai_slot_id = size() - 1;
    }
    
    // now move the seqregion to the selected slot
    // and increase its use count
    SeqRegionsWithCount* out_seqregions = &at(avai_slot_id);
    out_seqregions->setSeqRegions(std::move(mut_integrated_seqregions));
    out_seqregions->increaseCount();
    
    // inherit the local_ref_vec
    if (in_seqregions->getLocalRefVec().size())
    {
        out_seqregions->setLocalRefVec(in_seqregions->getLocalRefVec());
    }
    
    return out_seqregions;
}
        
template <const cmaple::StateType num_states>
auto cmaple::SeqRegionsMem::getMutIntegratedSeqRegions(
    SeqRegionsWithCount* const in_seqregions,
    const Alignment* aln,
    const cmaple::Index& local_ref_index,
    std::unique_ptr<SeqRegions>& mutations,
    const bool inverse)
    -> SeqRegionsWithCount*
{
    // first call the simple version of this function first
    SeqRegionsWithCount* out_seqregion_ptr =
        getMutIntegratedSeqRegions<num_states>(
        in_seqregions, aln, mutations, inverse);
    
    // if we integrate the mutations/local ref, record it
    if (mutations)
    {
        assert(local_ref_index.getMiniIndex() != UNDEFINED);
        out_seqregion_ptr->addLocalRefIndex(local_ref_index);
    }
    
    return out_seqregion_ptr;
}

}  // namespace cmaple
