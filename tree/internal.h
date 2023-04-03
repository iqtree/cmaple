#include "alignment/seqregions.h"

#pragma once

/** An internal node of the tree containing 3 minonodes*/
struct InternalNode
{
    // There partial_lh(s) for the three mininodes, which represent the lower; upper left; upper right likelihoods
    std::array<std::unique_ptr<SeqRegions>,3> partial_lh3_;
    
    // We need to keep track of the index of the neighbor miniNode
    std::array<Index,3> neighbor_index3_;
    
}; // 40 (includes 4 byte padding at the end)

