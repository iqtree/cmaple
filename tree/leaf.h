#include "../alignment/seqregions.h"

#pragma once

namespace cmaple
{
    /** A leaf node of the tree*/
    struct LeafNode
    {
        // vector is wasteful here.. (24 bytes)
        std::vector<cmaple::NumSeqsType> less_info_seqs_; // leafs only (contains seq_names)
        // better use 2x4 bytes and store sequence names in a separate (global) vector
        //uint32_t index_start;
        //uint32_t index_end;   // but this requires some sorting... so lets not go there for now
        // a quick improvement would be to implement our own vector
        // which only takes 16 bytes (8 bytes for a pointer to heap memory + 4byte for size + 4 byte for capacity)
        
        /// however we can quickly optimize the seq_name to just a char* and store sequence names as a really long contatenated global string
        // and we can do even better when just storing an index into a combined string or a vector<string> (which would need to be provided externally)
        cmaple::NumSeqsType seq_name_index_;
        
        // index to connect to its neighbor node
        cmaple::Index neighbor_index_;
        
        // The partial_lh (lower likelihood)
        std::unique_ptr<SeqRegions> partial_lh_;
        
        /** constructor */
        // LeafNode() {};
        
        /** constructor */
        LeafNode(cmaple::NumSeqsType new_seq_name_index):seq_name_index_(new_seq_name_index) {};
    }; // size: 40 bytes. no padding! :) -- we could bring this down to 32 bytes if need be
}
