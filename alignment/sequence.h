#include "seqregions.h"

#ifndef CMAPLE_SEQUENCE_H
#define CMAPLE_SEQUENCE_H

namespace cmaple
{
    class Alignment;
    class SeqRegions;

    /** Class present a sequence */
    class Sequence: public std::vector<Mutation> {
    public:        
        /**
         Name of the sequence
         */
        std::string seq_name;
        
        /**
         *  Sequence constructor
         */
        Sequence() = default;
        
        /**
         *  Sequence constructor
         */
        Sequence(std::string n_seq_name);
        
        /**
         *  Sequence constructor
         */
        Sequence(std::string n_seq_name, vector<Mutation> n_mutations);
        
        /// Move Ctor
        Sequence(Sequence&&) noexcept = default;
        /// Move assignment
        Sequence& operator=(Sequence&&) = default;
        
        /**
         *  Sequence destructor
         */
        ~Sequence() = default;
        
        /**
         Extract the lower likelihood vector (converting a vector of Mutations into a vector of Regions)
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        std::unique_ptr<SeqRegions> getLowerLhVector(const cmaple::PositionType sequence_length, const cmaple::StateType num_states, const cmaple::SeqRegion::SeqType seq_type);
    };
}
#endif
