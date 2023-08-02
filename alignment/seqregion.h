#include "mutation.h"

#pragma once
#include <array>
#include <memory> // for unique_ptr

namespace cmaple
{
    /** A region in a sequence */
    class SeqRegion: public Mutation {
    public:
        using LHType = std::array<cmaple::RealNumType, NUM_STATES>;
        using LHPtrType = std::unique_ptr<LHType>;
    private:
        /**
         *  Compute the relative likelihood for an ambiguous state
         */
        void computeLhAmbiguity(const LHType &entries);
        
        /**
         *  Convert Ambiguious state into typical states: nucleotides/amino-acids...; N; O; R
         *  @throw std::invalid\_argument if seq\_type is unknown/unsupported
         *  @throw std::logic\_error if invalid type of seqregion found
         */
        void convertAmbiguiousState(cmaple::SeqType seq_type, int max_num_states);
        
        /**
         *  Convert Ambiguious state (of DNA data) into typical states: A/C/G/T; N; O; R
         *  @throw std::logic\_error if invalid type of seqregion found
         */
        void convertAmbiguiousStateDNA(int max_num_states);
        
        /**
         *  Convert Ambiguious state (of AA data) into typical states: amino-acid; N; O; R
         *  @throw std::logic\_error if invalid type of seqregion found
         */
        void convertAmbiguiousStateAA(int max_num_states);
        
    public:
        /**
         Length of the path between the current phylo node and the node where the likelihood is calculated
         */
        cmaple::RealNumType plength_observation2node = -1;
        
        /**
         Distance separates the observation at root from the observation at the current node
         To take into account that the observation might have occurred on the other side of the phylogeny with respect to the root;
         for example if the entry is (type = 1, pos = 234, plength_observation = 0.0001, plength_from_root = 0.0002) then this means that a "C" was observed at genome position 234, and the observation is separated from the root by a distance of 0.0001, while a distance of 0.0002 separates the root from the current node (or position along a branch) considered.
         */
        cmaple::RealNumType plength_observation2root = -1;
        
        /**
         The relative partial likelihood
         */
        LHPtrType likelihood;
        cmaple::RealNumType getLH(int pos) const
        {
            return (*likelihood)[pos];
        }
        
        /**
         *  Region constructor
         */
        SeqRegion() = default;
        
        /**
         *  Region constructor
         */
        SeqRegion(cmaple::StateType n_type, cmaple::PositionType n_position, cmaple::RealNumType n_plength_observation = -1, cmaple::RealNumType n_plength_from_root = -1, LHPtrType n_likelihood = nullptr);
        
        /**
         *  Region constructor
         */
        SeqRegion(cmaple::StateType n_type, cmaple::PositionType n_position, cmaple::RealNumType n_plength_observation, cmaple::RealNumType n_plength_from_root, const LHType& n_likelihood);
        
        /**
         *  Region constructor
         *  @throw std::invalid\_argument if any of the following situations occur.
         *  - n\_type is invalid
         *  - seq\_type is unknown/unsupported
         */
        SeqRegion(cmaple::StateType n_type, cmaple::PositionType n_position, cmaple::SeqType seq_type, int max_num_states);
        
        /**
         *  Region constructor
         *  @throw std::invalid\_argument if any of the following situations occur.
         *  - the type of n_mutation  is invalid
         *  - seq\_type is unknown/unsupported
         */
        SeqRegion(Mutation* n_mutation, cmaple::SeqType seq_type, int max_num_states);
        
        /**
         *  Region constructor
         */
        //SeqRegion(SeqRegion* region, StateType num_states, bool copy_likelihood = true);
        
        /// Move CTor
        SeqRegion(SeqRegion&& region) noexcept = default;
        /// Move Assignment
        SeqRegion& operator=(SeqRegion&& region) noexcept = default;
        
        static inline SeqRegion clone(const SeqRegion& from)
        {
            if (from.likelihood)
                return SeqRegion(from.type, from.position, from.plength_observation2node, from.plength_observation2root, *from.likelihood);
            else
                return SeqRegion(from.type, from.position, from.plength_observation2node, from.plength_observation2root);
        }
        
        
        /**
         *  Region destructor
         */
        ~SeqRegion() = default;
        
        /**
         For testing only, export codes to re-contruct this seqregions
         */
        void writeConstructionCodes(const std::string regions_name, std::ofstream& out, const cmaple::StateType num_states) const;
        
        /**
         Compare two regions
         */
        bool operator==(const SeqRegion& seqregion_1) const;
    };
}
