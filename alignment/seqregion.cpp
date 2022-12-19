//
//  region.cpp
//  alignment
//
//  Created by NhanLT on 31/3/2022.
//

#include "seqregion.h"
#include <iomanip>


SeqRegion::SeqRegion(StateType n_type, PositionType n_position, RealNumType n_plength_observation, RealNumType n_plength_from_root, LHPtrType n_likelihood)
  : Mutation(n_type, n_position),
    plength_observation2node(n_plength_observation),
    plength_observation2root(n_plength_from_root)
{
    if (n_likelihood) likelihood = std::move(n_likelihood);
}

SeqRegion::SeqRegion(StateType n_type, PositionType n_position, RealNumType n_plength_observation, RealNumType n_plength_from_root, const LHType& n_likelihood)
: Mutation(n_type, n_position),
  plength_observation2node(n_plength_observation),
  plength_observation2root(n_plength_from_root),
  likelihood(std::make_unique<LHType>(n_likelihood))
{
}


SeqRegion::SeqRegion(StateType n_type, PositionType n_position, SeqType seq_type, int max_num_states)
  : Mutation(n_type, n_position)
{
    convertAmbiguiousState(seq_type, max_num_states);
}

SeqRegion::SeqRegion(Mutation* n_mutation, SeqType seq_type, int max_num_states)
  : Mutation(n_mutation->type, n_mutation->position + n_mutation->getLength() - 1)
{
    convertAmbiguiousState(seq_type, max_num_states);
}


void SeqRegion::convertAmbiguiousState(SeqType seq_type, int max_num_states)
{
    ASSERT(type >= 0 && type < TYPE_INVALID);
    switch (seq_type) {
        case SEQ_DNA:
            convertAmbiguiousStateDNA(max_num_states);
            break;
            
        default:
            outError("Sorry! Currently, we only support DNA data.");
            break;
    }
}

void SeqRegion::convertAmbiguiousStateDNA(int max_num_states)
{
    // do nothing if it is not an ambiguious state
    if (type < max_num_states || type == TYPE_N || type == TYPE_R)
        return;
    
    switch (type) {
        case TYPE_DEL: // convert '-' into type_N
            type = TYPE_N;
            break;
        case 1+4+3: // 'R' -> A or G, Purine
        {
            static constexpr LHType entries{0.5,0,0.5,0};
            computeLhAmbiguity(entries);
            break;
        }
        case 2+8+3: // 'Y' -> C or T, Pyrimidine
        {
          static constexpr LHType  entries{0, 0.5, 0, 0.5};
            computeLhAmbiguity(entries);
            break;
        }
        case 1+8+3: // 'W' -> A or T, Weak
        {
          static constexpr LHType  entries{0.5, 0, 0, 0.5};
            computeLhAmbiguity(entries);
            break;
        }
        case 2+4+3: // 'S' -> G or C, Strong
        {
          static constexpr LHType  entries{0, 0.5, 0.5, 0};
            computeLhAmbiguity(entries);
            break;
        }
        case 1+2+3:  // 'M' -> A or C, Amino
        {
          static constexpr LHType  entries{0.5, 0.5, 0, 0};
            computeLhAmbiguity(entries);
            break;
        }
        case 4+8+3: // 'K' -> G or T, Keto
        {
          static constexpr LHType  entries{0, 0, 0.5, 0.5};
            computeLhAmbiguity(entries);
            break;
        }
        case 2+4+8+3: // 'B' -> C or G or T
        {
          static constexpr LHType  entries{0, 1.0/3, 1.0/3, 1.0/3};
            computeLhAmbiguity(entries);
            break;
        }
        case 1+2+8+3: // 'H' -> A or C or T
        {
          static constexpr LHType  entries{ 1.0/3, 1.0/3, 0, 1.0/3};
            computeLhAmbiguity(entries);
            break;
        }
        case 1+4+8+3: // 'D' -> A or G or T
        {
          static constexpr LHType  entries{ 1.0/3, 0, 1.0/3, 1.0/3};
            computeLhAmbiguity(entries);
            break;
        }
        case 1+2+4+3: // 'V' -> A or G or C
        {
          static constexpr LHType  entries{ 1.0/3, 1.0/3, 1.0/3, 0};
            computeLhAmbiguity(entries);
            break;
        }
        default:
            outError("Invalid character for a genome entry. Please check and try again!");
    }
}

void SeqRegion::computeLhAmbiguity(const LHType &entries)
{
    // change type to 'O'
    type = TYPE_O;
    if (!likelihood) likelihood = std::make_unique<LHType>();
    (*likelihood) = entries;
}

void SeqRegion::writeConstructionCodes(const std::string regions_name, std::ofstream& out, const StateType num_states) const
{
    // export lilkehood
    std::string lh_str = "";
    if (likelihood)
    {
        out << "auto " << regions_name << "_new_lh" << position << " = std::make_unique<SeqRegion::LHType>();" << std::endl;
        out << "auto& " << regions_name << "_new_lh_value" << position << " = *" << regions_name << "_new_lh" << position << ";" << std::endl;
        for (StateType i = 0; i < num_states; ++i)
            out << regions_name << "_new_lh_value" << position << "[" << convertIntToString(i) << "] = " << std::setprecision(50) << (*likelihood)[i] << ";" << std::endl;
        
        lh_str = ",std::move(" + regions_name + "_new_lh" + convertPosTypeToString(position) + ")";
    }
    
    out << regions_name << "->emplace_back(" << convertIntToString(type) << "," << position << "," << std::setprecision(50) << plength_observation2node << "," << std::setprecision(50) << plength_observation2root << lh_str << ");" << std::endl;
}

bool SeqRegion::operator==(const SeqRegion& seqregion_1) const
{
    if (type != seqregion_1.type
        || position != seqregion_1.position
        || fabs(plength_observation2node - seqregion_1.plength_observation2node) > 1e-50
        || fabs (plength_observation2root - seqregion_1.plength_observation2root) > 1e-50)
        return false;
    
    if ((likelihood && !seqregion_1.likelihood)
        || (!likelihood && seqregion_1.likelihood))
        return false;
    
    if (likelihood && seqregion_1.likelihood)
    {
        for (StateType i = 0; i < seqregion_1.likelihood->size(); i++)
            if (fabs(likelihood->at(i) - seqregion_1.likelihood->at(i)) > 1e-50)
                return false;
    }
    
    return true;
}
