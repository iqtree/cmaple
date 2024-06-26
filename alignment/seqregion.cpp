//
//  region.cpp
//  alignment
//
//  Created by NhanLT on 31/3/2022.
//

#include "seqregion.h"
#include <iomanip>
using namespace cmaple;

cmaple::SeqRegion::SeqRegion(StateType n_type,
                             PositionType n_position,
                             RealNumType n_plength_observation,
                             RealNumType n_plength_from_root,
                             LHPtrType n_likelihood)
    : Mutation(n_type, n_position),
      plength_observation2node(n_plength_observation),
      plength_observation2root(n_plength_from_root) {
  if (n_likelihood) {
    likelihood = std::move(n_likelihood);
  }
}

cmaple::SeqRegion::SeqRegion(StateType n_type,
                             PositionType n_position,
                             RealNumType n_plength_observation,
                             RealNumType n_plength_from_root,
                             const LHType& n_likelihood)
    : Mutation(n_type, n_position),
      plength_observation2node(n_plength_observation),
      plength_observation2root(n_plength_from_root),
      likelihood(cmaple::make_unique<LHType>(n_likelihood)) {}

cmaple::SeqRegion::SeqRegion(StateType n_type,
                             PositionType n_position,
                             SeqType seq_type,
                             int max_num_states)
    : Mutation(n_type, n_position) {
  convertAmbiguiousState(seq_type, max_num_states);
}

cmaple::SeqRegion::SeqRegion(Mutation* n_mutation,
                             SeqType seq_type,
                             int max_num_states)
    : Mutation(n_mutation->type,
               n_mutation->position + n_mutation->getLength() - 1) {
  convertAmbiguiousState(seq_type, max_num_states);
}

void cmaple::SeqRegion::convertAmbiguiousState(SeqType seq_type,
                                               int max_num_states) {
    assert(type >= 0 && type < TYPE_INVALID);
    assert(max_num_states > 0);
    assert(seq_type == cmaple::SeqRegion::SEQ_DNA || seq_type == cmaple::SeqRegion::SEQ_PROTEIN);
    
  if (type >= TYPE_INVALID) {
    throw std::logic_error("Invalid type of seqregion");
  }

  switch (seq_type) {
    case cmaple::SeqRegion::SEQ_DNA:
      convertAmbiguiousStateDNA(max_num_states);
      break;
    case cmaple::SeqRegion::SEQ_PROTEIN:
      convertAmbiguiousStateAA(max_num_states);
      break;

    case cmaple::SeqRegion::SEQ_AUTO:
    case cmaple::SeqRegion::SEQ_UNKNOWN:
    default:
      throw std::invalid_argument(
          "Sorry! Currently, we only support DNA and Protein data.");
  }
}

void cmaple::SeqRegion::convertAmbiguiousStateAA(int max_num_states) {
  assert(type >= 0 && type < TYPE_INVALID);
  assert(max_num_states > 0);
    
  // do nothing if it is not an ambiguious state
  if (type < max_num_states || type == TYPE_N || type == TYPE_R) {
    return;
  }

  switch (type) {
    case TYPE_DEL:  // convert '-' into type_N
      type = TYPE_N;
      break;
    default:
      throw std::logic_error(
          "Invalid character for a genome entry. Please check and try again!");
  }
}

void cmaple::SeqRegion::convertAmbiguiousStateDNA(int max_num_states) {
  assert(type >= 0 && type < TYPE_INVALID);
  assert(max_num_states > 0);
  
  // do nothing if it is not an ambiguious state
  if (type < max_num_states || type == TYPE_N || type == TYPE_R) {
    return;
  }
  const StateType EIGHT = 8;
  switch (type) {
    case TYPE_DEL:  // convert '-' into type_N
      type = TYPE_N;
      break;
    case 1 + 4 +
        3:  // 'R' -> A or G, Purine
    {
      static constexpr LHType entries{0.5, 0, 0.5, 0};
      computeLhAmbiguity(entries);
      break;
    }
    case 2 + EIGHT +
        3:  // 'Y' -> C or T, Pyrimidine
    {
      static constexpr LHType entries{0, 0.5, 0, 0.5};
      computeLhAmbiguity(entries);
      break;
    }
    case 1 + EIGHT +
        3:  // 'W' -> A or T, Weak
    {
      static constexpr LHType entries{0.5, 0, 0, 0.5};
      computeLhAmbiguity(entries);
      break;
    }
    case 2 + 4 +
        3:  // 'S' -> G or C, Strong
    {
      static constexpr LHType entries{0, 0.5, 0.5, 0};
      computeLhAmbiguity(entries);
      break;
    }
    case 1 + 2 +
        3:  // 'M' -> A or C, Amino
    {
      static constexpr LHType entries{0.5, 0.5, 0, 0};
      computeLhAmbiguity(entries);
      break;
    }
    case 4 + EIGHT +
        3:  // 'K' -> G or T, Keto
    {
      static constexpr LHType entries{0, 0, 0.5, 0.5};
      computeLhAmbiguity(entries);
      break;
    }
    case 2 + 4 + EIGHT +
        3:  // 'B' -> C or G or T
    {
      static constexpr LHType entries{0, 1.0 / 3, 1.0 / 3, 1.0 / 3};
      computeLhAmbiguity(entries);
      break;
    }
    case 1 + 2 + EIGHT +
        3:  // 'H' -> A or C or T
    {
      static constexpr LHType entries{1.0 / 3, 1.0 / 3, 0, 1.0 / 3};
      computeLhAmbiguity(entries);
      break;
    }
    case 1 + 4 + EIGHT +
        3:  // 'D' -> A or G or T
    {
      static constexpr LHType entries{1.0 / 3, 0, 1.0 / 3, 1.0 / 3};
      computeLhAmbiguity(entries);
      break;
    }
    case 1 + 2 + 4 +
        3:  // 'V' -> A or G or C
    {
      static constexpr LHType entries{1.0 / 3, 1.0 / 3, 1.0 / 3, 0};
      computeLhAmbiguity(entries);
      break;
    }
    default:
      throw std::logic_error(
          "Invalid character for a genome entry. Please check and try again!");
  }
}

void cmaple::SeqRegion::computeLhAmbiguity(const LHType& entries) {
  // change type to 'O'
  type = TYPE_O;
  if (!likelihood) {
    likelihood = cmaple::make_unique<LHType>();
  }
  (*likelihood) = entries;
}

auto cmaple::SeqRegion::operator==(const SeqRegion& seqregion_1) const -> bool {
  if (type != seqregion_1.type || position != seqregion_1.position ||
      fabs(plength_observation2node - seqregion_1.plength_observation2node) >
          1e-50 ||
      fabs(plength_observation2root - seqregion_1.plength_observation2root) >
          1e-50) {
    return false;
  }

  if ((likelihood && !seqregion_1.likelihood) ||
      (!likelihood && seqregion_1.likelihood)) {
    return false;
  }

  if (likelihood && seqregion_1.likelihood) {
    for (StateType i = 0; i < seqregion_1.likelihood->size(); i++) {
      if (fabs(likelihood->at(i) - seqregion_1.likelihood->at(i)) > 1e-50) {
        return false;
      }
    }
  }

  return true;
}

auto cmaple::SeqRegion::parseSeqType(const std::string& n_seqtype_str)
    -> cmaple::SeqRegion::SeqType {
  // transform to uppercase
  std::string seqtype_str(n_seqtype_str);
  transform(seqtype_str.begin(), seqtype_str.end(), seqtype_str.begin(),
            ::toupper);
  if (seqtype_str == "DNA") {
    return cmaple::SeqRegion::SEQ_DNA;
  }
  if (seqtype_str == "AA") {
    return cmaple::SeqRegion::SEQ_PROTEIN;
  }
  if (seqtype_str == "AUTO") {
    return cmaple::SeqRegion::SEQ_AUTO;
  }

  // default
  return cmaple::SeqRegion::SEQ_UNKNOWN;
}
