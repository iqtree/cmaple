#include "model_aa.h"
#include "../libraries/nclextra/modelsblock.h"
using namespace std;
using namespace cmaple;

/*
    following are definitions for various protein models encoded in a string.
    This string contains the lower triangle of the rate matrix and the state
   frequencies at the end. It should follow the amino acid order: A   R   N   D
   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V Ala Arg Asn Asp
   Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val

   The actual data is pulled in via #include to avoid a thousand lines of 'code'
   here (which will also upset the linter)
*/
const char* cmaple::builtin_prot_models =
#include "model_aa.nexus"
;
/*R"(
#include <model_aa.nexus>
)";*/

cmaple::ModelAA::ModelAA(const cmaple::ModelBase::SubModel sub_model)
    : ModelBase(sub_model, 20) {
  try {
    init();
  } catch (std::logic_error e) {
    throw std::invalid_argument(e.what());
  }
}

cmaple::ModelAA::~ModelAA() {
  if (model_block != nullptr) {
    delete model_block;
  }
}

void cmaple::ModelAA::initMutationMat() {
  // init variable pointers
  initPointers();
  
  assert(row_index);
  assert(root_freqs);
  assert(root_log_freqs);
  assert(inverse_root_freqs);
  assert(mutation_mat);
  assert(transposed_mut_mat);
  assert(diagonal_mut_mat);
  assert(freqi_freqj_qij);
  assert(freq_j_transposed_ij);

  // Get the model_name
  string name_upper = getModelName();
    
  // init the normalized factor
  normalized_factor = 1.0;

  // read model params from string/file
  bool reversible;
  assert(num_states_ == 20);
  model_block = readModelsDefinition(builtin_prot_models);
  if (model_block == nullptr) {
    throw std::logic_error("model_block uninitialized");
  }
  NxsModel* nxs_model = model_block->findModel(name_upper);
  if (nxs_model != nullptr) {
    if (nxs_model->flag != NM_ATOMIC) {
      throw std::logic_error("Invalid protein model name " + name_upper);
    }

    reversible = readParametersString(nxs_model->description);

    // compute root_log_freqs and inverse_root_freqs
    for (StateType i = 0; i < num_states_; ++i) {
      inverse_root_freqs[i] = 1.0 / root_freqs[i];
      root_log_freqs[i] = log(root_freqs[i]);
    }

    // reversible models
    if (reversible) {
      // rescale the lower diagonal rates
      rescaleLowerDiagonalRates();

      // refill the upper diagonal entries by the lower ones
      RealNumType* mutation_mat_ptr = mutation_mat;
      StateType inverse_index = num_states_;
      for (StateType row = 0; row < num_states_; ++row,
                     mutation_mat_ptr += num_states_,
                     inverse_index += num_states_) {
        StateType tmp_index = inverse_index + row;
        for (StateType column = row + 1; column < num_states_;
             ++column, tmp_index += num_states_) {
          mutation_mat_ptr[column] = mutation_mat[tmp_index];
        }
      }

      // compute the diagonal entries
      RealNumType* mutation_mat_row = mutation_mat;
      for (StateType i = 0; i < num_states_;
           ++i, mutation_mat_row += num_states_) {
        RealNumType sum = 0;
        mutation_mat_row[i] = 0;
        for (StateType j = 0; j < num_states_; ++j) {
          sum += mutation_mat_row[j];
        }
        mutation_mat_row[i] = -sum;
      }
    }
    // non-reversible models
    // do nothing

    // normalize the QMatrix so that the expected number of substitutions per
    // site is 1
    normalizeQMatrix();

    // initialize transposed_mut_mat, diagonal_mut_mat, and freqi_freqj_qij
    RealNumType* mutation_mat_row = mutation_mat;
    RealNumType* freqi_freqj_qij_row = freqi_freqj_qij;
    for (StateType i = 0; i < num_states_; ++i, mutation_mat_row += num_states_,
                   freqi_freqj_qij_row += num_states_) {
      for (StateType j = 0; j < num_states_; ++j) {
        // update freqi_freqj_qij
        if (i != j) {
          freqi_freqj_qij_row[j] =
              root_freqs[i] * inverse_root_freqs[j] * mutation_mat_row[j];
              // root_freqs[i] / root_freqs[j] * mutation_mat_row[j];
        } else {
          freqi_freqj_qij_row[j] = mutation_mat_row[j];
        }

        // update the transposed mutation matrix
        transposed_mut_mat[row_index[j] + i] = mutation_mat_row[j];
      }

      // update diagonal
      diagonal_mut_mat[i] = mutation_mat_row[i];
    }

    // pre-compute freq_j_transposed_ij_row matrix to speedup
    RealNumType* transposed_mut_mat_row = transposed_mut_mat;
    RealNumType* freq_j_transposed_ij_row = freq_j_transposed_ij;

    for (StateType i = 0; i < num_states_; ++i,
                   transposed_mut_mat_row += num_states_,
                   freq_j_transposed_ij_row += num_states_) {
      setVecByProduct<20>(freq_j_transposed_ij_row, root_freqs,
                          transposed_mut_mat_row);
    }
  }
  // GTR20 or NONREV
  else if (name_upper.compare("NONREV") == 0 ||
           name_upper.compare("GTR20") == 0) {
    // init pseu_mutation_counts
    string model_rates = "1.0";
    for (StateType i = 0; i < row_index[num_states_] - 1; ++i) {
      model_rates += " 1.0";
    }
    convert_real_numbers(pseu_mutation_count, model_rates);

    updateMutationMat<20>();
  } else {
    throw std::logic_error("Model not found: " + name_upper);
  }
}

void cmaple::ModelAA::rescaleLowerDiagonalRates() {
  assert(mutation_mat);
  assert(num_states_ > 0);
    
  RealNumType max_rate = 0.0;

  RealNumType* mutation_mat_row = mutation_mat;
  for (StateType i = 0; i < num_states_; ++i, mutation_mat_row += num_states_) {
    for (StateType j = 0; j < i; ++j) {
      max_rate =
          max_rate > mutation_mat_row[j] ? max_rate : mutation_mat_row[j];
    }
  }

  const RealNumType AA_SCALE = 10.0;
  RealNumType scaler = AA_SCALE / max_rate;
    
  // record sum as the normalized factor
  normalized_factor /= scaler;

  /* SCALING HAS BEEN RE-INTRODUCED TO RESOLVE NUMERICAL  PROBLEMS */

  mutation_mat_row = mutation_mat;
  for (StateType i = 0; i < num_states_; ++i, mutation_mat_row += num_states_) {
    for (StateType j = 0; j < i; ++j) {
      mutation_mat_row[j] *= scaler;
    }
  }
}

void cmaple::ModelAA::rescaleAllRates() {
  assert(mutation_mat);
  assert(num_states_ > 0);
    
  RealNumType max_rate = 0.0;

  RealNumType* mutation_mat_row = mutation_mat;
  for (StateType i = 0; i < num_states_; ++i, mutation_mat_row += num_states_) {
    for (StateType j = 0; j < num_states_; ++j) {
      max_rate =
          max_rate > mutation_mat_row[j] ? max_rate : mutation_mat_row[j];
    }
  }

  const RealNumType AA_SCALE = 10.0;
  RealNumType scaler = AA_SCALE / max_rate;
    
  // record sum as the normalized factor
  normalized_factor /= scaler;

  /* SCALING HAS BEEN RE-INTRODUCED TO RESOLVE NUMERICAL  PROBLEMS */
  mutation_mat_row = mutation_mat;
  for (StateType i = 0; i < num_states_; ++i, mutation_mat_row += num_states_) {
    for (StateType j = 0; j < num_states_; ++j) {
      mutation_mat_row[j] *= scaler;
    }
  }
}

auto cmaple::ModelAA::updateMutationMatEmpirical() -> bool {    
  // only handle GTR20 or NONREV
  if (!fixed_params && (sub_model == GTR20 || sub_model == NONREV)) {
    return updateMutationMatEmpiricalTemplate<20>();
  }

  // no update -> return false;
  return false;
}

void cmaple::ModelAA::updatePesudoCount(const Alignment* aln,
                                        const SeqRegions& regions1,
                                        const SeqRegions& regions2) {
  assert(aln);
  assert(regions1.size() > 0);
  assert(regions2.size() > 0);
    
  // only handle GTR20 or NONREV
  if (!fixed_params && (sub_model == GTR20 || sub_model == NONREV)) {
    ModelBase::updatePesudoCount(aln, regions1, regions2);
  }
}

void cmaple::ModelAA::extractRootFreqs(const Alignment* aln) {
  assert(aln);
    
  // only extract root freqs for GTR20 or NONREV
  if (sub_model == GTR20 || sub_model == NONREV) {
    ModelBase::extractRootFreqs(aln);
  }
}
