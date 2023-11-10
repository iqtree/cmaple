#include "model_dna.h"
using namespace std;
using namespace cmaple;

cmaple::ModelDNA::ModelDNA(const cmaple::ModelBase::SubModel sub_model)
    : ModelBase(sub_model, 4) {
  try {
    init();
  } catch (std::logic_error e) {
    throw std::invalid_argument(e.what());
  }
}

void cmaple::ModelDNA::initMutationMatJC() {
  // update root_freqs
  const RealNumType log_freq = log(0.25);

  for (StateType i = 0; i < num_states_; ++i) {
    root_freqs[i] = 0.25;
    inverse_root_freqs[i] = 4;
    root_log_freqs[i] = log_freq;
  }

  // compute some fixed value for JC model
  const RealNumType jc_rate = 1.0 / 3.0;
  const RealNumType freq_j_jc_rate = 0.25 * jc_rate;

  RealNumType* mutation_mat_row = mutation_mat;
  RealNumType* freqi_freqj_qij_row = freqi_freqj_qij;
  RealNumType* freq_j_transposed_ij_row = freq_j_transposed_ij;
  for (StateType i = 0; i < num_states_; ++i, mutation_mat_row += num_states_,
                 freqi_freqj_qij_row += num_states_,
                 freq_j_transposed_ij_row += num_states_) {
    for (StateType j = 0; j < num_states_; ++j) {
      if (i == j) {
        mutation_mat_row[j] = -1;
        transposed_mut_mat[row_index[j] + i] = -1;
        // update freqi_freqj_qij
        freqi_freqj_qij_row[j] = -1;
        // update freq_j_transposed_ij_row
        freq_j_transposed_ij_row[j] =
            -0.25;  // 0.25 * -1 (freq(j) * transposed_ij)
      } else {
        mutation_mat_row[j] = jc_rate;
        transposed_mut_mat[row_index[j] + i] = jc_rate;
        // update freqi_freqj_qij
        freqi_freqj_qij_row[j] =
            root_freqs[i] * inverse_root_freqs[j] * jc_rate;
        // update freq_j_transposed_ij_row
        freq_j_transposed_ij_row[j] =
            freq_j_jc_rate;  // 0.25 * 0.333 (freq(j) * transposed_ij)
      }
    }

    // update diagonal entries
    diagonal_mut_mat[i] = -1;
  }
}

void cmaple::ModelDNA::initMutationMat() {
  // init variable pointers
  initPointers();

  // for JC model
  if (sub_model == JC) {
    initMutationMatJC();
    // for other models
  } else {
    // validate model
    if (ModelBase::detectSeqType(sub_model) != cmaple::SeqRegion::SEQ_DNA) {
      throw std::logic_error("Invalid or unsupported DNA model: " +
                             getModelName());
    }

    // init pseu_mutation_counts
    string model_rates =
        "0.0 1.0 5.0 2.0 2.0 0.0 1.0 40.0 5.0 2.0 0.0 20.0 2.0 3.0 1.0 0.0";
    convert_real_numbers(pseu_mutation_count, model_rates);

    updateMutationMat<4>();
  }
}

void cmaple::ModelDNA::extractRootFreqs(const Alignment* aln) {
  // Keep state freqs equally for some models (e.g., JC)
  if (sub_model != JC) {
    ModelBase::extractRootFreqs(aln);
  }
}

auto cmaple::ModelDNA::updateMutationMatEmpirical(const Alignment* aln)
    -> bool {
  // don't update JC model parameters
  if (!fixed_params && sub_model != JC) {
    return updateMutationMatEmpiricalTemplate<4>(aln);
  }

  // no update -> return false;
  return false;
}

void cmaple::ModelDNA::updatePesudoCount(const Alignment* aln,
                                         const SeqRegions& regions1,
                                         const SeqRegions& regions2) {
  if (!fixed_params && sub_model != JC) {
    ModelBase::updatePesudoCount(aln, regions1, regions2);
  }
}
