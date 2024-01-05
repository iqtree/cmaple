
#pragma once

#include <cassert>
#include <istream>
#include <string>
#include "../alignment/alignment.h"
#include "../utils/matrix.h"
#include "../utils/tools.h"

class ModelsBlock;  // do not pull in external headers!

namespace cmaple {
class SeqRegions;

/** Base class of evolutionary models */
class ModelBase {
 private:
  /**
   Initialize equal state frequencies
   */
  void initEqualStateFreqs();

  /**
   Get SeqType from num_states
   */
  cmaple::SeqRegion::SeqType getSeqType();

  /**
   Read root state frequencies from string/file
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  void readStateFreq(std::istream& in);

  /**
   Update the mutation rate matrix regarding the pseu_mutation_count

   @throw  std::logic\_error if the substitution model is unknown/unsupported
   */
  void updateMutMatbyMutCount();

 protected:
  // NHANLT: we can change to use unique_ptr(s) instead of normal pointers
   /*! \cond PRIVATE */
   /**
   Model definitions
   */
  ModelsBlock* model_block = nullptr;

  /**
   Initialize the model
   @throw std::logic\_error if the substitution model is unknown/unsupported
   */
  void init();

  /**
   Read model definitions from string/file
   */
  ModelsBlock* readModelsDefinition(const char* builtin_models);

  /**
   Read model parameters
   @return TRUE if the model is reversible
   @throw std::logic\_error if failing to read the parameters
   */
  bool readParametersString(std::string& model_str);

  /**
   Read model's rates from string/file
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  void readRates(std::istream& in, const bool is_reversible);

  /**
   Normalize the Q matrix so that the expected number of subtitution is 1
   @throw std::logic\_error if the Q matrix is empty
   */
  void normalizeQMatrix();

  /**
   Extract root freqs from the reference sequence
   @throw std::logic\_error if the reference genome is empty
   */
  virtual void extractRootFreqs(const Alignment* aln);

  /**
   Init pointers
   */
  void initPointers();

  /**
   Update the mutation rate matrix
   @throw  std::logic\_error if the substitution model is unknown/unsupported
   */
  template <cmaple::StateType num_states>
  void updateMutationMat();

  /**
   Update the mutation matrix periodically from the empirical count of mutations
   (template)
   @return TRUE if the mutation matrix is updated
   @throw  std::logic\_error if the substitution model is unknown/unsupported
   */
  template <cmaple::StateType num_states>
  bool updateMutationMatEmpiricalTemplate();
  /*! \endcond */
    
 public:
  /*!
   * List of substitution models. See [Substitution
   * models](http://www.iqtree.org/doc/Substitution-Models) for references of
   * those models.
   */
  enum SubModel {
    ///@{
    /// See [DNA
    /// models](http://www.iqtree.org/doc/Substitution-Models#dna-model)
    JC,
    GTR,
    UNREST,
    ///@}

    ///@{
    /// See [Protein
    /// models](http://www.iqtree.org/doc/Substitution-Models#protein-models)
    GTR20,
    NONREV,
    LG,
    WAG,
    JTT,
    Q_PFAM,
    Q_BIRD,
    Q_MAMMAL,
    Q_INSECT,
    Q_PLANT,
    Q_YEAST,
    JTTDCMUT,
    DCMUT,
    VT,
    PMB,
    BLOSUM62,
    DAYHOFF,
    MTREV,
    MTART,
    MTZOA,
    MTMET,
    MTVER,
    MTINV,
    MTMAM,
    FLAVI,
    HIVB,
    HIVW,
    FLU,
    RTREV,
    CPREV,
    NQ_PFAM,
    NQ_BIRD,
    NQ_MAMMAL,
    NQ_INSECT,
    NQ_PLANT,
    NQ_YEAST,
    ///@}

    DEFAULT,  ///<  Default - GTR for DNA, and LG for Protein data
    UNKNOWN,  ///<  Unknown model
  };

  /*! \cond PRIVATE */
  /**
   Constructor
   */
  ModelBase() = delete;

  /**
   Constructor
   */
  ModelBase(const SubModel sub_model, const cmaple::StateType num_states);

  /**
   Destructor
   */
   virtual ~ModelBase();

  /**
   * mapping between DNA model names and their enums
   */
  const static std::map<std::string, SubModel> dna_models_mapping;

  /**
   * mapping between Protein model names and their enums
   */
  const static std::map<std::string, SubModel> aa_models_mapping;

  // NHANLT: we can change to use unique_ptr(s) instead of normal pointers in
  // the following
  /**
   Substitution model
   */
  SubModel sub_model;

  /**
   Number of states
   */
  const cmaple::StateType num_states_ = 0;

  /**
   Pseudo mutation count
   */
  cmaple::RealNumType* pseu_mutation_count = nullptr;

  /**
   State frequencies
   */
  cmaple::RealNumType* root_freqs = nullptr;

  /**
   Mutation matrix
   */
  cmaple::RealNumType* mutation_mat = nullptr;

  /**
   log of state frequencies
   */
  cmaple::RealNumType* root_log_freqs = nullptr;

  /**
   diagonal of the mutation matrix
   */
  cmaple::RealNumType* diagonal_mut_mat = nullptr;

  /**
   the transposed matrix of the mutation matrix
   */
  cmaple::RealNumType* transposed_mut_mat = nullptr;

  /**
   the inversed values of state frequencies
   */
  cmaple::RealNumType* inverse_root_freqs = nullptr;

  /**
   freq(i) / freq(j) * Qij
   */
  cmaple::RealNumType* freqi_freqj_qij = nullptr;

  /**
   freq[j] * transposed[i][j]
   */
  cmaple::RealNumType* freq_j_transposed_ij = nullptr;

  /**
   the starting index of row i: i * num_states
   */
  cmaple::StateType* row_index = nullptr;

  /**
   TRUE to keep the model parameters unchanged
   */
  bool fixed_params = false;

  /**
   Get the model name
   */
  std::string getModelName() const;

  /**
   Export state frequencies at root
   */
  std::string exportRootFrequenciesStr();

  /**
   Export Q matrix
   */
  std::string exportQMatrixStr();

  /**
   Extract reference-related info (freqs, log_freqs)

   @throw std::logic\_error if the reference genome is empty
   */
  void extractRefInfo(const Alignment* aln);

  /**
   Init the mutation rate matrix from a model
   @throw std::logic\_error if the substitution model is unknown/unsupported
   */
  virtual void initMutationMat(){}

  /**
   Update the mutation matrix periodically from the empirical count of mutations
   @return TRUE if the mutation matrix is updated
   @throw std::logic\_error if any of the following situations occur.
   - the substitution model is unknown/unsupported
   - the reference genome is empty
   */
  virtual bool updateMutationMatEmpirical() {
    return false;
  }

  /**
   Update pseudocounts from new sample to improve the estimate of the
   substitution rates
   @param regions1 the genome list at the node where the appending happens;
   @param regions2 the genome list for the new sample.
   */
  virtual void updatePesudoCount(const Alignment* aln,
                                 const SeqRegions& regions1,
                                 const SeqRegions& regions2);

  /**
   Detect SeqType from a SubModel enum
   @param[in] sub_model SubModel enum
   */
  static cmaple::SeqRegion::SeqType detectSeqType(const SubModel sub_model);

  /**
   * Parse model from its name in a string
   * @param n_seqtype_str a sequence type in string
   */
  static cmaple::ModelBase::SubModel parseModel(const std::string& model_name);
  /*! \endcond */
};

/*! \cond PRIVATE */
template <StateType num_states>
void cmaple::ModelBase::updateMutationMat() {
  // update Mutation matrix regarding the pseudo muation count
  updateMutMatbyMutCount();

  // compute the total rate regarding the root freqs
  RealNumType total_rate = 0;
  total_rate -= dotProduct<num_states>(root_freqs, diagonal_mut_mat);

  // inverse total_rate
  total_rate = 1.0 / total_rate;

  // normalize the mutation_mat
  RealNumType* mutation_mat_row = mutation_mat;
  RealNumType* freqi_freqj_qij_row = freqi_freqj_qij;
  for (StateType i = 0; i < num_states_; ++i, mutation_mat_row += num_states_,
                 freqi_freqj_qij_row += num_states_) {
    for (StateType j = 0; j < num_states_; ++j) {
      mutation_mat_row[j] *= total_rate;
      // mutation_mat_row[j] /= total_rate;

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

  // pre-compute matrix to speedup
  RealNumType* transposed_mut_mat_row = transposed_mut_mat;
  RealNumType* freq_j_transposed_ij_row = freq_j_transposed_ij;

  for (StateType i = 0; i < num_states_; ++i,
                 transposed_mut_mat_row += num_states_,
                 freq_j_transposed_ij_row += num_states_) {
    setVecByProduct<num_states>(freq_j_transposed_ij_row, root_freqs,
                                transposed_mut_mat_row);
  }
}

template <StateType num_states>
auto cmaple::ModelBase::updateMutationMatEmpiricalTemplate()
    -> bool {
  bool update = false;

  if (!fixed_params) {
    // clone the current mutation matrix
    RealNumType* tmp_mut_mat = new RealNumType[row_index[num_states_]];
    memcpy(tmp_mut_mat, mutation_mat,
           row_index[num_states_] * sizeof(RealNumType));

    // update the mutation matrix regarding the pseu_mutation_count
    updateMutationMat<num_states>();

    // set update = true if the mutation matrix changes more than a threshold
    const RealNumType change_thresh = 1e-3;
    RealNumType sum_change = 0;
    RealNumType* tmp_mut_mat_ptr = tmp_mut_mat;
    RealNumType* mutation_mat_ptr = mutation_mat;
    for (StateType i = 0; i < num_states_;
         ++i, tmp_mut_mat_ptr += num_states_, mutation_mat_ptr += num_states_) {
      for (StateType j = 0; j < num_states_; ++j) {
        sum_change += fabs(tmp_mut_mat_ptr[j] - mutation_mat_ptr[j]);
      }

      sum_change -= fabs(tmp_mut_mat_ptr[i] - mutation_mat_ptr[i]);
    }

    update = sum_change > change_thresh;

    // delete tmp_diagonal_mutation_mat
    delete[] tmp_mut_mat;
  }

  // return update
  return update;
}
/*! \endcond */
}  // namespace cmaple
