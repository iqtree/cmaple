#include "../alignment/seqregions.h"

using namespace std;
#include "../libraries/ncl/ncl.h"
#include "../libraries/nclextra/modelsblock.h"
#include "../libraries/nclextra/myreader.h"
#include "modelbase.h"
using namespace cmaple;

const std::map<std::string, cmaple::ModelBase::SubModel>
    cmaple::ModelBase::dna_models_mapping = {
        {"JC", cmaple::ModelBase::JC},
        {"GTR", cmaple::ModelBase::GTR},
        {"UNREST", cmaple::ModelBase::UNREST}};

const std::map<std::string, cmaple::ModelBase::SubModel>
    cmaple::ModelBase::aa_models_mapping = {
        {"GTR20", cmaple::ModelBase::GTR20},
        {"NONREV", cmaple::ModelBase::NONREV},
        {"LG", cmaple::ModelBase::LG},
        {"WAG", cmaple::ModelBase::WAG},
        {"JTT", cmaple::ModelBase::JTT},
        {"Q.PFAM", cmaple::ModelBase::Q_PFAM},
        {"Q.BIRD", cmaple::ModelBase::Q_BIRD},
        {"Q.MAMMAL", cmaple::ModelBase::Q_MAMMAL},
        {"Q.INSECT", cmaple::ModelBase::Q_INSECT},
        {"Q.PLANT", cmaple::ModelBase::Q_PLANT},
        {"Q.YEAST", cmaple::ModelBase::Q_YEAST},
        {"JTTDCMUT", cmaple::ModelBase::JTTDCMUT},
        {"DCMUT", cmaple::ModelBase::DCMUT},
        {"VT", cmaple::ModelBase::VT},
        {"PMB", cmaple::ModelBase::PMB},
        {"BLOSUM62", cmaple::ModelBase::BLOSUM62},
        {"DAYHOFF", cmaple::ModelBase::DAYHOFF},
        {"MTREV", cmaple::ModelBase::MTREV},
        {"MTART", cmaple::ModelBase::MTART},
        {"MTZOA", cmaple::ModelBase::MTZOA},
        {"MTMET", cmaple::ModelBase::MTMET},
        {"MTVER", cmaple::ModelBase::MTVER},
        {"MTINV", cmaple::ModelBase::MTINV},
        {"MTMAM", cmaple::ModelBase::MTMAM},
        {"FLAVI", cmaple::ModelBase::FLAVI},
        {"HIVB", cmaple::ModelBase::HIVB},
        {"HIVW", cmaple::ModelBase::HIVW},
        {"FLU", cmaple::ModelBase::FLU},
        {"RTREV", cmaple::ModelBase::RTREV},
        {"CPREV", cmaple::ModelBase::CPREV},
        {"NQ.PFAM", cmaple::ModelBase::NQ_PFAM},
        {"NQ.BIRD", cmaple::ModelBase::NQ_BIRD},
        {"NQ.MAMMAL", cmaple::ModelBase::NQ_MAMMAL},
        {"NQ.INSECT", cmaple::ModelBase::NQ_INSECT},
        {"NQ.PLANT", cmaple::ModelBase::NQ_PLANT},
        {"NQ.YEAST", cmaple::ModelBase::NQ_YEAST}};

cmaple::ModelBase::ModelBase(const cmaple::ModelBase::SubModel n_sub_model,
                             const cmaple::StateType num_states)
    : num_states_(num_states), sub_model(n_sub_model) {}

cmaple::ModelBase::~ModelBase() {
  if (mutation_mat) {
    delete[] mutation_mat;
    mutation_mat = nullptr;
  }

  if (diagonal_mut_mat) {
    delete[] diagonal_mut_mat;
    diagonal_mut_mat = nullptr;
  }

  if (transposed_mut_mat) {
    delete[] transposed_mut_mat;
    transposed_mut_mat = nullptr;
  }

  if (freqi_freqj_qij) {
    delete[] freqi_freqj_qij;
    freqi_freqj_qij = nullptr;
  }

  if (freq_j_transposed_ij) {
    delete[] freq_j_transposed_ij;
    freq_j_transposed_ij = nullptr;
  }

  if (root_freqs) {
    delete[] root_freqs;
    root_freqs = nullptr;
  }

  if (root_log_freqs) {
    delete[] root_log_freqs;
    root_log_freqs = nullptr;
  }

  if (inverse_root_freqs) {
    delete[] inverse_root_freqs;
    inverse_root_freqs = nullptr;
  }

  if (row_index) {
    delete[] row_index;
    row_index = nullptr;
  }

  if (pseu_mutation_count) {
    delete[] pseu_mutation_count;
    pseu_mutation_count = nullptr;
  }
}

void cmaple::ModelBase::initEqualStateFreqs() {
  // init variables
  if (!root_freqs) {
    root_freqs = new RealNumType[num_states_];
  }
  if (!root_log_freqs) {
    root_log_freqs = new RealNumType[num_states_];
  }
  if (!inverse_root_freqs) {
    inverse_root_freqs = new RealNumType[num_states_];
  }

  RealNumType state_freq = 1.0 / num_states_;
  RealNumType log_state_freq = log(state_freq);
  RealNumType inverse_state_freq = num_states_;
  for (auto i = 0; i < num_states_; ++i) {
    root_freqs[i] = state_freq;
    root_log_freqs[i] = log_state_freq;
    inverse_root_freqs[i] = inverse_state_freq;
  }
}

void cmaple::ModelBase::init() {
  // init equal state_freqs
  initEqualStateFreqs();

  // init the mutation matrix
  initMutationMat();
}

void cmaple::ModelBase::extractRefInfo(const Alignment* aln) {
  assert(aln);
    
  if (!fixed_params) {
    // init variables
    if (!root_freqs) {
      root_freqs = new RealNumType[num_states_];
    }
    if (!root_log_freqs) {
      root_log_freqs = new RealNumType[num_states_];
    }
    if (!inverse_root_freqs) {
      inverse_root_freqs = new RealNumType[num_states_];
    }

    // extract root freqs from the ref_seq
    extractRootFreqs(aln);
  }
}

void cmaple::ModelBase::extractRootFreqs(const Alignment* aln) {
  assert(aln);
    
  // init variables
  const vector<StateType>& ref_seq = aln->ref_seq;
  if (!ref_seq.size()) {
    throw std::logic_error("The reference genome is empty!");
  }
  const std::vector<cmaple::StateType>::size_type seq_length = ref_seq.size();

  // init root_freqs
  switch (aln->getSeqType()) {
    case cmaple::SeqRegion::SEQ_PROTEIN:
      resetVec<20>(root_freqs);
      break;

    default:  // dna
      resetVec<4>(root_freqs);
      break;
  }

  // browse all sites in the ref one by one to count bases
  for (std::vector<cmaple::StateType>::size_type i = 0; i < seq_length; ++i) {
    ++root_freqs[ref_seq[i]];
  }

  // update root_freqs and root_log_freqs
  RealNumType inverse_seq_length = 1.0 / seq_length;
  for (StateType i = 0; i < num_states_; ++i) {
    // root_freqs[i] /= seq_length;
    root_freqs[i] *= inverse_seq_length;
    inverse_root_freqs[i] = 1.0 / root_freqs[i];
    root_log_freqs[i] = log(root_freqs[i]);
  }
}

std::string cmaple::ModelBase::exportRootFrequenciesStr() {
  // Handle cases when root_freqs is null
  if (!root_freqs) {
    return "\n \n";
  }

  string output{};
  string header{};

  // Get the seqtype
  const cmaple::SeqRegion::SeqType seqtype = getSeqType();

  for (StateType i = 0; i < num_states_; ++i) {
    header += cmaple::Alignment::convertState2Char(i, seqtype);
    header += "\t\t\t";
    output += convertDoubleToString(root_freqs[i]) + "\t";
  }

  return header + "\n" + output + "\n";
}

std::string cmaple::ModelBase::exportQMatrixStr() {
  // Handle cases when mutation_mat is null
  if (!mutation_mat) {
    return "\n";
  }

  string output{};

  // Get the seqtype
  const cmaple::SeqRegion::SeqType seqtype = getSeqType();

  // generate header
  output += "\t";
  for (StateType i = 0; i < num_states_; ++i) {
    output += cmaple::Alignment::convertState2Char(i, seqtype);
    output += "\t\t\t";
  }
  output += "\n";

  RealNumType* mut_mat_row = mutation_mat;
  for (StateType i = 0; i < num_states_; ++i, mut_mat_row += num_states_) {
    output += cmaple::Alignment::convertState2Char(i, seqtype);
    output += "\t";

    for (StateType j = 0; j < num_states_; ++j) {
      output += convertDoubleToString(mut_mat_row[j]) + "\t";
    }

    output += "\n";
  }

  return output;
}

ModelsBlock* cmaple::ModelBase::readModelsDefinition(
    const char* builtin_models) {
  ModelsBlock* models_block = new ModelsBlock;

  /*try
  {
      // loading internal model definitions
      stringstream in(builtin_mixmodels_definition);
      assert(in && "stringstream is OK");
      NxsReader nexus;
      nexus.Add(models_block);
      MyToken token(in);
      nexus.Execute(token);
  } catch (...) {
      assert(0 && "predefined mixture models not initialized");
  }*/

  try {
    // loading internal protei model definitions
    stringstream in(builtin_models);
    assert(in && "stringstream is OK");
    NxsReader nexus;
    nexus.Add(models_block);
    MyToken token(in);
    nexus.Execute(token);
  } catch (...) {
    assert(0 && "predefined protein models not initialized");
  }

  /*if (params.model_def_file) {
      cout << "Reading model definition file " << params.model_def_file << " ...
  "; MyReader nexus(params.model_def_file); nexus.Add(models_block); MyToken
  token(nexus.inf); nexus.Execute(token); int num_model = 0, num_freq = 0; for
  (ModelsBlock::iterator it = models_block->begin(); it != models_block->end();
  it++) if (it->second.flag & NM_FREQ) num_freq++; else num_model++; cout <<
  num_model << " models and " << num_freq << " frequency vectors loaded" <<
  endl;
  }*/
  return models_block;
}

bool cmaple::ModelBase::readParametersString(string& model_str) {
  // if detect if reading full matrix or half matrix by the first entry
  PositionType end_pos;
  RealNumType d = convert_real_number(model_str.c_str(), end_pos);
  const bool is_reversible = (d >= 0);
  try {
    stringstream in(model_str);
    readRates(in, is_reversible);
    readStateFreq(in);
  } catch (const char* str) {
    throw std::logic_error(str);
  }
  return is_reversible;
}

void cmaple::ModelBase::readStateFreq(istream& in) {
  StateType i;
  for (i = 0; i < num_states_; i++) {
    string tmp_value;
    in >> tmp_value;
    if (!tmp_value.length()) {
      throw "State frequencies could not be read";
    }
    root_freqs[i] = convert_real_number(tmp_value.c_str());
    if (root_freqs[i] < 0.0) {
      throw "Negative state frequencies found";
    }
  }

  RealNumType sum = 0.0;
  for (i = 0; i < num_states_; i++) {
    sum += root_freqs[i];
  }
  if (fabs(sum - 1.0) >= 1e-7) {
    if (cmaple::verbose_mode >= cmaple::VB_MED) {
      outWarning(
          "Normalizing state frequencies so that sum of them equals to 1");
    }
    sum = 1.0 / sum;
    for (i = 0; i < num_states_; i++) {
      // root_freqs[i] /= sum;
      root_freqs[i] *= sum;
    }
  }
}

void cmaple::ModelBase::normalizeQMatrix() {
  assert(root_freqs && mutation_mat);

  RealNumType sum = 0.0;
  RealNumType* mutation_mat_row = mutation_mat;
  for (StateType i = 0; i < num_states_; ++i, mutation_mat_row += num_states_) {
    sum -= mutation_mat_row[i] * root_freqs[i];
  }

  if (sum == 0.0) {
    throw std::logic_error("Empty Q matrix");
  }

  double delta = 1.0 / sum;

  mutation_mat_row = mutation_mat;
  for (StateType i = 0; i < num_states_; ++i, mutation_mat_row += num_states_) {
    for (StateType j = 0; j < num_states_; ++j) {
      mutation_mat_row[j] *= delta;
      // mutation_mat_row[j] /= sum;
    }
  }
}

void cmaple::ModelBase::initPointers() {
  // init row_index
  row_index = new StateType[num_states_ + 1];
  uint16_t start_index = 0;
  for (StateType i = 0; i < num_states_ + 1; i++, start_index += num_states_) {
    row_index[i] = start_index;
  }

  // init root_freqs, root_log_freqs, inverse_root_freqs if they have not yet
  // been initialized
  if (!root_freqs) {
    root_freqs = new RealNumType[num_states_];
  }
  if (!root_log_freqs) {
    root_log_freqs = new RealNumType[num_states_];
  }
  if (!inverse_root_freqs) {
    inverse_root_freqs = new RealNumType[num_states_];
  }

  // init mutation_mat, transposed_mut_mat, and diagonal_mut_mat
  uint16_t mat_size = row_index[num_states_];
  mutation_mat = new RealNumType[mat_size];
  transposed_mut_mat = new RealNumType[mat_size];
  diagonal_mut_mat = new RealNumType[num_states_];
  freqi_freqj_qij = new RealNumType[mat_size];
  freq_j_transposed_ij = new RealNumType[mat_size];
}

void cmaple::ModelBase::updateMutMatbyMutCount() {
  assert(pseu_mutation_count);
  assert(mutation_mat);
    
  RealNumType* pseu_mutation_count_row = pseu_mutation_count;
  RealNumType* mutation_mat_row = mutation_mat;

  // init UNREST model
  if (sub_model == UNREST || sub_model == NONREV) {
    for (StateType i = 0; i < num_states_; ++i,
                   pseu_mutation_count_row += num_states_,
                   mutation_mat_row += num_states_) {
      RealNumType sum_rate = 0;

      for (StateType j = 0; j < num_states_; ++j) {
        if (i != j) {
          RealNumType new_rate =
              // pseu_mutation_count_row[j] / root_freqs[i];
              pseu_mutation_count_row[j] * inverse_root_freqs[i];
          mutation_mat_row[j] = new_rate;
          sum_rate += new_rate;
        }
      }

      // update the diagonal entry
      mutation_mat_row[i] = -sum_rate;
      diagonal_mut_mat[i] = -sum_rate;
    }
  }
  // init GTR model
  else if (sub_model == GTR || sub_model == GTR20) {
    for (StateType i = 0; i < num_states_; ++i,
                   pseu_mutation_count_row += num_states_,
                   mutation_mat_row += num_states_) {
      RealNumType sum_rate = 0;

      for (StateType j = 0; j < num_states_; ++j) {
        if (i != j) {
          mutation_mat_row[j] = (pseu_mutation_count_row[j] +
                                 pseu_mutation_count[row_index[j] + i]) *
                                inverse_root_freqs[i];
                                 /*pseu_mutation_count[row_index[j] + i]) /
                                  root_freqs[i];*/
          sum_rate += mutation_mat_row[j];
        }
      }

      // update the diagonal entry
      mutation_mat_row[i] = -sum_rate;
      diagonal_mut_mat[i] = -sum_rate;
    }
  }
  // handle other model names
  else {
    throw std::logic_error("Unsupported model! Please check and try again!");
  }
}

void cmaple::ModelBase::updatePesudoCount(const Alignment* aln,
                                          const SeqRegions& regions1,
                                          const SeqRegions& regions2) {
  assert(aln);
  assert(regions1.size() > 0);
  assert(regions2.size() > 0);
    
  if (!fixed_params) {
    // init variables
    PositionType pos = 0;
    const SeqRegions& seq1_regions = regions1;
    const SeqRegions& seq2_regions = regions2;
    size_t iseq1 = 0;
    size_t iseq2 = 0;
    const std::vector<cmaple::StateType>& ref_seq = aln->ref_seq;
    const PositionType seq_length = (PositionType) ref_seq.size();

    while (pos < seq_length) {
      PositionType end_pos;

      // get the next shared segment in the two sequences
      SeqRegions::getNextSharedSegment(pos, seq1_regions, seq2_regions, iseq1,
                                       iseq2, end_pos);
      const auto* const seq1_region = &seq1_regions[iseq1];
      const auto* const seq2_region = &seq2_regions[iseq2];

      if (seq1_region->type != seq2_region->type &&
          (seq1_region->type < num_states_ || seq1_region->type == TYPE_R) &&
          (seq2_region->type < num_states_ || seq2_region->type == TYPE_R)) {
        if (seq1_region->type == TYPE_R) {
          pseu_mutation_count[row_index[ref_seq[(std::vector<cmaple::StateType>::size_type) end_pos]] +
                              seq2_region->type] += 1;
        } else if (seq2_region->type == TYPE_R) {
          pseu_mutation_count[row_index[seq1_region->type] +
                              ref_seq[(std::vector<cmaple::StateType>::size_type) end_pos]] += 1;
        } else {
          pseu_mutation_count[row_index[seq1_region->type] +
                              seq2_region->type] += 1;
        }
      }

      // update pos
      pos = end_pos + 1;
    }
  }
}

cmaple::SeqRegion::SeqType cmaple::ModelBase::getSeqType() {
  switch (num_states_) {
    case 20:
      return cmaple::SeqRegion::SEQ_PROTEIN;
      break;
    case 4:
      return cmaple::SeqRegion::SEQ_DNA;
      break;

    default:  // unkown
      return cmaple::SeqRegion::SEQ_UNKNOWN;
      break;
  }
}

cmaple::SeqRegion::SeqType cmaple::ModelBase::detectSeqType(
    const cmaple::ModelBase::SubModel sub_model) {
  // search in the list of dna models
  for (auto& it : dna_models_mapping)
    if (it.second == sub_model)
      return cmaple::SeqRegion::SEQ_DNA;

  // search in the list of protein models
  for (auto& it : aa_models_mapping)
    if (it.second == sub_model)
      return cmaple::SeqRegion::SEQ_PROTEIN;

  // not found
  return cmaple::SeqRegion::SEQ_UNKNOWN;
}

cmaple::ModelBase::SubModel cmaple::ModelBase::parseModel(
    const std::string& n_model_name) {
  // transform to uppercase
  string model_name(n_model_name);
  transform(model_name.begin(), model_name.end(), model_name.begin(),
            ::toupper);

  // parse model
  // handle "DEFAULT"
  if (n_model_name == "DEFAULT") {
    return DEFAULT;
  }
  // search in list of dna models
  auto it = dna_models_mapping.find(model_name);
  if (it != dna_models_mapping.end()) {
    return it->second;
  }
  // search in list of protein models
  it = aa_models_mapping.find(model_name);
  if (it != aa_models_mapping.end()) {
    return it->second;
  }

  // model not found
  return UNKNOWN;
}

std::string cmaple::ModelBase::getModelName() const {
  // Look for the model name from the list of DNA models
  for (auto& it : dna_models_mapping)
    if (it.second == sub_model)
      return it.first;

  // Look for the model name from the list of Protein models
  for (auto& it : aa_models_mapping)
    if (it.second == sub_model)
      return it.first;

  // if not found -> return ""
  return "";
}
