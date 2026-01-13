
#pragma once

#include "../alignment/seqregions.h"
#include "modelbase.h"

namespace cmaple {
/** Class represents a substitution model */
class Model {
 public:
  /*!
   A structure to store model parameters
   */
  struct ModelParams {
    /*!
     Name of the model in string
     */
    std::string model_name;

    /*!
     State frequencies in string
     */
    std::string state_freqs;

    /*!
     Mutation rates in string
     */
    std::string mut_rates;
  };

  /*! \brief Constructor from a model name
   * @param[in] sub_model a substitution model. Default: DEFAULT - GTR for DNA,
   * and LG for Protein data. List of supported models: <br>**DNA models**: JC,
   * GTR, UNREST;
   * <br>**Protein models**: GTR20, NONREV, LG, WAG, JTT, Q_PFAM, Q_BIRD,
   * Q_MAMMAL, Q_INSECT, Q_PLANT, Q_YEAST, JTTDCMUT, DCMUT, VT, PMB, BLOSUM62,
   * DAYHOFF, MTREV, MTART, MTZOA, MTMET, MTVER, MTINV, MTMAM, FLAVI, HIVB,
   * HIVW, FLU, RTREV, CPREV, NQ_PFAM, NQ_BIRD, NQ_MAMMAL, NQ_INSECT, NQ_PLANT,
   * NQ_YEAST; <br> <em> See [**Substitution
   * models**](http://www.iqtree.org/doc/Substitution-Models) for references of
   * those models.</em>
   * @param[in] seqtype Data type of sequences (optional): SEQ_DNA (nucleotide
   * data), SEQ_PROTEIN (amino acid data), or SEQ_AUTO (auto detection)
   * @param[in] rate_variation True to use rate variation (optional). Default: False;
   * @param[in] site_specific_rates True to use site-specific rates (optional). Default: False;
   * @param[in] seq_length Length of the genomes (optional). Default: 0;
   * @param[in] wt_pseudocount Pseudocount used for waiting times
   * when estimating site-specific rate matrices (optional). Default: 0.1;
   * @param[in] rates_filename Name of the file that contains rates (optional). Default: "";
   * @throw std::invalid\_argument if any of the following situations occur.
   * - sub\_model is unknown/unsupported
   * - sub_model is DEFAULT and seqtype is SEQ_AUTO
   */
  Model(
      const cmaple::ModelBase::SubModel sub_model = cmaple::ModelBase::DEFAULT,
      const cmaple::SeqRegion::SeqType seqtype = cmaple::SeqRegion::SEQ_AUTO,
      const bool rate_variation = false,
      const bool site_specific_rates = false,
      const cmaple::PositionType seq_length = 0,
      const cmaple::RealNumType wt_pseudocount = 0.1,
      const std::string& rates_filename = "");

  /*! \brief Destructor
   */
  ~Model();

  /*! \brief Keep the model parameters unchanged. This API is useful when using
   * one model for multiple trees - after the model parameters are estimated
   * according to a tree, users can keep them unchanged, then use the model for
   * other trees.
   * @param[in] fixed_model_params TRUE to keep the model parameters unchanged
   */
  void fixParameters(const bool& fixed_model_params);

  /*! \brief Export the substitution model and its parameters to a ModelParams
   * structure
   * @return A ModelParams structure.
   */
  ModelParams getParams();

  // TODO: allow users to specify model parameters

  // declare Tree as a friend class
  friend class Tree;

 private:
  /**
   A base instance of Model
   */
  ModelBase* model_base;

  /**
   Using rate variation.
   Currently only for DNA model.
   */
  bool rate_variation = false;
};
}  // namespace cmaple
