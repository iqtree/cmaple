
#pragma once

#include "model.h"
#include "modelbase.h"

namespace cmaple {
/** Class of DNA evolutionary models */
class ModelDNA : public ModelBase {
 private:
  /**
   Init the mutation rate matrix from JC model
   */
  void initMutationMatJC();

  /**
   extract root freqs from the reference sequence
   */
  virtual void extractRootFreqs(const Alignment* aln);

 public:
  /**
   Constructor
   @throw std::invalid\_argument if sub\_model is unknown/unsupported
   */
  ModelDNA(const cmaple::ModelBase::SubModel sub_model);

  /**
   Init the mutation rate matrix from a model
   @throw std::logic\_error if the substitution model is unknown/unsupported
   */
  virtual void initMutationMat();

  /**
   Update the mutation matrix periodically from the empirical count of mutations
   @return TRUE if the mutation matrix is updated
   @throw  std::logic\_error if any of the following situations occur.
   - the substitution model is unknown/unsupported
   - the reference genome is empty
   */
  virtual bool updateMutationMatEmpirical();

  /**
   Update pseudocounts from new sample to improve the estimate of the
   substitution rates
   @param regions1 the genome list at the node where the appending happens;
   @param regions2 the genome list for the new sample.
   */
  virtual void updatePesudoCount(const Alignment* aln,
                                 const SeqRegions& regions1,
                                 const SeqRegions& regions2);
};
}  // namespace cmaple
