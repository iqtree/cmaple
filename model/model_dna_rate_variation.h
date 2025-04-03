#pragma once

#include "model.h"
#include "modelbase.h"
#include "model_dna.h"


namespace cmaple {

class Tree;

/** Class of DNA evolutionary models with rate variation */
class ModelDNARateVariation : public ModelDNA {
public:
    ModelDNARateVariation(  const cmaple::ModelBase::SubModel sub_model, PositionType _genomeSize, 
                            bool _useSiteRates, cmaple::RealNumType _wtPseudocount);
    virtual ~ModelDNARateVariation();

    void estimateRates(cmaple::Tree* tree);

    void estimateRatePerSite(cmaple::Tree* tree);

    void estimateRatesPerSitePerEntry(cmaple::Tree* tree);

    virtual inline const cmaple::RealNumType *const getMutationMatrix(PositionType i) const override {
        return mutationMatrices + (i * matSize);
    }; 

    virtual inline const cmaple::RealNumType *const getMutationMatrixRow(StateType row, PositionType i) const override {
        return mutationMatrices + (i * matSize) + row_index[row];
    }; 

    virtual inline const cmaple::RealNumType *const getTransposedMutationMatrix(PositionType i) const override {
        return transposedMutationMatrices + (i * matSize);
    };

    virtual inline const cmaple::RealNumType *const getTransposedMutationMatrixRow(StateType row, PositionType i) const override {
        return transposedMutationMatrices + (i * matSize) + row_index[row];
    }; 

    virtual inline cmaple::RealNumType getMutationMatrixEntry(StateType row, StateType column, PositionType i) const override {
        return mutationMatrices[i * matSize + row_index[row] + column];
    }

    virtual inline cmaple::RealNumType getTransposedMutationMatrixEntry(StateType row, StateType column, PositionType i) const override {
        return transposedMutationMatrices[i * matSize + row_index[row] + column];
    }

    virtual inline cmaple::RealNumType getDiagonalMutationMatrixEntry(StateType j, PositionType i) const override {
        return diagonalMutationMatrices[i * num_states_ + j];
    }

    virtual inline cmaple::RealNumType getFreqiFreqjQij(StateType row, StateType column, PositionType i) const override {
        return freqiFreqjQijs[i * matSize + row_index[row] + column];
    }

    virtual inline const cmaple::RealNumType* const getFreqjTransposedijRow(StateType row, PositionType i) const override {
        return freqjTransposedijs + (i * matSize) + row_index[row];
    }

    const cmaple::RealNumType* const getOriginalRateMatrix() {
        return mutation_mat;
    }

    /**
   Update the mutation matrix periodically from the empirical count of mutations
   @return TRUE if the mutation matrix is updated
   @throw  std::logic\_error if any of the following situations occur.
   - the substitution model is unknown/unsupported
   - the reference genome is empty
   */
  virtual bool updateMutationMatEmpirical() override;

  void setAllMatricesToDefault();
  void setMatrixAtPosition(RealNumType* matrix, PositionType i);

  void printMatrix(const RealNumType* matrix, std::ostream* outStream);
  void printCountsAndWaitingTimes(const RealNumType* counts, const RealNumType* waitingTImes, std::ostream* outStream);

private:

    void updateCountsAndWaitingTimesAcrossRoot( PositionType start, PositionType end, 
                                                StateType parentState, StateType childState,
                                                RealNumType distToRoot, RealNumType distToObserved,
                                                RealNumType** waitingTimes, RealNumType** counts,
                                                RealNumType weight = 1.);

    cmaple::PositionType genomeSize;

    cmaple::RealNumType* mutationMatrices = nullptr;
    cmaple::RealNumType* diagonalMutationMatrices = nullptr;
    cmaple::RealNumType* transposedMutationMatrices = nullptr;
    cmaple::RealNumType* freqiFreqjQijs = nullptr;
    cmaple::RealNumType* freqjTransposedijs = nullptr;
    cmaple::RealNumType* rates = nullptr;
    uint16_t matSize;
    bool useSiteRates = false;
    bool ratesEstimated = false;

    cmaple::RealNumType waitingTimePseudoCount;

};
}