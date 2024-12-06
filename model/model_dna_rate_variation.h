#pragma once

#include "model.h"
#include "modelbase.h"
#include "model_dna.h"


namespace cmaple {

class Tree;

/** Class of DNA evolutionary models with rate variation */
class ModelDNARateVariation : public ModelDNA {
public:
    ModelDNARateVariation(const cmaple::ModelBase::SubModel sub_model, PositionType _genomeSize);
    ~ModelDNARateVariation();

    void estimateRatesPerSite(cmaple::Tree* tree);

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

    /**
   Update the mutation matrix periodically from the empirical count of mutations
   @return TRUE if the mutation matrix is updated
   @throw  std::logic\_error if any of the following situations occur.
   - the substitution model is unknown/unsupported
   - the reference genome is empty
   */
  virtual bool updateMutationMatEmpirical() override;

  void printMatrix(const RealNumType* matrix, std::ostream* outStream);

private:

    void updateCountsAndWaitingTimesAcrossRoot( PositionType start, PositionType end, 
                                                StateType parentState, StateType childState,
                                                RealNumType distToRoot, RealNumType distToObserved,
                                                RealNumType** waitingTimes, RealNumType** counts,
                                                RealNumType weight = 1.);

    cmaple::PositionType genomeSize;
    /**
     Position rate multiplier for Mutation matrix
    */
    cmaple::RealNumType* rates = nullptr;
    cmaple::RealNumType* mutationMatrices = nullptr;
    cmaple::RealNumType* diagonalMutationMatrices = nullptr;
    cmaple::RealNumType* transposedMutationMatrices = nullptr;
    uint16_t matSize;

};
}