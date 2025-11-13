#pragma once

#include "model.h"
#include "modelbase.h"
#include "model_dna.h"


namespace cmaple {

class Tree;

/** Class of DNA evolutionary models with rate variation */
class ModelDNARateVariation : public ModelDNA {
public:
    ModelDNARateVariation(  const cmaple::ModelBase::SubModel sub_model, PositionType _genome_size, 
                            bool _use_site_rates, cmaple::RealNumType _wt_pseudocount, std::string _rates_filename);
    virtual ~ModelDNARateVariation();

    void estimateRates(cmaple::Tree* tree);

    void estimateRatePerSite(cmaple::Tree* tree);

    void estimateRatesPerSitePerEntry(cmaple::Tree* tree);

    virtual inline const cmaple::RealNumType *const getMutationMatrix(PositionType i) const override {
        return mutation_matrices + (i * mat_size);
    }; 

    virtual inline const cmaple::RealNumType *const getMutationMatrixRow(StateType row, PositionType i) const override {
        return mutation_matrices + (i * mat_size) + row_index[row];
    }; 

    virtual inline const cmaple::RealNumType *const getTransposedMutationMatrix(PositionType i) const override {
        return transposed_mutation_matrices + (i * mat_size);
    };

    virtual inline const cmaple::RealNumType *const getTransposedMutationMatrixRow(StateType row, PositionType i) const override {
        return transposed_mutation_matrices + (i * mat_size) + row_index[row];
    }; 

    virtual inline cmaple::RealNumType getMutationMatrixEntry(StateType row, StateType column, PositionType i) const override {
        return mutation_matrices[i * mat_size + row_index[row] + column];
    }

    virtual inline cmaple::RealNumType getTransposedMutationMatrixEntry(StateType row, StateType column, PositionType i) const override {
        return transposed_mutation_matrices[i * mat_size + row_index[row] + column];
    }

    virtual inline cmaple::RealNumType getDiagonalMutationMatrixEntry(StateType j, PositionType i) const override {
        return diagonal_mutation_matrices[i * num_states_ + j];
    }

    virtual inline cmaple::RealNumType getFreqiFreqjQij(StateType row, StateType column, PositionType i) const override {
        return freqi_freqj_Qijs[i * mat_size + row_index[row] + column];
    }

    virtual inline const cmaple::RealNumType* const getFreqjTransposedijRow(StateType row, PositionType i) const override {
        return freqj_transposedijs + (i * mat_size) + row_index[row];
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

  void printMatrix(const RealNumType* matrix, std::ostream* out_stream);
  void printCountsAndWaitingTimes(const RealNumType* counts, const RealNumType* waiting_times, std::ostream* out_stream);

private:

    void updateCountsAndWaitingTimesAcrossRoot( PositionType start, PositionType end, 
                                                StateType parent_state, StateType child_state,
                                                RealNumType dist_to_root, RealNumType dist_to_observed,
                                                RealNumType* waiting_times, RealNumType* counts,
                                                RealNumType weight = 1.);
    
    void readRatesFile();

    cmaple::PositionType genome_size;

    cmaple::RealNumType* mutation_matrices = nullptr;
    cmaple::RealNumType* diagonal_mutation_matrices = nullptr;
    cmaple::RealNumType* transposed_mutation_matrices = nullptr;
    cmaple::RealNumType* freqi_freqj_Qijs = nullptr;
    cmaple::RealNumType* freqj_transposedijs = nullptr;
    cmaple::RealNumType* rates = nullptr;
    uint16_t mat_size;
    bool use_site_rates = false;
    bool rates_estimated = false;

    cmaple::RealNumType waiting_time_pseudocount;

    std::string rates_filename;

};
}