#include "gtest/gtest.h"
#include "../model/model_dna.h"

using namespace cmaple;

/*
 Test void initMutationMat(const cmaple::ModelBase::SubModel sub_model, const StateType num_states)
 */
TEST(Model, initMutationMat)
{
    ModelDNA model_JC = ModelDNA(cmaple::ModelBase::JC);
    
    EXPECT_EQ(model_JC.sub_model, cmaple::ModelBase::JC);
    EXPECT_EQ(model_JC.root_freqs[1], 0.25);
    EXPECT_EQ(model_JC.root_freqs[3], 0.25);
    EXPECT_EQ(model_JC.inverse_root_freqs[0], 4);
    EXPECT_EQ(model_JC.inverse_root_freqs[2], 4);
    EXPECT_EQ(model_JC.root_log_freqs[0], log(0.25));
    EXPECT_EQ(model_JC.root_log_freqs[3], model_JC.root_log_freqs[0]);
    EXPECT_EQ(model_JC.diagonal_mut_mat[1], -1);
    EXPECT_EQ(model_JC.diagonal_mut_mat[2], -1);
    EXPECT_EQ(model_JC.mutation_mat[5], -1);
    EXPECT_EQ(model_JC.mutation_mat[10], -1);
    EXPECT_EQ(model_JC.transposed_mut_mat[0], -1);
    EXPECT_EQ(model_JC.transposed_mut_mat[15], -1);
    EXPECT_EQ(model_JC.freqi_freqj_qij[5], -1);
    EXPECT_EQ(model_JC.freqi_freqj_qij[10], -1);
    EXPECT_EQ(model_JC.freq_j_transposed_ij[0], -0.25);
    EXPECT_EQ(model_JC.freq_j_transposed_ij[15], -0.25);
}

// NOT YET DONE
