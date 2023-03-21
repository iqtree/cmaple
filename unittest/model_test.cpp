#include "gtest/gtest.h"
#include "model/model_dna.h"

/*
 Test void initMutationMat(const std::string n_model_name, const StateType num_states)
 */
TEST(Model, initMutationMat)
{
    std::unique_ptr<Model> model1 = std::make_unique<ModelDNA>();
    model1->initMutationMat("JC", 4);
    EXPECT_EQ(model1->model_name, "JC");
    EXPECT_EQ(model1->root_freqs[1], 0.25);
    EXPECT_EQ(model1->root_freqs[3], 0.25);
    EXPECT_EQ(model1->inverse_root_freqs[0], 4);
    EXPECT_EQ(model1->inverse_root_freqs[2], 4);
    EXPECT_EQ(model1->root_log_freqs[0], log(0.25));
    EXPECT_EQ(model1->root_log_freqs[3], model1->root_log_freqs[0]);
    EXPECT_EQ(model1->diagonal_mut_mat[1], -1);
    EXPECT_EQ(model1->diagonal_mut_mat[2], -1);
    EXPECT_EQ(model1->mutation_mat[5], -1);
    EXPECT_EQ(model1->mutation_mat[10], -1);
    EXPECT_EQ(model1->transposed_mut_mat[0], -1);
    EXPECT_EQ(model1->transposed_mut_mat[15], -1);
    EXPECT_EQ(model1->freqi_freqj_qij[5], -1);
    EXPECT_EQ(model1->freqi_freqj_qij[10], -1);
    EXPECT_EQ(model1->freq_j_transposed_ij[0], -0.25);
    EXPECT_EQ(model1->freq_j_transposed_ij[15], -0.25);

    
    // Test num_states != 4
    std::unique_ptr<Model> model2 = std::make_unique<ModelDNA>();
    EXPECT_DEATH(model2->initMutationMat("JC",20), ".*");
}

// NOT YET DONE
