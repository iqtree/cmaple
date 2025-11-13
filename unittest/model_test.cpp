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
    EXPECT_EQ(model_JC.getRootFreq(1), 0.25);
    EXPECT_EQ(model_JC.getRootFreq(3), 0.25);
    EXPECT_EQ(model_JC.getInverseRootFreq(0), 4);
    EXPECT_EQ(model_JC.getInverseRootFreq(2), 4);
    EXPECT_EQ(model_JC.getRootLogFreq(0), log(0.25));
    EXPECT_EQ(model_JC.getRootLogFreq(3), model_JC.getRootLogFreq(0));
    EXPECT_EQ(model_JC.getDiagonalMutationMatrixEntry(1,0), -1);
    EXPECT_EQ(model_JC.getDiagonalMutationMatrixEntry(2,0), -1);
    EXPECT_EQ(model_JC.getMutationMatrix(0)[5], -1);
    EXPECT_EQ(model_JC.getMutationMatrix(0)[10], -1);
    EXPECT_EQ(model_JC.getTransposedMutationMatrix(0)[0], -1);
    EXPECT_EQ(model_JC.getTransposedMutationMatrix(0)[15], -1);
    EXPECT_EQ(model_JC.getFreqiFreqjQij(0, 5, 0), -1);
    EXPECT_EQ(model_JC.getFreqiFreqjQij(0, 10, 0), -1);
    EXPECT_EQ(model_JC.getFreqjTransposedijRow(0, 0)[0], -0.25);
    EXPECT_EQ(model_JC.getFreqjTransposedijRow(3, 0)[3], -0.25);
}

// NOT YET DONE
