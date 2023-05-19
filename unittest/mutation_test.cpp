#include "gtest/gtest.h"
#include "../alignment/mutation.h"

using namespace cmaple;

/*
 Test the default Mutation() constructor
 */
TEST(Mutation, constructor_default)
{
    Mutation m;
    EXPECT_EQ(m.type, TYPE_N);
    EXPECT_EQ(m.position, 0);
    EXPECT_EQ(m.getLength(), 1);
}

/*
 Test the Mutation(StateType n_type, PositionType n_position) constructor
 */
TEST(Mutation, constructor_1)
{
    // TODO: we may need apply several contraints here, e.g.,
    // type must be a valid type (e.g., ACGT, R, O, etc)
    // position must be positive;
    
    Mutation m1(0, 1);
    Mutation m2(TYPE_N, 267);
    Mutation m3(TYPE_O, 30211);
    Mutation m4(TYPE_R, 10);
    Mutation m5(TYPE_DEL, 1012);
    
    EXPECT_EQ(m1.type, 0);
    EXPECT_EQ(m1.position, 1);
    EXPECT_EQ(m1.getLength(), 1);
    
    EXPECT_EQ(m2.type, TYPE_N);
    EXPECT_EQ(m2.position, 267);
    EXPECT_EQ(m2.getLength(), 1);
    
    EXPECT_EQ(m3.type, TYPE_O);
    EXPECT_EQ(m3.position, 30211);
    EXPECT_EQ(m3.getLength(), 1);
    
    EXPECT_EQ(m4.type, TYPE_R);
    EXPECT_EQ(m4.position, 10);
    EXPECT_EQ(m4.getLength(), 1);
    
    EXPECT_EQ(m5.type, TYPE_DEL);
    EXPECT_EQ(m5.position, 1012);
    EXPECT_EQ(m5.getLength(), 1);
}

/*
 Test the Mutation(StateType n_type, PositionType n_position, LengthTypeLarge n_length) constructor
 */
TEST(Mutation, constructor_2)
{
    // TODO: we may need apply several contraints here, e.g.,
    // type must be a valid type (e.g., ACGT, R, O, etc)
    // position must be positive;
    // length must be positive.
    
    Mutation m1(1, 132, 1);
    Mutation m2(TYPE_N, 5434, 943);
    Mutation m3(TYPE_O, 9563, 1);
    Mutation m4(TYPE_R, 38432, 45);
    Mutation m5(TYPE_DEL, 153, 3);
    
    EXPECT_EQ(m1.type, 1);
    EXPECT_EQ(m1.position, 132);
    EXPECT_EQ(m1.getLength(), 1);
    
    EXPECT_EQ(m2.type, TYPE_N);
    EXPECT_EQ(m2.position, 5434);
    EXPECT_EQ(m2.getLength(), 943);
    
    EXPECT_EQ(m3.type, TYPE_O);
    EXPECT_EQ(m3.position, 9563);
    EXPECT_EQ(m3.getLength(), 1);
    
    EXPECT_EQ(m4.type, TYPE_R);
    EXPECT_EQ(m4.position, 38432);
    EXPECT_EQ(m4.getLength(), 45);
    
    EXPECT_EQ(m5.type, TYPE_DEL);
    EXPECT_EQ(m5.position, 153);
    EXPECT_EQ(m5.getLength(), 3);
    
    // Test constructor with invalid length for A/C/G/T
    EXPECT_EXIT(Mutation invalidMutation(0, 1000, 5), ::testing::ExitedWithCode(2), ".*");
    
    // Test constructor with overflow length
    LengthTypeLarge overflow_length = (std::numeric_limits<LengthType>::max)() + 1;
    EXPECT_EXIT(Mutation invalidMutation(TYPE_R, 1000, overflow_length), ::testing::ExitedWithCode(2), ".*");
}
