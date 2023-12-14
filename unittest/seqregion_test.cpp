#include "gtest/gtest.h"
#include "../alignment/seqregion.h"

using namespace cmaple;

/*
 Test the default Mutation() constructor
 */
TEST(SeqRegion, constructor_default)
{
    SeqRegion seqregion;
    EXPECT_EQ(seqregion.type, TYPE_N);
    EXPECT_EQ(seqregion.position, 0);
    EXPECT_EQ(seqregion.getLength(), 1);
    EXPECT_EQ(seqregion.plength_observation2node, -1);
    EXPECT_EQ(seqregion.plength_observation2root, -1);
    EXPECT_EQ(seqregion.likelihood, nullptr);
}

/*
 Test SeqRegion(StateType n_type, PositionType n_position, RealNumType n_plength_observation = -1, RealNumType n_plength_from_root = -1, LHPtrType n_likelihood = nullptr) constructor
 */
TEST(SeqRegion, constructor_1)
{
    // TODO: we may need apply several contraints here, e.g.,
    // only SeqRegion type O can have likelihood
    SeqRegion seqregion1(TYPE_R, 382);
    SeqRegion seqregion2(0, 7323, 0);
    SeqRegion seqregion3(3, 14432, 0.432, 0);
    SeqRegion seqregion4(TYPE_N, 654, -1, 0.321);
    
    auto new_lh = cmaple::make_unique<SeqRegion::LHType>(); // = new RealNumType[num_states];
    SeqRegion seqregion5(TYPE_O, 8432, 0, -1, std::move(new_lh));
    
    new_lh = cmaple::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value{0.1,0.3,0.2,0.4};
    (*new_lh) = new_lh_value;
    SeqRegion seqregion6(TYPE_O, 39324, -1, -1, std::move(new_lh));
    
    EXPECT_EQ(seqregion1.type, TYPE_R);
    EXPECT_EQ(seqregion1.position, 382);
    EXPECT_EQ(seqregion1.getLength(), 1);
    EXPECT_EQ(seqregion1.plength_observation2node, -1);
    EXPECT_EQ(seqregion1.plength_observation2root, -1);
    EXPECT_EQ(seqregion1.likelihood, nullptr);
    
    EXPECT_EQ(seqregion2.type, 0);
    EXPECT_EQ(seqregion2.position, 7323);
    EXPECT_EQ(seqregion2.getLength(), 1);
    EXPECT_EQ(seqregion2.plength_observation2node, 0);
    EXPECT_EQ(seqregion2.plength_observation2root, -1);
    EXPECT_EQ(seqregion2.likelihood, nullptr);
    
    EXPECT_EQ(seqregion3.type, 3);
    EXPECT_EQ(seqregion3.position, 14432);
    EXPECT_EQ(seqregion3.getLength(), 1);
    EXPECT_EQ(seqregion3.plength_observation2node, 0.432);
    EXPECT_EQ(seqregion3.plength_observation2root, 0);
    EXPECT_EQ(seqregion3.likelihood, nullptr);
    
    EXPECT_EQ(seqregion4.type, TYPE_N);
    EXPECT_EQ(seqregion4.position, 654);
    EXPECT_EQ(seqregion4.getLength(), 1);
    EXPECT_EQ(seqregion4.plength_observation2node, -1);
    EXPECT_EQ(seqregion4.plength_observation2root, 0.321);
    EXPECT_EQ(seqregion4.likelihood, nullptr);
    
    EXPECT_EQ(seqregion5.type, TYPE_O);
    EXPECT_EQ(seqregion5.position, 8432);
    EXPECT_EQ(seqregion5.getLength(), 1);
    EXPECT_EQ(seqregion5.plength_observation2node, 0);
    EXPECT_EQ(seqregion5.plength_observation2root, -1);
    EXPECT_EQ(seqregion5.getLH(0), 0);
    
    EXPECT_EQ(seqregion6.type, TYPE_O);
    EXPECT_EQ(seqregion6.position, 39324);
    EXPECT_EQ(seqregion6.getLength(), 1);
    EXPECT_EQ(seqregion6.plength_observation2node, -1);
    EXPECT_EQ(seqregion6.plength_observation2root, -1);
    EXPECT_EQ(seqregion6.getLH(0), 0.1);
    EXPECT_EQ(seqregion6.getLH(1), 0.3);
    EXPECT_EQ(seqregion6.getLH(2), 0.2);
    EXPECT_EQ(seqregion6.getLH(3), 0.4);
}

/*
 Test SeqRegion(StateType n_type, PositionType n_position, RealNumType n_plength_observation, RealNumType n_plength_from_root, const LHType& n_likelihood) constructor
 */
TEST(SeqRegion, constructor_2)
{
    // TODO: we may need apply several contraints here, e.g.,
    // only SeqRegion type O can have likelihood
    SeqRegion::LHType new_lh_value1;
    SeqRegion seqregion1(TYPE_O, 18432, 0.4324, -1, new_lh_value1);
    SeqRegion::LHType new_lh_value2{};
    SeqRegion seqregion2(TYPE_O, 32424, -1, -1, new_lh_value2);
    SeqRegion::LHType new_lh_value3{0.2,0.05,0.1,0.65};
    SeqRegion seqregion3(TYPE_O, 6453, 0, -1, new_lh_value3);
    
    EXPECT_EQ(seqregion1.type, TYPE_O);
    EXPECT_EQ(seqregion1.position, 18432);
    EXPECT_EQ(seqregion1.getLength(), 1);
    EXPECT_EQ(seqregion1.plength_observation2node, 0.4324);
    EXPECT_EQ(seqregion1.plength_observation2root, -1);
    EXPECT_TRUE(seqregion1.likelihood != nullptr);
    
    EXPECT_EQ(seqregion2.type, TYPE_O);
    EXPECT_EQ(seqregion2.position, 32424);
    EXPECT_EQ(seqregion2.getLength(), 1);
    EXPECT_EQ(seqregion2.plength_observation2node, -1);
    EXPECT_EQ(seqregion2.plength_observation2root, -1);
    EXPECT_EQ(seqregion2.getLH(0), 0);
    
    EXPECT_EQ(seqregion3.type, TYPE_O);
    EXPECT_EQ(seqregion3.position, 6453);
    EXPECT_EQ(seqregion3.getLength(), 1);
    EXPECT_EQ(seqregion3.plength_observation2node, 0);
    EXPECT_EQ(seqregion3.plength_observation2root, -1);
    EXPECT_EQ(seqregion3.getLH(0), 0.2);
    EXPECT_EQ(seqregion3.getLH(1), 0.05);
    EXPECT_EQ(seqregion3.getLH(2), 0.1);
    EXPECT_EQ(seqregion3.getLH(3), 0.65);
}

/*
 Test SeqRegion(StateType n_type, PositionType n_position, SeqType seq_type, int max_num_states) constructor
 Also test convertAmbiguiousState(), convertAmbiguiousStateDNA(), and computeLhAmbiguity()
 */
TEST(SeqRegion, constructor_3)
{
    // TODO: don't allow TYPE_INVALID or TYPE_O as an input value for n_type in this constructor
    
    SeqRegion seqregion1(3, 3213, cmaple::SeqRegion::SEQ_DNA, 4);
    SeqRegion seqregion2(TYPE_N, 2131, cmaple::SeqRegion::SEQ_DNA, 4);
    SeqRegion seqregion3(TYPE_R, 21, cmaple::SeqRegion::SEQ_DNA, 4);
    SeqRegion seqregion4(TYPE_DEL, 4543, cmaple::SeqRegion::SEQ_DNA, 4);
    
    SeqRegion seqregion5(1+4+3, 12, cmaple::SeqRegion::SEQ_DNA, 4); // {0.5,0,0.5,0}
    SeqRegion seqregion6(1+2+8+3, 3243, cmaple::SeqRegion::SEQ_DNA, 4); // { 1.0/3, 1.0/3, 0, 1.0/3}
    SeqRegion seqregion7(1+2+3, 553, cmaple::SeqRegion::SEQ_DNA, 4); // {0.5, 0.5, 0, 0}
    SeqRegion seqregion8(1+2+4+3, 49, cmaple::SeqRegion::SEQ_DNA, 4); // { 1.0/3, 1.0/3, 1.0/3, 0}
    
    EXPECT_EQ(seqregion1.type, 3);
    EXPECT_EQ(seqregion1.position, 3213);
    EXPECT_EQ(seqregion1.getLength(), 1);
    EXPECT_EQ(seqregion1.plength_observation2node, -1);
    EXPECT_EQ(seqregion1.plength_observation2root, -1);
    EXPECT_EQ(seqregion1.likelihood, nullptr);
    
    EXPECT_EQ(seqregion2.type, TYPE_N);
    EXPECT_EQ(seqregion2.position, 2131);
    EXPECT_EQ(seqregion2.getLength(), 1);
    EXPECT_EQ(seqregion2.plength_observation2node, -1);
    EXPECT_EQ(seqregion2.plength_observation2root, -1);
    EXPECT_EQ(seqregion2.likelihood, nullptr);
    
    EXPECT_EQ(seqregion3.type, TYPE_R);
    EXPECT_EQ(seqregion3.position, 21);
    EXPECT_EQ(seqregion3.getLength(), 1);
    EXPECT_EQ(seqregion3.plength_observation2node, -1);
    EXPECT_EQ(seqregion3.plength_observation2root, -1);
    EXPECT_EQ(seqregion3.likelihood, nullptr);
    
    EXPECT_EQ(seqregion4.type, TYPE_N);
    EXPECT_EQ(seqregion4.position, 4543);
    EXPECT_EQ(seqregion4.getLength(), 1);
    EXPECT_EQ(seqregion4.plength_observation2node, -1);
    EXPECT_EQ(seqregion4.plength_observation2root, -1);
    EXPECT_EQ(seqregion4.likelihood, nullptr);
    
    EXPECT_EQ(seqregion5.type, TYPE_O);
    EXPECT_EQ(seqregion5.position, 12);
    EXPECT_EQ(seqregion5.getLength(), 1);
    EXPECT_EQ(seqregion5.plength_observation2node, -1);
    EXPECT_EQ(seqregion5.plength_observation2root, -1);
    EXPECT_EQ(seqregion5.getLH(0), 0.5);
    EXPECT_EQ(seqregion5.getLH(1), 0);
    EXPECT_EQ(seqregion5.getLH(2), 0.5);
    EXPECT_EQ(seqregion5.getLH(3), 0);
    
    EXPECT_EQ(seqregion6.type, TYPE_O);
    EXPECT_EQ(seqregion6.position, 3243);
    EXPECT_EQ(seqregion6.getLength(), 1);
    EXPECT_EQ(seqregion6.plength_observation2node, -1);
    EXPECT_EQ(seqregion6.plength_observation2root, -1);
    EXPECT_EQ(seqregion6.getLH(0), 1.0 / 3);
    EXPECT_EQ(seqregion6.getLH(1), 1.0 / 3);
    EXPECT_EQ(seqregion6.getLH(2), 0);
    EXPECT_EQ(seqregion6.getLH(3), 1.0 / 3);
    
    EXPECT_EQ(seqregion7.type, TYPE_O);
    EXPECT_EQ(seqregion7.position, 553);
    EXPECT_EQ(seqregion7.getLength(), 1);
    EXPECT_EQ(seqregion7.plength_observation2node, -1);
    EXPECT_EQ(seqregion7.plength_observation2root, -1);
    EXPECT_EQ(seqregion7.getLH(0), 0.5);
    EXPECT_EQ(seqregion7.getLH(1), 0.5);
    EXPECT_EQ(seqregion7.getLH(2), 0);
    EXPECT_EQ(seqregion7.getLH(3), 0);
    
    EXPECT_EQ(seqregion8.type, TYPE_O);
    EXPECT_EQ(seqregion8.position, 49);
    EXPECT_EQ(seqregion8.getLength(), 1);
    EXPECT_EQ(seqregion8.plength_observation2node, -1);
    EXPECT_EQ(seqregion8.plength_observation2root, -1);
    EXPECT_EQ(seqregion8.getLH(0), 1.0 / 3);
    EXPECT_EQ(seqregion8.getLH(1), 1.0 / 3);
    EXPECT_EQ(seqregion8.getLH(2), 1.0 / 3);
    EXPECT_EQ(seqregion8.getLH(3), 0);
}

/*
 Test SeqRegion(Mutation* n_mutation, SeqType seq_type, int max_num_states)
 */
TEST(SeqRegion, constructor_4)
{
    Mutation m1(1, 132, 1);
    SeqRegion seqregion1(&m1, cmaple::SeqRegion::SEQ_DNA, 4);
    Mutation m2(TYPE_N, 5434, 943);
    SeqRegion seqregion2(&m2, cmaple::SeqRegion::SEQ_DNA, 4);
    Mutation m3(TYPE_R, 38432, 45);
    SeqRegion seqregion3(&m3, cmaple::SeqRegion::SEQ_DNA, 4);
    Mutation m4(TYPE_DEL, 153, 3);
    SeqRegion seqregion4(&m4, cmaple::SeqRegion::SEQ_DNA, 4);
    Mutation m5(2+8+3, 43823, 1); // {0, 0.5, 0, 0.5}
    SeqRegion seqregion5(&m5, cmaple::SeqRegion::SEQ_DNA, 4);
    Mutation m6(1+4+8+3, 543, 1); // { 1.0/3, 0, 1.0/3, 1.0/3}
    SeqRegion seqregion6(&m6, cmaple::SeqRegion::SEQ_DNA, 4);
    
    EXPECT_EQ(seqregion1.type, 1);
    EXPECT_EQ(seqregion1.position, 132);
    EXPECT_EQ(seqregion1.getLength(), 1);
    EXPECT_EQ(seqregion1.plength_observation2node, -1);
    EXPECT_EQ(seqregion1.plength_observation2root, -1);
    EXPECT_EQ(seqregion1.likelihood, nullptr);
    
    EXPECT_EQ(seqregion2.type, TYPE_N);
    EXPECT_EQ(seqregion2.position, 5434 + 943 - 1);
    EXPECT_EQ(seqregion2.getLength(), 1);
    EXPECT_EQ(seqregion2.plength_observation2node, -1);
    EXPECT_EQ(seqregion2.plength_observation2root, -1);
    EXPECT_EQ(seqregion2.likelihood, nullptr);
    
    EXPECT_EQ(seqregion3.type, TYPE_R);
    EXPECT_EQ(seqregion3.position, 38432 + 45 - 1);
    EXPECT_EQ(seqregion3.getLength(), 1);
    EXPECT_EQ(seqregion3.plength_observation2node, -1);
    EXPECT_EQ(seqregion3.plength_observation2root, -1);
    EXPECT_EQ(seqregion3.likelihood, nullptr);
    
    EXPECT_EQ(seqregion4.type, TYPE_N);
    EXPECT_EQ(seqregion4.position, 153 + 3 - 1);
    EXPECT_EQ(seqregion4.getLength(), 1);
    EXPECT_EQ(seqregion4.plength_observation2node, -1);
    EXPECT_EQ(seqregion4.plength_observation2root, -1);
    EXPECT_EQ(seqregion4.likelihood, nullptr);
    
    EXPECT_EQ(seqregion5.type, TYPE_O);
    EXPECT_EQ(seqregion5.position, 43823);
    EXPECT_EQ(seqregion5.getLength(), 1);
    EXPECT_EQ(seqregion5.plength_observation2node, -1);
    EXPECT_EQ(seqregion5.plength_observation2root, -1);
    EXPECT_EQ(seqregion5.getLH(0), 0);
    EXPECT_EQ(seqregion5.getLH(1), 0.5);
    EXPECT_EQ(seqregion5.getLH(2), 0);
    EXPECT_EQ(seqregion5.getLH(3), 0.5);
    
    EXPECT_EQ(seqregion6.type, TYPE_O);
    EXPECT_EQ(seqregion6.position, 543);
    EXPECT_EQ(seqregion6.getLength(), 1);
    EXPECT_EQ(seqregion6.plength_observation2node, -1);
    EXPECT_EQ(seqregion6.plength_observation2root, -1);
    EXPECT_EQ(seqregion6.getLH(0), 1.0 / 3);
    EXPECT_EQ(seqregion6.getLH(1), 0);
    EXPECT_EQ(seqregion6.getLH(2), 1.0 / 3);
    EXPECT_EQ(seqregion6.getLH(3), 1.0 / 3);
}

/*
 Test operators: =, ==
 */
TEST(SeqRegion, operators)
{
    Mutation m1(1+4+8+3, 543, 1); // { 1.0/3, 0, 1.0/3, 1.0/3}
    SeqRegion seqregion1(&m1, cmaple::SeqRegion::SEQ_DNA, 4);
    
    auto new_lh = cmaple::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value{1.0/3, 0, 1.0/3, 1.0/3};
    (*new_lh) = new_lh_value;
    SeqRegion seqregion2(TYPE_O, 543, -1, -1, std::move(new_lh));
    EXPECT_TRUE(seqregion1 == seqregion2);
    
    // Test Move CTor
    SeqRegion seqregion3(std::move(seqregion2));
    EXPECT_TRUE(seqregion1 == seqregion3);
    
    // Test Move Assignment
    SeqRegion seqregion4 = std::move(seqregion3);
    EXPECT_TRUE(seqregion1 == seqregion4);
    
    // Test == operator
    SeqRegion seqregion5 = SeqRegion::clone(seqregion4);
    seqregion5.type = 0;
    EXPECT_FALSE(seqregion1 == seqregion5);
    
    SeqRegion seqregion6 = SeqRegion::clone(seqregion4);
    seqregion6.position = 432;
    EXPECT_FALSE(seqregion1 == seqregion6);
    
    SeqRegion seqregion7 = SeqRegion::clone(seqregion4);
    seqregion7.plength_observation2node += 1e-60;
    EXPECT_TRUE(seqregion1 == seqregion7);
    
    SeqRegion seqregion8 = SeqRegion::clone(seqregion4);
    seqregion8.plength_observation2root += 0.1;
    EXPECT_FALSE(seqregion1 == seqregion8);
    
    SeqRegion seqregion9 = SeqRegion::clone(seqregion4);
    seqregion9.likelihood = nullptr;
    EXPECT_FALSE(seqregion1 == seqregion9);
    
    SeqRegion seqregion10 = SeqRegion::clone(seqregion4);
    (*seqregion10.likelihood)[1] += 1e-51;
    EXPECT_TRUE(seqregion1 == seqregion10);
    
    SeqRegion seqregion11 = SeqRegion::clone(seqregion4);
    (*seqregion11.likelihood)[1] += 2e-50;
    EXPECT_FALSE(seqregion1 == seqregion11);
}
