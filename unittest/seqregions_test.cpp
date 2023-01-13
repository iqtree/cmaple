#include "gtest/gtest.h"
#include "alignment/seqregions.h"

/*
 Test addNonConsecutiveRRegion()
 */
TEST(SeqRegions, addNonConsecutiveRRegion)
{
    RealNumType threshold_prob = 1e-8;
    SeqRegions seqregions;
    
    SeqRegions::addNonConsecutiveRRegion(&seqregions, 0, -1, -1, 24, threshold_prob);
    EXPECT_EQ(seqregions.size(), 1);
    EXPECT_EQ(seqregions.back().position, 24);
    EXPECT_EQ(seqregions.back().plength_observation2node, -1);
    EXPECT_EQ(seqregions.back().plength_observation2root, -1);
    
    SeqRegions::addNonConsecutiveRRegion(&seqregions, TYPE_R, -1, -1, 113, threshold_prob);
    EXPECT_EQ(seqregions.size(), 2);
    EXPECT_EQ(seqregions.back().position, 113);
    EXPECT_EQ(seqregions.back().plength_observation2node, -1);
    EXPECT_EQ(seqregions.back().plength_observation2root, -1);
    
    SeqRegions::addNonConsecutiveRRegion(&seqregions, TYPE_R, -1, -1, 159, threshold_prob);
    EXPECT_EQ(seqregions.size(), 2); // merge two consecotive R
    EXPECT_EQ(seqregions.back().position, 159);
    EXPECT_EQ(seqregions.back().plength_observation2node, -1);
    EXPECT_EQ(seqregions.back().plength_observation2root, -1);
    
    SeqRegions::addNonConsecutiveRRegion(&seqregions, 3, -1, 0, 160, threshold_prob);
    EXPECT_EQ(seqregions.size(), 3);
    EXPECT_EQ(seqregions.back().position, 160);
    EXPECT_EQ(seqregions.back().plength_observation2node, -1);
    EXPECT_EQ(seqregions.back().plength_observation2root, 0);
    
    SeqRegions::addNonConsecutiveRRegion(&seqregions, TYPE_R, -1, -1, 223, threshold_prob);
    EXPECT_EQ(seqregions.size(), 4);
    EXPECT_EQ(seqregions.back().position, 223);
    EXPECT_EQ(seqregions.back().plength_observation2node, -1);
    EXPECT_EQ(seqregions.back().plength_observation2root, -1);
    
    SeqRegions::addNonConsecutiveRRegion(&seqregions, TYPE_R, -1, 1e-10, 240, threshold_prob);
    EXPECT_EQ(seqregions.size(), 5);
    EXPECT_EQ(seqregions.back().position, 240);
    EXPECT_EQ(seqregions.back().plength_observation2node, -1);
    EXPECT_EQ(seqregions.back().plength_observation2root, 1e-10);
    
    SeqRegions::addNonConsecutiveRRegion(&seqregions, TYPE_R, -1, 1e-11, 264, threshold_prob);
    EXPECT_EQ(seqregions.size(), 5); // merge two consecotive R
    EXPECT_EQ(seqregions.back().position, 264);
    EXPECT_EQ(seqregions.back().plength_observation2node, -1);
    EXPECT_EQ(seqregions.back().plength_observation2root, 1e-11);
    
    SeqRegions::addNonConsecutiveRRegion(&seqregions, TYPE_R, 0, 1e-11, 289, threshold_prob);
    EXPECT_EQ(seqregions.size(), 6);
    EXPECT_EQ(seqregions.back().position, 289);
    EXPECT_EQ(seqregions.back().plength_observation2node, 0);
    EXPECT_EQ(seqregions.back().plength_observation2root, 1e-11);
    
    SeqRegions::addNonConsecutiveRRegion(&seqregions, TYPE_R, 1e-9, 0, 299, threshold_prob);
    EXPECT_EQ(seqregions.size(), 6); // merge two consecotive R
    EXPECT_EQ(seqregions.back().position, 299);
    EXPECT_EQ(seqregions.back().plength_observation2node, 1e-9);
    EXPECT_EQ(seqregions.back().plength_observation2root, 0);
    
    SeqRegions::addNonConsecutiveRRegion(&seqregions, TYPE_N, 1e-9, 0.2, 311, threshold_prob);
    EXPECT_EQ(seqregions.size(), 7);
    EXPECT_EQ(seqregions.back().position, 311);
    EXPECT_EQ(seqregions.back().plength_observation2node, 1e-9);
    EXPECT_EQ(seqregions.back().plength_observation2root, 0.2);
    
    SeqRegions::addNonConsecutiveRRegion(&seqregions, TYPE_R, 1e-9, 0, 324, threshold_prob);
    EXPECT_EQ(seqregions.size(), 8);
    EXPECT_EQ(seqregions.back().position, 324);
    EXPECT_EQ(seqregions.back().plength_observation2node, 1e-9);
    EXPECT_EQ(seqregions.back().plength_observation2root, 0);
    
    SeqRegions::addNonConsecutiveRRegion(&seqregions, TYPE_R, 1e-10, 1e-11, 324, threshold_prob);
    EXPECT_EQ(seqregions.size(), 8); // merge two consecotive R
    EXPECT_EQ(seqregions.back().position, 324);
    EXPECT_EQ(seqregions.back().plength_observation2node, 1e-10);
    EXPECT_EQ(seqregions.back().plength_observation2root, 1e-11);
    
    SeqRegions::addNonConsecutiveRRegion(&seqregions, TYPE_R, 1e-15, 1e-10, 356, threshold_prob);
    EXPECT_EQ(seqregions.size(), 8); // merge two consecotive R
    EXPECT_EQ(seqregions.back().position, 356);
    EXPECT_EQ(seqregions.back().plength_observation2node, 1e-15);
    EXPECT_EQ(seqregions.back().plength_observation2root, 1e-10);
    
    SeqRegions::addNonConsecutiveRRegion(&seqregions, TYPE_R, 1e-15, 1e-300, 382, threshold_prob);
    EXPECT_EQ(seqregions.size(), 8); // merge two consecotive R
    EXPECT_EQ(seqregions.back().position, 382);
    EXPECT_EQ(seqregions.back().plength_observation2node, 1e-15);
    EXPECT_EQ(seqregions.back().plength_observation2root, 1e-300);
    
    SeqRegions::addNonConsecutiveRRegion(&seqregions, TYPE_R, 1e-5, 1e-301, 545, threshold_prob);
    EXPECT_EQ(seqregions.size(), 9);
    EXPECT_EQ(seqregions.back().position, 545);
    EXPECT_EQ(seqregions.back().plength_observation2node, 1e-5);
    EXPECT_EQ(seqregions.back().plength_observation2root, 1e-301);
}

/*
    Generate testing data (seqregions1, seqregions2)
 */
void genTestData(SeqRegions& seqregions1, SeqRegions& seqregions2)
{
    seqregions1.emplace_back(TYPE_R, 382, -1, 0.321);
    seqregions1.emplace_back(TYPE_N, 654);
    seqregions1.emplace_back(0, 655, 0);
    seqregions1.emplace_back(TYPE_N, 1431);
    seqregions1.emplace_back(3, 1432, 0.432, 0);
    seqregions1.emplace_back(TYPE_N, 2431);
    auto new_lh = std::make_unique<SeqRegion::LHType>();
    seqregions1.emplace_back(TYPE_O, 2432, 0, -1, std::move(new_lh));
    seqregions1.emplace_back(TYPE_N, 3381);
    seqregions1.emplace_back(TYPE_R, 3500, 0, 0.1321);
    
    seqregions2.emplace_back(TYPE_N, 15);
    seqregions2.emplace_back(TYPE_R, 381);
    seqregions2.emplace_back(1, 382, -1, 0.321);
    seqregions2.emplace_back(TYPE_R, 984, -1, 0);
    new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value1{0.2,0.2,0.2,0.4};
    (*new_lh) = new_lh_value1;
    seqregions2.emplace_back(TYPE_O, 985, 1e-3, -1, std::move(new_lh));
    seqregions2.emplace_back(TYPE_N, 1210);
    seqregions2.emplace_back(2, 1211, -1, 0);
    seqregions2.emplace_back(TYPE_R, 2432, -1, 1e-50);
    seqregions2.emplace_back(TYPE_R, 3500, 0, 0.1321);
}

/*
 Test getNextSharedSegment(PositionType current_pos, const SeqRegions& seq1_region, const SeqRegions& seq2_region, size_t& i1, size_t& i2, PositionType &end_pos)
 */
TEST(SeqRegions, getNextSharedSegment)
{
    PositionType current_pos{0}, end_pos;
    size_t i1{0}, i2{0};
    SeqRegions seqregions1, seqregions2;
    genTestData(seqregions1, seqregions2);
    
    // test getNextSharedSegment()
    SeqRegions::getNextSharedSegment(current_pos, seqregions1, seqregions2, i1, i2, end_pos);
    EXPECT_EQ(i1, 0);
    EXPECT_EQ(i2, 0);
    EXPECT_EQ(end_pos, 15);
    current_pos = end_pos + 1;
    
    SeqRegions::getNextSharedSegment(current_pos, seqregions1, seqregions2, i1, i2, end_pos);
    EXPECT_EQ(i1, 0);
    EXPECT_EQ(i2, 1);
    EXPECT_EQ(end_pos, 381);
    current_pos = end_pos + 1;
    
    SeqRegions::getNextSharedSegment(current_pos, seqregions1, seqregions2, i1, i2, end_pos);
    EXPECT_EQ(i1, 0);
    EXPECT_EQ(i2, 2);
    EXPECT_EQ(end_pos, 382);
    current_pos = end_pos + 1;
    
    SeqRegions::getNextSharedSegment(current_pos, seqregions1, seqregions2, i1, i2, end_pos);
    EXPECT_EQ(i1, 1);
    EXPECT_EQ(i2, 3);
    EXPECT_EQ(end_pos, 654);
    current_pos = end_pos + 1;
    
    SeqRegions::getNextSharedSegment(current_pos, seqregions1, seqregions2, i1, i2, end_pos);
    EXPECT_EQ(i1, 2);
    EXPECT_EQ(i2, 3);
    EXPECT_EQ(end_pos, 655);
    current_pos = end_pos + 1;
    
    SeqRegions::getNextSharedSegment(current_pos, seqregions1, seqregions2, i1, i2, end_pos);
    EXPECT_EQ(i1, 3);
    EXPECT_EQ(i2, 3);
    EXPECT_EQ(end_pos, 984);
    current_pos = end_pos + 1;
    
    SeqRegions::getNextSharedSegment(current_pos, seqregions1, seqregions2, i1, i2, end_pos);
    EXPECT_EQ(i1, 3);
    EXPECT_EQ(i2, 4);
    EXPECT_EQ(end_pos, 985);
    current_pos = end_pos + 1;
    
    SeqRegions::getNextSharedSegment(current_pos, seqregions1, seqregions2, i1, i2, end_pos);
    EXPECT_EQ(i1, 3);
    EXPECT_EQ(i2, 5);
    EXPECT_EQ(end_pos, 1210);
    current_pos = end_pos + 1;
    
    SeqRegions::getNextSharedSegment(current_pos, seqregions1, seqregions2, i1, i2, end_pos);
    EXPECT_EQ(i1, 3);
    EXPECT_EQ(i2, 6);
    EXPECT_EQ(end_pos, 1211);
    current_pos = end_pos + 1;
    
    SeqRegions::getNextSharedSegment(current_pos, seqregions1, seqregions2, i1, i2, end_pos);
    EXPECT_EQ(i1, 3);
    EXPECT_EQ(i2, 7);
    EXPECT_EQ(end_pos, 1431);
    current_pos = end_pos + 1;
    
    SeqRegions::getNextSharedSegment(current_pos, seqregions1, seqregions2, i1, i2, end_pos);
    EXPECT_EQ(i1, 4);
    EXPECT_EQ(i2, 7);
    EXPECT_EQ(end_pos, 1432);
    current_pos = end_pos + 1;
    
    SeqRegions::getNextSharedSegment(current_pos, seqregions1, seqregions2, i1, i2, end_pos);
    EXPECT_EQ(i1, 5);
    EXPECT_EQ(i2, 7);
    EXPECT_EQ(end_pos, 2431);
    current_pos = end_pos + 1;
    
    SeqRegions::getNextSharedSegment(current_pos, seqregions1, seqregions2, i1, i2, end_pos);
    EXPECT_EQ(i1, 6);
    EXPECT_EQ(i2, 7);
    EXPECT_EQ(end_pos, 2432);
    current_pos = end_pos + 1;
    
    SeqRegions::getNextSharedSegment(current_pos, seqregions1, seqregions2, i1, i2, end_pos);
    EXPECT_EQ(i1, 7);
    EXPECT_EQ(i2, 8);
    EXPECT_EQ(end_pos, 3381);
    current_pos = end_pos + 1;
    
    SeqRegions::getNextSharedSegment(current_pos, seqregions1, seqregions2, i1, i2, end_pos);
    EXPECT_EQ(i1, 8);
    EXPECT_EQ(i2, 8);
    EXPECT_EQ(end_pos, 3500);
    current_pos = end_pos + 1;
}

/*
 Test countSharedSegments(const SeqRegions& seq2_regions, const size_t seq_length) const
 */
TEST(SeqRegions, countSharedSegments)
{
    SeqRegions seqregions1, seqregions2;
    
    // test empty data
    EXPECT_EQ(seqregions1.countSharedSegments(seqregions2, 0), 1);
    
    // generate testing data
    genTestData(seqregions1, seqregions2);
    EXPECT_EQ(seqregions1.countSharedSegments(seqregions2, 100), 3);
    EXPECT_EQ(seqregions1.countSharedSegments(seqregions2, 200), 3);
    EXPECT_EQ(seqregions1.countSharedSegments(seqregions2, 500), 5);
    EXPECT_EQ(seqregions1.countSharedSegments(seqregions2, 1000), 9);
    EXPECT_EQ(seqregions1.countSharedSegments(seqregions2, 1500), 13);
    EXPECT_EQ(seqregions1.countSharedSegments(seqregions2, 2000), 13);
    EXPECT_EQ(seqregions1.countSharedSegments(seqregions2, 2000), 13);
    EXPECT_EQ(seqregions1.countSharedSegments(seqregions2, 3500), 16);
}

/*
 Test compareWithSample(const SeqRegions& sequence2, PositionType seq_length, StateType num_states) const
 */
TEST(SeqRegions, compareWithSample)
{
    SeqRegions seqregions1, seqregions2;
    
    // don't need to test empty data -> we have an assert in compareWithSample to make sure seq_length > 0
    // generate testing data
    genTestData(seqregions1, seqregions2);
    
    // test data
    EXPECT_EQ(seqregions1.compareWithSample(seqregions2, 3500, 4), 0);
    
    SeqRegions seqregions3(&seqregions1);
    EXPECT_EQ(seqregions1.compareWithSample(seqregions3, 3500, 4), 1);
    
    seqregions1.emplace(seqregions1.begin(), TYPE_N, 54);
    EXPECT_EQ(seqregions1.compareWithSample(seqregions3, 3500, 4), -1);
    
    seqregions1.emplace(seqregions1.begin() + 4, TYPE_N, 999);
    seqregions1.emplace(seqregions1.begin() + 5, 1, 1000);
    EXPECT_EQ(seqregions1.compareWithSample(seqregions3, 3500, 4), 0);
    
    seqregions3.emplace(seqregions3.begin(), TYPE_N, 54);
    EXPECT_EQ(seqregions1.compareWithSample(seqregions3, 3500, 4), 1);
    
    seqregions3.emplace(seqregions3.begin() + 4, TYPE_N, 999);
    seqregions3.emplace(seqregions3.begin() + 5, 1, 1000);
    EXPECT_EQ(seqregions1.compareWithSample(seqregions3, 3500, 4), 1);
    
    auto new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value1{0.5,0,0,0.5};
    (*new_lh) = new_lh_value1;
    seqregions1.emplace(seqregions1.begin() + seqregions1.size() - 1, TYPE_O, 3382, 1e-3, -1, std::move(new_lh));
    EXPECT_EQ(seqregions1.compareWithSample(seqregions3, 3500, 4), -1);
    
    new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value2{0.5,0,0,0.5};
    (*new_lh) = new_lh_value2;
    seqregions3.emplace(seqregions3.begin() + seqregions3.size() - 1, TYPE_O, 3382, 1e-3, -1, std::move(new_lh));
    EXPECT_EQ(seqregions1.compareWithSample(seqregions3, 3500, 4), 1);
    
    (*seqregions3[seqregions3.size() - 2].likelihood)[2] = 0.5;
    (*seqregions3[seqregions3.size() - 2].likelihood)[3] = 0;
    EXPECT_EQ(seqregions1.compareWithSample(seqregions3, 3500, 4), 0);
    
    (*seqregions1[seqregions1.size() - 2].likelihood)[0] = 1.0 / 3;
    (*seqregions1[seqregions1.size() - 2].likelihood)[1] = 1.0 / 3;
    (*seqregions1[seqregions1.size() - 2].likelihood)[2] = 1.0 / 3;
    (*seqregions1[seqregions1.size() - 2].likelihood)[3] = 0;
    EXPECT_EQ(seqregions1.compareWithSample(seqregions3, 3500, 4), -1);
}
