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
void genTestData1(SeqRegions& seqregions1, SeqRegions& seqregions2)
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
    genTestData1(seqregions1, seqregions2);
    
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
    genTestData1(seqregions1, seqregions2);
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
    genTestData1(seqregions1, seqregions2);
    
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

/*
 Test areDiffFrom(const SeqRegions& regions2, PositionType seq_length, StateType num_states, const Params* params) const
 */
TEST(SeqRegions, areDiffFrom)
{
    SeqRegions seqregions1, seqregions2;
    
    // init testing data
    const PositionType seq_length = 3500;
    Params params = Params::getInstance();
    initDefaultValue(params);
    seqregions1.emplace_back(TYPE_R, 382, -1, 0.321);
    seqregions1.emplace_back(TYPE_N, 654);
    seqregions1.emplace_back(0, 655, 0);
    seqregions1.emplace_back(TYPE_N, 1431);
    seqregions1.emplace_back(3, 1432, 0.432, 0);
    seqregions1.emplace_back(TYPE_N, 2431);
    auto new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value{0.1,0.3,0.2,0.4};
    (*new_lh) = new_lh_value;
    seqregions1.emplace_back(TYPE_O, 2432, 0, -1, std::move(new_lh));
    seqregions1.emplace_back(TYPE_N, 3381);
    seqregions1.emplace_back(TYPE_R, seq_length, 0, 0.1321);
    
    // ---- seqregions2 is empty -----
    EXPECT_EQ(seqregions1.areDiffFrom(seqregions2, seq_length, 4, &params), true);
    // ---- seqregions2 is empty -----
    
    // ---- other tests -----
    seqregions2.emplace_back(TYPE_R, seq_length, -1, 0.321);
    EXPECT_EQ(seqregions1.areDiffFrom(seqregions2, seq_length, 4, &params), true);
    
    seqregions2.back().position = 382;
    seqregions2.emplace_back(TYPE_N, 654);
    seqregions2.emplace_back(0, 655, 0);
    seqregions2.emplace_back(TYPE_N, 1431);
    seqregions2.emplace_back(3, 1432, 0.432, 0);
    seqregions2.emplace_back(TYPE_N, 2431);
    new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value1{0.1,0.3,0.2,0.4};
    (*new_lh) = new_lh_value1;
    seqregions2.emplace_back(TYPE_O, 2432, 0, -1, std::move(new_lh));
    seqregions2.emplace_back(TYPE_N, 3381);
    seqregions2.emplace_back(TYPE_R, seq_length, 0, 0.1321);
    EXPECT_EQ(seqregions1.areDiffFrom(seqregions2, seq_length, 4, &params), false);
    
    seqregions2.data()[3].type = TYPE_R;
    EXPECT_EQ(seqregions1.areDiffFrom(seqregions2, seq_length, 4, &params), true);
    
    seqregions2.data()[3].type = TYPE_N;
    seqregions2.data()[4].plength_observation2root = -1;
    EXPECT_EQ(seqregions1.areDiffFrom(seqregions2, seq_length, 4, &params), true);
    
    seqregions2.data()[4].plength_observation2root = 1e-9;
    EXPECT_EQ(seqregions1.areDiffFrom(seqregions2, seq_length, 4, &params), false);
    
    seqregions1.back().plength_observation2node = 0.0113;
    EXPECT_EQ(seqregions1.areDiffFrom(seqregions2, seq_length, 4, &params), true);
    
    seqregions2.back().plength_observation2node = 0.0113 + 1e-10;
    EXPECT_EQ(seqregions1.areDiffFrom(seqregions2, seq_length, 4, &params), false);
    
    // test on type O
    RealNumType thresh_diff_update = params.thresh_diff_update / 2; // difference is too small less then thresh_diff_update
    seqregions2.data()[6].likelihood->data()[0] += thresh_diff_update;
    seqregions2.data()[6].likelihood->data()[2] -= thresh_diff_update;
    EXPECT_EQ(seqregions1.areDiffFrom(seqregions2, seq_length, 4, &params), false);
    
    seqregions2.data()[6].likelihood->data()[0] = thresh_diff_update; // reset lh so that all pairs of lh between seqregions1 and seqregions2 equal to each other
    seqregions2.data()[6].likelihood->data()[2] += seqregions2.data()[6].likelihood->data()[0];
    seqregions1.data()[6].likelihood->data()[0] = seqregions2.data()[6].likelihood->data()[0];
    seqregions1.data()[6].likelihood->data()[2] = seqregions2.data()[6].likelihood->data()[2];
    EXPECT_EQ(seqregions1.areDiffFrom(seqregions2, seq_length, 4, &params), false);
    
    seqregions2.data()[6].likelihood->data()[2] += seqregions2.data()[6].likelihood->data()[0];
    seqregions2.data()[6].likelihood->data()[0] = 0; // one lh = 0 but diff != 0
    EXPECT_EQ(seqregions1.areDiffFrom(seqregions2, seq_length, 4, &params), true);
    
    seqregions2.data()[6].likelihood->data()[2] -= thresh_diff_update / 10;
    seqregions2.data()[6].likelihood->data()[0] = thresh_diff_update / 10;
    EXPECT_EQ(seqregions1.areDiffFrom(seqregions2, seq_length, 4, &params), true);
    // ---- other tests -----
}

/*
 Test simplifyO(RealNumType* const partial_lh, StateType ref_state, StateType num_states, RealNumType threshold) const
 */
TEST(SeqRegions, simplifyO)
{
    // init testing data
    SeqRegions seqregions1;
    Params params = Params::getInstance();
    initDefaultValue(params);
    RealNumType threshold_prob = params.threshold_prob;
    
    // tests
    auto new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value{0.1,0.3,0.2,0.4};
    (*new_lh) = new_lh_value;
    EXPECT_EQ(seqregions1.simplifyO(new_lh->data(), 2, 4, threshold_prob), TYPE_O) ;
    
    SeqRegion::LHType new_lh_value1{1.0 - 3 * threshold_prob, threshold_prob, threshold_prob, threshold_prob};
    (*new_lh) = new_lh_value1;
    EXPECT_EQ(seqregions1.simplifyO(new_lh->data(), 2, 4, threshold_prob), 0) ;
    
    new_lh->data()[1] += new_lh->data()[2];
    new_lh->data()[2] = 0;
    EXPECT_EQ(seqregions1.simplifyO(new_lh->data(), 2, 4, threshold_prob), TYPE_O) ;
    
    new_lh->data()[2] = new_lh->data()[0];
    new_lh->data()[0] = new_lh->data()[3];
    new_lh->data()[1] = new_lh->data()[3];
    EXPECT_EQ(seqregions1.simplifyO(new_lh->data(), 2, 4, threshold_prob), TYPE_R) ;
}

/*
    Generate testing data (seqregions1, seqregions2, seqregions3)
 */
void genTestData2(SeqRegions& seqregions1, SeqRegions& seqregions2, SeqRegions& seqregions3)
{
    // seqregions1
    seqregions1.emplace_back(250,239,-1,-1);
    auto new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value1{2.6766774402933635499746492514283602304203668609262e-05,0.39996252651583585890904259940725751221179962158203,2.6766774402933635499746492514283602304203668609262e-05,0.59998393993535814594508792652050033211708068847656};
    (*new_lh) = new_lh_value1;
    seqregions1.emplace_back(251,240,0,0,std::move(new_lh));
    seqregions1.emplace_back(250,3035,-1,-1);
    seqregions1.emplace_back(1,3036,-1,-1);
    seqregions1.emplace_back(250,8780,-1,-1);
    new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value2{2.6766774402933635499746492514283602304203668609262e-05,0.39996252651583585890904259940725751221179962158203,2.6766774402933635499746492514283602304203668609262e-05,0.59998393993535814594508792652050033211708068847656};
    (*new_lh) = new_lh_value2;
    seqregions1.emplace_back(251,8781,0,0,std::move(new_lh));
    seqregions1.emplace_back(250,14406,-1,-1);
    seqregions1.emplace_back(1,14407,-1,-1);
    seqregions1.emplace_back(250,18754,-1,-1);
    new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value3{2.6766774402933635499746492514283602304203668609262e-05,2.6766774402933635499746492514283602304203668609262e-05,0.39996252651583585890904259940725751221179962158203,0.59998393993535814594508792652050033211708068847656};
    (*new_lh) = new_lh_value3;
    seqregions1.emplace_back(251,18755,0,0,std::move(new_lh));
    seqregions1.emplace_back(250,22466,-1,-1);
    new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value4{2.6766774402933635499746492514283602304203668609262e-05,2.6766774402933635499746492514283602304203668609262e-05,0.39996252651583585890904259940725751221179962158203,0.59998393993535814594508792652050033211708068847656};
    (*new_lh) = new_lh_value4;
    seqregions1.emplace_back(251,22467,0,0,std::move(new_lh));
    seqregions1.emplace_back(250,23401,-1,-1);
    new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value5{0.39996252651583585890904259940725751221179962158203,2.6766774402933635499746492514283602304203668609262e-05,0.59998393993535814594508792652050033211708068847656,2.6766774402933635499746492514283602304203668609262e-05};
    (*new_lh) = new_lh_value5;
    seqregions1.emplace_back(251,23402,0,0,std::move(new_lh));
    seqregions1.emplace_back(250,28142,-1,-1);
    new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value6{2.6766774402933635499746492514283602304203668609262e-05,0.39996252651583585890904259940725751221179962158203,2.6766774402933635499746492514283602304203668609262e-05,0.59998393993535814594508792652050033211708068847656};
    (*new_lh) = new_lh_value6;
    seqregions1.emplace_back(251,28143,0,0,std::move(new_lh));
    seqregions1.emplace_back(250,28876,-1,-1);
    new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value7{0.59998393993535814594508792652050033211708068847656,2.6766774402933635499746492514283602304203668609262e-05,0.39996252651583585890904259940725751221179962158203,2.6766774402933635499746492514283602304203668609262e-05};
    (*new_lh) = new_lh_value7;
    seqregions1.emplace_back(251,28877,0,0,std::move(new_lh));
    seqregions1.emplace_back(250,29708,-1,-1);
    new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value8{0.59998393993535814594508792652050033211708068847656,2.6766774402933635499746492514283602304203668609262e-05,2.6766774402933635499746492514283602304203668609262e-05,0.39996252651583585890904259940725751221179962158203};
    (*new_lh) = new_lh_value8;
    seqregions1.emplace_back(251,29709,0,0,std::move(new_lh));
    seqregions1.emplace_back(250,29740,-1,-1);
    new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value9{0.59998393993535814594508792652050033211708068847656,2.6766774402933635499746492514283602304203668609262e-05,0.39996252651583585890904259940725751221179962158203,2.6766774402933635499746492514283602304203668609262e-05};
    (*new_lh) = new_lh_value9;
    seqregions1.emplace_back(251,29741,0,0,std::move(new_lh));
    seqregions1.emplace_back(250,29890,-1,-1);
    
    // seqregions2
    seqregions2.emplace_back(250,239,-1,-1);
    seqregions2.emplace_back(1,240,-1,-1);
    seqregions2.emplace_back(250,1267,-1,-1);
    new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value10{0.80000178432426516383912940000300295650959014892578,8.9216213262436184458807272856795123061601771041751e-06,0.19998037243308225408000566858390811830759048461914,8.9216213262436184458807272856795123061601771041751e-06};
    (*new_lh) = new_lh_value10;
    seqregions2.emplace_back(251,1268,0,0,std::move(new_lh));
    seqregions2.emplace_back(250,3035,-1,-1);
    seqregions2.emplace_back(1,3036,-1,-1);
    seqregions2.emplace_back(250,8780,-1,-1);
    seqregions2.emplace_back(3,8781,-1,-1);
    seqregions2.emplace_back(250,14371,-1,-1);
    new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value11{0.80000178432426516383912940000300295650959014892578,0.19998037243308225408000566858390811830759048461914,8.9216213262436184458807272856795123061601771041751e-06,8.9216213262436184458807272856795123061601771041751e-06};
    (*new_lh) = new_lh_value11;
    seqregions2.emplace_back(251,14372,0,0,std::move(new_lh));
    seqregions2.emplace_back(250,14406,-1,-1);
    seqregions2.emplace_back(1,14407,-1,-1);
    seqregions2.emplace_back(250,21362,-1,-1);
    new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value12{8.9216213262436184458807272856795123061601771041751e-06,0.80000178432426516383912940000300295650959014892578,8.9216213262436184458807272856795123061601771041751e-06,0.19998037243308225408000566858390811830759048461914};
    (*new_lh) = new_lh_value12;
    seqregions2.emplace_back(251,21363,0,0,std::move(new_lh));
    seqregions2.emplace_back(250,23401,-1,-1);
    seqregions2.emplace_back(0,23402,-1,-1);
    seqregions2.emplace_back(250,28142,-1,-1);
    seqregions2.emplace_back(1,28143,-1,-1);
    seqregions2.emplace_back(250,28165,-1,-1);
    new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value13{0.19998037243308225408000566858390811830759048461914,8.9216213262436184458807272856795123061601771041751e-06,0.80000178432426516383912940000300295650959014892578,8.9216213262436184458807272856795123061601771041751e-06};
    (*new_lh) = new_lh_value13;
    seqregions2.emplace_back(251,28166,0,0,std::move(new_lh));
    seqregions2.emplace_back(250,28876,-1,-1);
    seqregions2.emplace_back(0,28877,-1,-1);
    seqregions2.emplace_back(250,29740,-1,-1);
    seqregions2.emplace_back(0,29741,-1,-1);
    seqregions2.emplace_back(250,29890,-1,-1);
    
    // seqregions3
    seqregions3.emplace_back(250,53,3.3454886086112878360986772063867533688608091324568e-05,-1);
    seqregions3.emplace_back(250,239,-1,-1);
    seqregions3.emplace_back(1,240,-1,-1);
    seqregions3.emplace_back(250,3035,-1,-1);
    seqregions3.emplace_back(1,3036,-1,-1);
    seqregions3.emplace_back(250,8780,-1,-1);
    seqregions3.emplace_back(3,8781,-1,-1);
    seqregions3.emplace_back(250,14406,-1,-1);
    seqregions3.emplace_back(1,14407,-1,-1);
    seqregions3.emplace_back(250,17745,-1,-1);
    new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value14{1.2436575683940661210420103588707597258578019250308e-10,1.1151877415789591228433182135137968771232408471406e-05,1.2436575683940661210420103588707597258578019250308e-10,0.99998884787385255989988763758447021245956420898438};
    (*new_lh) = new_lh_value14;
    seqregions3.emplace_back(251,17746,0,0,std::move(new_lh));
    seqregions3.emplace_back(250,17856,-1,-1);
    seqregions3.emplace_back(2,17857,-1,-1);
    seqregions3.emplace_back(250,18058,-1,-1);
    seqregions3.emplace_back(3,18059,-1,-1);
    seqregions3.emplace_back(250,21135,-1,-1);
    new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value15{0.99998884787385255989988763758447021245956420898438,1.2436575683940661210420103588707597258578019250308e-10,1.1151877415789591228433182135137968771232408471406e-05,1.2436575683940661210420103588707597258578019250308e-10};
    (*new_lh) = new_lh_value15;
    seqregions3.emplace_back(251,21136,0,0,std::move(new_lh));
    seqregions3.emplace_back(250,23401,-1,-1);
    seqregions3.emplace_back(0,23402,-1,-1);
    seqregions3.emplace_back(250,28142,-1,-1);
    seqregions3.emplace_back(1,28143,-1,-1);
    seqregions3.emplace_back(250,29835,-1,-1);
    seqregions3.emplace_back(250,29881,3.3454886086112878360986772063867533688608091324568e-05,-1);
    seqregions3.emplace_back(250,29890,6.6909772172225756721973544127735067377216182649136e-05,-1);
}

/*
    Generate output data (seqregions1, seqregions2, seqregions3)
 */
void genOutputData2(SeqRegions& seqregions2_total_lh, SeqRegions& seqregions3_total_lh, SeqRegions& seqregions4_total_lh)
{
    // seqregions2_total_lh
    seqregions2_total_lh.emplace_back(250,239,-1,-1);
    auto new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value{3.0088842169314083832640932536506284122879151254892e-05,0.27601702884096879220265918775112368166446685791016,1.973391907849880502474719523586799141412484459579e-05,0.72393314839778344360610162766533903777599334716797};
    (*new_lh) = new_lh_value;
    seqregions2_total_lh.emplace_back(251,240,0,0,std::move(new_lh));
    seqregions2_total_lh.emplace_back(250,3035,-1,-1);
    seqregions2_total_lh.emplace_back(1,3036,-1,-1);
    seqregions2_total_lh.emplace_back(250,8780,-1,-1);
    new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value1{3.0088842169314083832640932536506284122879151254892e-05,0.27601702884096879220265918775112368166446685791016,1.973391907849880502474719523586799141412484459579e-05,0.72393314839778344360610162766533903777599334716797};
    (*new_lh) = new_lh_value1;
    seqregions2_total_lh.emplace_back(251,8781,0,0,std::move(new_lh));
    seqregions2_total_lh.emplace_back(250,14406,-1,-1);
    seqregions2_total_lh.emplace_back(1,14407,-1,-1);
    seqregions2_total_lh.emplace_back(250,18754,-1,-1);
    new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value2{2.9531992963613746839288479173468715544004226103425e-05,1.8130087168915766599255889879316328006098046898842e-05,0.28941690008858006466496703978918958455324172973633,0.71053543783128736421872417849954217672348022460938};
    (*new_lh) = new_lh_value2;
    seqregions2_total_lh.emplace_back(251,18755,0,0,std::move(new_lh));
    seqregions2_total_lh.emplace_back(250,22466,-1,-1);
    new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value3{2.9531992963613746839288479173468715544004226103425e-05,1.8130087168915766599255889879316328006098046898842e-05,0.28941690008858006466496703978918958455324172973633,0.71053543783128736421872417849954217672348022460938};
    (*new_lh) = new_lh_value3;
    seqregions2_total_lh.emplace_back(251,22467,0,0,std::move(new_lh));
    seqregions2_total_lh.emplace_back(250,23401,-1,-1);
    new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value4{0.50404208207157885635041338900919072329998016357422,2.0708596602308600767008636700516888140555238351226e-05,0.49590100229928235631149391338112764060497283935547,3.6207032536410212360256100083688579616136848926544e-05};
    (*new_lh) = new_lh_value4;
    seqregions2_total_lh.emplace_back(251,23402,0,0,std::move(new_lh));
    seqregions2_total_lh.emplace_back(250,28142,-1,-1);
    new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value5{3.0088842169314083832640932536506284122879151254892e-05,0.27601702884096879220265918775112368166446685791016,1.973391907849880502474719523586799141412484459579e-05,0.72393314839778344360610162766533903777599334716797};
    (*new_lh) = new_lh_value5;
    seqregions2_total_lh.emplace_back(251,28143,0,0,std::move(new_lh));
    seqregions2_total_lh.emplace_back(250,28876,-1,-1);
    new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value6{0.69575734124132038171950398464105091989040374755859,1.9055543772650076496909254952782930558896623551846e-05,0.30419028638969553002269208263896871358156204223633,3.3316825211536308379896981213263984500372316688299e-05};
    (*new_lh) = new_lh_value6;
    seqregions2_total_lh.emplace_back(251,28877,0,0,std::move(new_lh));
    seqregions2_total_lh.emplace_back(250,29708,-1,-1);
    new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value7{0.58289030162667865297976277361158281564712524414062,1.5964318303107729663447775236839731860527535900474e-05,1.7054975866173025634371304692926685220299987122416e-05,0.41707667907915202398783094395184889435768127441406};
    (*new_lh) = new_lh_value7;
    seqregions2_total_lh.emplace_back(251,29709,0,0,std::move(new_lh));
    seqregions2_total_lh.emplace_back(250,29740,-1,-1);
    new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value8{0.69575734124132038171950398464105091989040374755859,1.9055543772650076496909254952782930558896623551846e-05,0.30419028638969553002269208263896871358156204223633,3.3316825211536308379896981213263984500372316688299e-05};
    (*new_lh) = new_lh_value8;
    seqregions2_total_lh.emplace_back(251,29741,0,0,std::move(new_lh));
    seqregions2_total_lh.emplace_back(250,29890,-1,-1);
    
    // seqregions3_total_lh
    seqregions3_total_lh.emplace_back(250,239,1.0000000000000000818030539140313095458623138256371e-05,0);
    seqregions3_total_lh.emplace_back(1,240,1.0000000000000000818030539140313095458623138256371e-05,0);
    seqregions3_total_lh.emplace_back(250,1267,1.0000000000000000818030539140313095458623138256371e-05,0);
    new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value9{0.85912905483272095086277886366588063538074493408203,6.5230177179574731834046756595579807935791905038059e-06,0.14085255687431794124897521669481648132205009460449,1.1865275243000726553043236433104823390749515965581e-05};
    (*new_lh) = new_lh_value9;
    seqregions3_total_lh.emplace_back(251,1268,0,0,std::move(new_lh));
    seqregions3_total_lh.emplace_back(250,3035,1.0000000000000000818030539140313095458623138256371e-05,0);
    seqregions3_total_lh.emplace_back(1,3036,1.0000000000000000818030539140313095458623138256371e-05,0);
    seqregions3_total_lh.emplace_back(250,8780,1.0000000000000000818030539140313095458623138256371e-05,0);
    seqregions3_total_lh.emplace_back(3,8781,1.0000000000000000818030539140313095458623138256371e-05,0);
    seqregions3_total_lh.emplace_back(250,14371,1.0000000000000000818030539140313095458623138256371e-05,0);
    new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value10{0.86693729064253011173946106282528489828109741210938,0.13304159211598842738055736845126375555992126464844,8.1954130681283745082791117320120122258231276646256e-06,1.292182841329131426106463509384525423229206353426e-05};
    (*new_lh) = new_lh_value10;
    seqregions3_total_lh.emplace_back(251,14372,0,0,std::move(new_lh));
    seqregions3_total_lh.emplace_back(250,14406,1.0000000000000000818030539140313095458623138256371e-05,0);
    seqregions3_total_lh.emplace_back(1,14407,1.0000000000000000818030539140313095458623138256371e-05,0);
    seqregions3_total_lh.emplace_back(250,21362,1.0000000000000000818030539140313095458623138256371e-05,0);
    new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value11{1.3542631186466222171496318060768260238546645268798e-05,0.69584313629962879499402106375782750546932220458984,1.0150959877123927614823994947101937214029021561146e-05,0.30413317010930768224952203127031680196523666381836};
    (*new_lh) = new_lh_value11;
    seqregions3_total_lh.emplace_back(251,21363,0,0,std::move(new_lh));
    seqregions3_total_lh.emplace_back(250,23401,1.0000000000000000818030539140313095458623138256371e-05,0);
    seqregions3_total_lh.emplace_back(0,23402,1.0000000000000000818030539140313095458623138256371e-05,0);
    seqregions3_total_lh.emplace_back(250,28142,1.0000000000000000818030539140313095458623138256371e-05,0);
    seqregions3_total_lh.emplace_back(1,28143,1.0000000000000000818030539140313095458623138256371e-05,0);
    seqregions3_total_lh.emplace_back(250,28165,1.0000000000000000818030539140313095458623138256371e-05,0);
    new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value12{0.27595781923293211113090706021466758102178573608398,8.3817209383128356075479820086471249851456377655268e-06,0.72401575181023847260775028189527802169322967529297,1.8047235891267821277644811672757896303664892911911e-05};
    (*new_lh) = new_lh_value12;
    seqregions3_total_lh.emplace_back(251,28166,0,0,std::move(new_lh));
    seqregions3_total_lh.emplace_back(250,28876,1.0000000000000000818030539140313095458623138256371e-05,0);
    seqregions3_total_lh.emplace_back(0,28877,1.0000000000000000818030539140313095458623138256371e-05,0);
    seqregions3_total_lh.emplace_back(250,29740,1.0000000000000000818030539140313095458623138256371e-05,0);
    seqregions3_total_lh.emplace_back(0,29741,1.0000000000000000818030539140313095458623138256371e-05,0);
    seqregions3_total_lh.emplace_back(250,29890,1.0000000000000000818030539140313095458623138256371e-05,0);
    
    // seqregions4_total_lh
    seqregions4_total_lh.emplace_back(250,53,0.037653454886086110131593329697352601215243339538574,0);
    seqregions4_total_lh.emplace_back(250,239,0.037620000000000000661692922676593298092484474182129,0);
    seqregions4_total_lh.emplace_back(1,240,0.037620000000000000661692922676593298092484474182129,0);
    seqregions4_total_lh.emplace_back(250,3035,0.037620000000000000661692922676593298092484474182129,0);
    seqregions4_total_lh.emplace_back(1,3036,0.037620000000000000661692922676593298092484474182129,0);
    seqregions4_total_lh.emplace_back(250,8780,0.037620000000000000661692922676593298092484474182129,0);
    seqregions4_total_lh.emplace_back(3,8781,0.037620000000000000661692922676593298092484474182129,0);
    seqregions4_total_lh.emplace_back(250,14406,0.037620000000000000661692922676593298092484474182129,0);
    seqregions4_total_lh.emplace_back(1,14407,0.037620000000000000661692922676593298092484474182129,0);
    seqregions4_total_lh.emplace_back(250,17745,0.037620000000000000661692922676593298092484474182129,0);
    new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value13{0.0027898145230272866808285403550371484016068279743195,0.029996251268427902986202226998102560173720121383667,0.014646426545977904096207389272876753238961100578308,0.95256750766256692575240094811306335031986236572266};
    (*new_lh) = new_lh_value13;
    seqregions4_total_lh.emplace_back(251,17746,0,0,std::move(new_lh));
    seqregions4_total_lh.emplace_back(250,17856,0.037620000000000000661692922676593298092484474182129,0);
    seqregions4_total_lh.emplace_back(2,17857,0.037620000000000000661692922676593298092484474182129,0);
    seqregions4_total_lh.emplace_back(250,18058,0.037620000000000000661692922676593298092484474182129,0);
    seqregions4_total_lh.emplace_back(3,18059,0.037620000000000000661692922676593298092484474182129,0);
    seqregions4_total_lh.emplace_back(250,21135,0.037620000000000000661692922676593298092484474182129,0);
    new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value14{0.98726625100220433139952547207940369844436645507812,0.0022458812800641840025372975020445664995349943637848,0.0074932173880506384305855149818853533361107110977173,0.002994650329681008363996719268129709234926849603653};
    (*new_lh) = new_lh_value14;
    seqregions4_total_lh.emplace_back(251,21136,0,0,std::move(new_lh));
    seqregions4_total_lh.emplace_back(250,23401,0.037620000000000000661692922676593298092484474182129,0);
    seqregions4_total_lh.emplace_back(0,23402,0.037620000000000000661692922676593298092484474182129,0);
    seqregions4_total_lh.emplace_back(250,28142,0.037620000000000000661692922676593298092484474182129,0);
    seqregions4_total_lh.emplace_back(1,28143,0.037620000000000000661692922676593298092484474182129,0);
    seqregions4_total_lh.emplace_back(250,29835,0.037620000000000000661692922676593298092484474182129,0);
    seqregions4_total_lh.emplace_back(250,29881,0.037653454886086110131593329697352601215243339538574,0);
    seqregions4_total_lh.emplace_back(250,29890,0.037686909772172226540387640625340281985700130462646,0);
}

/*
 Test computeAbsoluteLhAtRoot(const Alignment& aln, const Model& model, std::vector< std::vector<PositionType> > &cumulative_base)
 */
TEST(SeqRegions, computeAbsoluteLhAtRoot)
{
    Alignment aln;
    Model model;
    Params params = Params::getInstance();
    RealNumType *cumulative_rate = nullptr;
    std::vector< std::vector<PositionType> > cumulative_base;
    
    // Init params, aln, and model
    initDefaultValue(params);
    params.model_name = "GTR";
    std::string diff_file_path("../../example/test_5K.diff");
    char* diff_file_path_ptr = new char[diff_file_path.length() + 1];
    strcpy(diff_file_path_ptr, diff_file_path.c_str());
    aln.readDiff(diff_file_path_ptr, NULL);
    // extract related info (freqs, log_freqs) of the ref sequence
    model.extractRefInfo(aln.ref_seq, aln.num_states);
    // init the mutation matrix from a model name
    model.initMutationMat(params.model_name, aln.num_states);
    // compute cumulative rates of the ref sequence
    model.computeCumulativeRate(cumulative_rate, cumulative_base, aln);
    // dummy variables
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    
    // ----- Simple tests -----
    SeqRegions seqregions1;
    seqregions1.emplace_back(TYPE_R, seq_length - 1);
    EXPECT_EQ(seqregions1.computeAbsoluteLhAtRoot(num_states, model, cumulative_base), -40547.865582541948);
    
    seqregions1.emplace(seqregions1.begin(), 3, 3242);
    seqregions1.emplace(seqregions1.begin(), TYPE_R, 3241);
    EXPECT_EQ(seqregions1.computeAbsoluteLhAtRoot(num_states, model, cumulative_base), -40547.372963956892);
    
    auto new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value{0.1,0.3,0.2,0.4};
    (*new_lh) = new_lh_value;
    seqregions1.emplace(seqregions1.begin(), TYPE_O, 243, 0, -1, std::move(new_lh));
    seqregions1.emplace(seqregions1.begin(), TYPE_R, 242);
    EXPECT_EQ(seqregions1.computeAbsoluteLhAtRoot(num_states, model, cumulative_base), -40547.053844652924);
    // ----- Simple tests -----
    
    // Generate complex seqregions
    SeqRegions seqregions2, seqregions3, seqregions4;
    genTestData2(seqregions2, seqregions3, seqregions4);
    
    // ----- Test 1 on a more complex seqregions -----
    EXPECT_EQ(seqregions2.computeAbsoluteLhAtRoot(num_states, model, cumulative_base), -40547.644608026116);
    // ----- Test 1 on a more complex seqregions -----
    
    // ----- Test 2 on a more complex seqregions -----
    EXPECT_EQ(seqregions3.computeAbsoluteLhAtRoot(num_states, model, cumulative_base), -40548.188644939095);
    // ----- Test 2 on a more complex seqregions -----
    
    // ----- Test 3 on a more complex seqregions -----
    EXPECT_EQ(seqregions4.computeAbsoluteLhAtRoot(num_states, model, cumulative_base), -40548.424295613549);
    // ----- Test 3 on a more complex seqregions -----
}

/*
 Test computeTotalLhAtRoot(StateType num_states, const Model& model, RealNumType blength = -1)
 */
TEST(SeqRegions, computeTotalLhAtRoot)
{
    Alignment aln;
    Model model;
    Params params = Params::getInstance();
    RealNumType *cumulative_rate = nullptr;
    std::vector< std::vector<PositionType> > cumulative_base;
    
    // Init params, aln, and model
    initDefaultValue(params);
    params.model_name = "GTR";
    std::string diff_file_path("../../example/test_5K.diff");
    char* diff_file_path_ptr = new char[diff_file_path.length() + 1];
    strcpy(diff_file_path_ptr, diff_file_path.c_str());
    aln.readDiff(diff_file_path_ptr, NULL);
    // extract related info (freqs, log_freqs) of the ref sequence
    model.extractRefInfo(aln.ref_seq, aln.num_states);
    // init the mutation matrix from a model name
    model.initMutationMat(params.model_name, aln.num_states);
    // compute cumulative rates of the ref sequence
    model.computeCumulativeRate(cumulative_rate, cumulative_base, aln);
    // dummy variables
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    
    // ----- Simple tests -----
    SeqRegions seqregions1;
    seqregions1.emplace_back(TYPE_R, seq_length - 1);
    SeqRegions* seqregions_total_lh = seqregions1.computeTotalLhAtRoot(num_states, model);
    EXPECT_EQ(seqregions_total_lh->size(), 1);
    EXPECT_EQ(seqregions_total_lh->front().plength_observation2root, -1);
    EXPECT_EQ(seqregions_total_lh->front().plength_observation2node, -1);
    
    seqregions1.emplace(seqregions1.begin(), 3, 3242);
    seqregions1.emplace(seqregions1.begin(), TYPE_R, 3241);
    delete seqregions_total_lh;
    seqregions_total_lh = seqregions1.computeTotalLhAtRoot(num_states, model);
    EXPECT_EQ(seqregions_total_lh->size(), 3);
    EXPECT_EQ(seqregions_total_lh->front().type, seqregions_total_lh->back().type); // TYPE_R
    EXPECT_EQ(seqregions_total_lh->data()[1].type, 3);
    EXPECT_EQ(seqregions_total_lh->data()[1].position, 3242);
    
    auto new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value{0.1,0.3,0.2,0.4};
    (*new_lh) = new_lh_value;
    seqregions1.emplace(seqregions1.begin(), TYPE_O, 243, 0, -1, std::move(new_lh));
    seqregions1.emplace(seqregions1.begin(), TYPE_N, 242);
    delete seqregions_total_lh;
    seqregions_total_lh = seqregions1.computeTotalLhAtRoot(num_states, model);
    EXPECT_EQ(seqregions_total_lh->size(), 5);
    EXPECT_EQ(seqregions_total_lh->front().type, TYPE_N);
    EXPECT_EQ(seqregions_total_lh->front().plength_observation2root, -1);
    EXPECT_TRUE(seqregions_total_lh->data()[1].likelihood != nullptr);
    EXPECT_EQ(seqregions_total_lh->data()[1].plength_observation2root, -1);
    EXPECT_EQ(seqregions_total_lh->data()[3].position, 3242);
    EXPECT_EQ(seqregions_total_lh->data()[2].type, seqregions_total_lh->back().type); // TYPE_R
    // ----- Simple tests -----
    
    // Generate complex seqregions
    SeqRegions seqregions2, seqregions3, seqregions4;
    genTestData2(seqregions2, seqregions3, seqregions4);
    SeqRegions seqregions2_total_lh, seqregions3_total_lh, seqregions4_total_lh;
    genOutputData2(seqregions2_total_lh, seqregions3_total_lh, seqregions4_total_lh);
    
    // ----- Test 1 on a more complex seqregions -----
    delete seqregions_total_lh;
    seqregions_total_lh = seqregions2.computeTotalLhAtRoot(num_states, model);
    EXPECT_EQ(*seqregions_total_lh, seqregions2_total_lh);
    // ----- Test 1 on a more complex seqregions -----
    
    // ----- Test 2 on a more complex seqregions -----
    delete seqregions_total_lh;
    seqregions_total_lh = seqregions3.computeTotalLhAtRoot(num_states, model, 1e-5);
    EXPECT_EQ(*seqregions_total_lh, seqregions3_total_lh);
    // ----- Test 2 on a more complex seqregions -----
    
    // ----- Test 3 on a more complex seqregions -----
    delete seqregions_total_lh;
    seqregions_total_lh = seqregions4.computeTotalLhAtRoot(num_states, model, 0.03762);
    EXPECT_EQ(*seqregions_total_lh, seqregions4_total_lh);
    // ----- Test 3 on a more complex seqregions -----
    
}
