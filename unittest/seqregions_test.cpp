#include "gtest/gtest.h"
#include "alignment/seqregions.h"

/*
 Test addNonConsecutiveRRegion()
 */
TEST(SeqRegions, addNonConsecutiveRRegion)
{
    RealNumType threshold_prob = 1e-8;
    SeqRegions seqregions;
    
    SeqRegions::addNonConsecutiveRRegion(seqregions, 0, -1, -1, 24, threshold_prob);
    EXPECT_EQ(seqregions.size(), 1);
    EXPECT_EQ(seqregions.back().position, 24);
    EXPECT_EQ(seqregions.back().plength_observation2node, -1);
    EXPECT_EQ(seqregions.back().plength_observation2root, -1);
    
    SeqRegions::addNonConsecutiveRRegion(seqregions, TYPE_R, -1, -1, 113, threshold_prob);
    EXPECT_EQ(seqregions.size(), 2);
    EXPECT_EQ(seqregions.back().position, 113);
    EXPECT_EQ(seqregions.back().plength_observation2node, -1);
    EXPECT_EQ(seqregions.back().plength_observation2root, -1);
    
    SeqRegions::addNonConsecutiveRRegion(seqregions, TYPE_R, -1, -1, 159, threshold_prob);
    EXPECT_EQ(seqregions.size(), 2); // merge two consecotive R
    EXPECT_EQ(seqregions.back().position, 159);
    EXPECT_EQ(seqregions.back().plength_observation2node, -1);
    EXPECT_EQ(seqregions.back().plength_observation2root, -1);
    
    SeqRegions::addNonConsecutiveRRegion(seqregions, 3, -1, 0, 160, threshold_prob);
    EXPECT_EQ(seqregions.size(), 3);
    EXPECT_EQ(seqregions.back().position, 160);
    EXPECT_EQ(seqregions.back().plength_observation2node, -1);
    EXPECT_EQ(seqregions.back().plength_observation2root, 0);
    
    SeqRegions::addNonConsecutiveRRegion(seqregions, TYPE_R, -1, -1, 223, threshold_prob);
    EXPECT_EQ(seqregions.size(), 4);
    EXPECT_EQ(seqregions.back().position, 223);
    EXPECT_EQ(seqregions.back().plength_observation2node, -1);
    EXPECT_EQ(seqregions.back().plength_observation2root, -1);
    
    SeqRegions::addNonConsecutiveRRegion(seqregions, TYPE_R, -1, 1e-10, 240, threshold_prob);
    EXPECT_EQ(seqregions.size(), 5);
    EXPECT_EQ(seqregions.back().position, 240);
    EXPECT_EQ(seqregions.back().plength_observation2node, -1);
    EXPECT_EQ(seqregions.back().plength_observation2root, 1e-10);
    
    SeqRegions::addNonConsecutiveRRegion(seqregions, TYPE_R, -1, 1e-11, 264, threshold_prob);
    EXPECT_EQ(seqregions.size(), 5); // merge two consecotive R
    EXPECT_EQ(seqregions.back().position, 264);
    EXPECT_EQ(seqregions.back().plength_observation2node, -1);
    EXPECT_EQ(seqregions.back().plength_observation2root, 1e-11);
    
    SeqRegions::addNonConsecutiveRRegion(seqregions, TYPE_R, 0, 1e-11, 289, threshold_prob);
    EXPECT_EQ(seqregions.size(), 6);
    EXPECT_EQ(seqregions.back().position, 289);
    EXPECT_EQ(seqregions.back().plength_observation2node, 0);
    EXPECT_EQ(seqregions.back().plength_observation2root, 1e-11);
    
    SeqRegions::addNonConsecutiveRRegion(seqregions, TYPE_R, 1e-9, 0, 299, threshold_prob);
    EXPECT_EQ(seqregions.size(), 6); // merge two consecotive R
    EXPECT_EQ(seqregions.back().position, 299);
    EXPECT_EQ(seqregions.back().plength_observation2node, 1e-9);
    EXPECT_EQ(seqregions.back().plength_observation2root, 0);
    
    SeqRegions::addNonConsecutiveRRegion(seqregions, TYPE_N, 1e-9, 0.2, 311, threshold_prob);
    EXPECT_EQ(seqregions.size(), 7);
    EXPECT_EQ(seqregions.back().position, 311);
    EXPECT_EQ(seqregions.back().plength_observation2node, 1e-9);
    EXPECT_EQ(seqregions.back().plength_observation2root, 0.2);
    
    SeqRegions::addNonConsecutiveRRegion(seqregions, TYPE_R, 1e-9, 0, 324, threshold_prob);
    EXPECT_EQ(seqregions.size(), 8);
    EXPECT_EQ(seqregions.back().position, 324);
    EXPECT_EQ(seqregions.back().plength_observation2node, 1e-9);
    EXPECT_EQ(seqregions.back().plength_observation2root, 0);
    
    SeqRegions::addNonConsecutiveRRegion(seqregions, TYPE_R, 1e-10, 1e-11, 324, threshold_prob);
    EXPECT_EQ(seqregions.size(), 8); // merge two consecotive R
    EXPECT_EQ(seqregions.back().position, 324);
    EXPECT_EQ(seqregions.back().plength_observation2node, 1e-10);
    EXPECT_EQ(seqregions.back().plength_observation2root, 1e-11);
    
    SeqRegions::addNonConsecutiveRRegion(seqregions, TYPE_R, 1e-15, 1e-10, 356, threshold_prob);
    EXPECT_EQ(seqregions.size(), 8); // merge two consecotive R
    EXPECT_EQ(seqregions.back().position, 356);
    EXPECT_EQ(seqregions.back().plength_observation2node, 1e-15);
    EXPECT_EQ(seqregions.back().plength_observation2root, 1e-10);
    
    SeqRegions::addNonConsecutiveRRegion(seqregions, TYPE_R, 1e-15, 1e-300, 382, threshold_prob);
    EXPECT_EQ(seqregions.size(), 8); // merge two consecotive R
    EXPECT_EQ(seqregions.back().position, 382);
    EXPECT_EQ(seqregions.back().plength_observation2node, 1e-15);
    EXPECT_EQ(seqregions.back().plength_observation2root, 1e-300);
    
    SeqRegions::addNonConsecutiveRRegion(seqregions, TYPE_R, 1e-5, 1e-301, 545, threshold_prob);
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
    EXPECT_EQ(SeqRegions::simplifyO(new_lh->data(), 2, 4, threshold_prob), TYPE_O) ;
    
    SeqRegion::LHType new_lh_value1{1.0 - 3 * threshold_prob, threshold_prob, threshold_prob, threshold_prob};
    (*new_lh) = new_lh_value1;
    EXPECT_EQ(SeqRegions::simplifyO(new_lh->data(), 2, 4, threshold_prob), 0) ;
    
    new_lh->data()[1] += new_lh->data()[2];
    new_lh->data()[2] = 0;
    EXPECT_EQ(SeqRegions::simplifyO(new_lh->data(), 2, 4, threshold_prob), TYPE_O) ;
    
    new_lh->data()[2] = new_lh->data()[0];
    new_lh->data()[0] = new_lh->data()[3];
    new_lh->data()[1] = new_lh->data()[3];
    EXPECT_EQ(SeqRegions::simplifyO(new_lh->data(), 2, 4, threshold_prob), TYPE_R) ;
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
    Initialize Alignment, Model, and Parameters
 */
void initAlnModelParams(Params& params, Alignment& aln, Model& model, const std::string model_name = "GTR")
{
    // Init params, aln, and model
    initDefaultValue(params);
    params.model_name = model_name;
    std::string diff_file_path("../../example/test_5K.diff");
    char* diff_file_path_ptr = new char[diff_file_path.length() + 1];
    strcpy(diff_file_path_ptr, diff_file_path.c_str());
    aln.readDiff(diff_file_path_ptr, NULL);
    // extract related info (freqs, log_freqs) of the ref sequence
    model.extractRefInfo(aln.ref_seq, aln.num_states);
    // init the mutation matrix from a model name
    model.initMutationMat(params.model_name, aln.num_states);
    // compute cumulative rates of the ref sequence
    model.computeCumulativeRate(aln);
}

/*
 Test computeAbsoluteLhAtRoot(const Alignment& aln, const Model& model)
 */
TEST(SeqRegions, computeAbsoluteLhAtRoot)
{
    Alignment aln;
    Model model;
    Params params = Params::getInstance();
    
    // Init params, aln, and model
    initAlnModelParams(params, aln, model);
    // dummy variables
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    
    // ----- Simple tests -----
    SeqRegions seqregions1;
    seqregions1.emplace_back(TYPE_R, seq_length - 1);
    EXPECT_EQ(seqregions1.computeAbsoluteLhAtRoot(num_states, model), -40547.865582541948);
    
    seqregions1.emplace(seqregions1.begin(), 3, 3242);
    seqregions1.emplace(seqregions1.begin(), TYPE_R, 3241);
    EXPECT_EQ(seqregions1.computeAbsoluteLhAtRoot(num_states, model), -40547.372963956892);
    
    auto new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value{0.1,0.3,0.2,0.4};
    (*new_lh) = new_lh_value;
    seqregions1.emplace(seqregions1.begin(), TYPE_O, 243, 0, -1, std::move(new_lh));
    seqregions1.emplace(seqregions1.begin(), TYPE_R, 242);
    EXPECT_EQ(seqregions1.computeAbsoluteLhAtRoot(num_states, model), -40547.053844652924);
    // ----- Simple tests -----
    
    // Generate complex seqregions
    SeqRegions seqregions2, seqregions3, seqregions4;
    genTestData2(seqregions2, seqregions3, seqregions4);
    
    // ----- Test 1 on a more complex seqregions -----
    EXPECT_EQ(seqregions2.computeAbsoluteLhAtRoot(num_states, model), -40547.644608026116);
    // ----- Test 1 on a more complex seqregions -----
    
    // ----- Test 2 on a more complex seqregions -----
    EXPECT_EQ(seqregions3.computeAbsoluteLhAtRoot(num_states, model), -40548.188644939095);
    // ----- Test 2 on a more complex seqregions -----
    
    // ----- Test 3 on a more complex seqregions -----
    EXPECT_EQ(seqregions4.computeAbsoluteLhAtRoot(num_states, model), -40548.424295613549);
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
    
    // Init params, aln, and model
    initAlnModelParams(params, aln, model);
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

/*
 Test merge_N_O(const RealNumType lower_plength, const SeqRegion& reg_o, const Model& model,
 const PositionType end_pos, const StateType num_states, SeqRegions& merged_target)
 */
TEST(SeqRegions, merge_N_O)
{
    Alignment aln;
    Model model;
    Params params = Params::getInstance();
    
    // Init params, aln, and model
    initAlnModelParams(params, aln, model);
    // dummy variables
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    
    // ----- Test 1 -----
    SeqRegions merged_regions;
    const PositionType end_pos = 8623;
    RealNumType lower_plength = -1;
    auto new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value{0.27595781923293211113090706021466758102178573608398,8.3817209383128356075479820086471249851456377655268e-06,0.72401575181023847260775028189527802169322967529297,1.8047235891267821277644811672757896303664892911911e-05};
    (*new_lh) = new_lh_value;
    SeqRegion seqregion1(TYPE_O, 243, -1, -1, std::move(new_lh));
    merge_N_O(lower_plength, seqregion1, model,
              end_pos, num_states, merged_regions);
    EXPECT_EQ(merged_regions.size(), 1);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge{0.3675361745020691,0.0000068532680661245309,0.63243117236202639,0.000025799867838435163};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    seqregion1.plength_observation2root = 0.311; // plength_observation2root does NOT affect the likelihood of new merged region
    merge_N_O(lower_plength, seqregion1, model,
              end_pos, num_states, merged_regions);
    EXPECT_EQ(merged_regions.size(), 2);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    seqregion1.plength_observation2node = 0; // plength_observation2node = 0 does NOT affect the likelihood of new merged region
    merge_N_O(lower_plength, seqregion1, model,
              end_pos, num_states, merged_regions);
    EXPECT_EQ(merged_regions.size(), 3);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge);
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    seqregion1.plength_observation2node = 0.1252; // plength_observation2node > 0 affects the likelihood of new merged region
    merge_N_O(lower_plength, seqregion1, model,
              end_pos, num_states, merged_regions);
    EXPECT_EQ(merged_regions.size(), 4);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge1{0.37599425541307963,0.0099624992407331657,0.55990605200203047,0.054137193344156891};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge1);
    // ----- Test 4 -----
    
    // ----- Test 5 -----
    lower_plength = 0; // lower_plength = 0; plength_observation2node > 0
    merge_N_O(lower_plength, seqregion1, model,
              end_pos, num_states, merged_regions);
    EXPECT_EQ(merged_regions.size(), 5);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge1);
    // ----- Test 5 -----
    
    // ----- Test 6 -----
    lower_plength = 0; // lower_plength = 0; plength_observation2node = 0
    seqregion1.plength_observation2node = 0;
    merge_N_O(lower_plength, seqregion1, model,
              end_pos, num_states, merged_regions);
    EXPECT_EQ(merged_regions.size(), 6);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge);
    // ----- Test 6 -----
    
    // ----- Test 7 -----
    lower_plength = 0; // lower_plength = 0; plength_observation2node = -1
    seqregion1.plength_observation2node = -1;
    merge_N_O(lower_plength, seqregion1, model,
              end_pos, num_states, merged_regions);
    EXPECT_EQ(merged_regions.size(), 7);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge);
    // ----- Test 7 -----
    
    // ----- Test 8 -----
    lower_plength = 0.02331; // lower_plength > 0; plength_observation2node = -1
    seqregion1.plength_observation2node = -1;
    merge_N_O(lower_plength, seqregion1, model,
              end_pos, num_states, merged_regions);
    EXPECT_EQ(merged_regions.size(), 8);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge2{0.36911091784101202,0.0018604164279931914,0.61892829252503356,0.010100373205961303};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge2);
    // ----- Test 8 -----
    
    // ----- Test 9 -----
    lower_plength = 0.02331; // lower_plength > 0; plength_observation2node = 0
    seqregion1.plength_observation2node = 0;
    merge_N_O(lower_plength, seqregion1, model,
              end_pos, num_states, merged_regions);
    EXPECT_EQ(merged_regions.size(), 9);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge2);
    // ----- Test 9 -----
    
    // ----- Test 10 -----
    lower_plength = 0.02331; // lower_plength > 0; plength_observation2node > 0
    seqregion1.plength_observation2node = 0.1e-6;
    merge_N_O(lower_plength, seqregion1, model,
              end_pos, num_states, merged_regions);
    EXPECT_EQ(merged_regions.size(), 10);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge3{0.36911092459666772,0.0018604243797870988,0.61892823459762114,0.010100416425924146};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge3);
    // ----- Test 10 -----
}

/*
 Test merge_N_RACGT(const SeqRegion& reg_racgt, const RealNumType lower_plength, const PositionType end_pos,
 const RealNumType threshold_prob, SeqRegions& merged_regions)
 */
TEST(SeqRegions, merge_N_RACGT)
{
    Alignment aln;
    Model model;
    Params params = Params::getInstance();
    
    // Init params, aln, and model
    initAlnModelParams(params, aln, model);
    // dummy variables
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    
    // ----- Test 1 -----
    SeqRegions merged_regions;
    const PositionType end_pos = 6543;
    RealNumType lower_plength = -1;
    SeqRegion seqregion1(TYPE_R, 3442, -1, -1);
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 1);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    lower_plength = -1;
    seqregion1.plength_observation2node = -1;
    seqregion1.plength_observation2root = 0;
    seqregion1.type = 1; // C
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 2);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    lower_plength = -1;
    seqregion1.plength_observation2node = -1;
    seqregion1.plength_observation2root = 0.3231;
    seqregion1.type = TYPE_R;
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 3);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    lower_plength = -1;
    seqregion1.plength_observation2node = 0;
    seqregion1.plength_observation2root = -1;
    seqregion1.type = 3;
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 4);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 4 -----
    
    // ----- Test 5 -----
    lower_plength = -1;
    seqregion1.plength_observation2node = 0;
    seqregion1.plength_observation2root = 0;
    seqregion1.type = 2;
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 5);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 5 -----
    
    // ----- Test 6 -----
    lower_plength = -1;
    seqregion1.plength_observation2node = 0;
    seqregion1.plength_observation2root = 1e-3;
    seqregion1.type = TYPE_R;
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 6);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 6 -----
    
    // ----- Test 7 -----
    lower_plength = -1;
    seqregion1.plength_observation2node = 1e-10;
    seqregion1.plength_observation2root = -1;
    seqregion1.type = 1;
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 7);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 1e-10);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 7 -----
    
    // ----- Test 8 -----
    lower_plength = -1;
    seqregion1.plength_observation2node = 1e-10;
    seqregion1.plength_observation2root = 0;
    seqregion1.type = 2;
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 8);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 1e-10);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 8 -----
    
    // ----- Test 9 -----
    lower_plength = -1;
    seqregion1.plength_observation2node = 1e-10;
    seqregion1.plength_observation2root = 0.1;
    seqregion1.type = 0;
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 9);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 1e-10);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 9 -----
    
    // ----- Test 10 -----
    lower_plength = 0;
    seqregion1.plength_observation2node = -1;
    seqregion1.plength_observation2root = -1;
    seqregion1.type = 1;
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 10);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 10 -----
    
    // ----- Test 11 -----
    lower_plength = 0;
    seqregion1.plength_observation2node = -1;
    seqregion1.plength_observation2root = 0;
    seqregion1.type = TYPE_R;
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 11);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 11 -----
    
    // ----- Test 12 -----
    lower_plength = 0;
    seqregion1.plength_observation2node = -1;
    seqregion1.plength_observation2root = 0.3231;
    seqregion1.type = 3;
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 12);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 12 -----
    
    // ----- Test 13 -----
    lower_plength = 0;
    seqregion1.plength_observation2node = 0;
    seqregion1.plength_observation2root = -1;
    seqregion1.type = 1;
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 13);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 13 -----
    
    // ----- Test 14 -----
    lower_plength = 0;
    seqregion1.plength_observation2node = 0;
    seqregion1.plength_observation2root = 0;
    seqregion1.type = 0;
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 14);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 14 -----
    
    // ----- Test 15 -----
    lower_plength = 0;
    seqregion1.plength_observation2node = 0;
    seqregion1.plength_observation2root = 1e-3;
    seqregion1.type = 0;
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 15);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 15 -----
    
    // ----- Test 16 -----
    lower_plength = 0;
    seqregion1.plength_observation2node = 1e-10;
    seqregion1.plength_observation2root = -1;
    seqregion1.type = 1;
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 16);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 1e-10);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 16 -----
    
    // ----- Test 17 -----
    lower_plength = 0;
    seqregion1.plength_observation2node = 1e-10;
    seqregion1.plength_observation2root = 0;
    seqregion1.type = 1;
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 17);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 1e-10);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 17 -----
    
    // ----- Test 18 -----
    lower_plength = 0;
    seqregion1.plength_observation2node = 1e-10;
    seqregion1.plength_observation2root = 0.1;
    seqregion1.type = TYPE_R;
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 18);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 1e-10);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 18 -----
    
    // ----- Test 19 -----
    lower_plength = 1e-5;
    seqregion1.plength_observation2node = -1;
    seqregion1.plength_observation2root = -1;
    seqregion1.type = 3;
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 19);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 1e-5);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 19 -----
    
    // ----- Test 20 -----
    lower_plength = 1e-5;
    seqregion1.plength_observation2node = -1;
    seqregion1.plength_observation2root = 0;
    seqregion1.type = TYPE_R;
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 20);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 1e-5);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 20 -----
    
    // ----- Test 21 -----
    lower_plength = 0.212;
    seqregion1.plength_observation2node = -1;
    seqregion1.plength_observation2root = 0.3231;
    seqregion1.type = 0;
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 21);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0.212);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 21 -----
    
    // ----- Test 22 -----
    lower_plength = 0.011;
    seqregion1.plength_observation2node = 0;
    seqregion1.plength_observation2root = -1;
    seqregion1.type = 3;
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 22);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0.011);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 22 -----
    
    // ----- Test 23 -----
    lower_plength = 1e-3;
    seqregion1.plength_observation2node = 0;
    seqregion1.plength_observation2root = 0;
    seqregion1.type = 3;
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 23);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 1e-3);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 23 -----
    
    // ----- Test 24 -----
    lower_plength = 0.2312;
    seqregion1.plength_observation2node = 0;
    seqregion1.plength_observation2root = 1e-3;
    seqregion1.type = 2;
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 24);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0.2312);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 24 -----
    
    // ----- Test 25 -----
    lower_plength = 0.04123;
    seqregion1.plength_observation2node = 1e-10;
    seqregion1.plength_observation2root = -1;
    seqregion1.type = 0;
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 25);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0.04123 + 1e-10);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 25 -----
    
    // ----- Test 26 -----
    lower_plength = 0.001243;
    seqregion1.plength_observation2node = 1e-3;
    seqregion1.plength_observation2root = 0;
    seqregion1.type = TYPE_R;
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 26);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0.002243);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 26 -----
    
    // ----- Test 27 -----
    lower_plength = 0.0001211;
    seqregion1.plength_observation2node = 1e-1;
    seqregion1.plength_observation2root = 0.1;
    seqregion1.type = 1;
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 27);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0.1001211);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 27 -----
}

/*
 Test merge_O_N(const SeqRegion& reg_o, const RealNumType upper_plength, const PositionType end_pos, const Model& model, const StateType num_states, SeqRegions& merged_regions)
 */
TEST(SeqRegions, merge_O_N)
{
    Alignment aln;
    Model model;
    Params params = Params::getInstance();
    
    // Init params, aln, and model
    initAlnModelParams(params, aln, model);
    // dummy variables
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    
    // ----- Test 1 -----
    SeqRegions merged_regions;
    const PositionType end_pos = 29180;
    RealNumType upper_plength = -1;
    auto new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value{7.4346402191731947664294494204639818235591519623995e-06,0.66666418845326025355291221785591915249824523925781,7.4346402191731947664294494204639818235591519623995e-06,0.33332094226630137878686355179524980485439300537109};
    (*new_lh) = new_lh_value;
    SeqRegion seqregion1(TYPE_O, 243, -1, -1, std::move(new_lh));
    merge_O_N(seqregion1, upper_plength, end_pos, model, num_states, merged_regions);
    EXPECT_EQ(merged_regions.size(), 1);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge{0.0000074346402191731948,0.66666418845326025,0.0000074346402191731948,0.33332094226630138};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    seqregion1.plength_observation2root = 0.311; // plength_observation2root does NOT affect the likelihood of new merged region
    merge_O_N(seqregion1, upper_plength, end_pos, model, num_states, merged_regions);
    EXPECT_EQ(merged_regions.size(), 2);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    seqregion1.plength_observation2node = 0; // plength_observation2node = 0 does NOT affect the likelihood of new merged region
    merge_O_N(seqregion1, upper_plength, end_pos, model, num_states, merged_regions);
    EXPECT_EQ(merged_regions.size(), 3);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge);
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    seqregion1.plength_observation2node = 0.1252; // plength_observation2node > 0 affects the likelihood of new merged region
    merge_O_N(seqregion1, upper_plength, end_pos, model, num_states, merged_regions);
    EXPECT_EQ(merged_regions.size(), 4);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge1{0.011218665480774669,0.5673625983091064,0.024370520108437949,0.39704821610168095};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge1);
    // ----- Test 4 -----
    
    // ----- Test 5 -----
    upper_plength = 0; // lower_plength = 0; plength_observation2node > 0
    merge_O_N(seqregion1, upper_plength, end_pos, model, num_states, merged_regions);
    EXPECT_EQ(merged_regions.size(), 5);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge1);
    // ----- Test 5 -----
    
    // ----- Test 6 -----
    upper_plength = 0; // lower_plength = 0; plength_observation2node = 0
    seqregion1.plength_observation2node = 0;
    merge_O_N(seqregion1, upper_plength, end_pos, model, num_states, merged_regions);
    EXPECT_EQ(merged_regions.size(), 6);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge);
    // ----- Test 6 -----
    
    // ----- Test 7 -----
    upper_plength = 0; // lower_plength = 0; plength_observation2node = -1
    seqregion1.plength_observation2node = -1;
    merge_O_N(seqregion1, upper_plength, end_pos, model, num_states, merged_regions);
    EXPECT_EQ(merged_regions.size(), 7);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge);
    // ----- Test 7 -----
    
    // ----- Test 8 -----
    upper_plength = 0.02331; // lower_plength > 0; plength_observation2node = -1
    seqregion1.plength_observation2node = -1;
    merge_O_N(seqregion1, upper_plength, end_pos, model, num_states, merged_regions);
    EXPECT_EQ(merged_regions.size(), 8);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge2{0.0020947652384088584,0.64817600901028716,0.00454340526533243,0.34518582048597146};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge2);
    // ----- Test 8 -----
    
    // ----- Test 9 -----
    upper_plength = 0.02331; // lower_plength > 0; plength_observation2node = 0
    seqregion1.plength_observation2node = 0;
    merge_O_N(seqregion1, upper_plength, end_pos, model, num_states, merged_regions);
    EXPECT_EQ(merged_regions.size(), 9);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge2);
    // ----- Test 9 -----
    
    // ----- Test 10 -----
    upper_plength = 0.02331; // lower_plength > 0; plength_observation2node > 0
    seqregion1.plength_observation2node = 0.1e-6;
    merge_O_N(seqregion1, upper_plength, end_pos, model, num_states, merged_regions);
    EXPECT_EQ(merged_regions.size(), 10);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge3{0.0020947741930660798,0.6481759296959182,0.0045434247246658715,0.34518587138634999};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge3);
    // ----- Test 10 -----
}

/*
 Test merge_RACGT_N(const SeqRegion& reg_n, const RealNumType upper_plength, const PositionType end_pos,
 const RealNumType threshold_prob, SeqRegions& merged_regions)
 */
TEST(SeqRegions, merge_RACGT_N)
{
    Alignment aln;
    Model model;
    Params params = Params::getInstance();
    
    // Init params, aln, and model
    initAlnModelParams(params, aln, model);
    // dummy variables
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    
    // ----- Test 1 -----
    SeqRegions merged_regions;
    const PositionType end_pos = 6543;
    RealNumType lower_plength = -1;
    SeqRegion seqregion1(TYPE_R, 3442, -1, -1);
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 1);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    lower_plength = -1;
    seqregion1.plength_observation2node = -1;
    seqregion1.plength_observation2root = 0;
    seqregion1.type = 1; // C
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 2);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    lower_plength = -1;
    seqregion1.plength_observation2node = -1;
    seqregion1.plength_observation2root = 0.3231;
    seqregion1.type = TYPE_R;
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 3);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0.3231);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    lower_plength = -1;
    seqregion1.plength_observation2node = 0;
    seqregion1.plength_observation2root = -1;
    seqregion1.type = 3;
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 4);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 4 -----
    
    // ----- Test 5 -----
    lower_plength = -1;
    seqregion1.plength_observation2node = 0;
    seqregion1.plength_observation2root = 0;
    seqregion1.type = 2;
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 5);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 5 -----
    
    // ----- Test 6 -----
    lower_plength = -1;
    seqregion1.plength_observation2node = 0;
    seqregion1.plength_observation2root = 1e-3;
    seqregion1.type = TYPE_R;
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 6);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 1e-3);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 6 -----
    
    // ----- Test 7 -----
    lower_plength = -1;
    seqregion1.plength_observation2node = 1e-10;
    seqregion1.plength_observation2root = -1;
    seqregion1.type = 1;
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 7);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 1e-10);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 7 -----
    
    // ----- Test 8 -----
    lower_plength = -1;
    seqregion1.plength_observation2node = 1e-10;
    seqregion1.plength_observation2root = 0;
    seqregion1.type = 2;
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 8);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 1e-10);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 8 -----
    
    // ----- Test 9 -----
    lower_plength = -1;
    seqregion1.plength_observation2node = 1e-10;
    seqregion1.plength_observation2root = 0.1;
    seqregion1.type = 0;
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 9);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0.1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 1e-10);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 9 -----
    
    // ----- Test 10 -----
    lower_plength = 0;
    seqregion1.plength_observation2node = -1;
    seqregion1.plength_observation2root = -1;
    seqregion1.type = 1;
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 10);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 10 -----
    
    // ----- Test 11 -----
    lower_plength = 0;
    seqregion1.plength_observation2node = -1;
    seqregion1.plength_observation2root = 0;
    seqregion1.type = TYPE_R;
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 11);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 11 -----
    
    // ----- Test 12 -----
    lower_plength = 0;
    seqregion1.plength_observation2node = -1;
    seqregion1.plength_observation2root = 0.3231;
    seqregion1.type = 3;
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 12);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0.3231);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 12 -----
    
    // ----- Test 13 -----
    lower_plength = 0;
    seqregion1.plength_observation2node = 0;
    seqregion1.plength_observation2root = -1;
    seqregion1.type = 1;
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 13);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 13 -----
    
    // ----- Test 14 -----
    lower_plength = 0;
    seqregion1.plength_observation2node = 0;
    seqregion1.plength_observation2root = 0;
    seqregion1.type = 0;
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 14);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 14 -----
    
    // ----- Test 15 -----
    lower_plength = 0;
    seqregion1.plength_observation2node = 0;
    seqregion1.plength_observation2root = 1e-3;
    seqregion1.type = 0;
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 15);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 1e-3);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 15 -----
    
    // ----- Test 16 -----
    lower_plength = 0;
    seqregion1.plength_observation2node = 1e-10;
    seqregion1.plength_observation2root = -1;
    seqregion1.type = 1;
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 16);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 1e-10);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 16 -----
    
    // ----- Test 17 -----
    lower_plength = 0;
    seqregion1.plength_observation2node = 1e-10;
    seqregion1.plength_observation2root = 0;
    seqregion1.type = 1;
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 17);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 1e-10);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 17 -----
    
    // ----- Test 18 -----
    lower_plength = 0;
    seqregion1.plength_observation2node = 1e-10;
    seqregion1.plength_observation2root = 0.1;
    seqregion1.type = TYPE_R;
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 18);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0.1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 1e-10);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 18 -----
    
    // ----- Test 19 -----
    lower_plength = 1e-5;
    seqregion1.plength_observation2node = -1;
    seqregion1.plength_observation2root = -1;
    seqregion1.type = 3;
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 19);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 1e-5);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 19 -----
    
    // ----- Test 20 -----
    lower_plength = 1e-5;
    seqregion1.plength_observation2node = -1;
    seqregion1.plength_observation2root = 0;
    seqregion1.type = TYPE_R;
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 20);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 1e-5);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 20 -----
    
    // ----- Test 21 -----
    lower_plength = 0.212;
    seqregion1.plength_observation2node = -1;
    seqregion1.plength_observation2root = 0.3231;
    seqregion1.type = 0;
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 21);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0.5351);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 21 -----
    
    // ----- Test 22 -----
    lower_plength = 0.011;
    seqregion1.plength_observation2node = 0;
    seqregion1.plength_observation2root = -1;
    seqregion1.type = 3;
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 22);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0.011);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 22 -----
    
    // ----- Test 23 -----
    lower_plength = 1e-3;
    seqregion1.plength_observation2node = 0;
    seqregion1.plength_observation2root = 0;
    seqregion1.type = 3;
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 23);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 1e-3);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 23 -----
    
    // ----- Test 24 -----
    lower_plength = 0.2312;
    seqregion1.plength_observation2node = 0;
    seqregion1.plength_observation2root = 1e-3;
    seqregion1.type = 2;
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 24);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0.2322);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 24 -----
    
    // ----- Test 25 -----
    lower_plength = 0.04123;
    seqregion1.plength_observation2node = 1e-10;
    seqregion1.plength_observation2root = -1;
    seqregion1.type = 0;
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 25);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0.04123 + 1e-10);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 25 -----
    
    // ----- Test 26 -----
    lower_plength = 0.001243;
    seqregion1.plength_observation2node = 1e-3;
    seqregion1.plength_observation2root = 0;
    seqregion1.type = TYPE_R;
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 26);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0.001243);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 1e-3);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 26 -----
    
    // ----- Test 27 -----
    lower_plength = 0.0001211;
    seqregion1.plength_observation2node = 1e-1;
    seqregion1.plength_observation2root = 0.1;
    seqregion1.type = 1;
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params.threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 27);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0.1001211);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0.1);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 27 -----
}

/*
 Test merge_Zero_Distance(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType total_blength_1, const RealNumType total_blength_2, const PositionType end_pos, const RealNumType threshold_prob, const StateType num_states, SeqRegions* &merged_regions)
 */
TEST(SeqRegions, merge_Zero_Distance)
{
    Alignment aln;
    Model model;
    Params params = Params::getInstance();
    
    // Init params, aln, and model
    initAlnModelParams(params, aln, model);
    // dummy variables
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    
    // ----- Test 1 -----
    SeqRegions* merged_regions_ptr = new SeqRegions();
    const RealNumType threshold_prob = params.threshold_prob;
    const PositionType end_pos = 6543;
    RealNumType total_blength_1 = -1;
    RealNumType total_blength_2 = -1;
    SeqRegion seqregion1(TYPE_R, 3442, -1, -1);
    SeqRegion seqregion2(0, 234, -1, 0);
    EXPECT_TRUE(merge_Zero_Distance(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos, threshold_prob, num_states, merged_regions_ptr));
    EXPECT_EQ(merged_regions_ptr, nullptr);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    merged_regions_ptr = new SeqRegions();
    total_blength_1 = -1;
    total_blength_2 = 0;
    seqregion1.type = TYPE_R;
    seqregion2.type = 0;
    EXPECT_TRUE(merge_Zero_Distance(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos, threshold_prob, num_states, merged_regions_ptr));
    EXPECT_EQ(merged_regions_ptr, nullptr);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    merged_regions_ptr = new SeqRegions();
    total_blength_1 = -1;
    total_blength_2 = 1e-3;
    seqregion1.type = 3;
    seqregion2.type = 0;
    EXPECT_TRUE(merge_Zero_Distance(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos, threshold_prob, num_states, merged_regions_ptr));
    EXPECT_TRUE(merged_regions_ptr != nullptr);
    EXPECT_EQ(merged_regions_ptr->size(), 1);
    EXPECT_EQ(merged_regions_ptr->back().type, seqregion1.type);
    EXPECT_EQ(merged_regions_ptr->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions_ptr->back().plength_observation2node, -1);
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    total_blength_1 = 0;
    total_blength_2 = -1;
    seqregion1.type = 1;
    seqregion2.type = TYPE_R;
    EXPECT_TRUE(merge_Zero_Distance(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos, threshold_prob, num_states, merged_regions_ptr));
    EXPECT_EQ(merged_regions_ptr, nullptr);
    // ----- Test 4 -----
    
    // ----- Test 5 -----
    merged_regions_ptr = new SeqRegions();
    total_blength_1 = 0;
    total_blength_2 = 0;
    seqregion1.type = 0;
    seqregion2.type = 0;
    EXPECT_TRUE(merge_Zero_Distance(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos, threshold_prob, num_states, merged_regions_ptr));
    EXPECT_EQ(merged_regions_ptr, nullptr);
    // ----- Test 5 -----
    
    // ----- Test 6 -----
    merged_regions_ptr = new SeqRegions();
    total_blength_1 = 0;
    total_blength_2 = 1e-10;
    seqregion1.type = TYPE_R;
    seqregion2.type = TYPE_R;
    EXPECT_TRUE(merge_Zero_Distance(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos, threshold_prob, num_states, merged_regions_ptr));
    EXPECT_TRUE(merged_regions_ptr != nullptr);
    EXPECT_EQ(merged_regions_ptr->size(), 1);
    EXPECT_EQ(merged_regions_ptr->back().type, seqregion2.type);
    EXPECT_EQ(merged_regions_ptr->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions_ptr->back().plength_observation2node, -1);
    // ----- Test 6 -----
    
    // ----- Test 7 -----
    total_blength_1 = 0.0131;
    total_blength_2 = -1;
    seqregion1.type = 1;
    seqregion2.type = 2;
    EXPECT_TRUE(merge_Zero_Distance(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos, threshold_prob, num_states, merged_regions_ptr));
    EXPECT_TRUE(merged_regions_ptr != nullptr);
    EXPECT_EQ(merged_regions_ptr->size(), 2);
    EXPECT_EQ(merged_regions_ptr->back().type, seqregion2.type);
    EXPECT_EQ(merged_regions_ptr->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions_ptr->back().plength_observation2node, -1);
    // ----- Test 7 -----
    
    // ----- Test 8 -----
    total_blength_1 = 1e-4;
    total_blength_2 = 0;
    seqregion1.type = 1;
    seqregion2.type = TYPE_O;
    EXPECT_FALSE(merge_Zero_Distance(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos, threshold_prob, num_states, merged_regions_ptr));
    EXPECT_TRUE(merged_regions_ptr != nullptr);
    EXPECT_EQ(merged_regions_ptr->size(), 2);
    // ----- Test 8 -----
    
    // ----- Test 9 -----
    total_blength_1 = 0.1321;
    total_blength_2 = 1e-5;
    seqregion1.type = TYPE_O;
    seqregion2.type = TYPE_R;
    EXPECT_FALSE(merge_Zero_Distance(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos, threshold_prob, num_states, merged_regions_ptr));
    EXPECT_TRUE(merged_regions_ptr != nullptr);
    EXPECT_EQ(merged_regions_ptr->size(), 2);
    delete merged_regions_ptr;
    // ----- Test 9 -----
}

/*
 Test merge_O_ORACGT(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType total_blength_1, const RealNumType total_blength_2, const PositionType end_pos, const RealNumType threshold_prob, const Model& model, const Alignment& aln, SeqRegions& merged_regions)
 */
TEST(SeqRegions, merge_O_ORACGT)
{
    Alignment aln;
    Model model;
    Params params = Params::getInstance();
    
    // Init params, aln, and model
    initAlnModelParams(params, aln, model);
    // dummy variables
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    
    // ----- Test 1 -----
    SeqRegions merged_regions;
    auto new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value{0.27595781923293211113090706021466758102178573608398,8.3817209383128356075479820086471249851456377655268e-06,0.72401575181023847260775028189527802169322967529297,1.8047235891267821277644811672757896303664892911911e-05};
    (*new_lh) = new_lh_value;
    SeqRegion seqregion1(TYPE_O, 243, -1, -1, std::move(new_lh));
    new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value1{0.27595781923293211113090706021466758102178573608398,8.3817209383128356075479820086471249851456377655268e-06,0.72401575181023847260775028189527802169322967529297,1.8047235891267821277644811672757896303664892911911e-05};
    (*new_lh) = new_lh_value1;
    SeqRegion seqregion2(TYPE_O, 243, -1, 1e-3, std::move(new_lh));
    RealNumType total_blength_1 = -1;
    RealNumType total_blength_2 = -1;
    const RealNumType threshold_prob = params.threshold_prob;
    const PositionType end_pos = 3213;
    merge_O_ORACGT(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 1);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge{0.12684687976595482,1.1702018350525203E-10,0.87315311957450514,5.425200212298543E-10};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    total_blength_1 = -1;
    total_blength_2 = 1e-5;
    merge_O_ORACGT(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 2);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge1{0.12684809781766063,1.3059895886789165E-10,0.87315190141833243,6.3340794416577515E-10};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge1);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    total_blength_1 = 0.01;
    total_blength_2 = -1;
    merge_O_ORACGT(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 3);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge2{0.12842468075362262,1.1709109654051516E-8,0.87157516057518813,1.4696207956334386E-7};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge2);
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    total_blength_1 = 0.011321;
    total_blength_2 = 1e-3;
    merge_O_ORACGT(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 4);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge3{0.12875794255340739,1.6716816801441426E-7,0.87123893271963804,0.0000029575587865296985};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge3);
    // ----- Test 4 -----
    
    // ----- Test 5 -----
    total_blength_1 = -1;
    total_blength_2 = -1;
    SeqRegion seqregion3(TYPE_R, 4324, 0.121, 0);
    merge_O_ORACGT(seqregion1, seqregion3, total_blength_1, total_blength_2, end_pos, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 5);
    EXPECT_EQ(merged_regions.back().type, seqregion3.type);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    // ----- Test 5 -----
    
    // ----- Test 6 -----
    total_blength_1 = -1;
    total_blength_2 = 0; // either total_blength_2 = -1 or = 0 leads to the same result
    seqregion3.type = 2;
    merge_O_ORACGT(seqregion1, seqregion3, total_blength_1, total_blength_2, end_pos, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 6);
    EXPECT_EQ(merged_regions.back().type, seqregion3.type);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    // ----- Test 6 -----
    
    // ----- Test 7 -----
    total_blength_1 = -1;
    total_blength_2 = 0;
    seqregion3.type = 0;
    merge_O_ORACGT(seqregion1, seqregion3, total_blength_1, total_blength_2, end_pos, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 7);
    EXPECT_EQ(merged_regions.back().type, seqregion3.type);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    // ----- Test 7 -----
    
    // ----- Test 8 -----
    total_blength_1 = 0.231;
    total_blength_2 = 0;
    seqregion3.type = 0;
    merge_O_ORACGT(seqregion1, seqregion3, total_blength_1, total_blength_2, end_pos, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 8);
    EXPECT_EQ(merged_regions.back().type, seqregion3.type);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    // ----- Test 8 -----
    
    // ----- Test 9 -----
    total_blength_1 = 0.231;
    total_blength_2 = 0.3121;
    seqregion3.type = 2;
    merge_O_ORACGT(seqregion1, seqregion3, total_blength_1, total_blength_2, end_pos, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 9);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge4{0.045597223147683601,0.0013901861006516344,0.92036965105659874,0.032642939695066098};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge4);
    // ----- Test 9 -----
}

/*
 Test merge_RACGT_O(const SeqRegion& seq2_region, const RealNumType total_blength_2, const PositionType end_pos, SeqRegion::LHType& new_lh, const RealNumType threshold_prob, const Model& model, const Alignment& aln, SeqRegions& merged_regions)
 */
TEST(SeqRegions, merge_RACGT_O)
{
    Alignment aln;
    Model model;
    Params params = Params::getInstance();
    
    // Init params, aln, and model
    initAlnModelParams(params, aln, model);
    // dummy variables
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    
    // ----- Test 1 -----
    SeqRegions merged_regions;
    auto new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value{0.49999442406129068761089229155913926661014556884766,5.5759387092817078802929955938516570768115343526006e-06,0.49999442406129068761089229155913926661014556884766,5.5759387092817078802929955938516570768115343526006e-06};
    (*new_lh) = new_lh_value;
    SeqRegion seqregion1(TYPE_O, 243, -1, -1, std::move(new_lh));
    RealNumType total_blength_1 = -1;
    const PositionType end_pos = 3213;
    const RealNumType threshold_prob = params.threshold_prob;
    new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value1{0.99997490659690790870683940738672390580177307128906,2.5092097243268992381179730011275808010395849123597e-05,6.5292438309787796439726148030194621818544931102224e-10,6.5292438309787796439726148030194621818544931102224e-10};
    (*new_lh) = new_lh_value1;
    merge_RACGT_O(seqregion1, total_blength_1, end_pos, *new_lh, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 1);
    EXPECT_EQ(merged_regions.back().type, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    SeqRegion::LHType new_lh_value2{0.24998327274350143345493791002809302881360054016113,0.74997769674260927885711680573876947164535522460938,1.9515256944661845116837164959555650511902058497071e-05,1.9515256944661845116837164959555650511902058497071e-05};
    (*new_lh) = new_lh_value2;
    merge_RACGT_O(seqregion1, total_blength_1, end_pos, *new_lh, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 2);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge{0.99988848806554697,0.000033453518158479732,0.000078057545796795158,8.704978900055221E-10};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    total_blength_1 = 0; // either total_blength_1 = 0 or = -1 does not affect the merged region
    (*new_lh) = new_lh_value2;
    merge_RACGT_O(seqregion1, total_blength_1, end_pos, *new_lh, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 3);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge);
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    total_blength_1 = 1e-4;
    (*new_lh) = new_lh_value2;
    merge_RACGT_O(seqregion1, total_blength_1, end_pos, *new_lh, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 4);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge1{0.99983015675026043,0.000091790366085780112,0.000078048395545374228,4.4881083285432782E-9};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge1);
    // ----- Test 4 -----
    
    // ----- Test 5 -----
    seqregion1.plength_observation2root = 0.231; // neither plength_observation2root nor plength_observation2node affect the merged region
    total_blength_1 = 1e-4;
    (*new_lh) = new_lh_value2;
    merge_RACGT_O(seqregion1, total_blength_1, end_pos, *new_lh, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 5);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge1);
    // ----- Test 5 -----
    
    // ----- Test 6 -----
    seqregion1.plength_observation2node = 1e-7; // neither plength_observation2root nor plength_observation2node affect the merged region
    total_blength_1 = 1e-4;
    (*new_lh) = new_lh_value2;
    merge_RACGT_O(seqregion1, total_blength_1, end_pos, *new_lh, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 6);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge1);
    // ----- Test 6 -----
}

/*
 Test merge_RACGT_RACGT(const SeqRegion& seq2_region, const RealNumType total_blength_2, const PositionType end_pos, SeqRegion::LHType& new_lh, const Model& model, const Alignment& aln, SeqRegions& merged_regions)
 */
TEST(SeqRegions, merge_RACGT_RACGT)
{
    Alignment aln;
    Model model;
    Params params = Params::getInstance();
    
    // Init params, aln, and model
    initAlnModelParams(params, aln, model);
    // dummy variables
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    
    // ----- Test 1 -----
    SeqRegions merged_regions;
    SeqRegion seqregion1(TYPE_R, 4231);
    RealNumType total_blength_1 = -1;
    const PositionType end_pos = 3213;
    const RealNumType threshold_prob = params.threshold_prob;
    auto new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value1{5.5759387092817078802929955938516570768115343526006e-06,0.49999442406129068761089229155913926661014556884766,5.5759387092817078802929955938516570768115343526006e-06,0.49999442406129068761089229155913926661014556884766};
    (*new_lh) = new_lh_value1;
    merge_RACGT_RACGT(seqregion1, total_blength_1, end_pos, *new_lh, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 1);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge{0,0,0,1};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    SeqRegion::LHType new_lh_value2{8.3640013382402129011325767060647251582850003615022e-06,8.3640013382402129011325767060647251582850003615022e-06,0.24998327199732350845096107150311581790447235107422,0.74999999999999988897769753748434595763683319091797};
    (*new_lh) = new_lh_value2;
    merge_RACGT_RACGT(seqregion1, total_blength_1, end_pos, *new_lh, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 2);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    total_blength_1 = 0; // either total_blength_1 = 0 or = -1 does not affect the merged region
    (*new_lh) = new_lh_value2;
    merge_RACGT_RACGT(seqregion1, total_blength_1, end_pos, *new_lh, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 3);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge);
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    total_blength_1 = 1e-4;
    (*new_lh) = new_lh_value2;
    merge_RACGT_RACGT(seqregion1, total_blength_1, end_pos, *new_lh, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 4);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge1{8.8777603721537831E-11,1.5545501848967752E-9,0.000021239862924016271,0.99997875849374817};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge1);
    // ----- Test 4 -----
    
    // ----- Test 5 -----
    seqregion1.plength_observation2root = 0.231; // neither plength_observation2root nor plength_observation2node affect the merged region
    total_blength_1 = 1e-4;
    (*new_lh) = new_lh_value2;
    merge_RACGT_RACGT(seqregion1, total_blength_1, end_pos, *new_lh, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 5);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge1);
    // ----- Test 5 -----
    
    // ----- Test 6 -----
    seqregion1.plength_observation2node = 1e-7; // neither plength_observation2root nor plength_observation2node affect the merged region
    total_blength_1 = 1e-4;
    (*new_lh) = new_lh_value2;
    merge_RACGT_RACGT(seqregion1, total_blength_1, end_pos, *new_lh, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 6);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge1);
    // ----- Test 6 -----
    
    // ----- Test 7 -----
    SeqRegion seqregion2(1, 5345);
    total_blength_1 = -1;
    (*new_lh) = new_lh_value1;
    merge_RACGT_RACGT(seqregion2, total_blength_1, end_pos, *new_lh, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 7);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge2{0,1,0,0};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge2);
    // ----- Test 7 -----
    
    // ----- Test 8 -----
    (*new_lh) = new_lh_value2;
    merge_RACGT_RACGT(seqregion2, total_blength_1, end_pos, *new_lh, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 8);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge2);
    // ----- Test 8 -----
    
    // ----- Test 9 -----
    total_blength_1 = 0; // either total_blength_1 = 0 or = -1 does not affect the merged region
    (*new_lh) = new_lh_value2;
    merge_RACGT_RACGT(seqregion2, total_blength_1, end_pos, *new_lh, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 9);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge2);
    // ----- Test 9 -----
    
    // ----- Test 10 -----
    total_blength_1 = 1e-4;
    (*new_lh) = new_lh_value2;
    merge_RACGT_RACGT(seqregion2, total_blength_1, end_pos, *new_lh, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 10);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge3{7.0898429515507956E-7,0.11874116559913409,0.032309155233313958,0.84894897018325677};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge3);
    // ----- Test 10 -----
    
    // ----- Test 11 -----
    seqregion1.plength_observation2root = 0.231; // neither plength_observation2root nor plength_observation2node affect the merged region
    total_blength_1 = 1e-4;
    (*new_lh) = new_lh_value2;
    merge_RACGT_RACGT(seqregion2, total_blength_1, end_pos, *new_lh, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 11);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge3);
    // ----- Test 11 -----
    
    // ----- Test 12 -----
    seqregion1.plength_observation2node = 1e-7; // neither plength_observation2root nor plength_observation2node affect the merged region
    total_blength_1 = 1e-4;
    (*new_lh) = new_lh_value2;
    merge_RACGT_RACGT(seqregion2, total_blength_1, end_pos, *new_lh, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 12);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge3);
    // ----- Test 12 -----
}

/*
 Test merge_RACGT_ORACGT(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType total_blength_1, const RealNumType total_blength_2, const RealNumType upper_plength, const PositionType end_pos, const RealNumType threshold_prob, const Model& model, const Alignment& aln, SeqRegions& merged_regions)
 */
TEST(SeqRegions, merge_RACGT_ORACGT)
{
    // NOTE: if plength_observation2root > 0 then must be plength_observation2node != -1;
    Alignment aln;
    Model model;
    Params params = Params::getInstance();
    
    // Init params, aln, and model
    initAlnModelParams(params, aln, model);
    // dummy variables
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    
    // ----- Test 1 -----
    SeqRegions merged_regions;
    const RealNumType threshold_prob = params.threshold_prob;
    RealNumType upper_plength = -1;
    RealNumType total_blength_1 = -1;
    RealNumType total_blength_2 = -1;
    const PositionType end_pos = 5432;
    SeqRegion seqregion1(TYPE_R, 2131);
    auto new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value1{5.5759387092817078802929955938516570768115343526006e-06,0.49999442406129068761089229155913926661014556884766,5.5759387092817078802929955938516570768115343526006e-06,0.49999442406129068761089229155913926661014556884766};
    (*new_lh) = new_lh_value1;
    SeqRegion seqregion2(TYPE_O, 543, -1, -1, std::move(new_lh));
    merge_RACGT_ORACGT(seqregion1, seqregion2, total_blength_1, total_blength_2, upper_plength, end_pos, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 1);
    EXPECT_EQ(merged_regions.back().type, TYPE_R);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    upper_plength = -1;
    total_blength_1 = 1e-3;
    total_blength_2 = -1;
    seqregion1.type = TYPE_R;
    seqregion1.plength_observation2root = -1;
    seqregion1.plength_observation2node = -1;
    seqregion2.type = 0;
    seqregion2.plength_observation2root = -1;
    seqregion2.plength_observation2node = -1;
    merge_RACGT_ORACGT(seqregion1, seqregion2, total_blength_1, total_blength_2, upper_plength, end_pos, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 2);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge{1,0,0,0};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    upper_plength = 0;
    total_blength_1 = 1e-3;
    total_blength_2 = -1;
    seqregion1.type = 2;
    seqregion1.plength_observation2root = 1e-5;
    seqregion1.plength_observation2node = 2e-4;
    seqregion2.type = TYPE_O;
    seqregion2.plength_observation2root = -1;
    seqregion2.plength_observation2node = 0.0233;
    merge_RACGT_ORACGT(seqregion1, seqregion2, total_blength_1, total_blength_2, upper_plength, end_pos, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 3);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge1{0.0000043308850656905248,0.11650995967890213,0.067956494829030142,0.81552921460700212};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge1);
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    upper_plength = -1;
    total_blength_1 = 0;
    total_blength_2 = 1e-5;
    seqregion1.type = TYPE_R;
    seqregion1.plength_observation2root = 0;
    seqregion1.plength_observation2node = 0.043231;
    seqregion2.type = 3;
    seqregion2.plength_observation2root = 0.12;
    seqregion2.plength_observation2node = 3e-4;
    merge_RACGT_ORACGT(seqregion1, seqregion2, total_blength_1, total_blength_2, upper_plength, end_pos, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 4);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge2{3.7896075858290732E-7,0.0000019907504571924796,0.00022095280891100915,0.99977667747987319};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge2);
    // ----- Test 4 -----
    
    // ----- Test 5 -----
    upper_plength = -1;
    total_blength_1 = 1e-3;
    total_blength_2 = 1e-5;
    seqregion1.type = 0;
    seqregion1.plength_observation2root = 0;
    seqregion1.plength_observation2node = 0;
    seqregion2.type = TYPE_R;
    seqregion2.plength_observation2root = -1;
    seqregion2.plength_observation2node = 1e-7;
    merge_RACGT_ORACGT(seqregion1, seqregion2, total_blength_1, total_blength_2, upper_plength, end_pos, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 5);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge3{1,0,0,0};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge3);
    // ----- Test 5 -----
    
    // ----- Test 6 -----
    upper_plength = 0.311;
    total_blength_1 = 1e-5;
    total_blength_2 = 0.2311;
    seqregion1.type = TYPE_R;
    seqregion1.plength_observation2root = -1;
    seqregion1.plength_observation2node = -1;
    seqregion2.type = 2;
    seqregion2.plength_observation2root = -1;
    seqregion2.plength_observation2node = 0;
    merge_RACGT_ORACGT(seqregion1, seqregion2, total_blength_1, total_blength_2, upper_plength, end_pos, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 6);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge4{1.8321558474427985E-7,2.6859491508084263E-8,0.99999903717347061,7.5275145311095654E-7};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge4);
    // ----- Test 6 -----
    
    // ----- Test 7 -----
    upper_plength = -1;
    total_blength_1 = 0;
    total_blength_2 = 1e-4;
    seqregion1.type = 2;
    seqregion1.plength_observation2root = -1;
    seqregion1.plength_observation2node = -1;
    seqregion2.type = TYPE_R;
    seqregion2.plength_observation2root = -1;
    seqregion2.plength_observation2node = 0;
    merge_RACGT_ORACGT(seqregion1, seqregion2, total_blength_1, total_blength_2, upper_plength, end_pos, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 7);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge5{0,0,0.99999999999999988,0};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge5);
    // ----- Test 7 -----
    
    // ----- Test 8 -----
    upper_plength = 0;
    total_blength_1 = 0.321;
    total_blength_2 = 0;
    seqregion1.type = TYPE_R;
    seqregion1.plength_observation2root = 0;
    seqregion1.plength_observation2node = 1e-3;
    seqregion2.type = 1;
    seqregion2.plength_observation2root = 0;
    seqregion2.plength_observation2node = 0;
    merge_RACGT_ORACGT(seqregion1, seqregion2, total_blength_1, total_blength_2, upper_plength, end_pos, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 8);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge6{0,1,0,0};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge6);
    // ----- Test 8 -----
    
    // ----- Test 9 -----
    upper_plength = 0.3321;
    total_blength_1 = -1;
    total_blength_2 = -1;
    seqregion1.type = 1;
    seqregion1.plength_observation2root = 0;
    seqregion1.plength_observation2node = 1e-4;
    seqregion2.type = TYPE_O;
    seqregion2.plength_observation2root = 0;
    seqregion2.plength_observation2node = 1e-4;
    merge_RACGT_ORACGT(seqregion1, seqregion2, total_blength_1, total_blength_2, upper_plength, end_pos, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 9);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge7{3.8512118129099356E-7,0.50512531494378898,3.852643996135802E-7,0.49487391467062997};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge7);
    // ----- Test 9 -----
    
    // ----- Test 10 -----
    upper_plength = 0;
    total_blength_1 = 0;
    total_blength_2 = 1e-6;
    seqregion1.type = 3;
    seqregion1.plength_observation2root = -1;
    seqregion1.plength_observation2node = 0.11;
    seqregion2.type = TYPE_R;
    seqregion2.plength_observation2root = 0;
    seqregion2.plength_observation2node = 1e-5;
    merge_RACGT_ORACGT(seqregion1, seqregion2, total_blength_1, total_blength_2, upper_plength, end_pos, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 10);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge8{0,0,0,1};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge8);
    // ----- Test 10 -----
    
    // ----- Test 11 -----
    upper_plength = 0.123;
    total_blength_1 = -1;
    total_blength_2 = -1;
    seqregion1.type = TYPE_R;
    seqregion1.plength_observation2root = 0.432;
    seqregion1.plength_observation2node = 3e-6;
    seqregion2.type = TYPE_O;
    seqregion2.plength_observation2root = -1;
    seqregion2.plength_observation2node = -1;
    merge_RACGT_ORACGT(seqregion1, seqregion2, total_blength_1, total_blength_2, upper_plength, end_pos, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 11);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge9{0.0000046465940568638611,0.1249996558599968,0.000011794893992852322,0.87498390265195358};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge9);
    // ----- Test 11 -----
    
    // ----- Test 12 -----
    upper_plength = -1;
    total_blength_1 = 1e-5;
    total_blength_2 = 0;
    seqregion1.type = 0;
    seqregion1.plength_observation2root = 0;
    seqregion1.plength_observation2node = 0;
    seqregion2.type = TYPE_O;
    seqregion2.plength_observation2root = 0;
    seqregion2.plength_observation2node = 1e-7;
    merge_RACGT_ORACGT(seqregion1, seqregion2, total_blength_1, total_blength_2, upper_plength, end_pos, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 12);
    EXPECT_EQ(merged_regions.back().type, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    // ----- Test 12 -----
    
    // ----- Test 13 -----
    upper_plength = 1e-4;
    total_blength_1 = 1e-7;
    total_blength_2 = 0;
    seqregion1.type = TYPE_R;
    seqregion1.plength_observation2root = -1;
    seqregion1.plength_observation2node = -1;
    seqregion2.type = 3;
    seqregion2.plength_observation2root = 0;
    seqregion2.plength_observation2node = 0.0032;
    merge_RACGT_ORACGT(seqregion1, seqregion2, total_blength_1, total_blength_2, upper_plength, end_pos, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 13);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge10{0,0,0,1};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge10);
    // ----- Test 13 -----
    
    // ----- Test 14 -----
    upper_plength = 1e-5;
    total_blength_1 = 12e-6;
    total_blength_2 = -1;
    seqregion1.type = 1;
    seqregion1.plength_observation2root = -1;
    seqregion1.plength_observation2node = 2e-4;
    seqregion2.type = 3;
    seqregion2.plength_observation2root = 0;
    seqregion2.plength_observation2node = 0;
    merge_RACGT_ORACGT(seqregion1, seqregion2, total_blength_1, total_blength_2, upper_plength, end_pos, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 14);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge11{0,0,0,1};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge11);
    // ----- Test 14 -----
    
    // ----- Test 15 -----
    upper_plength = 3e-2;
    total_blength_1 = 2e-5;
    total_blength_2 = -1;
    seqregion1.type = 1;
    seqregion1.plength_observation2root = -1;
    seqregion1.plength_observation2node = -1;
    seqregion2.type = 0;
    seqregion2.plength_observation2root = -1;
    seqregion2.plength_observation2node = -1;
    merge_RACGT_ORACGT(seqregion1, seqregion2, total_blength_1, total_blength_2, upper_plength, end_pos, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 15);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge12{1,0,0,0};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge12);
    // ----- Test 15 -----
    
    // ----- Test 16 -----
    upper_plength = 4e-3;
    total_blength_1 = 0;
    total_blength_2 = 1e-4;
    seqregion1.type = TYPE_R;
    seqregion1.plength_observation2root = 1e-5;
    seqregion1.plength_observation2node = 2e-4;
    seqregion2.type = TYPE_O;
    seqregion2.plength_observation2root = 3e-2;
    seqregion2.plength_observation2node = 5e-5;
    merge_RACGT_ORACGT(seqregion1, seqregion2, total_blength_1, total_blength_2, upper_plength, end_pos, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 16);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge13{0.000010173897266335654,0.12178798877322355,0.026555202184610452,0.85164663514489958};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge13);
    // ----- Test 16 -----
    
    // ----- Test 17 -----
    upper_plength = 3e-3;
    total_blength_1 = 0;
    total_blength_2 = 231e-6;
    seqregion1.type = 2;
    seqregion1.plength_observation2root = -1;
    seqregion1.plength_observation2node = 3e-7;
    seqregion2.type = 1;
    seqregion2.plength_observation2root = -1;
    seqregion2.plength_observation2node = 0;
    merge_RACGT_ORACGT(seqregion1, seqregion2, total_blength_1, total_blength_2, upper_plength, end_pos, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 17);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge14{0,0,1,0};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge14);
    // ----- Test 17 -----
    
    // ----- Test 18 -----
    upper_plength = -1;
    total_blength_1 = 2e-4;
    total_blength_2 = 3e-8;
    seqregion1.type = 0;
    seqregion1.plength_observation2root = -1;
    seqregion1.plength_observation2node = -1;
    seqregion2.type = TYPE_O;
    seqregion2.plength_observation2root = 2e-5;
    seqregion2.plength_observation2node = 121e-6;
    merge_RACGT_ORACGT(seqregion1, seqregion2, total_blength_1, total_blength_2, upper_plength, end_pos, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 18);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge15{0.28592397557801535,0.30602769829658166,0.000011398356098642252,0.40803692776930434};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge15);
    // ----- Test 18 -----
    
    // ----- Test 19 -----
    upper_plength = 4e-2;
    total_blength_1 = 1e-5;
    total_blength_2 = 0;
    seqregion1.type = 1;
    seqregion1.plength_observation2root = 0;
    seqregion1.plength_observation2node = 3e-6;
    seqregion2.type = TYPE_R;
    seqregion2.plength_observation2root = -1;
    seqregion2.plength_observation2node = 0;
    merge_RACGT_ORACGT(seqregion1, seqregion2, total_blength_1, total_blength_2, upper_plength, end_pos, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 19);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge16{0,0,1,0};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge16);
    // ----- Test 19 -----
    
    // ----- Test 20 -----
    upper_plength = -1;
    total_blength_1 = -1;
    total_blength_2 = 4e-3;
    seqregion1.type = 0;
    seqregion1.plength_observation2root = -1;
    seqregion1.plength_observation2node = -1;
    seqregion2.type = TYPE_O;
    seqregion2.plength_observation2root = -1;
    seqregion2.plength_observation2node = -1;
    merge_RACGT_ORACGT(seqregion1, seqregion2, total_blength_1, total_blength_2, upper_plength, end_pos, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 20);
    EXPECT_EQ(merged_regions.back().type, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    // ----- Test 20 -----
    
    // ----- Test 21 -----
    upper_plength = 54e-8;
    total_blength_1 = -1;
    total_blength_2 = 1e-5;
    seqregion1.type = 3;
    seqregion1.plength_observation2root = 1e-3;
    seqregion1.plength_observation2node = 2e-7;
    seqregion2.type = TYPE_R;
    seqregion2.plength_observation2root = 0;
    seqregion2.plength_observation2node = 0.0142;
    merge_RACGT_ORACGT(seqregion1, seqregion2, total_blength_1, total_blength_2, upper_plength, end_pos, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 21);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge17{3.7529836729309572E-7,0.0000019715097300387874,0.99011618892809094,0.0098814642638117549};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge17);
    // ----- Test 21 -----
    
    // ----- Test 22 -----
    upper_plength = -1;
    total_blength_1 = 21e-3;
    total_blength_2 = 1e-8;
    seqregion1.type = TYPE_R;
    seqregion1.plength_observation2root = 1e-5;
    seqregion1.plength_observation2node = 2e-7;
    seqregion2.type = 1;
    seqregion2.plength_observation2root = -1;
    seqregion2.plength_observation2node = -1;
    merge_RACGT_ORACGT(seqregion1, seqregion2, total_blength_1, total_blength_2, upper_plength, end_pos, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 22);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge18{1.9880125157816739E-9,0.99902052127670382,0.00097942098640511682,5.5748878587250403E-8};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge18);
    // ----- Test 22 -----
    
    // ----- Test 23 -----
    upper_plength = 13e-8;
    total_blength_1 = 0;
    total_blength_2 = 156e-9;
    seqregion1.type = 0;
    seqregion1.plength_observation2root = -1;
    seqregion1.plength_observation2node = 0.431;
    seqregion2.type = TYPE_R;
    seqregion2.plength_observation2root = 0;
    seqregion2.plength_observation2node = 0.0042;
    merge_RACGT_ORACGT(seqregion1, seqregion2, total_blength_1, total_blength_2, upper_plength, end_pos, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 23);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge19{1,0,0,0};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge19);
    // ----- Test 23 -----
    
    // ----- Test 24 -----
    upper_plength = -1;
    total_blength_1 = -1;
    total_blength_2 = 41e-6;
    seqregion1.type = TYPE_R;
    seqregion1.plength_observation2root = 5e-5;
    seqregion1.plength_observation2node = 134e-8;
    seqregion2.type = TYPE_O;
    seqregion2.plength_observation2root = 342e-7;
    seqregion2.plength_observation2node = 13e-4;
    merge_RACGT_ORACGT(seqregion1, seqregion2, total_blength_1, total_blength_2, upper_plength, end_pos, threshold_prob, model, aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 24);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge20{0.0000033509369075667414,0.059613997536883151,0.52309211091516972,0.41729054061103954};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge20);
    // ----- Test 24 -----
}

/*
    Generate testing data (seqregions1, seqregions2) for testing mergeUpperLower() and mergeTwoLowers()
 */
void genTestData3_0(SeqRegions& seqregions1, SeqRegions& seqregions2, int test_case)
{
    seqregions1.clear();
    seqregions2.clear();
    
    switch (test_case)
    {
        case 1:
        {
            seqregions1.emplace_back(250,239,-1,-1);
            seqregions1.emplace_back(1,240,-1,-1);
            seqregions1.emplace_back(250,3035,-1,-1);
            seqregions1.emplace_back(1,3036,-1,-1);
            seqregions1.emplace_back(250,14406,-1,-1);
            seqregions1.emplace_back(1,14407,-1,-1);
            seqregions1.emplace_back(250,23401,-1,-1);
            seqregions1.emplace_back(0,23402,-1,-1);
            seqregions1.emplace_back(250,28142,-1,-1);
            seqregions1.emplace_back(1,28143,-1,-1);
            seqregions1.emplace_back(250,29890,-1,-1);
            seqregions2.emplace_back(250,239,-1,-1);
            seqregions2.emplace_back(1,240,-1,-1);
            seqregions2.emplace_back(250,3035,-1,-1);
            seqregions2.emplace_back(1,3036,-1,-1);
            seqregions2.emplace_back(250,11188,-1,-1);
            seqregions2.emplace_back(3,11189,-1,-1);
            seqregions2.emplace_back(250,14406,-1,-1);
            seqregions2.emplace_back(1,14407,-1,-1);
            seqregions2.emplace_back(250,23401,-1,-1);
            seqregions2.emplace_back(0,23402,-1,-1);
            seqregions2.emplace_back(250,28142,-1,-1);
            seqregions2.emplace_back(1,28143,-1,-1);
            seqregions2.emplace_back(250,29890,-1,-1);
            break;
        }
        case 2:
        {
            seqregions1.emplace_back(250,239,-1,-1);
            seqregions1.emplace_back(1,240,-1,-1);
            seqregions1.emplace_back(250,3035,-1,-1);
            seqregions1.emplace_back(1,3036,-1,-1);
            seqregions1.emplace_back(250,14406,-1,-1);
            seqregions1.emplace_back(1,14407,-1,-1);
            seqregions1.emplace_back(250,23401,-1,-1);
            seqregions1.emplace_back(0,23402,-1,-1);
            seqregions1.emplace_back(250,28142,-1,-1);
            seqregions1.emplace_back(1,28143,-1,-1);
            seqregions1.emplace_back(250,29890,-1,-1);
            seqregions2.emplace_back(250,239,-1,-1);
            seqregions2.emplace_back(1,240,-1,-1);
            seqregions2.emplace_back(250,3035,-1,-1);
            seqregions2.emplace_back(1,3036,-1,-1);
            seqregions2.emplace_back(250,14406,-1,-1);
            seqregions2.emplace_back(1,14407,-1,-1);
            seqregions2.emplace_back(250,28165,-1,-1);
            seqregions2.emplace_back(0,28166,-1,-1);
            seqregions2.emplace_back(250,28873,-1,-1);
            seqregions2.emplace_back(0,28874,-1,-1);
            seqregions2.emplace_back(250,28876,-1,-1);
            seqregions2.emplace_back(0,28877,-1,-1);
            seqregions2.emplace_back(250,29740,-1,-1);
            seqregions2.emplace_back(0,29741,-1,-1);
            seqregions2.emplace_back(250,29890,-1,-1);
            break;
        }
        case 3:
        {
            seqregions1.emplace_back(250,239,-1,-1);
            seqregions1.emplace_back(1,240,-1,-1);
            seqregions1.emplace_back(250,3035,-1,-1);
            seqregions1.emplace_back(1,3036,-1,-1);
            seqregions1.emplace_back(250,14406,-1,-1);
            seqregions1.emplace_back(1,14407,-1,-1);
            seqregions1.emplace_back(250,23401,-1,-1);
            seqregions1.emplace_back(0,23402,-1,-1);
            seqregions1.emplace_back(250,28142,-1,-1);
            seqregions1.emplace_back(1,28143,-1,-1);
            seqregions1.emplace_back(250,29890,-1,-1);
            seqregions2.emplace_back(250,3035,-1,-1);
            seqregions2.emplace_back(1,3036,-1,-1);
            seqregions2.emplace_back(250,8780,-1,-1);
            seqregions2.emplace_back(3,8781,-1,-1);
            seqregions2.emplace_back(250,14406,-1,-1);
            seqregions2.emplace_back(1,14407,-1,-1);
            seqregions2.emplace_back(250,18754,-1,-1);
            seqregions2.emplace_back(3,18755,-1,-1);
            seqregions2.emplace_back(250,22466,-1,-1);
            seqregions2.emplace_back(3,22467,-1,-1);
            seqregions2.emplace_back(250,28876,-1,-1);
            seqregions2.emplace_back(0,28877,-1,-1);
            seqregions2.emplace_back(250,29708,-1,-1);
            seqregions2.emplace_back(0,29709,-1,-1);
            seqregions2.emplace_back(250,29740,-1,-1);
            seqregions2.emplace_back(0,29741,-1,-1);
            seqregions2.emplace_back(250,29890,-1,-1);
            break;
        }
        case 4:
        {
            seqregions1.emplace_back(250,239,-1,-1);
            seqregions1.emplace_back(1,240,-1,-1);
            seqregions1.emplace_back(250,3035,-1,-1);
            seqregions1.emplace_back(1,3036,-1,-1);
            seqregions1.emplace_back(250,8780,-1,-1);
            seqregions1.emplace_back(3,8781,-1,-1);
            seqregions1.emplace_back(250,14406,-1,-1);
            seqregions1.emplace_back(1,14407,-1,-1);
            seqregions1.emplace_back(250,23401,-1,-1);
            seqregions1.emplace_back(0,23402,-1,-1);
            seqregions1.emplace_back(250,28142,-1,-1);
            seqregions1.emplace_back(1,28143,-1,-1);
            seqregions1.emplace_back(250,29890,-1,-1);
            seqregions2.emplace_back(250,239,-1,-1);
            seqregions2.emplace_back(1,240,-1,-1);
            seqregions2.emplace_back(250,3035,-1,-1);
            seqregions2.emplace_back(1,3036,-1,-1);
            seqregions2.emplace_back(250,8780,-1,-1);
            seqregions2.emplace_back(3,8781,-1,-1);
            seqregions2.emplace_back(250,14406,-1,-1);
            seqregions2.emplace_back(1,14407,-1,-1);
            seqregions2.emplace_back(250,17856,-1,-1);
            seqregions2.emplace_back(2,17857,-1,-1);
            seqregions2.emplace_back(250,23401,-1,-1);
            seqregions2.emplace_back(0,23402,-1,-1);
            seqregions2.emplace_back(250,28142,-1,-1);
            seqregions2.emplace_back(1,28143,-1,-1);
            seqregions2.emplace_back(250,29890,-1,-1);
            break;
        }
        case 5:
        {
            seqregions1.emplace_back(250,239,-1,-1);
            seqregions1.emplace_back(1,240,-1,-1);
            seqregions1.emplace_back(250,3035,-1,-1);
            seqregions1.emplace_back(1,3036,-1,-1);
            seqregions1.emplace_back(250,8780,-1,-1);
            seqregions1.emplace_back(3,8781,-1,-1);
            seqregions1.emplace_back(250,14406,-1,-1);
            seqregions1.emplace_back(1,14407,-1,-1);
            seqregions1.emplace_back(250,23401,-1,-1);
            seqregions1.emplace_back(0,23402,-1,-1);
            seqregions1.emplace_back(250,26086,-1,-1);
            seqregions1.emplace_back(3,26087,-1,-1);
            seqregions1.emplace_back(250,28142,-1,-1);
            seqregions1.emplace_back(1,28143,-1,-1);
            seqregions1.emplace_back(250,29890,-1,-1);
            seqregions2.emplace_back(252,22,-1,-1);
            seqregions2.emplace_back(250,24,3.3454886086112878360986772063867533688608091324568e-05,-1);
            seqregions2.emplace_back(250,239,-1,-1);
            seqregions2.emplace_back(1,240,-1,-1);
            seqregions2.emplace_back(250,3035,-1,-1);
            auto new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value{2.4873891076105684410204367675674072546421200513578e-10,0.99998513013917567082700088576530106365680694580078,2.4873891076105684410204367675674072546421200513578e-10,1.4869363346385217151198084029051216248262790031731e-05};
            (*new_lh) = new_lh_value;
            seqregions2.emplace_back(251,3036,0,0,std::move(new_lh));
            seqregions2.emplace_back(250,8780,-1,-1);
            seqregions2.emplace_back(3,8781,-1,-1);
            seqregions2.emplace_back(250,14406,-1,-1);
            seqregions2.emplace_back(1,14407,-1,-1);
            seqregions2.emplace_back(250,26086,-1,-1);
            seqregions2.emplace_back(3,26087,-1,-1);
            seqregions2.emplace_back(250,28142,-1,-1);
            seqregions2.emplace_back(1,28143,-1,-1);
            seqregions2.emplace_back(250,29846,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value1{5.5761563391267229329063825904455597992637194693089e-05,6.2185305443590081052239281053303776580776229820913e-10,6.2185305443590081052239281053303776580776229820913e-10,0.99994423719290270735626791065442375838756561279297};
            (*new_lh) = new_lh_value1;
            seqregions2.emplace_back(251,29847,0,0,std::move(new_lh));
            seqregions2.emplace_back(250,29857,-1,-1);
            seqregions2.emplace_back(250,29858,8.3637215215282185738071563108064765401650220155716e-05,-1);
            seqregions2.emplace_back(252,29890,-1,-1);
            break;
        }
        case 6:
        {
            seqregions1.emplace_back(250,239,-1,-1);
            seqregions1.emplace_back(1,240,-1,-1);
            seqregions1.emplace_back(250,3035,-1,-1);
            seqregions1.emplace_back(1,3036,-1,-1);
            seqregions1.emplace_back(250,8780,-1,-1);
            seqregions1.emplace_back(3,8781,-1,-1);
            seqregions1.emplace_back(250,11228,-1,-1);
            seqregions1.emplace_back(3,11229,-1,-1);
            seqregions1.emplace_back(250,14406,-1,-1);
            seqregions1.emplace_back(1,14407,-1,-1);
            seqregions1.emplace_back(250,20913,-1,-1);
            auto new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value{1.1151877415789591228433182135137968771232408471406e-05,1.2436575683940661210420103588707597258578019250308e-10,0.99998884787385255989988763758447021245956420898438,1.2436575683940661210420103588707597258578019250308e-10};
            (*new_lh) = new_lh_value;
            seqregions1.emplace_back(251,20914,0,0,std::move(new_lh));
            seqregions1.emplace_back(250,23401,-1,-1);
            seqregions1.emplace_back(0,23402,-1,-1);
            seqregions1.emplace_back(250,27598,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value1{1.2436575683940661210420103588707597258578019250308e-10,0.99998884787385255989988763758447021245956420898438,1.2436575683940661210420103588707597258578019250308e-10,1.1151877415789591228433182135137968771232408471406e-05};
            (*new_lh) = new_lh_value1;
            seqregions1.emplace_back(251,27599,0,0,std::move(new_lh));
            seqregions1.emplace_back(250,28142,-1,-1);
            seqregions1.emplace_back(1,28143,-1,-1);
            seqregions1.emplace_back(250,28165,-1,-1);
            seqregions1.emplace_back(0,28166,-1,-1);
            seqregions1.emplace_back(250,28876,-1,-1);
            seqregions1.emplace_back(0,28877,-1,-1);
            seqregions1.emplace_back(250,29740,-1,-1);
            seqregions1.emplace_back(0,29741,-1,-1);
            seqregions1.emplace_back(250,29890,-1,-1);
            seqregions2.emplace_back(250,239,-1,-1);
            seqregions2.emplace_back(1,240,-1,-1);
            seqregions2.emplace_back(250,3035,-1,-1);
            seqregions2.emplace_back(1,3036,-1,-1);
            seqregions2.emplace_back(250,8780,-1,-1);
            seqregions2.emplace_back(3,8781,-1,-1);
            seqregions2.emplace_back(250,11228,-1,-1);
            seqregions2.emplace_back(3,11229,-1,-1);
            seqregions2.emplace_back(250,14406,-1,-1);
            seqregions2.emplace_back(1,14407,-1,-1);
            seqregions2.emplace_back(250,23401,-1,-1);
            seqregions2.emplace_back(0,23402,-1,-1);
            seqregions2.emplace_back(250,28142,-1,-1);
            seqregions2.emplace_back(1,28143,-1,-1);
            seqregions2.emplace_back(250,28165,-1,-1);
            seqregions2.emplace_back(0,28166,-1,-1);
            seqregions2.emplace_back(250,28876,-1,-1);
            seqregions2.emplace_back(0,28877,-1,-1);
            seqregions2.emplace_back(250,29736,-1,-1);
            seqregions2.emplace_back(3,29737,-1,-1);
            seqregions2.emplace_back(250,29740,-1,-1);
            seqregions2.emplace_back(0,29741,-1,-1);
            seqregions2.emplace_back(250,29890,-1,-1);
            break;
        }
        case 7:
        {
            seqregions1.emplace_back(250,11,6.6909772172225756721973544127735067377216182649136e-05,-1);
            seqregions1.emplace_back(250,239,-1,-1);
            seqregions1.emplace_back(1,240,-1,-1);
            seqregions1.emplace_back(250,3035,-1,-1);
            seqregions1.emplace_back(1,3036,-1,-1);
            seqregions1.emplace_back(250,3177,-1,-1);
            seqregions1.emplace_back(0,3178,-1,-1);
            seqregions1.emplace_back(250,6980,-1,-1);
            seqregions1.emplace_back(3,6981,-1,-1);
            seqregions1.emplace_back(250,8780,-1,-1);
            seqregions1.emplace_back(3,8781,-1,-1);
            seqregions1.emplace_back(250,14406,-1,-1);
            seqregions1.emplace_back(1,14407,-1,-1);
            seqregions1.emplace_back(250,23401,-1,-1);
            seqregions1.emplace_back(0,23402,-1,-1);
            seqregions1.emplace_back(250,28142,-1,-1);
            seqregions1.emplace_back(1,28143,-1,-1);
            seqregions1.emplace_back(250,29866,-1,-1);
            seqregions1.emplace_back(250,29867,6.6909772172225756721973544127735067377216182649136e-05,-1);
            seqregions1.emplace_back(250,29890,-1,-1);
            seqregions2.emplace_back(252,28,-1,-1);
            seqregions2.emplace_back(250,239,-1,-1);
            seqregions2.emplace_back(1,240,-1,-1);
            seqregions2.emplace_back(250,3035,-1,-1);
            seqregions2.emplace_back(1,3036,-1,-1);
            seqregions2.emplace_back(250,3177,-1,-1);
            seqregions2.emplace_back(0,3178,-1,-1);
            seqregions2.emplace_back(250,6980,-1,-1);
            seqregions2.emplace_back(3,6981,-1,-1);
            seqregions2.emplace_back(250,8780,-1,-1);
            seqregions2.emplace_back(3,8781,-1,-1);
            seqregions2.emplace_back(250,14406,-1,-1);
            seqregions2.emplace_back(1,14407,-1,-1);
            seqregions2.emplace_back(250,23401,-1,-1);
            seqregions2.emplace_back(0,23402,-1,-1);
            seqregions2.emplace_back(250,28142,-1,-1);
            seqregions2.emplace_back(1,28143,-1,-1);
            seqregions2.emplace_back(250,29865,-1,-1);
            seqregions2.emplace_back(252,29867,-1,-1);
            seqregions2.emplace_back(250,29890,-1,-1);
            break;
        }
        case 8:
        {
            seqregions1.emplace_back(250,239,-1,-1);
            seqregions1.emplace_back(1,240,-1,-1);
            seqregions1.emplace_back(250,1267,-1,-1);
            auto new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value{0.80000178432426516383912940000300295650959014892578,8.9216213262436184458807272856795123061601771041751e-06,0.19998037243308225408000566858390811830759048461914,8.9216213262436184458807272856795123061601771041751e-06};
            (*new_lh) = new_lh_value;
            seqregions1.emplace_back(251,1268,0,0,std::move(new_lh));
            seqregions1.emplace_back(250,3035,-1,-1);
            seqregions1.emplace_back(1,3036,-1,-1);
            seqregions1.emplace_back(250,8780,-1,-1);
            seqregions1.emplace_back(3,8781,-1,-1);
            seqregions1.emplace_back(250,14371,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value1{0.80000178432426516383912940000300295650959014892578,0.19998037243308225408000566858390811830759048461914,8.9216213262436184458807272856795123061601771041751e-06,8.9216213262436184458807272856795123061601771041751e-06};
            (*new_lh) = new_lh_value1;
            seqregions1.emplace_back(251,14372,0,0,std::move(new_lh));
            seqregions1.emplace_back(250,14406,-1,-1);
            seqregions1.emplace_back(1,14407,-1,-1);
            seqregions1.emplace_back(250,21362,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value2{8.9216213262436184458807272856795123061601771041751e-06,0.80000178432426516383912940000300295650959014892578,8.9216213262436184458807272856795123061601771041751e-06,0.19998037243308225408000566858390811830759048461914};
            (*new_lh) = new_lh_value2;
            seqregions1.emplace_back(251,21363,0,0,std::move(new_lh));
            seqregions1.emplace_back(250,23401,-1,-1);
            seqregions1.emplace_back(0,23402,-1,-1);
            seqregions1.emplace_back(250,28142,-1,-1);
            seqregions1.emplace_back(1,28143,-1,-1);
            seqregions1.emplace_back(250,28165,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value3{0.19998037243308225408000566858390811830759048461914,8.9216213262436184458807272856795123061601771041751e-06,0.80000178432426516383912940000300295650959014892578,8.9216213262436184458807272856795123061601771041751e-06};
            (*new_lh) = new_lh_value3;
            seqregions1.emplace_back(251,28166,0,0,std::move(new_lh));
            seqregions1.emplace_back(250,28876,-1,-1);
            seqregions1.emplace_back(0,28877,-1,-1);
            seqregions1.emplace_back(250,29740,-1,-1);
            seqregions1.emplace_back(0,29741,-1,-1);
            seqregions1.emplace_back(250,29890,-1,-1);
            seqregions2.emplace_back(250,239,-1,-1);
            seqregions2.emplace_back(1,240,-1,-1);
            seqregions2.emplace_back(250,3035,-1,-1);
            seqregions2.emplace_back(1,3036,-1,-1);
            seqregions2.emplace_back(250,8780,-1,-1);
            seqregions2.emplace_back(3,8781,-1,-1);
            seqregions2.emplace_back(250,11228,-1,-1);
            seqregions2.emplace_back(3,11229,-1,-1);
            seqregions2.emplace_back(250,14406,-1,-1);
            seqregions2.emplace_back(1,14407,-1,-1);
            seqregions2.emplace_back(250,17857,-1,-1);
            seqregions2.emplace_back(1,17858,-1,-1);
            seqregions2.emplace_back(250,23401,-1,-1);
            seqregions2.emplace_back(0,23402,-1,-1);
            seqregions2.emplace_back(250,28142,-1,-1);
            seqregions2.emplace_back(1,28143,-1,-1);
            seqregions2.emplace_back(250,28165,-1,-1);
            seqregions2.emplace_back(0,28166,-1,-1);
            seqregions2.emplace_back(250,28876,-1,-1);
            seqregions2.emplace_back(0,28877,-1,-1);
            seqregions2.emplace_back(250,29740,-1,-1);
            seqregions2.emplace_back(0,29741,-1,-1);
            seqregions2.emplace_back(250,29890,-1,-1);
            break;
        }
        case 9:
        {
            seqregions1.emplace_back(250,239,-1,-1);
            seqregions1.emplace_back(1,240,-1,-1);
            seqregions1.emplace_back(250,1267,-1,-1);
            auto new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value{0.80000178432426516383912940000300295650959014892578,8.9216213262436184458807272856795123061601771041751e-06,0.19998037243308225408000566858390811830759048461914,8.9216213262436184458807272856795123061601771041751e-06};
            (*new_lh) = new_lh_value;
            seqregions1.emplace_back(251,1268,0,0,std::move(new_lh));
            seqregions1.emplace_back(250,3035,-1,-1);
            seqregions1.emplace_back(1,3036,-1,-1);
            seqregions1.emplace_back(250,8780,-1,-1);
            seqregions1.emplace_back(3,8781,-1,-1);
            seqregions1.emplace_back(250,14371,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value1{0.80000178432426516383912940000300295650959014892578,0.19998037243308225408000566858390811830759048461914,8.9216213262436184458807272856795123061601771041751e-06,8.9216213262436184458807272856795123061601771041751e-06};
            (*new_lh) = new_lh_value1;
            seqregions1.emplace_back(251,14372,0,0,std::move(new_lh));
            seqregions1.emplace_back(250,14406,-1,-1);
            seqregions1.emplace_back(1,14407,-1,-1);
            seqregions1.emplace_back(250,21362,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value3{8.9216213262436184458807272856795123061601771041751e-06,0.80000178432426516383912940000300295650959014892578,8.9216213262436184458807272856795123061601771041751e-06,0.19998037243308225408000566858390811830759048461914};
            (*new_lh) = new_lh_value3;
            seqregions1.emplace_back(251,21363,0,0,std::move(new_lh));
            seqregions1.emplace_back(250,23401,-1,-1);
            seqregions1.emplace_back(0,23402,-1,-1);
            seqregions1.emplace_back(250,28142,-1,-1);
            seqregions1.emplace_back(1,28143,-1,-1);
            seqregions1.emplace_back(250,28165,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value4{0.19998037243308225408000566858390811830759048461914,8.9216213262436184458807272856795123061601771041751e-06,0.80000178432426516383912940000300295650959014892578,8.9216213262436184458807272856795123061601771041751e-06};
            (*new_lh) = new_lh_value4;
            seqregions1.emplace_back(251,28166,0,0,std::move(new_lh));
            seqregions1.emplace_back(250,28876,-1,-1);
            seqregions1.emplace_back(0,28877,-1,-1);
            seqregions1.emplace_back(250,29740,-1,-1);
            seqregions1.emplace_back(0,29741,-1,-1);
            seqregions1.emplace_back(250,29890,-1,-1);
            seqregions2.emplace_back(250,239,-1,-1);
            seqregions2.emplace_back(1,240,-1,-1);
            seqregions2.emplace_back(250,3035,-1,-1);
            seqregions2.emplace_back(1,3036,-1,-1);
            seqregions2.emplace_back(250,8780,-1,-1);
            seqregions2.emplace_back(3,8781,-1,-1);
            seqregions2.emplace_back(250,14406,-1,-1);
            seqregions2.emplace_back(1,14407,-1,-1);
            seqregions2.emplace_back(250,16573,-1,-1);
            seqregions2.emplace_back(3,16574,-1,-1);
            seqregions2.emplace_back(250,17857,-1,-1);
            seqregions2.emplace_back(1,17858,-1,-1);
            seqregions2.emplace_back(250,23401,-1,-1);
            seqregions2.emplace_back(0,23402,-1,-1);
            seqregions2.emplace_back(250,26311,-1,-1);
            seqregions2.emplace_back(3,26312,-1,-1);
            seqregions2.emplace_back(250,28142,-1,-1);
            seqregions2.emplace_back(1,28143,-1,-1);
            seqregions2.emplace_back(250,28165,-1,-1);
            seqregions2.emplace_back(0,28166,-1,-1);
            seqregions2.emplace_back(250,28876,-1,-1);
            seqregions2.emplace_back(0,28877,-1,-1);
            seqregions2.emplace_back(250,29740,-1,-1);
            seqregions2.emplace_back(0,29741,-1,-1);
            seqregions2.emplace_back(250,29890,-1,-1);
            break;
        }
        case 10:
        {
            seqregions1.emplace_back(250,12,3.3454886086112878360986772063867533688608091324568e-05,-1);
            seqregions1.emplace_back(250,239,-1,-1);
            seqregions1.emplace_back(1,240,-1,-1);
            seqregions1.emplace_back(250,3035,-1,-1);
            seqregions1.emplace_back(1,3036,-1,-1);
            seqregions1.emplace_back(250,8780,-1,-1);
            seqregions1.emplace_back(3,8781,-1,-1);
            seqregions1.emplace_back(250,14406,-1,-1);
            seqregions1.emplace_back(1,14407,-1,-1);
            seqregions1.emplace_back(250,17745,-1,-1);
            seqregions1.emplace_back(3,17746,-1,-1);
            seqregions1.emplace_back(250,17856,-1,-1);
            seqregions1.emplace_back(2,17857,-1,-1);
            seqregions1.emplace_back(250,18058,-1,-1);
            seqregions1.emplace_back(3,18059,-1,-1);
            seqregions1.emplace_back(250,23401,-1,-1);
            seqregions1.emplace_back(0,23402,-1,-1);
            seqregions1.emplace_back(250,28142,-1,-1);
            seqregions1.emplace_back(1,28143,-1,-1);
            seqregions1.emplace_back(250,29874,-1,-1);
            seqregions1.emplace_back(250,29881,3.3454886086112878360986772063867533688608091324568e-05,-1);
            seqregions1.emplace_back(250,29890,6.6909772172225756721973544127735067377216182649136e-05,-1);
            seqregions2.emplace_back(250,239,-1,-1);
            seqregions2.emplace_back(1,240,-1,-1);
            seqregions2.emplace_back(250,3035,-1,-1);
            seqregions2.emplace_back(1,3036,-1,-1);
            seqregions2.emplace_back(250,8780,-1,-1);
            seqregions2.emplace_back(3,8781,-1,-1);
            seqregions2.emplace_back(250,14406,-1,-1);
            seqregions2.emplace_back(1,14407,-1,-1);
            seqregions2.emplace_back(250,17745,-1,-1);
            seqregions2.emplace_back(3,17746,-1,-1);
            seqregions2.emplace_back(250,17856,-1,-1);
            seqregions2.emplace_back(2,17857,-1,-1);
            seqregions2.emplace_back(250,18058,-1,-1);
            seqregions2.emplace_back(3,18059,-1,-1);
            seqregions2.emplace_back(250,23401,-1,-1);
            seqregions2.emplace_back(0,23402,-1,-1);
            seqregions2.emplace_back(250,28142,-1,-1);
            seqregions2.emplace_back(1,28143,-1,-1);
            seqregions2.emplace_back(250,29890,-1,-1);
            break;
        }
        default:
            break;
    }
}

/*
    Generate output data (merged_seqregions) for testing mergeUpperLower()
 */
void genOutputData3_1(SeqRegions& merged_regions, int test_case)
{
    merged_regions.clear();
    
    switch (test_case)
    {
        case 1:
        {
            merged_regions.emplace_back(250,239,-1,-1);
            merged_regions.emplace_back(1,240,-1,-1);
            merged_regions.emplace_back(250,3035,-1,-1);
            merged_regions.emplace_back(1,3036,-1,-1);
            merged_regions.emplace_back(250,11188,-1,-1);
            merged_regions.emplace_back(3,11189,-1,-1);
            merged_regions.emplace_back(250,14406,-1,-1);
            merged_regions.emplace_back(1,14407,-1,-1);
            merged_regions.emplace_back(250,23401,-1,-1);
            merged_regions.emplace_back(0,23402,-1,-1);
            merged_regions.emplace_back(250,28142,-1,-1);
            merged_regions.emplace_back(1,28143,-1,-1);
            merged_regions.emplace_back(250,29890,-1,-1);
            break;
        }
        case 2:
        {
            merged_regions.emplace_back(250,239,-1,-1);
            merged_regions.emplace_back(1,240,-1,-1);
            merged_regions.emplace_back(250,3035,-1,-1);
            merged_regions.emplace_back(1,3036,-1,-1);
            merged_regions.emplace_back(250,14406,-1,-1);
            merged_regions.emplace_back(1,14407,-1,-1);
            merged_regions.emplace_back(250,23401,-1,-1);
            auto new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value{0.12497560486006949187487435892762732692062854766846,9.7580559722090576469559833339140197949745925143361e-06,0.87500487902798607109389195102266967296600341796875,9.7580559722090576469559833339140197949745925143361e-06};
            (*new_lh) = new_lh_value;
            merged_regions.emplace_back(251,23402,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,28142,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value1{9.7580559722090576469559833339140197949745925143361e-06,0.12497560486006949187487435892762732692062854766846,9.7580559722090576469559833339140197949745925143361e-06,0.87500487902798607109389195102266967296600341796875};
            (*new_lh) = new_lh_value1;
            merged_regions.emplace_back(251,28143,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,28165,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value2{0.87500487902798607109389195102266967296600341796875,9.7580559722090576469559833339140197949745925143361e-06,0.12497560486006949187487435892762732692062854766846,9.7580559722090576469559833339140197949745925143361e-06};
            (*new_lh) = new_lh_value2;
            merged_regions.emplace_back(251,28166,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,28873,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            (*new_lh) = new_lh_value2;
            merged_regions.emplace_back(251,28874,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,28876,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            (*new_lh) = new_lh_value2;
            merged_regions.emplace_back(251,28877,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,29740,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            (*new_lh) = new_lh_value2;
            merged_regions.emplace_back(251,29741,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,29890,-1,-1);
            break;
        }
        case 3:
        {
            merged_regions.emplace_back(250,239,-1,-1);
            auto new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value3{2.6766774402933635499746492514283602304203668609262e-05,0.39996252651583585890904259940725751221179962158203,2.6766774402933635499746492514283602304203668609262e-05,0.59998393993535814594508792652050033211708068847656};
            (*new_lh) = new_lh_value3;
            merged_regions.emplace_back(251,240,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,3035,-1,-1);
            merged_regions.emplace_back(1,3036,-1,-1);
            merged_regions.emplace_back(250,8780,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            (*new_lh) = new_lh_value3;
            merged_regions.emplace_back(251,8781,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,14406,-1,-1);
            merged_regions.emplace_back(1,14407,-1,-1);
            merged_regions.emplace_back(250,18754,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value4{2.6766774402933635499746492514283602304203668609262e-05,2.6766774402933635499746492514283602304203668609262e-05,0.39996252651583585890904259940725751221179962158203,0.59998393993535814594508792652050033211708068847656};
            (*new_lh) = new_lh_value4;
            merged_regions.emplace_back(251,18755,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,22466,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            (*new_lh) = new_lh_value4;
            merged_regions.emplace_back(251,22467,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,23401,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value5{0.39996252651583585890904259940725751221179962158203,2.6766774402933635499746492514283602304203668609262e-05,0.59998393993535814594508792652050033211708068847656,2.6766774402933635499746492514283602304203668609262e-05};
            (*new_lh) = new_lh_value5;
            merged_regions.emplace_back(251,23402,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,28142,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            (*new_lh) = new_lh_value3;
            merged_regions.emplace_back(251,28143,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,28876,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value6{0.59998393993535814594508792652050033211708068847656,2.6766774402933635499746492514283602304203668609262e-05,0.39996252651583585890904259940725751221179962158203,2.6766774402933635499746492514283602304203668609262e-05};
            (*new_lh) = new_lh_value6;
            merged_regions.emplace_back(251,28877,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,29708,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value7{0.59998393993535814594508792652050033211708068847656,2.6766774402933635499746492514283602304203668609262e-05,2.6766774402933635499746492514283602304203668609262e-05,0.39996252651583585890904259940725751221179962158203};
            (*new_lh) = new_lh_value7;
            merged_regions.emplace_back(251,29709,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,29740,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            (*new_lh) = new_lh_value6;
            merged_regions.emplace_back(251,29741,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,29890,-1,-1);
            break;
        }
        case 4:
        {
            merged_regions.emplace_back(250,239,-1,-1);
            merged_regions.emplace_back(1,240,-1,-1);
            merged_regions.emplace_back(250,3035,-1,-1);
            merged_regions.emplace_back(1,3036,-1,-1);
            merged_regions.emplace_back(250,8780,-1,-1);
            merged_regions.emplace_back(3,8781,-1,-1);
            merged_regions.emplace_back(250,14406,-1,-1);
            merged_regions.emplace_back(1,14407,-1,-1);
            merged_regions.emplace_back(250,17856,-1,-1);
            merged_regions.emplace_back(2,17857,-1,-1);
            merged_regions.emplace_back(250,23401,-1,-1);
            merged_regions.emplace_back(0,23402,-1,-1);
            merged_regions.emplace_back(250,28142,-1,-1);
            merged_regions.emplace_back(1,28143,-1,-1);
            merged_regions.emplace_back(250,29890,-1,-1);
            break;
        }
        case 5:
        {
            merged_regions.emplace_back(250,239,-1,-1);
            merged_regions.emplace_back(1,240,-1,-1);
            merged_regions.emplace_back(250,3035,-1,-1);
            merged_regions.emplace_back(1,3036,-1,-1);
            merged_regions.emplace_back(250,8780,-1,-1);
            merged_regions.emplace_back(3,8781,-1,-1);
            merged_regions.emplace_back(250,14406,-1,-1);
            merged_regions.emplace_back(1,14407,-1,-1);
            merged_regions.emplace_back(250,23401,-1,-1);
            merged_regions.emplace_back(0,23402,-1,-1);
            merged_regions.emplace_back(250,26086,-1,-1);
            merged_regions.emplace_back(3,26087,-1,-1);
            merged_regions.emplace_back(250,28142,-1,-1);
            merged_regions.emplace_back(1,28143,-1,-1);
            merged_regions.emplace_back(250,29890,-1,-1);
            break;
        }
        case 6:
        {
            merged_regions.emplace_back(250,239,-1,-1);
            merged_regions.emplace_back(1,240,-1,-1);
            merged_regions.emplace_back(250,3035,-1,-1);
            merged_regions.emplace_back(1,3036,-1,-1);
            merged_regions.emplace_back(250,8780,-1,-1);
            merged_regions.emplace_back(3,8781,-1,-1);
            merged_regions.emplace_back(250,11228,-1,-1);
            merged_regions.emplace_back(3,11229,-1,-1);
            merged_regions.emplace_back(250,14406,-1,-1);
            merged_regions.emplace_back(1,14407,-1,-1);
            merged_regions.emplace_back(250,23401,-1,-1);
            merged_regions.emplace_back(0,23402,-1,-1);
            merged_regions.emplace_back(250,28142,-1,-1);
            merged_regions.emplace_back(1,28143,-1,-1);
            merged_regions.emplace_back(250,28165,-1,-1);
            merged_regions.emplace_back(0,28166,-1,-1);
            merged_regions.emplace_back(250,28876,-1,-1);
            merged_regions.emplace_back(0,28877,-1,-1);
            merged_regions.emplace_back(250,29736,-1,-1);
            merged_regions.emplace_back(3,29737,-1,-1);
            merged_regions.emplace_back(250,29740,-1,-1);
            merged_regions.emplace_back(0,29741,-1,-1);
            merged_regions.emplace_back(250,29890,-1,-1);
            break;
        }
        case 7:
        {
            merged_regions.emplace_back(250,11,6.6909772172225756721973544127735067377216182649136e-05,-1);
            merged_regions.emplace_back(250,239,-1,-1);
            merged_regions.emplace_back(1,240,-1,-1);
            merged_regions.emplace_back(250,3035,-1,-1);
            merged_regions.emplace_back(1,3036,-1,-1);
            merged_regions.emplace_back(250,3177,-1,-1);
            merged_regions.emplace_back(0,3178,-1,-1);
            merged_regions.emplace_back(250,6980,-1,-1);
            merged_regions.emplace_back(3,6981,-1,-1);
            merged_regions.emplace_back(250,8780,-1,-1);
            merged_regions.emplace_back(3,8781,-1,-1);
            merged_regions.emplace_back(250,14406,-1,-1);
            merged_regions.emplace_back(1,14407,-1,-1);
            merged_regions.emplace_back(250,23401,-1,-1);
            merged_regions.emplace_back(0,23402,-1,-1);
            merged_regions.emplace_back(250,28142,-1,-1);
            merged_regions.emplace_back(1,28143,-1,-1);
            merged_regions.emplace_back(250,29866,-1,-1);
            merged_regions.emplace_back(250,29867,6.6909772172225756721973544127735067377216182649136e-05,-1);
            merged_regions.emplace_back(250,29890,-1,-1);
            break;
        }
        case 8:
        {
            merged_regions.emplace_back(250,239,-1,-1);
            merged_regions.emplace_back(1,240,-1,-1);
            merged_regions.emplace_back(250,1267,-1,-1);
            auto new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value4{0.99999721161091070786852696983260102570056915283203,2.7982277946960444848486499948455024505689081593118e-10,2.7878294439446026602288583595701254580490058287978e-06,2.7982277946960444848486499948455024505689081593118e-10};
            (*new_lh) = new_lh_value4;
            merged_regions.emplace_back(251,1268,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,3035,-1,-1);
            merged_regions.emplace_back(1,3036,-1,-1);
            merged_regions.emplace_back(250,8780,-1,-1);
            merged_regions.emplace_back(3,8781,-1,-1);
            merged_regions.emplace_back(250,11228,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value5{5.5759387092817078802929955938516570768115343526006e-06,5.5759387092817078802929955938516570768115343526006e-06,0.49999442406129068761089229155913926661014556884766,0.49999442406129068761089229155913926661014556884766};
            (*new_lh) = new_lh_value5;
            merged_regions.emplace_back(251,11229,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,14371,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value6{0.99999721161091070786852696983260102570056915283203,2.7878294439446026602288583595701254580490058287978e-06,2.7982277946960444848486499948455024505689081593118e-10,2.7982277946960444848486499948455024505689081593118e-10};
            (*new_lh) = new_lh_value6;
            merged_regions.emplace_back(251,14372,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,14406,-1,-1);
            merged_regions.emplace_back(1,14407,-1,-1);
            merged_regions.emplace_back(250,17857,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value7{5.5759387092817078802929955938516570768115343526006e-06,0.49999442406129068761089229155913926661014556884766,5.5759387092817078802929955938516570768115343526006e-06,0.49999442406129068761089229155913926661014556884766};
            (*new_lh) = new_lh_value7;
            merged_regions.emplace_back(251,17858,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,21362,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value8{2.7982277946960444848486499948455024505689081593118e-10,0.99999721161091070786852696983260102570056915283203,2.7982277946960444848486499948455024505689081593118e-10,2.7878294439446026602288583595701254580490058287978e-06};
            (*new_lh) = new_lh_value8;
            merged_regions.emplace_back(251,21363,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,23401,-1,-1);
            merged_regions.emplace_back(0,23402,-1,-1);
            merged_regions.emplace_back(250,28142,-1,-1);
            merged_regions.emplace_back(1,28143,-1,-1);
            merged_regions.emplace_back(250,28165,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value9{0.99995538913255155311077260193997062742710113525391,1.1193098382510832E-9,4.4608628828755082513634472318742041352379601448774e-05,1.1193098382510834E-9};
            (*new_lh) = new_lh_value9;
            merged_regions.emplace_back(251,28166,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,28876,-1,-1);
            merged_regions.emplace_back(0,28877,-1,-1);
            merged_regions.emplace_back(250,29740,-1,-1);
            merged_regions.emplace_back(0,29741,-1,-1);
            merged_regions.emplace_back(250,29890,-1,-1);
            break;
        }
        case 9:
        {
            merged_regions.emplace_back(250,239,-1,-1);
            merged_regions.emplace_back(1,240,-1,-1);
            merged_regions.emplace_back(250,1267,-1,-1);
            auto new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value5{0.99999442305078967141440671184682287275791168212891,5.5966272240860645628430121287970669397004996881151e-10,5.5758298847802999098784669518291678969035274349153e-06,5.5966272240860645628430121287970669397004996881151e-10};
            (*new_lh) = new_lh_value5;
            merged_regions.emplace_back(251,1268,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,3035,-1,-1);
            merged_regions.emplace_back(1,3036,-1,-1);
            merged_regions.emplace_back(250,8780,-1,-1);
            merged_regions.emplace_back(3,8781,-1,-1);
            merged_regions.emplace_back(250,14371,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value6{0.99999442305078967141440671184682287275791168212891,5.5758298847802999098784669518291678969035274349153e-06,5.5966272240860645628430121287970669397004996881151e-10,5.5966272240860645628430121287970669397004996881151e-10};
            (*new_lh) = new_lh_value6;
            merged_regions.emplace_back(251,14372,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,14406,-1,-1);
            merged_regions.emplace_back(1,14407,-1,-1);
            merged_regions.emplace_back(250,16573,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value7{7.4346402191731947664294494204639818235591519623995e-06,0.66666418845326025355291221785591915249824523925781,7.4346402191731947664294494204639818235591519623995e-06,0.33332094226630137878686355179524980485439300537109};
            (*new_lh) = new_lh_value7;
            merged_regions.emplace_back(251,16574,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,17857,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value8{7.4346402191731947664294494204639818235591519623995e-06,0.33332094226630137878686355179524980485439300537109,7.4346402191731947664294494204639818235591519623995e-06,0.66666418845326025355291221785591915249824523925781};
            (*new_lh) = new_lh_value8;
            merged_regions.emplace_back(251,17858,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,21362,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value9{5.5966272240860645628430121287970669397004996881151e-10,0.99999442305078967141440671184682287275791168212891,5.5966272240860645628430121287970669397004996881151e-10,5.5758298847802999098784669518291678969035274349153e-06};
            (*new_lh) = new_lh_value9;
            merged_regions.emplace_back(251,21363,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,23401,-1,-1);
            merged_regions.emplace_back(0,23402,-1,-1);
            merged_regions.emplace_back(250,26311,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value10{7.4346402191731947664294494204639818235591519623995e-06,0.66666418845326025355291221785591915249824523925781,7.4346402191731947664294494204639818235591519623995e-06,0.33332094226630137878686355179524980485439300537109};
            (*new_lh) = new_lh_value10;
            merged_regions.emplace_back(251,26312,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,28142,-1,-1);
            merged_regions.emplace_back(1,28143,-1,-1);
            merged_regions.emplace_back(250,28165,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value11{0.99991077926048133051040167629253119230270385742188,2.238594701945873E-9,8.9216262329251821180768622365775399885023944079876e-05,2.2385947019458738E-9};
            (*new_lh) = new_lh_value11;
            merged_regions.emplace_back(251,28166,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,28876,-1,-1);
            merged_regions.emplace_back(0,28877,-1,-1);
            merged_regions.emplace_back(250,29740,-1,-1);
            merged_regions.emplace_back(0,29741,-1,-1);
            merged_regions.emplace_back(250,29890,-1,-1);
            break;
        }
        case 10:
        {
            merged_regions.emplace_back(250,239,-1,-1);
            merged_regions.emplace_back(1,240,-1,-1);
            merged_regions.emplace_back(250,3035,-1,-1);
            merged_regions.emplace_back(1,3036,-1,-1);
            merged_regions.emplace_back(250,8780,-1,-1);
            merged_regions.emplace_back(3,8781,-1,-1);
            merged_regions.emplace_back(250,14406,-1,-1);
            merged_regions.emplace_back(1,14407,-1,-1);
            merged_regions.emplace_back(250,17745,-1,-1);
            merged_regions.emplace_back(3,17746,-1,-1);
            merged_regions.emplace_back(250,17856,-1,-1);
            merged_regions.emplace_back(2,17857,-1,-1);
            merged_regions.emplace_back(250,18058,-1,-1);
            merged_regions.emplace_back(3,18059,-1,-1);
            merged_regions.emplace_back(250,23401,-1,-1);
            merged_regions.emplace_back(0,23402,-1,-1);
            merged_regions.emplace_back(250,28142,-1,-1);
            merged_regions.emplace_back(1,28143,-1,-1);
            merged_regions.emplace_back(250,29890,-1,-1);
            break;
        }
        default:
            break;
    }
}

/*
    Generate output data (merged_seqregions) for testing mergeTwoLowers()
 */
void genOutputData3_2(SeqRegions& merged_regions, int test_case)
{
    merged_regions.clear();
    
    switch (test_case)
    {
        case 1:
        {
            merged_regions.emplace_back(250,239,-1,-1);
            merged_regions.emplace_back(1,240,-1,-1);
            merged_regions.emplace_back(250,3035,-1,-1);
            merged_regions.emplace_back(1,3036,-1,-1);
            merged_regions.emplace_back(250,11188,-1,-1);
            merged_regions.emplace_back(3,11189,-1,-1);
            merged_regions.emplace_back(250,14406,-1,-1);
            merged_regions.emplace_back(1,14407,-1,-1);
            merged_regions.emplace_back(250,23401,-1,-1);
            merged_regions.emplace_back(0,23402,-1,-1);
            merged_regions.emplace_back(250,28142,-1,-1);
            merged_regions.emplace_back(1,28143,-1,-1);
            merged_regions.emplace_back(250,29890,-1,-1);
            break;
        }
        case 2:
        {
            merged_regions.emplace_back(250,239,-1,-1);
            merged_regions.emplace_back(1,240,-1,-1);
            merged_regions.emplace_back(250,3035,-1,-1);
            merged_regions.emplace_back(1,3036,-1,-1);
            merged_regions.emplace_back(250,14406,-1,-1);
            merged_regions.emplace_back(1,14407,-1,-1);
            merged_regions.emplace_back(250,23401,-1,-1);
            auto new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value{0.12497560486006949187487435892762732692062854766846,9.7580559722090576469559833339140197949745925143361e-06,0.87500487902798607109389195102266967296600341796875,9.7580559722090576469559833339140197949745925143361e-06};
            (*new_lh) = new_lh_value;
            merged_regions.emplace_back(251,23402,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,28142,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value1{9.7580559722090576469559833339140197949745925143361e-06,0.12497560486006949187487435892762732692062854766846,9.7580559722090576469559833339140197949745925143361e-06,0.87500487902798607109389195102266967296600341796875};
            (*new_lh) = new_lh_value1;
            merged_regions.emplace_back(251,28143,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,28165,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value2{0.87500487902798607109389195102266967296600341796875,9.7580559722090576469559833339140197949745925143361e-06,0.12497560486006949187487435892762732692062854766846,9.7580559722090576469559833339140197949745925143361e-06};
            (*new_lh) = new_lh_value2;
            merged_regions.emplace_back(251,28166,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,28873,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value3{0.87500487902798607109389195102266967296600341796875,9.7580559722090576469559833339140197949745925143361e-06,0.12497560486006949187487435892762732692062854766846,9.7580559722090576469559833339140197949745925143361e-06};
            (*new_lh) = new_lh_value3;
            merged_regions.emplace_back(251,28874,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,28876,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value4{0.87500487902798607109389195102266967296600341796875,9.7580559722090576469559833339140197949745925143361e-06,0.12497560486006949187487435892762732692062854766846,9.7580559722090576469559833339140197949745925143361e-06};
            (*new_lh) = new_lh_value4;
            merged_regions.emplace_back(251,28877,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,29740,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value5{0.87500487902798607109389195102266967296600341796875,9.7580559722090576469559833339140197949745925143361e-06,0.12497560486006949187487435892762732692062854766846,9.7580559722090576469559833339140197949745925143361e-06};
            (*new_lh) = new_lh_value5;
            merged_regions.emplace_back(251,29741,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,29890,-1,-1);
            break;
        }
        case 3:
        {
            merged_regions.emplace_back(250,239,-1,-1);
            auto new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value{2.6766774402933635499746492514283602304203668609262e-05,0.39996252651583585890904259940725751221179962158203,2.6766774402933635499746492514283602304203668609262e-05,0.59998393993535814594508792652050033211708068847656};
            (*new_lh) = new_lh_value;
            merged_regions.emplace_back(251,240,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,3035,-1,-1);
            merged_regions.emplace_back(1,3036,-1,-1);
            merged_regions.emplace_back(250,8780,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value1{2.6766774402933635499746492514283602304203668609262e-05,0.39996252651583585890904259940725751221179962158203,2.6766774402933635499746492514283602304203668609262e-05,0.59998393993535814594508792652050033211708068847656};
            (*new_lh) = new_lh_value1;
            merged_regions.emplace_back(251,8781,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,14406,-1,-1);
            merged_regions.emplace_back(1,14407,-1,-1);
            merged_regions.emplace_back(250,18754,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value2{2.6766774402933635499746492514283602304203668609262e-05,2.6766774402933635499746492514283602304203668609262e-05,0.39996252651583585890904259940725751221179962158203,0.59998393993535814594508792652050033211708068847656};
            (*new_lh) = new_lh_value2;
            merged_regions.emplace_back(251,18755,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,22466,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value3{2.6766774402933635499746492514283602304203668609262e-05,2.6766774402933635499746492514283602304203668609262e-05,0.39996252651583585890904259940725751221179962158203,0.59998393993535814594508792652050033211708068847656};
            (*new_lh) = new_lh_value3;
            merged_regions.emplace_back(251,22467,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,23401,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value4{0.39996252651583585890904259940725751221179962158203,2.6766774402933635499746492514283602304203668609262e-05,0.59998393993535814594508792652050033211708068847656,2.6766774402933635499746492514283602304203668609262e-05};
            (*new_lh) = new_lh_value4;
            merged_regions.emplace_back(251,23402,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,28142,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value5{2.6766774402933635499746492514283602304203668609262e-05,0.39996252651583585890904259940725751221179962158203,2.6766774402933635499746492514283602304203668609262e-05,0.59998393993535814594508792652050033211708068847656};
            (*new_lh) = new_lh_value5;
            merged_regions.emplace_back(251,28143,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,28876,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value6{0.59998393993535814594508792652050033211708068847656,2.6766774402933635499746492514283602304203668609262e-05,0.39996252651583585890904259940725751221179962158203,2.6766774402933635499746492514283602304203668609262e-05};
            (*new_lh) = new_lh_value6;
            merged_regions.emplace_back(251,28877,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,29708,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value7{0.59998393993535814594508792652050033211708068847656,2.6766774402933635499746492514283602304203668609262e-05,2.6766774402933635499746492514283602304203668609262e-05,0.39996252651583585890904259940725751221179962158203};
            (*new_lh) = new_lh_value7;
            merged_regions.emplace_back(251,29709,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,29740,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value8{0.59998393993535814594508792652050033211708068847656,2.6766774402933635499746492514283602304203668609262e-05,0.39996252651583585890904259940725751221179962158203,2.6766774402933635499746492514283602304203668609262e-05};
            (*new_lh) = new_lh_value8;
            merged_regions.emplace_back(251,29741,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,29890,-1,-1);
            break;
        }
        case 4:
        {
            merged_regions.emplace_back(250,239,-1,-1);
            merged_regions.emplace_back(1,240,-1,-1);
            merged_regions.emplace_back(250,3035,-1,-1);
            merged_regions.emplace_back(1,3036,-1,-1);
            merged_regions.emplace_back(250,8780,-1,-1);
            merged_regions.emplace_back(3,8781,-1,-1);
            merged_regions.emplace_back(250,14406,-1,-1);
            merged_regions.emplace_back(1,14407,-1,-1);
            merged_regions.emplace_back(250,17856,-1,-1);
            merged_regions.emplace_back(2,17857,-1,-1);
            merged_regions.emplace_back(250,23401,-1,-1);
            merged_regions.emplace_back(0,23402,-1,-1);
            merged_regions.emplace_back(250,28142,-1,-1);
            merged_regions.emplace_back(1,28143,-1,-1);
            merged_regions.emplace_back(250,29890,-1,-1);
            break;
        }
        case 5:
        {
            merged_regions.emplace_back(250,239,-1,-1);
            merged_regions.emplace_back(1,240,-1,-1);
            merged_regions.emplace_back(250,3035,-1,-1);
            merged_regions.emplace_back(1,3036,-1,-1);
            merged_regions.emplace_back(250,8780,-1,-1);
            merged_regions.emplace_back(3,8781,-1,-1);
            merged_regions.emplace_back(250,14406,-1,-1);
            merged_regions.emplace_back(1,14407,-1,-1);
            merged_regions.emplace_back(250,23401,-1,-1);
            merged_regions.emplace_back(0,23402,-1,-1);
            merged_regions.emplace_back(250,26086,-1,-1);
            merged_regions.emplace_back(3,26087,-1,-1);
            merged_regions.emplace_back(250,28142,-1,-1);
            merged_regions.emplace_back(1,28143,-1,-1);
            merged_regions.emplace_back(250,29890,-1,-1);
            break;
        }
        case 6:
        {
            merged_regions.emplace_back(250,239,-1,-1);
            merged_regions.emplace_back(1,240,-1,-1);
            merged_regions.emplace_back(250,3035,-1,-1);
            merged_regions.emplace_back(1,3036,-1,-1);
            merged_regions.emplace_back(250,8780,-1,-1);
            merged_regions.emplace_back(3,8781,-1,-1);
            merged_regions.emplace_back(250,11228,-1,-1);
            merged_regions.emplace_back(3,11229,-1,-1);
            merged_regions.emplace_back(250,14406,-1,-1);
            merged_regions.emplace_back(1,14407,-1,-1);
            merged_regions.emplace_back(250,23401,-1,-1);
            merged_regions.emplace_back(0,23402,-1,-1);
            merged_regions.emplace_back(250,28142,-1,-1);
            merged_regions.emplace_back(1,28143,-1,-1);
            merged_regions.emplace_back(250,28165,-1,-1);
            merged_regions.emplace_back(0,28166,-1,-1);
            merged_regions.emplace_back(250,28876,-1,-1);
            merged_regions.emplace_back(0,28877,-1,-1);
            merged_regions.emplace_back(250,29736,-1,-1);
            merged_regions.emplace_back(3,29737,-1,-1);
            merged_regions.emplace_back(250,29740,-1,-1);
            merged_regions.emplace_back(0,29741,-1,-1);
            merged_regions.emplace_back(250,29890,-1,-1);
            break;
        }
        case 7:
        {
            merged_regions.emplace_back(250,11,6.6909772172225756721973544127735067377216182649136e-05,-1);
            merged_regions.emplace_back(250,239,-1,-1);
            merged_regions.emplace_back(1,240,-1,-1);
            merged_regions.emplace_back(250,3035,-1,-1);
            merged_regions.emplace_back(1,3036,-1,-1);
            merged_regions.emplace_back(250,3177,-1,-1);
            merged_regions.emplace_back(0,3178,-1,-1);
            merged_regions.emplace_back(250,6980,-1,-1);
            merged_regions.emplace_back(3,6981,-1,-1);
            merged_regions.emplace_back(250,8780,-1,-1);
            merged_regions.emplace_back(3,8781,-1,-1);
            merged_regions.emplace_back(250,14406,-1,-1);
            merged_regions.emplace_back(1,14407,-1,-1);
            merged_regions.emplace_back(250,23401,-1,-1);
            merged_regions.emplace_back(0,23402,-1,-1);
            merged_regions.emplace_back(250,28142,-1,-1);
            merged_regions.emplace_back(1,28143,-1,-1);
            merged_regions.emplace_back(250,29866,-1,-1);
            merged_regions.emplace_back(250,29867,6.6909772172225756721973544127735067377216182649136e-05,-1);
            merged_regions.emplace_back(250,29890,-1,-1);
            break;
        }
        case 8:
        {
            merged_regions.emplace_back(250,239,-1,-1);
            merged_regions.emplace_back(1,240,-1,-1);
            merged_regions.emplace_back(250,1267,-1,-1);
            auto new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value{0.44442861163664077,1.2436159615250004E-10,0.0000012389946240511395,1.2436159615250004E-10};
            (*new_lh) = new_lh_value;
            merged_regions.emplace_back(251,1268,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,3035,-1,-1);
            merged_regions.emplace_back(1,3036,-1,-1);
            merged_regions.emplace_back(250,8780,-1,-1);
            merged_regions.emplace_back(3,8781,-1,-1);
            merged_regions.emplace_back(250,11228,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value1{5.5759387092817078802929955938516570768115343526006e-06,5.5759387092817078802929955938516570768115343526006e-06,0.49999442406129068761089229155913926661014556884766,0.49999442406129068761089229155913926661014556884766};
            (*new_lh) = new_lh_value1;
            merged_regions.emplace_back(251,11229,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,14371,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value2{0.44442861163664077,0.0000012389946240511395,1.2436159615250004E-10,1.2436159615250004E-10};
            (*new_lh) = new_lh_value2;
            merged_regions.emplace_back(251,14372,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,14406,-1,-1);
            merged_regions.emplace_back(1,14407,-1,-1);
            merged_regions.emplace_back(250,17857,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value3{5.5759387092817078802929955938516570768115343526006e-06,0.49999442406129068761089229155913926661014556884766,5.5759387092817078802929955938516570768115343526006e-06,0.49999442406129068761089229155913926661014556884766};
            (*new_lh) = new_lh_value3;
            merged_regions.emplace_back(251,17858,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,21362,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value4{1.2436159615250004E-10,0.44442861163664077,1.2436159615250004E-10,0.0000012389946240511395};
            (*new_lh) = new_lh_value4;
            merged_regions.emplace_back(251,21363,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,23401,-1,-1);
            merged_regions.emplace_back(0,23402,-1,-1);
            merged_regions.emplace_back(250,28142,-1,-1);
            merged_regions.emplace_back(1,28143,-1,-1);
            merged_regions.emplace_back(250,28165,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value5{0.16664870042553384,1.8653985162265207E-10,0.0000074343016727235307,1.8653985162265212E-10};
            (*new_lh) = new_lh_value5;
            merged_regions.emplace_back(251,28166,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,28876,-1,-1);
            merged_regions.emplace_back(0,28877,-1,-1);
            merged_regions.emplace_back(250,29740,-1,-1);
            merged_regions.emplace_back(0,29741,-1,-1);
            merged_regions.emplace_back(250,29890,-1,-1);
            break;
        }
        case 9:
        {
            merged_regions.emplace_back(250,239,-1,-1);
            merged_regions.emplace_back(1,240,-1,-1);
            merged_regions.emplace_back(250,1267,-1,-1);
            auto new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value{0.44441980007814447,2.4872658233953994E-10,0.000002478023022472343,2.4872658233953994E-10};
            (*new_lh) = new_lh_value;
            merged_regions.emplace_back(251,1268,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,3035,-1,-1);
            merged_regions.emplace_back(1,3036,-1,-1);
            merged_regions.emplace_back(250,8780,-1,-1);
            merged_regions.emplace_back(3,8781,-1,-1);
            merged_regions.emplace_back(250,14371,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value1{0.44441980007814447,0.000002478023022472343,2.4872658233953994E-10,2.4872658233953994E-10};
            (*new_lh) = new_lh_value1;
            merged_regions.emplace_back(251,14372,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,14406,-1,-1);
            merged_regions.emplace_back(1,14407,-1,-1);
            merged_regions.emplace_back(250,16573,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value2{7.4346402191731947664294494204639818235591519623995e-06,0.66666418845326025355291221785591915249824523925781,7.4346402191731947664294494204639818235591519623995e-06,0.33332094226630137878686355179524980485439300537109};
            (*new_lh) = new_lh_value2;
            merged_regions.emplace_back(251,16574,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,17857,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value3{7.4346402191731947664294494204639818235591519623995e-06,0.33332094226630137878686355179524980485439300537109,7.4346402191731947664294494204639818235591519623995e-06,0.66666418845326025355291221785591915249824523925781};
            (*new_lh) = new_lh_value3;
            merged_regions.emplace_back(251,17858,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,21362,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value4{2.4872658233953994E-10,0.44441980007814447,2.4872658233953994E-10,0.000002478023022472343};
            (*new_lh) = new_lh_value4;
            merged_regions.emplace_back(251,21363,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,23401,-1,-1);
            merged_regions.emplace_back(0,23402,-1,-1);
            merged_regions.emplace_back(250,26311,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value5{7.4346402191731947664294494204639818235591519623995e-06,0.66666418845326025355291221785591915249824523925781,7.4346402191731947664294494204639818235591519623995e-06,0.33332094226630137878686355179524980485439300537109};
            (*new_lh) = new_lh_value5;
            merged_regions.emplace_back(251,26312,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,28142,-1,-1);
            merged_regions.emplace_back(1,28143,-1,-1);
            merged_regions.emplace_back(250,28165,-1,-1);
            new_lh = std::make_unique<SeqRegion::LHType>();
            SeqRegion::LHType new_lh_value6{0.16664281519091029,3.7307900958885143E-10,0.000014868575700676886,3.7307900958885153E-10};
            (*new_lh) = new_lh_value6;
            merged_regions.emplace_back(251,28166,0,0,std::move(new_lh));
            merged_regions.emplace_back(250,28876,-1,-1);
            merged_regions.emplace_back(0,28877,-1,-1);
            merged_regions.emplace_back(250,29740,-1,-1);
            merged_regions.emplace_back(0,29741,-1,-1);
            merged_regions.emplace_back(250,29890,-1,-1);
            break;
        }
        case 10:
        {
            merged_regions.emplace_back(250,239,-1,-1);
            merged_regions.emplace_back(1,240,-1,-1);
            merged_regions.emplace_back(250,3035,-1,-1);
            merged_regions.emplace_back(1,3036,-1,-1);
            merged_regions.emplace_back(250,8780,-1,-1);
            merged_regions.emplace_back(3,8781,-1,-1);
            merged_regions.emplace_back(250,14406,-1,-1);
            merged_regions.emplace_back(1,14407,-1,-1);
            merged_regions.emplace_back(250,17745,-1,-1);
            merged_regions.emplace_back(3,17746,-1,-1);
            merged_regions.emplace_back(250,17856,-1,-1);
            merged_regions.emplace_back(2,17857,-1,-1);
            merged_regions.emplace_back(250,18058,-1,-1);
            merged_regions.emplace_back(3,18059,-1,-1);
            merged_regions.emplace_back(250,23401,-1,-1);
            merged_regions.emplace_back(0,23402,-1,-1);
            merged_regions.emplace_back(250,28142,-1,-1);
            merged_regions.emplace_back(1,28143,-1,-1);
            merged_regions.emplace_back(250,29890,-1,-1);
            break;
        }
        default:
            break;
    }
}

/*
    Generate testing data (seqregions1, seqregions2, merged_seqregions) for testing mergeUpperLower()
 */
void genTestData3(SeqRegions& seqregions1, SeqRegions& seqregions2, SeqRegions& merged_regions, int test_case)
{
    genTestData3_0(seqregions1, seqregions2, test_case);
    genOutputData3_1(merged_regions, test_case);
}
    
/*
    Generate testing data (seqregions1, seqregions2, merged_seqregions) for mergeTwoLowers()
 */
void genTestData4(SeqRegions& seqregions1, SeqRegions& seqregions2, SeqRegions& merged_regions, int test_case)
{
    genTestData3_0(seqregions1, seqregions2, test_case);
    genOutputData3_2(merged_regions, test_case);
}

/*
 Test mergeUpperLower(SeqRegions* &merged_regions, RealNumType upper_plength, const SeqRegions& lower_regions, RealNumType lower_plength, const
 Alignment& aln, const Model& model, RealNumType threshold) const
 */
TEST(SeqRegions, mergeUpperLower)
{
    Alignment aln;
    Model model;
    Params params = Params::getInstance();
    
    // Init params, aln, and model
    initAlnModelParams(params, aln, model, "JC");
    // dummy variables
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    const RealNumType threshold_prob = params.threshold_prob;
    
    // ----- Test 1 -----
    SeqRegions seqregions1, seqregions2, output_regions;
    SeqRegions* merged_regions_ptr = nullptr;
    
    genTestData3(seqregions1, seqregions2, output_regions, 1);
    seqregions1.mergeUpperLower(merged_regions_ptr, 3.3454886086112878360986772063867533688608091324568e-05, seqregions2, -1, aln, model, threshold_prob);
    EXPECT_EQ(*merged_regions_ptr, output_regions);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    genTestData3(seqregions1, seqregions2, output_regions, 2);
    seqregions1.mergeUpperLower(merged_regions_ptr, 0.00023418420260279014175064382641267002327367663383484, seqregions2, 3.3454886086112878360986772063867533688608091324568e-05, aln, model, threshold_prob);
    EXPECT_EQ(*merged_regions_ptr, output_regions);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    genTestData3(seqregions1, seqregions2, output_regions, 3);
    seqregions1.mergeUpperLower(merged_regions_ptr, 0.00020072931651667725661339347631439977703848853707314, seqregions2, 0.00013381954434445151344394708825547013475443236529827, aln, model, threshold_prob);
    EXPECT_EQ(*merged_regions_ptr, output_regions);
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    genTestData3(seqregions1, seqregions2, output_regions, 4);
    seqregions1.mergeUpperLower(merged_regions_ptr, 3.3454886086112878360986772063867533688608091324568e-05, seqregions2, -1, aln, model, threshold_prob);
    EXPECT_EQ(*merged_regions_ptr, output_regions);
    // ----- Test 4 -----
    
    // ----- Test 5 -----
    genTestData3(seqregions1, seqregions2, output_regions, 5);
    seqregions1.mergeUpperLower(merged_regions_ptr, -1, seqregions2, 5.0182329129169314153348369078599944259622134268284e-05, aln, model, threshold_prob);
    EXPECT_EQ(*merged_regions_ptr, output_regions);
    // ----- Test 5 -----
    
    // ----- Test 6 -----
    genTestData3(seqregions1, seqregions2, output_regions, 6);
    seqregions1.mergeUpperLower(merged_regions_ptr, 3.3454886086112878360986772063867533688608091324568e-05, seqregions2, -1, aln, model, threshold_prob);
    EXPECT_EQ(*merged_regions_ptr, output_regions);
    // ----- Test 6 -----
    
    // ----- Test 7 -----
    genTestData3(seqregions1, seqregions2, output_regions, 7);
    seqregions1.mergeUpperLower(merged_regions_ptr, -1, seqregions2, -1, aln, model, threshold_prob);
    EXPECT_EQ(*merged_regions_ptr, output_regions);
    // ----- Test 7 -----
    
    // ----- Test 8 -----
    genTestData3(seqregions1, seqregions2, output_regions, 8);
    seqregions1.mergeUpperLower(merged_regions_ptr, 3.3454886086112878360986772063867533688608091324568e-05, seqregions2, 3.3454886086112878360986772063867533688608091324568e-05, aln, model, threshold_prob);
    EXPECT_EQ(*merged_regions_ptr, output_regions);
    // ----- Test 8 -----
    
    // ----- Test 9 -----
    genTestData3(seqregions1, seqregions2, output_regions, 9);
    seqregions1.mergeUpperLower(merged_regions_ptr, 3.3454886086112878360986772063867533688608091324568e-05, seqregions2, 6.6909772172225756721973544127735067377216182649136e-05, aln, model, threshold_prob);
    EXPECT_EQ(*merged_regions_ptr, output_regions);
    // ----- Test 9 -----
    
    // ----- Test 10 -----
    genTestData3(seqregions1, seqregions2, output_regions, 10);
    seqregions1.mergeUpperLower(merged_regions_ptr, -1, seqregions2, -1, aln, model, threshold_prob);
    EXPECT_EQ(*merged_regions_ptr, output_regions);
    // ----- Test 10 -----
}

/*
 Test merge_N_O_TwoLowers(const SeqRegion& seq2_region, const PositionType end_pos, const RealNumType plength2, SeqRegions& merged_regions)
 */
TEST(SeqRegions, merge_N_O_TwoLowers)
{
    Alignment aln;
    Model model;
    Params params = Params::getInstance();
    
    // Init params, aln, and model
    initAlnModelParams(params, aln, model);
    // dummy variables
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    
    // ----- Test 1 -----
    SeqRegions merged_regions;
    const PositionType end_pos = 8623;
    RealNumType plength2 = -1;
    auto new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value{0.27595781923293211113090706021466758102178573608398,8.3817209383128356075479820086471249851456377655268e-06,0.72401575181023847260775028189527802169322967529297,1.8047235891267821277644811672757896303664892911911e-05};
    (*new_lh) = new_lh_value;
    SeqRegion seqregion1(TYPE_O, 243, -1, -1, std::move(new_lh));
    merge_N_O_TwoLowers(seqregion1, end_pos, plength2, merged_regions);
    EXPECT_EQ(merged_regions.size(), 1);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    seqregion1.plength_observation2root = 0.311;
    merge_N_O_TwoLowers(seqregion1, end_pos, plength2, merged_regions);
    EXPECT_EQ(merged_regions.size(), 2);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0.311);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    seqregion1.plength_observation2node = 0;
    merge_N_O_TwoLowers(seqregion1, end_pos, plength2, merged_regions);
    EXPECT_EQ(merged_regions.size(), 3);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0.311);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value);
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    seqregion1.plength_observation2node = 0.1252;
    merge_N_O_TwoLowers(seqregion1, end_pos, plength2, merged_regions);
    EXPECT_EQ(merged_regions.size(), 4);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0.311);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0.1252);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value);
    // ----- Test 4 -----
    
    // ----- Test 5 -----
    plength2 = 0; // plength2 = 0; plength_observation2node > 0
    merge_N_O_TwoLowers(seqregion1, end_pos, plength2, merged_regions);
    EXPECT_EQ(merged_regions.size(), 5);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0.311);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0.1252);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value);
    // ----- Test 5 -----
    
    // ----- Test 6 -----
    plength2 = 0; // plength2 = 0; plength_observation2node = 0
    seqregion1.plength_observation2node = 0;
    merge_N_O_TwoLowers(seqregion1, end_pos, plength2, merged_regions);
    EXPECT_EQ(merged_regions.size(), 6);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0.311);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value);
    // ----- Test 6 -----
    
    // ----- Test 7 -----
    plength2 = 0; // plength2 = 0; plength_observation2node = -1
    seqregion1.plength_observation2node = -1;
    merge_N_O_TwoLowers(seqregion1, end_pos, plength2, merged_regions);
    EXPECT_EQ(merged_regions.size(), 7);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0.311);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value);
    // ----- Test 7 -----
    
    // ----- Test 8 -----
    plength2 = 0.02331; // plength2 > 0; plength_observation2node = -1
    seqregion1.plength_observation2node = -1;
    merge_N_O_TwoLowers(seqregion1, end_pos, plength2, merged_regions);
    EXPECT_EQ(merged_regions.size(), 8);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0.311);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0.02331);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value);
    // ----- Test 8 -----
    
    // ----- Test 9 -----
    plength2 = 0.02331; // plength2 > 0; plength_observation2node = 0
    seqregion1.plength_observation2node = 0;
    merge_N_O_TwoLowers(seqregion1, end_pos, plength2, merged_regions);
    EXPECT_EQ(merged_regions.size(), 9);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0.311);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0.02331);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value);
    // ----- Test 9 -----
    
    // ----- Test 10 -----
    plength2 = 0.02331; // plength2 > 0; plength_observation2node > 0
    seqregion1.plength_observation2node = 0.1e-6;
    merge_N_O_TwoLowers(seqregion1, end_pos, plength2, merged_regions);
    EXPECT_EQ(merged_regions.size(), 10);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0.311);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0.0233101);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value);
    // ----- Test 10 -----
}

/*
 Test merge_N_RACGT_TwoLowers(const SeqRegion& seq2_region, const PositionType end_pos, const RealNumType plength2, const RealNumType threshold_prob, SeqRegions& merged_regions)
 */
TEST(SeqRegions, merge_N_RACGT_TwoLowers)
{
    Alignment aln;
    Model model;
    Params params = Params::getInstance();
    
    // Init params, aln, and model
    initAlnModelParams(params, aln, model);
    // dummy variables
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    const RealNumType threshold_prob = params.threshold_prob;
    
    // ----- Test 1 -----
    SeqRegions merged_regions;
    const PositionType end_pos = 6543;
    RealNumType plength2 = -1;
    SeqRegion seqregion1(TYPE_R, 3442, -1, -1);
    merge_N_RACGT_TwoLowers(seqregion1, end_pos, plength2, threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 1);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    plength2 = -1;
    seqregion1.plength_observation2node = -1;
    seqregion1.plength_observation2root = 0;
    seqregion1.type = 1; // C
    merge_N_RACGT_TwoLowers(seqregion1, end_pos, plength2, threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 2);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    plength2 = -1;
    seqregion1.plength_observation2node = -1;
    seqregion1.plength_observation2root = 0.3231;
    seqregion1.type = TYPE_R;
    merge_N_RACGT_TwoLowers(seqregion1, end_pos, plength2, threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 3);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    plength2 = -1;
    seqregion1.plength_observation2node = 0;
    seqregion1.plength_observation2root = -1;
    seqregion1.type = 3;
    merge_N_RACGT_TwoLowers(seqregion1, end_pos, plength2, threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 4);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 4 -----
    
    // ----- Test 5 -----
    plength2 = -1;
    seqregion1.plength_observation2node = 0;
    seqregion1.plength_observation2root = 0;
    seqregion1.type = 2;
    merge_N_RACGT_TwoLowers(seqregion1, end_pos, plength2, threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 5);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 5 -----
    
    // ----- Test 6 -----
    plength2 = -1;
    seqregion1.plength_observation2node = 0;
    seqregion1.plength_observation2root = 1e-3;
    seqregion1.type = TYPE_R;
    merge_N_RACGT_TwoLowers(seqregion1, end_pos, plength2, threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 6);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 6 -----
    
    // ----- Test 7 -----
    plength2 = -1;
    seqregion1.plength_observation2node = 1e-10;
    seqregion1.plength_observation2root = -1;
    seqregion1.type = 1;
    merge_N_RACGT_TwoLowers(seqregion1, end_pos, plength2, threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 7);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 1e-10);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 7 -----
    
    // ----- Test 8 -----
    plength2 = -1;
    seqregion1.plength_observation2node = 1e-10;
    seqregion1.plength_observation2root = 0;
    seqregion1.type = 2;
    merge_N_RACGT_TwoLowers(seqregion1, end_pos, plength2, threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 8);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 1e-10);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 8 -----
    
    // ----- Test 9 -----
    plength2 = -1;
    seqregion1.plength_observation2node = 1e-10;
    seqregion1.plength_observation2root = 0.1;
    seqregion1.type = 0;
    merge_N_RACGT_TwoLowers(seqregion1, end_pos, plength2, threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 9);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 1e-10);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 9 -----
    
    // ----- Test 10 -----
    plength2 = 0;
    seqregion1.plength_observation2node = -1;
    seqregion1.plength_observation2root = -1;
    seqregion1.type = 1;
    merge_N_RACGT_TwoLowers(seqregion1, end_pos, plength2, threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 10);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 10 -----
    
    // ----- Test 11 -----
    plength2 = 0;
    seqregion1.plength_observation2node = -1;
    seqregion1.plength_observation2root = 0;
    seqregion1.type = TYPE_R;
    merge_N_RACGT_TwoLowers(seqregion1, end_pos, plength2, threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 11);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 11 -----
    
    // ----- Test 12 -----
    plength2 = 0;
    seqregion1.plength_observation2node = -1;
    seqregion1.plength_observation2root = 0.3231;
    seqregion1.type = 3;
    merge_N_RACGT_TwoLowers(seqregion1, end_pos, plength2, threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 12);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 12 -----
    
    // ----- Test 13 -----
    plength2 = 0;
    seqregion1.plength_observation2node = 0;
    seqregion1.plength_observation2root = -1;
    seqregion1.type = 1;
    merge_N_RACGT_TwoLowers(seqregion1, end_pos, plength2, threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 13);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 13 -----
    
    // ----- Test 14 -----
    plength2 = 0;
    seqregion1.plength_observation2node = 0;
    seqregion1.plength_observation2root = 0;
    seqregion1.type = 0;
    merge_N_RACGT_TwoLowers(seqregion1, end_pos, plength2, threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 14);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 14 -----
    
    // ----- Test 15 -----
    plength2 = 0;
    seqregion1.plength_observation2node = 0;
    seqregion1.plength_observation2root = 1e-3;
    seqregion1.type = 0;
    merge_N_RACGT_TwoLowers(seqregion1, end_pos, plength2, threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 15);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 15 -----
    
    // ----- Test 16 -----
    plength2 = 0;
    seqregion1.plength_observation2node = 1e-10;
    seqregion1.plength_observation2root = -1;
    seqregion1.type = 1;
    merge_N_RACGT_TwoLowers(seqregion1, end_pos, plength2, threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 16);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 1e-10);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 16 -----
    
    // ----- Test 17 -----
    plength2 = 0;
    seqregion1.plength_observation2node = 1e-10;
    seqregion1.plength_observation2root = 0;
    seqregion1.type = 1;
    merge_N_RACGT_TwoLowers(seqregion1, end_pos, plength2, threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 17);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 1e-10);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 17 -----
    
    // ----- Test 18 -----
    plength2 = 0;
    seqregion1.plength_observation2node = 1e-10;
    seqregion1.plength_observation2root = 0.1;
    seqregion1.type = TYPE_R;
    merge_N_RACGT_TwoLowers(seqregion1, end_pos, plength2, threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 18);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 1e-10);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 18 -----
    
    // ----- Test 19 -----
    plength2 = 1e-5;
    seqregion1.plength_observation2node = -1;
    seqregion1.plength_observation2root = -1;
    seqregion1.type = 3;
    merge_N_RACGT_TwoLowers(seqregion1, end_pos, plength2, threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 19);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 1e-5);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 19 -----
    
    // ----- Test 20 -----
    plength2 = 1e-5;
    seqregion1.plength_observation2node = -1;
    seqregion1.plength_observation2root = 0;
    seqregion1.type = TYPE_R;
    merge_N_RACGT_TwoLowers(seqregion1, end_pos, plength2, threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 20);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 1e-5);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 20 -----
    
    // ----- Test 21 -----
    plength2 = 0.212;
    seqregion1.plength_observation2node = -1;
    seqregion1.plength_observation2root = 0.3231;
    seqregion1.type = 0;
    merge_N_RACGT_TwoLowers(seqregion1, end_pos, plength2, threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 21);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0.212);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 21 -----
    
    // ----- Test 22 -----
    plength2 = 0.011;
    seqregion1.plength_observation2node = 0;
    seqregion1.plength_observation2root = -1;
    seqregion1.type = 3;
    merge_N_RACGT_TwoLowers(seqregion1, end_pos, plength2, threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 22);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0.011);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 22 -----
    
    // ----- Test 23 -----
    plength2 = 1e-3;
    seqregion1.plength_observation2node = 0;
    seqregion1.plength_observation2root = 0;
    seqregion1.type = 3;
    merge_N_RACGT_TwoLowers(seqregion1, end_pos, plength2, threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 23);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 1e-3);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 23 -----
    
    // ----- Test 24 -----
    plength2 = 0.2312;
    seqregion1.plength_observation2node = 0;
    seqregion1.plength_observation2root = 1e-3;
    seqregion1.type = 2;
    merge_N_RACGT_TwoLowers(seqregion1, end_pos, plength2, threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 24);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0.2312);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 24 -----
    
    // ----- Test 25 -----
    plength2 = 0.04123;
    seqregion1.plength_observation2node = 1e-10;
    seqregion1.plength_observation2root = -1;
    seqregion1.type = 0;
    merge_N_RACGT_TwoLowers(seqregion1, end_pos, plength2, threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 25);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0.04123 + 1e-10);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 25 -----
    
    // ----- Test 26 -----
    plength2 = 0.001243;
    seqregion1.plength_observation2node = 1e-3;
    seqregion1.plength_observation2root = 0;
    seqregion1.type = TYPE_R;
    merge_N_RACGT_TwoLowers(seqregion1, end_pos, plength2, threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 26);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0.002243);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 26 -----
    
    // ----- Test 27 -----
    plength2 = 0.0001211;
    seqregion1.plength_observation2node = 1e-1;
    seqregion1.plength_observation2root = 0.1;
    seqregion1.type = 1;
    merge_N_RACGT_TwoLowers(seqregion1, end_pos, plength2, threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 27);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0.1001211);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 27 -----
}

/*
 Test merge_identicalRACGT_TwoLowers(const SeqRegion& seq1_region, const PositionType end_pos, RealNumType total_blength_1, RealNumType total_blength_2, const PositionType pos, const RealNumType threshold_prob, const Model& model, RealNumType &log_lh, SeqRegions& merged_regions, const bool return_log_lh)
 */
TEST(SeqRegions, merge_identicalRACGT_TwoLowers)
{
    Alignment aln;
    Model model;
    Params params = Params::getInstance();
    
    // Init params, aln, and model
    initAlnModelParams(params, aln, model);
    // dummy variables
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    const RealNumType threshold_prob = params.threshold_prob;
    
    // ----- Test 1 -----
    SeqRegions merged_regions;
    const PositionType end_pos = 6543;
    PositionType pos = 6489;
    RealNumType total_blength_1 = -1;
    RealNumType total_blength_2 = -1;
    RealNumType log_lh = 0;
    SeqRegion seqregion1(TYPE_R, 3442, -1, -1);
    merge_identicalRACGT_TwoLowers(seqregion1, end_pos, total_blength_1, total_blength_2, pos, threshold_prob, model, log_lh, merged_regions, true);
    EXPECT_EQ(merged_regions.size(), 1);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    EXPECT_EQ(log_lh, 0);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    total_blength_1 = 0;
    total_blength_2 = 0;
    log_lh = 0;
    seqregion1.type = 0;
    merge_identicalRACGT_TwoLowers(seqregion1, end_pos, total_blength_1, total_blength_2, pos, threshold_prob, model, log_lh, merged_regions, true);
    EXPECT_EQ(merged_regions.size(), 2);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    EXPECT_EQ(log_lh, 0);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    total_blength_1 = -1;
    total_blength_2 = 13e-4;
    log_lh = 0;
    seqregion1.type = 3;
    merge_identicalRACGT_TwoLowers(seqregion1, end_pos, total_blength_1, total_blength_2, pos, threshold_prob, model, log_lh, merged_regions, true);
    EXPECT_EQ(merged_regions.size(), 3);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    EXPECT_EQ(log_lh, -0.0016388829346472356);
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    total_blength_1 = 21e-6;
    total_blength_2 = -1;
    log_lh = 0;
    seqregion1.type = TYPE_R;
    merge_identicalRACGT_TwoLowers(seqregion1, end_pos, total_blength_1, total_blength_2, pos, threshold_prob, model, log_lh, merged_regions, true);
    EXPECT_EQ(merged_regions.size(), 4);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    EXPECT_EQ(log_lh, -0.00099527007095828778);
    // ----- Test 4 -----
    
    // ----- Test 5 -----
    total_blength_1 = 1e-5;
    total_blength_2 = 0.04231;
    log_lh = 0;
    seqregion1.type = 1;
    merge_identicalRACGT_TwoLowers(seqregion1, end_pos, total_blength_1, total_blength_2, pos, threshold_prob, model, log_lh, merged_regions, true);
    EXPECT_EQ(merged_regions.size(), 5);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    EXPECT_EQ(log_lh, -0.067217084471974248);
    // ----- Test 5 -----
}

/*
 Test merge_O_O_TwoLowers(const SeqRegion& seq2_region, RealNumType total_blength_2, const PositionType end_pos, const Alignment& aln, const Model& model, const RealNumType threshold_prob, RealNumType &log_lh, SeqRegion::LHType& new_lh, SeqRegions* merged_regions, const bool return_log_lh)
 */
TEST(SeqRegions, merge_O_O_TwoLowers)
{
    Alignment aln;
    Model model;
    Params params = Params::getInstance();
    
    // Init params, aln, and model
    initAlnModelParams(params, aln, model);
    // dummy variables
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    const RealNumType threshold_prob = params.threshold_prob;
    
    // ----- Test 1 -----
    SeqRegions* merged_regions = new SeqRegions();
    RealNumType log_lh = 0;
    auto new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value{0.12497560486006949187487435892762732692062854766846,9.7580559722090576469559833339140197949745925143361e-06,0.87500487902798607109389195102266967296600341796875,9.7580559722090576469559833339140197949745925143361e-06};
    (*new_lh) = new_lh_value;
    SeqRegion seqregion1(TYPE_O, 243, -1, -1, std::move(new_lh));
    auto new_lh1 = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value1{0.59998393993535814594508792652050033211708068847656,2.6766774402933635499746492514283602304203668609262e-05,2.6766774402933635499746492514283602304203668609262e-05,0.39996252651583585890904259940725751221179962158203};
    (*new_lh1) = new_lh_value1;
    RealNumType total_blength_2 = -1;
    const PositionType end_pos = 3213;
    EXPECT_TRUE(merge_O_O_TwoLowers(seqregion1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, *new_lh1, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 1);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge{0.99963572952385693,3.4820599267114684E-9,0.00031223631362821728,0.000052030680454886103};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge);
    EXPECT_EQ(log_lh, -2.5901247759055703);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    log_lh = 0;
    total_blength_2 = 1e-5;
    (*new_lh1) = new_lh_value1;
    EXPECT_TRUE(merge_O_O_TwoLowers(seqregion1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, *new_lh1, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 2);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge1{0.99961708509029623,3.8289366877191551E-9,0.00031222411049282615,0.000070686970274545192};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge1);
    EXPECT_EQ(log_lh, -2.5900955748485903);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    log_lh = 0;
    total_blength_2 = 0;
    (*new_lh1) = new_lh_value1;
    EXPECT_TRUE(merge_O_O_TwoLowers(seqregion1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, *new_lh1, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 3);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge2{0.99963572952385693,3.4820599267114684E-9,0.00031223631362821728,0.000052030680454886103};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge2);
    EXPECT_EQ(log_lh, -2.5901247759055703);
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    log_lh = 0;
    total_blength_2 = -1;
    SeqRegion::LHType new_lh_value2{1.1151877415789591228433182135137968771232408471406e-05,1.2436575683940661210420103588707597258578019250308e-10,0.99998884787385255989988763758447021245956420898438,1.2436575683940661210420103588707597258578019250308e-10};
    (*new_lh1) = new_lh_value2;
    EXPECT_TRUE(merge_O_O_TwoLowers(seqregion1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, *new_lh1, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 4);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge3{0.0000015928207737367849,1.3869404003893418E-15,0.99999840717922339,1.3869404003893418E-15};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge3);
    EXPECT_EQ(log_lh, -0.1335353759743724);
    // ----- Test 4 -----
    
    // ----- Test 5 -----
    log_lh = 0;
    total_blength_2 = 0;
    SeqRegion::LHType new_lh_value3{0.80000178432426516383912940000300295650959014892578,8.9216213262436184458807272856795123061601771041751e-06,0.19998037243308225408000566858390811830759048461914,8.9216213262436184458807272856795123061601771041751e-06};
    (*new_lh1) = new_lh_value3;
    EXPECT_TRUE(merge_O_O_TwoLowers(seqregion1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, *new_lh1, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 5);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge4{0.36361313457027911,3.1661424484350952E-10,0.63638686479649231,3.1661424484350952E-10};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge4);
    EXPECT_EQ(log_lh, -1.2911132491064328);
    // ----- Test 5 -----
    
    // ----- Test 6 -----
    log_lh = 0;
    total_blength_2 = 121e-5;
    SeqRegion::LHType new_lh_value4{0.99999721161091070786852696983260102570056915283203,2.7982277946960444848486499948455024505689081593118e-10,2.7878294439446026602288583595701254580490058287978e-06,2.7982277946960444848486499948455024505689081593118e-10};
    (*new_lh1) = new_lh_value4;
    EXPECT_TRUE(merge_O_O_TwoLowers(seqregion1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, *new_lh1, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 6);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge5{0.99998052979430418,2.8492225625088802E-13,0.000019470204442126495,9.6862198410969555E-13};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge5);
    EXPECT_EQ(log_lh, -2.0783443388603096);
    // ----- Test 6 -----
    
    // ----- Test 7 -----
    log_lh = 0;
    total_blength_2 = 121e-5;
    SeqRegion::LHType new_lh_value5{0.99999442305078967141440671184682287275791168212891,5.5758298847802999098784669518291678969035274349153e-06,5.5966272240860645628430121287970669397004996881151e-10,5.5966272240860645628430121287970669397004996881151e-10};
    *seqregion1.likelihood = new_lh_value5;
    EXPECT_TRUE(merge_O_O_TwoLowers(seqregion1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, *new_lh1, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 7);
    EXPECT_EQ(merged_regions->back().type, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -0.00043445905798766815652);
    // ----- Test 7 -----
    
    // ----- Test 8 -----
    log_lh = 0;
    total_blength_2 = -1;
    SeqRegion::LHType new_lh_value6{0.80000044606912734668213715849560685455799102783203,2.2303456366633140507259747825630213924341660458595e-06,2.2303456366633140507259747825630213924341660458595e-06,0.1999950932395993530299449503218056634068489074707};
    *seqregion1.likelihood = new_lh_value6;
    EXPECT_TRUE(merge_O_O_TwoLowers(seqregion1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, *new_lh1, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 8);
    EXPECT_EQ(merged_regions->back().type, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    SeqRegion::LHType new_lh_value_merge7{0.12684809781766063,1.3059895886789165E-10,0.87315190141833243,6.3340794416577515E-10};
    EXPECT_EQ(log_lh, -0.22314300087916599802);
    // ----- Test 8 -----
}

/*
 Test  merge_O_RACGT_TwoLowers(const SeqRegion& seq2_region, RealNumType total_blength_2, const PositionType end_pos, const Alignment& aln, const Model& model, const RealNumType threshold_prob, RealNumType &log_lh, SeqRegion::LHType& new_lh, RealNumType& sum_lh, SeqRegions* merged_regions, const bool return_log_lh)
 */
TEST(SeqRegions, merge_O_RACGT_TwoLowers)
{
    Alignment aln;
    Model model;
    Params params = Params::getInstance();
    
    // Init params, aln, and model
    initAlnModelParams(params, aln, model);
    // dummy variables
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    const RealNumType threshold_prob = params.threshold_prob;
    
    // ----- Test 1 -----
    SeqRegions* merged_regions = new SeqRegions();
    RealNumType log_lh = 0;
    RealNumType sum_lh = 0;
    SeqRegion seqregion1(TYPE_R, 243, -1, -1);
    auto new_lh1 = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value1{0.59998393993535814594508792652050033211708068847656,2.6766774402933635499746492514283602304203668609262e-05,2.6766774402933635499746492514283602304203668609262e-05,0.39996252651583585890904259940725751221179962158203};
    (*new_lh1) = new_lh_value1;
    RealNumType total_blength_2 = -1;
    const PositionType end_pos = 243;
    EXPECT_TRUE(merge_O_RACGT_TwoLowers(seqregion1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, *new_lh1, sum_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 1);
    EXPECT_EQ(merged_regions->back().type, TYPE_R);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -10.528349200671567);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    log_lh = 0;
    total_blength_2 = 1e-5;
    seqregion1.type = 1;
    (*new_lh1) = new_lh_value1;
    EXPECT_TRUE(merge_O_RACGT_TwoLowers(seqregion1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, *new_lh1, sum_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 2);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge1{0.011816164293335018,0.88299798641181038,8.0375755307632964E-7,0.10518504553730158};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge1);
    EXPECT_EQ(log_lh, -10.403932725076588);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    log_lh = 0;
    total_blength_2 = 0;
    seqregion1.type = 0;
    (*new_lh1) = new_lh_value1;
    EXPECT_TRUE(merge_O_RACGT_TwoLowers(seqregion1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, *new_lh1, sum_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 3);
    EXPECT_EQ(merged_regions->back().type, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -0.510852390898630326354634689778);
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    log_lh = 0;
    total_blength_2 = -1;
    seqregion1.type = 2;
    SeqRegion::LHType new_lh_value2{1.1151877415789591228433182135137968771232408471406e-05,1.2436575683940661210420103588707597258578019250308e-10,0.99998884787385255989988763758447021245956420898438,1.2436575683940661210420103588707597258578019250308e-10};
    (*new_lh1) = new_lh_value2;
    EXPECT_TRUE(merge_O_RACGT_TwoLowers(seqregion1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, *new_lh1, sum_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 4);
    EXPECT_EQ(merged_regions->back().type, 2);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -1.1152188332861237853011436571559755748239695094526e-05);
    // ----- Test 4 -----
    
    // ----- Test 5 -----
    log_lh = 0;
    total_blength_2 = 0;
    seqregion1.type = TYPE_R;
    SeqRegion::LHType new_lh_value3{0.80000178432426516383912940000300295650959014892578,8.9216213262436184458807272856795123061601771041751e-06,0.19998037243308225408000566858390811830759048461914,8.9216213262436184458807272856795123061601771041751e-06};
    (*new_lh1) = new_lh_value3;
    EXPECT_TRUE(merge_O_RACGT_TwoLowers(seqregion1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, *new_lh1, sum_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 5);
    EXPECT_EQ(merged_regions->back().type, TYPE_R);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -11.627032864857458349661101237870752811431884765625);
    // ----- Test 5 -----
    
    // ----- Test 6 -----
    log_lh = 0;
    total_blength_2 = 121e-5;
    seqregion1.type = 0;
    SeqRegion::LHType new_lh_value4{0.99999721161091070786852696983260102570056915283203,2.7982277946960444848486499948455024505689081593118e-10,2.7878294439446026602288583595701254580490058287978e-06,2.7982277946960444848486499948455024505689081593118e-10};
    (*new_lh1) = new_lh_value4;
    EXPECT_TRUE(merge_O_RACGT_TwoLowers(seqregion1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, *new_lh1, sum_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 6);
    EXPECT_EQ(merged_regions->back().type, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -0.00038188154788715934556328490678822618065169081091881);
    // ----- Test 6 -----
    
    // ----- Test 7 -----
    log_lh = 0;
    total_blength_2 = 121e-5;
    seqregion1.type = 3;
    EXPECT_TRUE(merge_O_RACGT_TwoLowers(seqregion1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, *new_lh1, sum_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 7);
    EXPECT_EQ(merged_regions->back().type, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -0.00028553809857495332426985390483764604141470044851303);
    // ----- Test 7 -----
    
    // ----- Test 8 -----
    log_lh = 0;
    total_blength_2 = -1;
    seqregion1.type = TYPE_R;
    EXPECT_TRUE(merge_O_RACGT_TwoLowers(seqregion1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, *new_lh1, sum_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 8);
    EXPECT_EQ(merged_regions->back().type, TYPE_R);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    SeqRegion::LHType new_lh_value_merge7{0.12684809781766063,1.3059895886789165E-10,0.87315190141833243,6.3340794416577515E-10};
    EXPECT_EQ(log_lh, -37.4289599328077855489027570001780986785888671875);
    // ----- Test 8 -----
}

/*
 Test merge_O_ORACGT_TwoLowers(const SeqRegion& seq1_region, const SeqRegion& seq2_region, RealNumType total_blength_1, RealNumType total_blength_2, const PositionType end_pos, const Alignment& aln, const Model& model, const RealNumType threshold_prob, RealNumType &log_lh, SeqRegions* merged_regions, const bool return_log_lh)
 */
TEST(SeqRegions, merge_O_ORACGT_TwoLowers)
{
    Alignment aln;
    Model model;
    Params params = Params::getInstance();
    
    // Init params, aln, and model
    initAlnModelParams(params, aln, model);
    // dummy variables
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    const RealNumType threshold_prob = params.threshold_prob;
    
    // ----- Test 1 -----
    SeqRegions* merged_regions = new SeqRegions();
    RealNumType log_lh = 0;
    auto new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value{0.12497560486006949187487435892762732692062854766846,9.7580559722090576469559833339140197949745925143361e-06,0.87500487902798607109389195102266967296600341796875,9.7580559722090576469559833339140197949745925143361e-06};
    (*new_lh) = new_lh_value;
    SeqRegion seqregion1(TYPE_O, 412, -1, -1, std::move(new_lh));
    auto new_lh1 = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value1{0.59998393993535814594508792652050033211708068847656,2.6766774402933635499746492514283602304203668609262e-05,2.6766774402933635499746492514283602304203668609262e-05,0.39996252651583585890904259940725751221179962158203};
    (*new_lh1) = new_lh_value1;
    SeqRegion seqregion2(TYPE_O, 412, -1, -1, std::move(new_lh1));
    RealNumType total_blength_1 = -1;
    RealNumType total_blength_2 = -1;
    const PositionType end_pos = 412;
    EXPECT_TRUE(merge_O_ORACGT_TwoLowers(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 1);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge{0.99963572952385693,3.4820599267114684E-9,0.00031223631362821728,0.000052030680454886103};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge);
    EXPECT_EQ(log_lh, -2.5901247759055703);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    log_lh = 0;
    total_blength_1 = -1;
    total_blength_2 = 32e-5;
    seqregion2.type = TYPE_R;
    EXPECT_TRUE(merge_O_ORACGT_TwoLowers(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 2);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge1{0.9993205977222771,2.4282817638871146E-9,0.00067939799763484592,1.851806273068163E-9};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge1);
    EXPECT_EQ(log_lh, -2.0790653485332990513256845588330179452896118164062);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    log_lh = 0;
    total_blength_1 = 143e-7;
    total_blength_2 = 0;
    SeqRegion::LHType new_lh_value2{2.7879382638950842936132173272012479969816922675818e-06,0.49999721206173608489820026079542003571987152099609,2.7879382638950842936132173272012479969816922675818e-06,0.49999721206173608489820026079542003571987152099609};
    (*seqregion1.likelihood) = new_lh_value2;
    seqregion2.type = 0;
    EXPECT_TRUE(merge_O_ORACGT_TwoLowers(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 3);
    EXPECT_EQ(merged_regions->back().type, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -12.484754332344007110577877028845250606536865234375);
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    log_lh = 0;
    total_blength_1 = 0;
    total_blength_2 = 12e-10;
    seqregion2.type = 3;
    EXPECT_TRUE(merge_O_ORACGT_TwoLowers(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 4);
    EXPECT_EQ(merged_regions->back().type, 3);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -0.69315275629224570863584631297271698713302612304688);
    // ----- Test 4 -----
    
    // ----- Test 5 -----
    log_lh = 0;
    total_blength_1 = 1e-5;
    total_blength_2 = 68e-4;
    SeqRegion::LHType new_lh_value3{8.3640013382402129011325767060647251582850003615022e-06,0.24998327199732350845096107150311581790447235107422,8.3640013382402129011325767060647251582850003615022e-06,0.74999999999999988897769753748434595763683319091797};
    (*seqregion1.likelihood) = new_lh_value3;
    seqregion2.type = TYPE_R;
    EXPECT_TRUE(merge_O_ORACGT_TwoLowers(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 5);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge5{0.0000090841953188726001,0.00016521455452187341,2.7570486522508116E-8,0.00037798930109716154};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge5);
    EXPECT_EQ(log_lh, 0.00055724163839606082743172166260592348407953977584839);
    // ----- Test 5 -----
    
    // ----- Test 6 -----
    log_lh = 0;
    total_blength_1 = 0;
    total_blength_2 = 15e-2;
    seqregion2.type = TYPE_O;
    SeqRegion::LHType new_lh_value6{0.39998126426090857554740409796067979186773300170898,1.3382670779607485874503763900733588343427982181311e-05,0.59999197039753215943136410714942030608654022216797,1.3382670779607485874503763900733588343427982181311e-05};
    (*seqregion2.likelihood) = new_lh_value6;
    EXPECT_TRUE(merge_O_ORACGT_TwoLowers(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 6);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge6{0.000099915947241675934,0.10965211286069129,0.00013202204674923939,0.890115949145318};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge6);
    EXPECT_EQ(log_lh, -3.40271550630699692874259199015796184539794921875);
    // ----- Test 6 -----
    
    // ----- Test 7 -----
    log_lh = 0;
    total_blength_1 = 0;
    total_blength_2 = 0;
    SeqRegion::LHType new_lh_value4{7.4346402191731947664294494204639818235591519623995e-06,0.33332094226630137878686355179524980485439300537109,7.4346402191731947664294494204639818235591519623995e-06,0.66666418845326025355291221785591915249824523925781};
    (*seqregion1.likelihood) = new_lh_value4;
    seqregion2.type = 1;
    EXPECT_TRUE(merge_O_ORACGT_TwoLowers(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 7);
    EXPECT_EQ(merged_regions->back().type, 1);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -1.0986494625601461727626428910298272967338562011719);
    // ----- Test 7 -----
    
    // ----- Test 8 -----
    log_lh = 0;
    total_blength_1 = 0;
    total_blength_2 = -1;
    seqregion2.type = TYPE_O;
    SeqRegion::LHType new_lh_value7{2.4874076016223502201974116939081453636628538106379e-10,7.4343638397288813108904244331132105116921593435109e-06,2.4874076016223502201974116939081453636628538106379e-10,0.99999256513867862405930964087019674479961395263672};
    (*seqregion2.likelihood) = new_lh_value7;
    EXPECT_TRUE(merge_O_ORACGT_TwoLowers(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 8);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge8{2.7739677142300221E-15,0.0000037170713771481666,2.7739677142300221E-15,0.9999962829286173};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge8);
    EXPECT_EQ(log_lh, -0.40547254324585230156330339923442807048559188842773);
    // ----- Test 8 -----
    
    // ----- Test 9 -----
    log_lh = 0;
    total_blength_1 = 437e-6;
    total_blength_2 = -1;
    SeqRegion::LHType new_lh_value5{0.49999442406129068761089229155913926661014556884766,5.5759387092817078802929955938516570768115343526006e-06,0.49999442406129068761089229155913926661014556884766,5.5759387092817078802929955938516570768115343526006e-06};
    (*seqregion1.likelihood) = new_lh_value5;
    seqregion2.type = TYPE_R;
    EXPECT_TRUE(merge_O_ORACGT_TwoLowers(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 9);
    EXPECT_EQ(merged_regions->back().type, TYPE_R);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -0.69321920665096659064374762238003313541412353515625);
    // ----- Test 9 -----
    
    // ----- Test 10 -----
    log_lh = 0;
    total_blength_1 = -1;
    total_blength_2 = 0;
    seqregion2.type = 2;
    EXPECT_TRUE(merge_O_ORACGT_TwoLowers(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 10);
    EXPECT_EQ(merged_regions->back().type, 2);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -0.69315833249954661443581471758079715073108673095703);
    // ----- Test 10 -----
}

/*
 Test merge_RACGT_O_TwoLowers(const SeqRegion& seq2_region, RealNumType total_blength_2, const PositionType end_pos, const Alignment& aln, const Model& model, const RealNumType threshold_prob, SeqRegion::LHType& new_lh, RealNumType &log_lh, SeqRegions* merged_regions, const bool return_log_lh)
 */
TEST(SeqRegions, merge_RACGT_O_TwoLowers)
{
    Alignment aln;
    Model model;
    Params params = Params::getInstance();
    
    // Init params, aln, and model
    initAlnModelParams(params, aln, model);
    // dummy variables
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    const RealNumType threshold_prob = params.threshold_prob;
    
    // ----- Test 1 -----
    SeqRegions* merged_regions = new SeqRegions();
    RealNumType log_lh = 0;
    auto new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value{0.12497560486006949187487435892762732692062854766846,9.7580559722090576469559833339140197949745925143361e-06,0.87500487902798607109389195102266967296600341796875,9.7580559722090576469559833339140197949745925143361e-06};
    (*new_lh) = new_lh_value;
    auto new_lh1 = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value1{0.59998393993535814594508792652050033211708068847656,2.6766774402933635499746492514283602304203668609262e-05,2.6766774402933635499746492514283602304203668609262e-05,0.39996252651583585890904259940725751221179962158203};
    (*new_lh1) = new_lh_value1;
    SeqRegion seqregion2(TYPE_O, 412, -1, -1, std::move(new_lh1));
    RealNumType total_blength_2 = -1;
    const PositionType end_pos = 412;
    EXPECT_TRUE(merge_RACGT_O_TwoLowers(seqregion2, total_blength_2, end_pos, aln, model, threshold_prob, *new_lh, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 1);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge{0.999635729523856930711644963594,3.48205992671146837419631824906e-09,0.000312236313628217279116106031012,5.20306804548861033901142880698e-05};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge);
    EXPECT_EQ(log_lh, -2.5901247759055703312469631782732903957366943359375);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    log_lh = 0;
    total_blength_2 = 1e-3;
    (*seqregion2.likelihood) = new_lh_value1;
    SeqRegion::LHType new_lh_value2{6.9955421312460765020272434505291649434188805400936e-11,6.9955421312460765020272434505291649434188805400936e-11,0.99999686351370775660996059741592034697532653808594,3.1363463814122085285697027340345854895531374495476e-06};
    (*new_lh) = new_lh_value2;
    EXPECT_TRUE(merge_RACGT_O_TwoLowers(seqregion2, total_blength_2, end_pos, aln, model, threshold_prob, *new_lh, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 2);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge2{9.02597705818903812311447803357e-08,9.66903454131082698149291275817e-11,0.997304647644765784875175995694,0.00269526199877322038961358074971};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge2);
    EXPECT_EQ(log_lh, -7.6737265825075411385114421136677265167236328125);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    log_lh = 0;
    total_blength_2 = 0;
    SeqRegion::LHType new_lh_value3{4.9563776809356656518181297177427779843128519132733e-06,4.9563776809356656518181297177427779843128519132733e-06,0.11109844481259316395505010177657823078334331512451,0.88889164243204510373885796070680953562259674072266};
    (*seqregion2.likelihood) = new_lh_value3;
    (*new_lh) = new_lh_value;
    EXPECT_TRUE(merge_RACGT_O_TwoLowers(seqregion2, total_blength_2, end_pos, aln, model, threshold_prob, *new_lh, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 3);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge3{6.37132369003884032904877487979e-06,4.9747095246365812296955045472e-10,0.999904410245811336999111063051,8.92179330276898546660951927478e-05};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge3);
    EXPECT_EQ(log_lh, -2.3307688028059012630421875655883923172950744628906);
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    log_lh = 0;
    total_blength_2 = 213e-7;
    (*seqregion2.likelihood) = new_lh_value1;
    SeqRegion::LHType new_lh_value4{2.6136903006998423677005767562508964374501374550164e-06,0.062492812351673081294745060176865081302821636199951,2.6136903006998423677005767562508964374501374550164e-06,0.93750196026772558699491355582722462713718414306641};
    (*new_lh) = new_lh_value4;
    EXPECT_TRUE(merge_RACGT_O_TwoLowers(seqregion2, total_blength_2, end_pos, aln, model, threshold_prob, *new_lh, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 4);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge4{4.18220729036683253488498879236e-06,6.6470825462486196933439668022e-06,2.51442378530120592451234655152e-10,0.999989170458721043921457294346};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge4);
    EXPECT_EQ(log_lh, -0.98093450214178112833707245954428799450397491455078);
    // ----- Test 4 -----
    
    // ----- Test 5 -----
    log_lh = 0;
    total_blength_2 = -1;
    SeqRegion::LHType new_lh_value5{9.9498152920206932387357935649403392619483099679201e-10,8.9219993723549547142565030455330088443588465452194e-05,9.9498152920206932387357935649403392619483099679201e-10,0.99991077801631333965559633725206367671489715576172};
    (*seqregion2.likelihood) = new_lh_value5;
    (*new_lh) = new_lh_value;
    EXPECT_TRUE(merge_RACGT_O_TwoLowers(seqregion2, total_blength_2, end_pos, aln, model, threshold_prob, *new_lh, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 5);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge5{1.27418556907770044148197133294e-05,8.92108976770603789461719368425e-05,8.92108976770603653936447807737e-05,0.999808836348955121131609757867};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge5);
    EXPECT_EQ(log_lh, -11.537315404593238454822312633041292428970336914062);
    // ----- Test 5 -----
    
    // ----- Test 6 -----
    log_lh = 0;
    total_blength_2 = 0;
    (*seqregion2.likelihood) = new_lh_value1;
    SeqRegion::LHType new_lh_value6{5.5759387092817078802929955938516570768115343526006e-06,0.49999442406129068761089229155913926661014556884766,5.5759387092817078802929955938516570768115343526006e-06,0.49999442406129068761089229155913926661014556884766};
    (*new_lh) = new_lh_value6;
    EXPECT_TRUE(merge_RACGT_O_TwoLowers(seqregion2, total_blength_2, end_pos, aln, model, threshold_prob, *new_lh, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 6);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge6{1.67277228426044177343658819757e-05,6.6917607757947737928509723826e-05,7.46265281119078188993421229841e-10,0.999916353923134160197605524445};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge6);
    EXPECT_EQ(log_lh, -1.6094591028973106450195018624071963131427764892578);
    // ----- Test 6 -----
    
    // ----- Test 7 -----
    log_lh = 0;
    total_blength_2 = 3163e-7;
    SeqRegion::LHType new_lh_value7{1.2436575683940661210420103588707597258578019250308e-10,1.1151877415789591228433182135137968771232408471406e-05,1.2436575683940661210420103588707597258578019250308e-10,0.99998884787385255989988763758447021245956420898438};
    (*seqregion2.likelihood) = new_lh_value7;
    (*new_lh) = new_lh_value;
    EXPECT_TRUE(merge_RACGT_O_TwoLowers(seqregion2, total_blength_2, end_pos, aln, model, threshold_prob, *new_lh, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 7);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge7{0.0166262913969323911089759349125,2.33062694270459372499967876102e-05,0.93180971485701069578766464474,0.0515406874766299802348434866417};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge7);
    EXPECT_EQ(log_lh, -8.5724436071309479956426002900116145610809326171875);
    // ----- Test 7 -----
    
    // ----- Test 8 -----
    log_lh = 0;
    total_blength_2 = 163e-10;
    (*seqregion2.likelihood) = new_lh_value1;
    SeqRegion::LHType new_lh_value8{0.74997769674260927885711680573876947164535522460938,0.24998327274350143345493791002809302881360054016113,1.9515256944661845116837164959555650511902058497071e-05,1.9515256944661845116837164959555650511902058497071e-05};
    (*new_lh) = new_lh_value8;
    EXPECT_TRUE(merge_RACGT_O_TwoLowers(seqregion2, total_blength_2, end_pos, aln, model, threshold_prob, *new_lh, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 8);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge8{0.999967777775584876209791218571,1.48753723892170687314673652168e-05,1.16113808009755283430647568722e-09,1.73456908878727982935719770241e-05};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge8);
    EXPECT_EQ(log_lh, -0.79853198337463082712162076859385706484317779541016);
    // ----- Test 8 -----
    
    // ----- Test 9 -----
    log_lh = 0;
    total_blength_2 = 0;
    SeqRegion::LHType new_lh_value9{0.99997490659690790870683940738672390580177307128906,2.5092097243268992381179730011275808010395849123597e-05,6.5292438309787796439726148030194621818544931102224e-10,6.5292438309787796439726148030194621818544931102224e-10};
    (*seqregion2.likelihood) = new_lh_value9;
    (*new_lh) = new_lh_value;
    EXPECT_TRUE(merge_RACGT_O_TwoLowers(seqregion2, total_blength_2, end_pos, aln, model, threshold_prob, *new_lh, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 9);
    EXPECT_EQ(merged_regions->back().type, TYPE_R);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -2.0796618090328400363375749293481931090354919433594);
    // ----- Test 9 -----
}

/*
 Test merge_RACGT_RACGT_TwoLowers(const SeqRegion& seq2_region, RealNumType total_blength_2, const PositionType end_pos, const Alignment& aln, const Model& model, const RealNumType threshold_prob, SeqRegion::LHType& new_lh, RealNumType& sum_lh, RealNumType &log_lh, SeqRegions* merged_regions, const bool return_log_lh)
 */
TEST(SeqRegions, merge_RACGT_RACGT_TwoLowers)
{
    Alignment aln;
    Model model;
    Params params = Params::getInstance();
    
    // Init params, aln, and model
    initAlnModelParams(params, aln, model);
    // dummy variables
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    const RealNumType threshold_prob = params.threshold_prob;
    
    // ----- Test 1 -----
    SeqRegions* merged_regions = new SeqRegions();
    RealNumType log_lh = 0;
    RealNumType sum_lh = 0;
    auto new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value{9.7580559722090576469559833339140197949745925143361e-06,9.7580559722090576469559833339140197949745925143361e-06,0.87500487902798607109389195102266967296600341796875,0.12497560486006949187487435892762732692062854766846};
    (*new_lh) = new_lh_value;
    SeqRegion seqregion2(TYPE_R, 412, -1, -1);
    RealNumType total_blength_2 = -1;
    const PositionType end_pos = 412;
    EXPECT_TRUE(merge_RACGT_RACGT_TwoLowers(seqregion2, total_blength_2, end_pos, aln, model, threshold_prob, *new_lh, sum_lh, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 1);
    EXPECT_EQ(merged_regions->back().type, seqregion2.type);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -11.537417360554178102916011994238942861557006835938);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    log_lh = 0;
    sum_lh = 0;
    total_blength_2 = 3.1363463814122085285697027340345854895531374495476e-06;
    seqregion2.type = 0;
    SeqRegion::LHType new_lh_value2{6.9955421312460765020272434505291649434188805400936e-11,3.1363463814122085285697027340345854895531374495476e-06,0.99999686351370775660996059741592034697532653808594,6.9955421312460765020272434505291649434188805400936e-11};
    (*new_lh) = new_lh_value2;
    EXPECT_TRUE(merge_RACGT_RACGT_TwoLowers(seqregion2, total_blength_2, end_pos, aln, model, threshold_prob, *new_lh, sum_lh, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 2);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge2{7.35069993527212464743542108536e-05,1.00511327029477612984332841189e-06,0.999925487870280349511631357018,1.70965639624293546874861604245e-11};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge2);
    EXPECT_EQ(log_lh, -13.865034049156509610156717826612293720245361328125);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    log_lh = 0;
    sum_lh = 0;
    total_blength_2 = 0;
    seqregion2.type = TYPE_R;
    (*new_lh) = new_lh_value;
    EXPECT_TRUE(merge_RACGT_RACGT_TwoLowers(seqregion2, total_blength_2, end_pos, aln, model, threshold_prob, *new_lh, sum_lh, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 3);
    EXPECT_EQ(merged_regions->back().type, seqregion2.type);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -11.537417360554178102916011994238942861557006835938);
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    log_lh = 0;
    sum_lh = 0;
    total_blength_2 = 0.062492812351673081294745060176865081302821636199951;
    seqregion2.type = 1;
    SeqRegion::LHType new_lh_value4{0.062492812351673081294745060176865081302821636199951,0.93750196026772558699491355582722462713718414306641,2.6136903006998423677005767562508964374501374550164e-06,2.6136903006998423677005767562508964374501374550164e-06};
    (*new_lh) = new_lh_value4;
    EXPECT_TRUE(merge_RACGT_RACGT_TwoLowers(seqregion2, total_blength_2, end_pos, aln, model, threshold_prob, *new_lh, sum_lh, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 4);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge4{0.000276015672152719990402325311862,0.9997238125720525614426037464,1.76015102091620514116593592541e-08,1.54154284403481275470771598261e-07};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge4);
    EXPECT_EQ(log_lh, -0.16879625050215649184615074318571714684367179870605);
    // ----- Test 4 -----
    
    // ----- Test 5 -----
    log_lh = 0;
    sum_lh = 0;
    total_blength_2 = -1;
    seqregion2.type = 0;
    (*new_lh) = new_lh_value;
    EXPECT_TRUE(merge_RACGT_RACGT_TwoLowers(seqregion2, total_blength_2, end_pos, aln, model, threshold_prob, *new_lh, sum_lh, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 5);
    EXPECT_EQ(merged_regions->back().type, seqregion2.type);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -11.537417360554178102916011994238942861557006835938);
    // ----- Test 5 -----
    
    // ----- Test 6 -----
    log_lh = 0;
    sum_lh = 0;
    total_blength_2 = 0.49999442406129068761089229155913926661014556884766;
    seqregion2.type = 3;
    SeqRegion::LHType new_lh_value6{0.49999442406129068761089229155913926661014556884766,0.49999442406129068761089229155913926661014556884766,5.5759387092817078802929955938516570768115343526006e-06,5.5759387092817078802929955938516570768115343526006e-06};
    (*new_lh) = new_lh_value6;
    EXPECT_TRUE(merge_RACGT_RACGT_TwoLowers(seqregion2, total_blength_2, end_pos, aln, model, threshold_prob, *new_lh, sum_lh, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 6);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge6{0.0540225020296241503769962832848,0.945967079514273279094993540639,4.82257467713597815899819951091e-06,5.59588142544069931037494305959e-06};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge6);
    EXPECT_EQ(log_lh, -0.9987216763314137324414332397282123565673828125);
    // ----- Test 6 -----
    
    // ----- Test 7 -----
    log_lh = 0;
    sum_lh = 0;
    total_blength_2 = 1.2436575683940661210420103588707597258578019250308e-10;
    seqregion2.type = 0;
    (*new_lh) = new_lh_value;
    EXPECT_TRUE(merge_RACGT_RACGT_TwoLowers(seqregion2, total_blength_2, end_pos, aln, model, threshold_prob, *new_lh, sum_lh, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 7);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge7{0.99999649823344849419726187989,1.20936802836750174748186692148e-11,3.38363623651822093845354007258e-06,1.18118221219596752754965953997e-07};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge7);
    EXPECT_EQ(log_lh, -11.537413858823567736067161604296416044235229492188);
    // ----- Test 7 -----
    
    // ----- Test 8 -----
    log_lh = 0;
    sum_lh = 0;
    total_blength_2 = 1.9515256944661845116837164959555650511902058497071e-05;
    seqregion2.type = TYPE_R;
    SeqRegion::LHType new_lh_value8{0.74997769674260927885711680573876947164535522460938,1.9515256944661845116837164959555650511902058497071e-05,1.9515256944661845116837164959555650511902058497071e-05,0.24998327274350143345493791002809302881360054016113};
    (*new_lh) = new_lh_value8;
    EXPECT_TRUE(merge_RACGT_RACGT_TwoLowers(seqregion2, total_blength_2, end_pos, aln, model, threshold_prob, *new_lh, sum_lh, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 8);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge8{0.999999517409988936122999803047,4.93812459495446616345023798386e-11,1.54077796167036529766850756536e-10,4.82386551824058610632602311225e-07};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge8);
    EXPECT_EQ(log_lh, -0.28771792989182387589863765242625959217548370361328);
    // ----- Test 8 -----
    
    // ----- Test 9 -----
    log_lh = 0;
    sum_lh = 0;
    total_blength_2 = 0;
    seqregion2.type = 1;
    (*new_lh) = new_lh_value;
    EXPECT_TRUE(merge_RACGT_RACGT_TwoLowers(seqregion2, total_blength_2, end_pos, aln, model, threshold_prob, *new_lh, sum_lh, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 9);
    EXPECT_EQ(merged_regions->back().type, seqregion2.type);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -11.537417360554178102916011994238942861557006835938);
    // ----- Test 9 -----
    
    // ----- Test 10 -----
    log_lh = 0;
    sum_lh = 0;
    total_blength_2 = -1;
    seqregion2.type = 2;
    (*new_lh) = new_lh_value;
    EXPECT_TRUE(merge_RACGT_RACGT_TwoLowers(seqregion2, total_blength_2, end_pos, aln, model, threshold_prob, *new_lh, sum_lh, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 10);
    EXPECT_EQ(merged_regions->back().type, seqregion2.type);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -0.13352581660808454655509081021591555327177047729492);
    // ----- Test 10 -----
    
    // ----- Test 11 -----
    log_lh = 0;
    sum_lh = 0;
    total_blength_2 = 6.5292438309787796439726148030194621818544931102224e-10;
    seqregion2.type = TYPE_R;
    SeqRegion::LHType new_lh_value11{2.5092097243268992381179730011275808010395849123597e-05,0.99997490659690790870683940738672390580177307128906,6.5292438309787796439726148030194621818544931102224e-10,6.5292438309787796439726148030194621818544931102224e-10};
    (*new_lh) = new_lh_value11;
    EXPECT_TRUE(merge_RACGT_RACGT_TwoLowers(seqregion2, total_blength_2, end_pos, aln, model, threshold_prob, *new_lh, sum_lh, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 11);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge11{0.999997469693585716576933464239,2.53030640795433316882884731969e-06,5.15495543614262290103202928001e-15,1.25992117026322114735542991217e-15};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge11);
    EXPECT_EQ(log_lh, -10.592955082179518200291568064130842685699462890625);
    // ----- Test 11 -----
}


/*
 Test merge_RACGT_ORACGT_TwoLowers(const SeqRegion& seq1_region, const SeqRegion& seq2_region, RealNumType total_blength_1, RealNumType total_blength_2, const PositionType end_pos, const Alignment& aln, const Model& model, const RealNumType threshold_prob, RealNumType &log_lh, SeqRegions* merged_regions, const bool return_log_lh)
 */
TEST(SeqRegions, merge_RACGT_ORACGT_TwoLowers)
{
    Alignment aln;
    Model model;
    Params params = Params::getInstance();
    
    // Init params, aln, and model
    initAlnModelParams(params, aln, model);
    // dummy variables
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    const RealNumType threshold_prob = params.threshold_prob;
    
    // ----- Test 1 -----
    SeqRegions* merged_regions = new SeqRegions();
    RealNumType log_lh = 0;
    SeqRegion seqregion1(TYPE_R, 43223);
    auto new_lh1 = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value1{2.6766774402933635499746492514283602304203668609262e-05,0.39996252651583585890904259940725751221179962158203,2.6766774402933635499746492514283602304203668609262e-05,0.59998393993535814594508792652050033211708068847656};
    (*new_lh1) = new_lh_value1;
    SeqRegion seqregion2(TYPE_O, 43223, -1, -1, std::move(new_lh1));
    RealNumType total_blength_1 = 1326e-5;
    RealNumType total_blength_2 = -1;
    const PositionType end_pos = 43223;
    EXPECT_TRUE(merge_RACGT_ORACGT_TwoLowers(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 1);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    EXPECT_EQ(log_lh, -6.7833577173188217557253665290772914886474609375);
    SeqRegion::LHType new_lh_value_merge1{0.0235298053213485250378944613203,0.455403998389659003809271098362,9.50936657505749543661116574e-05,0.520971102623241977269685776264};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge1);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    log_lh = 0;
    total_blength_1 = -1;
    total_blength_2 = 32e-5;
    seqregion1.type = 2;
    seqregion2.type = TYPE_R;
    EXPECT_TRUE(merge_RACGT_ORACGT_TwoLowers(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 2);
    EXPECT_EQ(merged_regions->back().type, 2);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -9.23984298142250537466679816134274005889892578125);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    log_lh = 0;
    total_blength_1 = 143e-7;
    total_blength_2 = 0;
    seqregion1.type = TYPE_R;
    seqregion2.type = 1;
    EXPECT_TRUE(merge_RACGT_ORACGT_TwoLowers(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 3);
    EXPECT_EQ(merged_regions->back().type, 1);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -13.485791369259754191034517134539783000946044921875);
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    log_lh = 0;
    total_blength_1 = 0;
    total_blength_2 = 12e-10;
    seqregion1.type = 0;
    seqregion2.type = 3;
    EXPECT_TRUE(merge_RACGT_ORACGT_TwoLowers(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 4);
    EXPECT_EQ(merged_regions->back().type, TYPE_R);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -23.07170390644744628616535919718444347381591796875);
    // ----- Test 4 -----
    
    // ----- Test 5 -----
    log_lh = 0;
    total_blength_1 = 1e-5;
    total_blength_2 = 68e-4;
    seqregion1.type = 1;
    seqregion2.type = TYPE_R;
    EXPECT_TRUE(merge_RACGT_ORACGT_TwoLowers(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 5);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge5{0.000899932704475660738205333721851,0.999091155743592085336501895654,2.83758240671473494537067026877e-06,6.0739695257571569122004295771e-06};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge5);
    EXPECT_EQ(log_lh, -7.3204796410357877434194051602389663457870483398438);
    // ----- Test 5 -----
    
    // ----- Test 6 -----
    log_lh = 0;
    total_blength_1 = 0;
    total_blength_2 = 15e-2;
    seqregion1.type = 2;
    seqregion2.type = TYPE_O;
    SeqRegion::LHType new_lh_value6{1.3382670779607485874503763900733588343427982181311e-05,1.3382670779607485874503763900733588343427982181311e-05,0.39998126426090857554740409796067979186773300170898,0.59999197039753215943136410714942030608654022216797};
    (*seqregion2.likelihood) = new_lh_value6;
    EXPECT_TRUE(merge_RACGT_ORACGT_TwoLowers(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 6);
    EXPECT_EQ(merged_regions->back().type, 2);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -0.92777401003210746566196576168294996023178100585938);
 
    // ----- Test 6 -----
    
    // ----- Test 7 -----
    log_lh = 0;
    total_blength_1 = 0.13121;
    total_blength_2 = 0;
    seqregion1.type = TYPE_R;
    seqregion2.type = 1;
    EXPECT_TRUE(merge_RACGT_ORACGT_TwoLowers(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 7);
    EXPECT_EQ(merged_regions->back().type, 1);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -4.3614965344258544988065295910928398370742797851562);
    // ----- Test 7 -----
    
    // ----- Test 8 -----
    log_lh = 0;
    total_blength_1 = 0;
    total_blength_2 = -1;
    seqregion1.type = 3;
    seqregion2.type = TYPE_O;
    SeqRegion::LHType new_lh_value7{0.99999256513867862405930964087019674479961395263672,7.4343638397288813108904244331132105116921593435109e-06,2.4874076016223502201974116939081453636628538106379e-10,2.4874076016223502201974116939081453636628538106379e-10};
    (*seqregion2.likelihood) = new_lh_value7;
    EXPECT_TRUE(merge_RACGT_ORACGT_TwoLowers(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 8);
    EXPECT_EQ(merged_regions->back().type, 3);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -22.114609885656182797220026259310543537139892578125);
    // ----- Test 8 -----
    
    // ----- Test 9 -----
    log_lh = 0;
    total_blength_1 = 437e-6;
    total_blength_2 = -1;
    seqregion1.type = 3;
    seqregion2.type = TYPE_R;
    EXPECT_TRUE(merge_RACGT_ORACGT_TwoLowers(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 9);
    EXPECT_EQ(merged_regions->back().type, TYPE_R);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -10.266336989163672654967740527354180812835693359375);
    // ----- Test 9 -----
    
    // ----- Test 10 -----
    log_lh = 0;
    total_blength_1 = 0;
    total_blength_2 = 42.165e-8;
    seqregion1.type = 2;
    seqregion2.type = 3;
    EXPECT_TRUE(merge_RACGT_ORACGT_TwoLowers(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 10);
    EXPECT_EQ(merged_regions->back().type, 2);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -15.129806325448386772336561989504843950271606445312);
    // ----- Test 10 -----
}

/*
 Test merge_notN_notN_TwoLowers(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType plength1, const RealNumType plength2, const PositionType end_pos, const PositionType pos, const Alignment& aln, const Model& model, const RealNumType threshold_prob, RealNumType &log_lh, SeqRegions* merged_regions, const bool return_log_lh)
 */
TEST(SeqRegions, merge_notN_notN_TwoLowers)
{
    Alignment aln;
    Model model;
    Params params = Params::getInstance();
    
    // Init params, aln, and model
    initAlnModelParams(params, aln, model);
    // dummy variables
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    const RealNumType threshold_prob = params.threshold_prob;
    
    // ----- Test 1 -----
    SeqRegions* merged_regions = new SeqRegions();
    RealNumType log_lh = 0;
    auto new_lh = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value{0.42855939517854246822992081433767452836036682128906,0.57141073458160607234646022334345616400241851806641,1.4935119925781397759922616841343767646321794018149e-05,1.4935119925781397759922616841343767646321794018149e-05};
    (*new_lh) = new_lh_value;
    SeqRegion seqregion1(TYPE_O, 1412, -1, -1, std::move(new_lh));
    auto new_lh1 = std::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value1{0.61110448521316429459915298139094375073909759521484,1.1926616304335032169438128579752600444408017210662e-05,1.1926616304335032169438128579752600444408017210662e-05,0.38887166155422708824218602785549592226743698120117};
    (*new_lh1) = new_lh_value1;
    SeqRegion seqregion2(TYPE_O, 1412, -1, -1, std::move(new_lh1));
    RealNumType plength1 = -1;
    RealNumType plength2 = -1;
    const PositionType end_pos = 1412;
    const PositionType pos = 1355;
    EXPECT_TRUE(merge_notN_notN_TwoLowers(seqregion1, seqregion2, plength1, plength2, end_pos, pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 1);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    EXPECT_EQ(log_lh, -1.3397650685348243548844493489013984799385070800781);
    SeqRegion::LHType new_lh_value_merge1{0.999951803463153265916218970233,2.60206546527806729897022708364e-05,6.80109025377646455925279108979e-10,2.21752020848110232805524416611e-05};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge1);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    log_lh = 0;
    plength1 = -1;
    plength2 = 1.3382670779607485874503763900733588343427982181311e-05;
    seqregion1.type = TYPE_R;
    seqregion2.type = TYPE_R;
    EXPECT_TRUE(merge_notN_notN_TwoLowers(seqregion1, seqregion2, plength1, plength2, end_pos, pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 2);
    EXPECT_EQ(merged_regions->back().type, TYPE_R);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -0.0007633997783012927721910112488501454208744689822197);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    log_lh = 0;
    plength1 = 3.4820599267114684E-9;
    plength2 = 3.283e-9;
    seqregion1.type = 2;
    seqregion2.type = 0;
    EXPECT_TRUE(merge_notN_notN_TwoLowers(seqregion1, seqregion2, plength1, plength2, end_pos, pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 3);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    EXPECT_EQ(log_lh, -20.199112063767703517669360735453665256500244140625);
    SeqRegion::LHType new_lh_value_merge3{0.410245869594259182644435668408,6.40012426019357435068582000991e-11,0.589754130146332378181739386491,1.95407241773354552292205525266e-10};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge3);
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    log_lh = 0;
    plength1 = 0;
    plength2 = 0.00031223631362821728;
    seqregion1.type = TYPE_R;
    seqregion2.type = 3;
    EXPECT_TRUE(merge_notN_notN_TwoLowers(seqregion1, seqregion2, plength1, plength2, end_pos, pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 4);
    EXPECT_EQ(merged_regions->back().type, TYPE_R);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -8.5224663158967093323781227809377014636993408203125);
    // ----- Test 4 -----
    
    // ----- Test 5 -----
    log_lh = 0;
    plength1 = 1e-5;
    plength2 = 0.000052030680454886103;
    seqregion1.type = 3;
    SeqRegion::LHType new_lh_value3{4.9383548301647924830444502664050787643645890057087e-05,0.71422131017848822231997019116533920168876647949219,4.9383548301647924830444502664050787643645890057087e-05,0.28567992272490849714472460618708282709121704101562};
    (*seqregion2.likelihood) = new_lh_value3;
    seqregion2.type = TYPE_O;
    EXPECT_TRUE(merge_notN_notN_TwoLowers(seqregion1, seqregion2, plength1, plength2, end_pos, pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 5);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    EXPECT_EQ(log_lh, -1.2528228994081356262313420302234590053558349609375);
    SeqRegion::LHType new_lh_value_merge5{1.47064689584769381918621458095e-10,3.48425531256729660231630241185e-05,1.38799271102831098152142401229e-09,0.99996515591181689419641998029};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge5);
    // ----- Test 5 -----
    
    // ----- Test 6 -----
    log_lh = 0;
    plength1 = 0;
    plength2 = 8.3640013382402129011325767060647251582850003615022e-06;
    seqregion1.type = 1;
    seqregion2.type = TYPE_O;
    SeqRegion::LHType new_lh_value6{6.6911562987199814600764238847752096717158565297723e-06,0.39999063238118182095348629445652477443218231201172,6.6911562987199814600764238847752096717158565297723e-06,0.59999598530622066938633452082285657525062561035156};
    (*seqregion2.likelihood) = new_lh_value6;
    EXPECT_TRUE(merge_notN_notN_TwoLowers(seqregion1, seqregion2, plength1, plength2, end_pos, pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 6);
    EXPECT_EQ(merged_regions->back().type, 1);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -0.91630994861674031071174795215483754873275756835938);
    // ----- Test 6 -----
    
    // ----- Test 7 -----
    log_lh = 0;
    plength1 = 0;
    plength2 = 0;
    SeqRegion::LHType new_lh_value4{4.430579981880772584663374230178360321796837695274e-10,0.99999581732016740165391865957644768059253692626953,4.430579981880772584663374230178360321796837695274e-10,4.181793716486326507387246559366289488934853579849e-06};
    (*seqregion1.likelihood) = new_lh_value4;
    seqregion1.type = TYPE_O;
    seqregion2.type = 1;
    EXPECT_TRUE(merge_notN_notN_TwoLowers(seqregion1, seqregion2, plength1, plength2, end_pos, pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 7);
    EXPECT_EQ(merged_regions->back().type, 1);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -4.1826885800280286031383250588966404848179081454873e-06);
    // ----- Test 7 -----
    
    // ----- Test 8 -----
    log_lh = 0;
    plength1 = 0;
    plength2 = -1;
    seqregion1.type = TYPE_R;
    seqregion2.type = TYPE_O;
    SeqRegion::LHType new_lh_value7{1.1815113248226830777267349123884135342343881802663e-09,0.9999702585027225865133004845120012760162353515625,1.1815113248226830777267349123884135342343881802663e-09,2.9739134254674524059092188821296076639555394649506e-05};
    (*seqregion2.likelihood) = new_lh_value7;
    EXPECT_TRUE(merge_notN_notN_TwoLowers(seqregion1, seqregion2, plength1, plength2, end_pos, pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 8);
    EXPECT_EQ(merged_regions->back().type, TYPE_R);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -20.556471434224643957122680149041116237640380859375);
    // ----- Test 8 -----
    
    // ----- Test 9 -----
    log_lh = 0;
    plength1 = 3.4820599267114684E-9;
    plength2 = -1;
    SeqRegion::LHType new_lh_value5{6.6911562987199814600764238847752096717158565297723e-06,0.39999063238118182095348629445652477443218231201172,6.6911562987199814600764238847752096717158565297723e-06,0.59999598530622066938633452082285657525062561035156};
    (*seqregion1.likelihood) = new_lh_value5;
    seqregion1.type = TYPE_O;
    seqregion2.type = TYPE_R;
    EXPECT_TRUE(merge_notN_notN_TwoLowers(seqregion1, seqregion2, plength1, plength2, end_pos, pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 8); // consecutive R regions -> merged
    EXPECT_EQ(merged_regions->back().type, TYPE_R);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -11.914505989884778713872037769760936498641967773438);
    // ----- Test 9 -----
    
    // ----- Test 10 -----
    log_lh = 0;
    plength1 = -1;
    plength2 = 231.65e-8;
    seqregion1.type = 2;
    seqregion2.type = 2;
    EXPECT_TRUE(merge_notN_notN_TwoLowers(seqregion1, seqregion2, plength1, plength2, end_pos, pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 9);
    EXPECT_EQ(merged_regions->back().type, 2);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -2.389727728309945835126872418219257099281094269827e-06);
    // ----- Test 10 -----
    
    // ----- Test 11 -----
    log_lh = 0;
    plength1 = 5.5759387092817078802929955938516570768115343526006e-06;
    plength2 = 0;
    SeqRegion::LHType new_lh_value11{4.181793716486327354420193813666628557257354259491e-06,0.99999581732016751267622112209210172295570373535156,4.430579981880772584663374230178360321796837695274e-10,4.430579981880772584663374230178360321796837695274e-10};
    (*seqregion1.likelihood) = new_lh_value11;
    seqregion1.type = TYPE_O;
    seqregion2.type = TYPE_R;
    EXPECT_TRUE(merge_notN_notN_TwoLowers(seqregion1, seqregion2, plength1, plength2, end_pos, pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 10);
    EXPECT_EQ(merged_regions->back().type, TYPE_R);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -14.492793795180125115962255222257226705551147460938);
    // ----- Test 11 -----
    
    // ----- Test 12 -----
    log_lh = 0;
    plength1 = -1;
    plength2 = 3.4820599267114684E-5;
    seqregion1.type = 3;
    seqregion2.type = 1;
    EXPECT_TRUE(merge_notN_notN_TwoLowers(seqregion1, seqregion2, plength1, plength2, end_pos, pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 11);
    EXPECT_EQ(merged_regions->back().type, 3);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -10.491958394947619837012098287232220172882080078125);
    // ----- Test 12 -----
    
    // ----- Test 13 -----
    log_lh = 0;
    plength1 = 143e-7;
    plength2 = 2.7879382638950842936132173272012479969816922675818e-06;
    SeqRegion::LHType new_lh_value13{0.33332094226630137878686355179524980485439300537109,7.4346402191731947664294494204639818235591519623995e-06,0.66666418845326025355291221785591915249824523925781,7.4346402191731947664294494204639818235591519623995e-06};
    (*seqregion1.likelihood) = new_lh_value13;
    seqregion1.type = TYPE_O;
    seqregion2.type = 0;
    EXPECT_TRUE(merge_notN_notN_TwoLowers(seqregion1, seqregion2, plength1, plength2, end_pos, pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 12);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    EXPECT_EQ(log_lh, 0.28767119359422399504921941115753725171089172363281);
    SeqRegion::LHType new_lh_value_merge13{0.249993404024937859730925993063,1.7944514993369046267887038812e-12,4.22949001176853426491232144577e-07,1.78315422997564743574952224173e-12};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge13);
    // ----- Test 13 -----
    
    // ----- Test 14 -----
    log_lh = 0;
    plength1 = 0;
    plength2 = 0.31223631362821728e-10;
    seqregion1.type = 1;
    seqregion2.type = 3;
    EXPECT_TRUE(merge_notN_notN_TwoLowers(seqregion1, seqregion2, plength1, plength2, end_pos, pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 13);
    EXPECT_EQ(merged_regions->back().type, 1);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -23.857798413868291476092053926549851894378662109375);
    // ----- Test 14 -----
    
    // ----- Test 15 -----
    log_lh = 0;
    plength1 = 8.364001338240212901132e-5;
    plength2 = 0.2168e-4;
    seqregion1.type = TYPE_R;
    SeqRegion::LHType new_lh_value15{9.9493714778138372440088605107603655919312757305306e-10,0.99996654175612176285170562550774775445461273193359,9.9493714778138372440088605107603655919312757305306e-10,3.345625400403152548828300538730218249838799238205e-05};
    (*seqregion2.likelihood) = new_lh_value15;
    seqregion2.type = TYPE_O;
    EXPECT_TRUE(merge_notN_notN_TwoLowers(seqregion1, seqregion2, plength1, plength2, end_pos, pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 14);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    EXPECT_EQ(log_lh, -11.502067070054362574182960088364779949188232421875);
    SeqRegion::LHType new_lh_value_merge15{2.13259105287923615202651250744e-06,0.804503715022567122971963726741,0.195330717008421000935314282287,0.000163435377959033966777449564667};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge15);
    // ----- Test 15 -----
    
    // ----- Test 16 -----
    log_lh = 0;
    plength1 = 0;
    plength2 = 1.09085755475e-2;
    seqregion1.type = 2;
    seqregion2.type = TYPE_O;
    SeqRegion::LHType new_lh_value16{1.5333235876614054963434918832376752106938511133194e-05,0.99996236262065618660699328756891191005706787109375,1.1152071733605872690944446623539931806590175256133e-05,1.1152071733605872690944446623539931806590175256133e-05};
    (*seqregion2.likelihood) = new_lh_value16;
    EXPECT_TRUE(merge_notN_notN_TwoLowers(seqregion1, seqregion2, plength1, plength2, end_pos, pos, aln, model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 15);
    EXPECT_EQ(merged_regions->back().type, TYPE_R);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -6.903698068549179112096680910326540470123291015625);
    // ----- Test 16 -----
}

/*
 Test mergeTwoLowers(SeqRegions* &merged_regions, RealNumType plength1, const SeqRegions* const regions2, RealNumType plength2, const Alignment& aln, const Model& model, RealNumType threshold_prob, bool return_log_lh)
 */
TEST(SeqRegions, mergeTwoLowers)
{
    Alignment aln;
    Model model;
    Params params = Params::getInstance();
    
    // Init params, aln, and model
    initAlnModelParams(params, aln, model, "JC");
    // dummy variables
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    const RealNumType threshold_prob = params.threshold_prob;
    
    // ----- Test 1 -----
    SeqRegions seqregions1, seqregions2, output_regions;
    SeqRegions* merged_regions_ptr = nullptr;
    
    genTestData4(seqregions1, seqregions2, output_regions, 1);
    seqregions1.mergeTwoLowers(merged_regions_ptr, 3.3454886086112878360986772063867533688608091324568e-05, seqregions2, -1, aln, model, threshold_prob);
    EXPECT_EQ(*merged_regions_ptr, output_regions);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    genTestData4(seqregions1, seqregions2, output_regions, 2);
    seqregions1.mergeTwoLowers(merged_regions_ptr, 0.00023418420260279014175064382641267002327367663383484, seqregions2, 3.3454886086112878360986772063867533688608091324568e-05, aln, model, threshold_prob);
    EXPECT_EQ(*merged_regions_ptr, output_regions);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    genTestData4(seqregions1, seqregions2, output_regions, 3);
    seqregions1.mergeTwoLowers(merged_regions_ptr, 0.00020072931651667725661339347631439977703848853707314, seqregions2, 0.00013381954434445151344394708825547013475443236529827, aln, model, threshold_prob);
    EXPECT_EQ(*merged_regions_ptr, output_regions);
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    genTestData4(seqregions1, seqregions2, output_regions, 4);
    seqregions1.mergeTwoLowers(merged_regions_ptr, 3.3454886086112878360986772063867533688608091324568e-05, seqregions2, -1, aln, model, threshold_prob);
    EXPECT_EQ(*merged_regions_ptr, output_regions);
    // ----- Test 4 -----
    
    // ----- Test 5 -----
    genTestData4(seqregions1, seqregions2, output_regions, 5);
    seqregions1.mergeTwoLowers(merged_regions_ptr, -1, seqregions2, 5.0182329129169314153348369078599944259622134268284e-05, aln, model, threshold_prob);
    EXPECT_EQ(*merged_regions_ptr, output_regions);
    // ----- Test 5 -----
    
    // ----- Test 6 -----
    genTestData4(seqregions1, seqregions2, output_regions, 6);
    seqregions1.mergeTwoLowers(merged_regions_ptr, 3.3454886086112878360986772063867533688608091324568e-05, seqregions2, -1, aln, model, threshold_prob);
    EXPECT_EQ(*merged_regions_ptr, output_regions);
    // ----- Test 6 -----
    
    // ----- Test 7 -----
    genTestData4(seqregions1, seqregions2, output_regions, 7);
    seqregions1.mergeTwoLowers(merged_regions_ptr, -1, seqregions2, -1, aln, model, threshold_prob);
    EXPECT_EQ(*merged_regions_ptr, output_regions);
    // ----- Test 7 -----
    
    // ----- Test 8 -----
    genTestData4(seqregions1, seqregions2, output_regions, 8);
    seqregions1.mergeTwoLowers(merged_regions_ptr, 3.3454886086112878360986772063867533688608091324568e-05, seqregions2, 3.3454886086112878360986772063867533688608091324568e-05, aln, model, threshold_prob);
    EXPECT_EQ(*merged_regions_ptr, output_regions);
    // ----- Test 8 -----
    
    // ----- Test 9 -----
    genTestData4(seqregions1, seqregions2, output_regions, 9);
    seqregions1.mergeTwoLowers(merged_regions_ptr, 3.3454886086112878360986772063867533688608091324568e-05, seqregions2, 6.6909772172225756721973544127735067377216182649136e-05, aln, model, threshold_prob);
    EXPECT_EQ(*merged_regions_ptr, output_regions);
    // ----- Test 9 -----
    
    // ----- Test 10 -----
    genTestData4(seqregions1, seqregions2, output_regions, 10);
    seqregions1.mergeTwoLowers(merged_regions_ptr, -1, seqregions2, -1, aln, model, threshold_prob);
    EXPECT_EQ(*merged_regions_ptr, output_regions);
    // ----- Test 10 -----
}
