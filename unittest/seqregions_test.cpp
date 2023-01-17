#include "gtest/gtest.h"
#include "alignment/seqregions.h"
#include "alignment/seqregions.cpp"

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
void initAlnModelParams(Params& params, Alignment& aln, Model& model, RealNumType *cumulative_rate, std::vector< std::vector<PositionType> >& cumulative_base)
{
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
    initAlnModelParams(params, aln, model, cumulative_rate, cumulative_base);
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
    initAlnModelParams(params, aln, model, cumulative_rate, cumulative_base);
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
    RealNumType *cumulative_rate = nullptr;
    std::vector< std::vector<PositionType> > cumulative_base;
    
    // Init params, aln, and model
    initAlnModelParams(params, aln, model, cumulative_rate, cumulative_base);
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
    RealNumType *cumulative_rate = nullptr;
    std::vector< std::vector<PositionType> > cumulative_base;
    
    // Init params, aln, and model
    initAlnModelParams(params, aln, model, cumulative_rate, cumulative_base);
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
    RealNumType *cumulative_rate = nullptr;
    std::vector< std::vector<PositionType> > cumulative_base;
    
    // Init params, aln, and model
    initAlnModelParams(params, aln, model, cumulative_rate, cumulative_base);
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
    RealNumType *cumulative_rate = nullptr;
    std::vector< std::vector<PositionType> > cumulative_base;
    
    // Init params, aln, and model
    initAlnModelParams(params, aln, model, cumulative_rate, cumulative_base);
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
    RealNumType *cumulative_rate = nullptr;
    std::vector< std::vector<PositionType> > cumulative_base;
    
    // Init params, aln, and model
    initAlnModelParams(params, aln, model, cumulative_rate, cumulative_base);
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
    RealNumType *cumulative_rate = nullptr;
    std::vector< std::vector<PositionType> > cumulative_base;
    
    // Init params, aln, and model
    initAlnModelParams(params, aln, model, cumulative_rate, cumulative_base);
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
    RealNumType *cumulative_rate = nullptr;
    std::vector< std::vector<PositionType> > cumulative_base;
    
    // Init params, aln, and model
    initAlnModelParams(params, aln, model, cumulative_rate, cumulative_base);
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
    RealNumType *cumulative_rate = nullptr;
    std::vector< std::vector<PositionType> > cumulative_base;
    
    // Init params, aln, and model
    initAlnModelParams(params, aln, model, cumulative_rate, cumulative_base);
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
