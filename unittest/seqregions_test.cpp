#include "gtest/gtest.h"
#include "../alignment/seqregions.h"
#include "../model/model_dna.h"
#include "../tree/tree.h"

using namespace cmaple;

/*
 Test constructor
 */
TEST(SeqRegions, constructor)
{
    SeqRegions seqregions1;
    seqregions1.resize(8);
    seqregions1.emplace_back(TYPE_R, 3500, 0, 0.1321);
    
    SeqRegions seqregions2(std::move(seqregions1));
    EXPECT_EQ(seqregions2.size(), 9);
    EXPECT_EQ(seqregions2.back().position, 3500);
    EXPECT_EQ(seqregions2.back().plength_observation2node, 0);
    EXPECT_EQ(seqregions2.back().plength_observation2root, 0.1321);
    
    // Test constructor with null param
    EXPECT_THROW(SeqRegions invalidSeqRegions(NULL), std::invalid_argument);
}

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
    
    SeqRegions::addNonConsecutiveRRegion(seqregions, TYPE_R, -1, -1, 159, threshold_prob);
    EXPECT_EQ(seqregions.size(), 2); // merge two consecotive R
    EXPECT_EQ(seqregions.back().position, 159);
    
    SeqRegions::addNonConsecutiveRRegion(seqregions, 3, -1, 0, 160, threshold_prob);
    EXPECT_EQ(seqregions.size(), 3);
    EXPECT_EQ(seqregions.back().position, 160);
    EXPECT_EQ(seqregions.back().plength_observation2root, 0);
    
    SeqRegions::addNonConsecutiveRRegion(seqregions, TYPE_R, -1, -1, 223, threshold_prob);
    EXPECT_EQ(seqregions.size(), 4);
    EXPECT_EQ(seqregions.back().position, 223);
    
    SeqRegions::addNonConsecutiveRRegion(seqregions, TYPE_R, -1, 1e-10, 240, threshold_prob);
    EXPECT_EQ(seqregions.size(), 5);
    EXPECT_EQ(seqregions.back().position, 240);
    EXPECT_EQ(seqregions.back().plength_observation2root, 1e-10);
    
    SeqRegions::addNonConsecutiveRRegion(seqregions, TYPE_R, -1, 1e-11, 264, threshold_prob);
    EXPECT_EQ(seqregions.size(), 5); // merge two consecotive R
    EXPECT_EQ(seqregions.back().position, 264);
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
    
    /*SeqRegions::addNonConsecutiveRRegion(seqregions, TYPE_R, 1e-9, 0, 324, threshold_prob);
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
    EXPECT_EQ(seqregions.back().plength_observation2root, 1e-301);*/
}

/*
    Load the alignment with 5K seqs
 */
cmaple::Alignment loadAln5K()
{
    // detect the path to the example directory
    std::string example_dir = "../../example/";
    if (!fileExists(example_dir + "example.maple"))
        example_dir = "../example/";
    
    // load aln 5K
    return Alignment(example_dir + "test_5K.maple");
}

/*
    Generate testing data (seqregions1, seqregions2)
 */
void genTestData1(std::unique_ptr<SeqRegions>& seqregions1, std::unique_ptr<SeqRegions>& seqregions2)
{
    Alignment aln = loadAln5K();
    Model model(cmaple::ModelBase::GTR);
    std::unique_ptr<Params> params = ParamsBuilder().build();
    Tree tree(&aln, &model);
    
    std::unique_ptr<SeqRegions> seqregions_1 = aln.data[0]
        .getLowerLhVector(aln.ref_seq.size(), aln.num_states, aln.getSeqType());
    std::unique_ptr<SeqRegions> seqregions_2 = aln.data[10]
        .getLowerLhVector(aln.ref_seq.size(), aln.num_states, aln.getSeqType());
    std::unique_ptr<SeqRegions> seqregions_3 = aln.data[100]
        .getLowerLhVector(aln.ref_seq.size(), aln.num_states, aln.getSeqType());
    
    seqregions_1->mergeTwoLowers<4>(seqregions1, 1e-5, *seqregions_2, 123e-3,
                tree.aln, tree.model, tree.cumulative_rate, params->threshold_prob);
    
    seqregions_2->mergeTwoLowers<4>(seqregions2, 214e-4, *seqregions_3, 13e-8,
                tree.aln, tree.model, tree.cumulative_rate, params->threshold_prob);
}

/*
 Test getNextSharedSegment(PositionType current_pos, const SeqRegions& seq1_region,
 const SeqRegions& seq2_region, size_t& i1, size_t& i2, PositionType &end_pos)
 */
TEST(SeqRegions, getNextSharedSegment)
{
    PositionType current_pos{0};
    PositionType end_pos{0};
    size_t i1{0};
    size_t i2{0};
    std::unique_ptr<SeqRegions> seqregions1 = nullptr;
    std::unique_ptr<SeqRegions> seqregions2 = nullptr;
    genTestData1(seqregions1, seqregions2);
    
    // test getNextSharedSegment()
    SeqRegions::getNextSharedSegment(current_pos, *seqregions1, *seqregions2, i1, i2, end_pos);
    EXPECT_EQ(i1, 0);
    EXPECT_EQ(i2, 0);
    EXPECT_EQ(end_pos, 239);
    current_pos = end_pos + 1;
    
    SeqRegions::getNextSharedSegment(current_pos, *seqregions1, *seqregions2, i1, i2, end_pos);
    EXPECT_EQ(i1, 1);
    EXPECT_EQ(i2, 1);
    EXPECT_EQ(end_pos, 240);
    current_pos = end_pos + 1;
    
    SeqRegions::getNextSharedSegment(current_pos, *seqregions1, *seqregions2, i1, i2, end_pos);
    EXPECT_EQ(i1, 2);
    EXPECT_EQ(i2, 2);
    EXPECT_EQ(end_pos, 377);
    current_pos = end_pos + 1;
    
    /*SeqRegions::getNextSharedSegment(current_pos, *seqregions1, *seqregions2, i1, i2, end_pos);
    EXPECT_EQ(i1, 2);
    EXPECT_EQ(i2, 3);
    EXPECT_EQ(end_pos, 378);
    current_pos = end_pos + 1;
    
    SeqRegions::getNextSharedSegment(current_pos, *seqregions1, *seqregions2, i1, i2, end_pos);
    EXPECT_EQ(i1, 2);
    EXPECT_EQ(i2, 4);
    EXPECT_EQ(end_pos, 1705);
    current_pos = end_pos + 1;
    
    SeqRegions::getNextSharedSegment(current_pos, *seqregions1, *seqregions2, i1, i2, end_pos);
    EXPECT_EQ(i1, 3);
    EXPECT_EQ(i2, 5);
    EXPECT_EQ(end_pos, 1706);
    current_pos = end_pos + 1;
    
    SeqRegions::getNextSharedSegment(current_pos, *seqregions1, *seqregions2, i1, i2, end_pos);
    EXPECT_EQ(i1, 4);
    EXPECT_EQ(i2, 6);
    EXPECT_EQ(end_pos, 3035);
    current_pos = end_pos + 1;
    
    SeqRegions::getNextSharedSegment(current_pos, *seqregions1, *seqregions2, i1, i2, end_pos);
    EXPECT_EQ(i1, 5);
    EXPECT_EQ(i2, 7);
    EXPECT_EQ(end_pos, 3036);
    current_pos = end_pos + 1;
    
    SeqRegions::getNextSharedSegment(current_pos, *seqregions1, *seqregions2, i1, i2, end_pos);
    EXPECT_EQ(i1, 6);
    EXPECT_EQ(i2, 8);
    EXPECT_EQ(end_pos, 14406);
    current_pos = end_pos + 1;
    
    SeqRegions::getNextSharedSegment(current_pos, *seqregions1, *seqregions2, i1, i2, end_pos);
    EXPECT_EQ(i1, 7);
    EXPECT_EQ(i2, 9);
    EXPECT_EQ(end_pos, 14407);
    current_pos = end_pos + 1;
    
    SeqRegions::getNextSharedSegment(current_pos, *seqregions1, *seqregions2, i1, i2, end_pos);
    EXPECT_EQ(i1, 8);
    EXPECT_EQ(i2, 10);
    EXPECT_EQ(end_pos, 23401);
    current_pos = end_pos + 1;
    
    SeqRegions::getNextSharedSegment(current_pos, *seqregions1, *seqregions2, i1, i2, end_pos);
    EXPECT_EQ(i1, 9);
    EXPECT_EQ(i2, 11);
    EXPECT_EQ(end_pos, 23402);
    current_pos = end_pos + 1;
    
    SeqRegions::getNextSharedSegment(current_pos, *seqregions1, *seqregions2, i1, i2, end_pos);
    EXPECT_EQ(i1, 10);
    EXPECT_EQ(i2, 12);
    EXPECT_EQ(end_pos, 26445);
    current_pos = end_pos + 1;
    
    SeqRegions::getNextSharedSegment(current_pos, *seqregions1, *seqregions2, i1, i2, end_pos);
    EXPECT_EQ(i1, 10);
    EXPECT_EQ(i2, 13);
    EXPECT_EQ(end_pos, 26446);
    current_pos = end_pos + 1;
    
    SeqRegions::getNextSharedSegment(current_pos, *seqregions1, *seqregions2, i1, i2, end_pos);
    EXPECT_EQ(i1, 10);
    EXPECT_EQ(i2, 14);
    EXPECT_EQ(end_pos, 28142);
    current_pos = end_pos + 1;*/
}

/*
 Test countSharedSegments(const SeqRegions& seq2_regions, const size_t seq_length) const
 */
TEST(SeqRegions, countSharedSegments)
{
    std::unique_ptr<SeqRegions> seqregions1 = nullptr;
    std::unique_ptr<SeqRegions> seqregions2 = nullptr;
    
    // test empty data
    EXPECT_EQ(seqregions1->countSharedSegments(*seqregions2, 0), 1);
    
    // generate testing data
    genTestData1(seqregions1, seqregions2);
    EXPECT_EQ(seqregions1->countSharedSegments(*seqregions2, 100), 2);
    EXPECT_EQ(seqregions1->countSharedSegments(*seqregions2, 200), 2);
    EXPECT_EQ(seqregions1->countSharedSegments(*seqregions2, 500), 6);
    EXPECT_EQ(seqregions1->countSharedSegments(*seqregions2, 1000), 6);
    EXPECT_EQ(seqregions1->countSharedSegments(*seqregions2, 1500), 6);
    EXPECT_EQ(seqregions1->countSharedSegments(*seqregions2, 2000), 8);
    EXPECT_EQ(seqregions1->countSharedSegments(*seqregions2, 3000), 8);
    EXPECT_EQ(seqregions1->countSharedSegments(*seqregions2, 3500), 10);
}

/*
 Test compareWithSample(const SeqRegions& sequence2, PositionType seq_length, StateType num_states) const
 */
TEST(SeqRegions, compareWithSample)
{
    Alignment aln = loadAln5K();
    std::unique_ptr<SeqRegions> seqregions1 = aln.data[0]
        .getLowerLhVector(aln.ref_seq.size(), aln.num_states, aln.getSeqType());
    std::unique_ptr<SeqRegions> seqregions2 = aln.data[10]
        .getLowerLhVector(aln.ref_seq.size(), aln.num_states, aln.getSeqType());
    std::vector<int> expected_results{1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
    std::vector<int> results(20);
    
    // compare the first 20 sequences with the first one
    for (int i = 0; i < results.size(); ++i)
        results[i] = seqregions1->compareWithSample(*aln.data[i]
                            .getLowerLhVector(aln.ref_seq.size(), aln.num_states,
                            aln.getSeqType()), aln.ref_seq.size(), &aln);
    EXPECT_EQ(expected_results, results);
    
    // ----- Test compareWithSample() invalid sequence length
#ifdef DEBUG
    EXPECT_DEATH(seqregions1->compareWithSample(*seqregions2, 0, &aln), ".*");
#endif
}

/*
 Test areDiffFrom(const SeqRegions& regions2, PositionType seq_length,
 StateType num_states, const Params* params) const
 */
TEST(SeqRegions, areDiffFrom)
{
    Alignment aln = loadAln5K();
    Model model(cmaple::ModelBase::GTR);
    Tree tree(&aln, &model);
    std::unique_ptr<Params>& params = tree.params;
    const PositionType seq_length = aln.ref_seq.size();
    
    std::unique_ptr<SeqRegions> seqregions1 = nullptr;
    std::unique_ptr<SeqRegions> seqregions2 = nullptr;
    
    std::unique_ptr<SeqRegions> seqregions_1 = aln.data[0]
        .getLowerLhVector(seq_length, aln.num_states, aln.getSeqType());
    std::unique_ptr<SeqRegions> seqregions_2 = aln.data[10]
        .getLowerLhVector(seq_length, aln.num_states, aln.getSeqType());
    std::unique_ptr<SeqRegions> seqregions_3 = aln.data[100]
        .getLowerLhVector(seq_length, aln.num_states, aln.getSeqType());
    
    seqregions_1->mergeTwoLowers<4>(seqregions1, 1e-5, *seqregions_2, 123e-3, tree.aln,
                        tree.model, tree.cumulative_rate, params->threshold_prob);
    
    seqregions_2->mergeTwoLowers<4>(seqregions2, 214e-4, *seqregions_3, 13e-8,
                    tree.aln, tree.model, tree.cumulative_rate, params->threshold_prob);
    
    // ---- seqregions4 is empty -----
    std::unique_ptr<SeqRegions> seqregions4 = cmaple::make_unique<SeqRegions>();
    EXPECT_EQ(seqregions1->areDiffFrom(seqregions4, seq_length, aln.num_states, *params), true);
    // ---- seqregions4 is empty -----
    
    // ---- other tests -----
    EXPECT_EQ(seqregions1->areDiffFrom(seqregions2, seq_length, aln.num_states, *params), true);
    
    std::unique_ptr<SeqRegions> seqregions5 = cmaple::make_unique<SeqRegions>(seqregions1);
    EXPECT_EQ(seqregions1->areDiffFrom(seqregions5, seq_length, aln.num_states, *params), false);
    
    seqregions5->data()[3].type = TYPE_R;
    EXPECT_EQ(seqregions1->areDiffFrom(seqregions5, seq_length, aln.num_states, *params), true);
    
    seqregions5->data()[3].type = TYPE_O;
    seqregions5->data()[0].plength_observation2root = 0;
    EXPECT_EQ(seqregions1->areDiffFrom(seqregions5, seq_length, aln.num_states, *params), true);
    
    seqregions1->data()[0].plength_observation2root = 0;
    seqregions5->data()[0].plength_observation2root = 1e-9;
    EXPECT_EQ(seqregions1->areDiffFrom(seqregions5, seq_length, aln.num_states, *params), false);
    
    seqregions1->data()[0].plength_observation2node = 0.0113;
    EXPECT_EQ(seqregions1->areDiffFrom(seqregions5, seq_length, aln.num_states, *params), true);
    
    seqregions5->data()[0].plength_observation2node = 0.0113 + 1e-10;
    EXPECT_EQ(seqregions1->areDiffFrom(seqregions5, seq_length, aln.num_states, *params), false);
    
    // test on type O
    // difference is too small less then thresh_diff_update
    RealNumType thresh_diff_update = params->thresh_diff_update / 2;
    seqregions5->data()[3].likelihood->data()[0] += thresh_diff_update;
    seqregions5->data()[3].likelihood->data()[2] -= thresh_diff_update;
    EXPECT_EQ(seqregions1->areDiffFrom(seqregions5, seq_length, aln.num_states, *params), false);
    
    // reset lh so that all pairs of lh between seqregions1 and seqregions2 equal to each other
    seqregions5->data()[3].likelihood->data()[0] = thresh_diff_update;
    seqregions5->data()[3].likelihood->data()[2] += seqregions5->data()[3].likelihood->data()[0];
    seqregions1->data()[3].likelihood->data()[0] = seqregions5->data()[3].likelihood->data()[0];
    seqregions1->data()[3].likelihood->data()[2] = seqregions5->data()[3].likelihood->data()[2];
    EXPECT_EQ(seqregions1->areDiffFrom(seqregions5, seq_length, aln.num_states, *params), false);
    
    seqregions5->data()[3].likelihood->data()[2] += seqregions5->data()[3].likelihood->data()[0];
    seqregions5->data()[3].likelihood->data()[0] = 0; // one lh = 0 but diff != 0
    EXPECT_EQ(seqregions1->areDiffFrom(seqregions5, seq_length, aln.num_states, *params), true);
    
    seqregions5->data()[3].likelihood->data()[2] -= thresh_diff_update / 10;
    seqregions5->data()[3].likelihood->data()[0] = thresh_diff_update / 10;
    EXPECT_EQ(seqregions1->areDiffFrom(seqregions5, seq_length, aln.num_states, *params), true);
    // ---- other tests -----
}

/*
 Test simplifyO(RealNumType* const partial_lh, StateType ref_state,
 StateType num_states, RealNumType threshold) const
 */
TEST(SeqRegions, simplifyO)
{
    // init testing data
    SeqRegions seqregions1;
    std::unique_ptr<Params> params = ParamsBuilder().build();
    RealNumType threshold_prob = params->threshold_prob;
    
    // tests
    auto new_lh = cmaple::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value{0.1,0.3,0.2,0.4};
    (*new_lh) = new_lh_value;
    EXPECT_EQ(SeqRegions::simplifyO(new_lh->data(), 2, 4, threshold_prob), TYPE_O);
    
    SeqRegion::LHType new_lh_value1{1.0 - 3 * threshold_prob, threshold_prob,
        threshold_prob, threshold_prob};
    (*new_lh) = new_lh_value1;
    EXPECT_EQ(SeqRegions::simplifyO(new_lh->data(), 2, 4, threshold_prob), 0);
    
    new_lh->data()[1] += new_lh->data()[2];
    new_lh->data()[2] = 0;
    EXPECT_EQ(SeqRegions::simplifyO(new_lh->data(), 2, 4, threshold_prob), TYPE_O);
    
    new_lh->data()[2] = new_lh->data()[0];
    new_lh->data()[0] = new_lh->data()[3];
    new_lh->data()[1] = new_lh->data()[3];
    EXPECT_EQ(SeqRegions::simplifyO(new_lh->data(), 2, 4, threshold_prob), TYPE_R);
    
    // Test null lh
#ifdef DEBUG
    EXPECT_DEATH(SeqRegions::simplifyO(NULL, 2, 4, threshold_prob), ".*");
#endif
}

/*
 Test computeAbsoluteLhAtRoot(const Alignment& aln, const ModelBase* model)
 */
TEST(SeqRegions, computeAbsoluteLhAtRoot)
{
    Alignment aln = loadAln5K();
    Model model(cmaple::ModelBase::GTR);
    std::unique_ptr<Params> params = ParamsBuilder().build();
    Tree tree(&aln, &model);
    
    std::unique_ptr<SeqRegions> seqregions1 = nullptr;
    std::unique_ptr<SeqRegions> seqregions2 = nullptr;
    std::unique_ptr<SeqRegions> seqregions3 = nullptr;
    
    // dummy variables
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    std::unique_ptr<SeqRegions> seqregions_1 = aln.data[0]
        .getLowerLhVector(aln.ref_seq.size(), aln.num_states, aln.getSeqType());
    std::unique_ptr<SeqRegions> seqregions_2 = aln.data[10]
        .getLowerLhVector(aln.ref_seq.size(), aln.num_states, aln.getSeqType());
    std::unique_ptr<SeqRegions> seqregions_3 = aln.data[100]
        .getLowerLhVector(aln.ref_seq.size(), aln.num_states, aln.getSeqType());
    
    seqregions_1->mergeTwoLowers<4>(seqregions1, 1e-5, *seqregions_2, 123e-3, \
            tree.aln, tree.model, tree.cumulative_rate, params->threshold_prob);
    
    seqregions_2->mergeTwoLowers<4>(seqregions2, 214e-4, *seqregions_3, 13e-8,
            tree.aln, tree.model, tree.cumulative_rate, params->threshold_prob);
    
    seqregions_3->mergeTwoLowers<4>(seqregions3, 5421e-8, *seqregions_1, 1073e-6,
            tree.aln, tree.model, tree.cumulative_rate, params->threshold_prob);
    
    // ----- Test 1 on a more complex seqregions -----
    EXPECT_EQ(seqregions_1->computeAbsoluteLhAtRoot<4>(tree.model,
                tree.cumulative_base), -40549.6785849070511176250874996185302734375);
    // ----- Test 1 on a more complex seqregions -----
    
    // ----- Test 2 on a more complex seqregions -----
    EXPECT_EQ(seqregions_2->computeAbsoluteLhAtRoot<4>(tree.model,
                tree.cumulative_base), -40549.19068355686613358557224273681640625);
    // ----- Test 2 on a more complex seqregions -----
    
    // ----- Test 3 on a more complex seqregions -----
    EXPECT_EQ(seqregions_3->computeAbsoluteLhAtRoot<4>(tree.model,
                tree.cumulative_base), -40548.41030627492000348865985870361328125);
    // ----- Test 3 on a more complex seqregions -----
}

/*
 Test computeTotalLhAtRoot(StateType num_states, const ModelBase* model, RealNumType blength = -1)
 */
TEST(SeqRegions, computeTotalLhAtRoot)
{
    Alignment aln = loadAln5K();
    Model model(cmaple::ModelBase::GTR);
    Tree tree(&aln, &model);
    std::unique_ptr<Params>& params = tree.params;
    
    std::unique_ptr<SeqRegions> seqregions1 = nullptr;
    std::unique_ptr<SeqRegions> seqregions2 = nullptr;
    std::unique_ptr<SeqRegions> seqregions3 = nullptr;
    std::unique_ptr<SeqRegions> seqregions_total_lh = nullptr;
    
    // Generate complex seqregions
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    std::unique_ptr<SeqRegions> seqregions_1 = aln.data[20]
        .getLowerLhVector(aln.ref_seq.size(), aln.num_states, aln.getSeqType());
    std::unique_ptr<SeqRegions> seqregions_2 = aln.data[200]
        .getLowerLhVector(aln.ref_seq.size(), aln.num_states, aln.getSeqType());
    std::unique_ptr<SeqRegions> seqregions_3 = aln.data[2000]
        .getLowerLhVector(aln.ref_seq.size(), aln.num_states, aln.getSeqType());
    
    seqregions_1->mergeTwoLowers<4>(seqregions1, 1073e-6, *seqregions_2, 13e-8,
                    tree.aln, tree.model, tree.cumulative_rate, params->threshold_prob);
    
    seqregions_2->mergeTwoLowers<4>(seqregions2, 123e-3, *seqregions_3, 5421e-8,
                    tree.aln, tree.model, tree.cumulative_rate, params->threshold_prob);
    
    seqregions_3->mergeTwoLowers<4>(seqregions3, 214e-4, *seqregions_1, 1e-5,
                tree.aln, tree.model, tree.cumulative_rate, params->threshold_prob);
    
    // ----- Test 1 on a more complex seqregions -----
    seqregions1->computeTotalLhAtRoot<4>(seqregions_total_lh, tree.model);
    EXPECT_EQ(seqregions_total_lh->size(), 22);
    SeqRegion::LHType lh_value1_10{7.2185298165439574223086912192976286051226963991212e-10,
        0.00012093455079658274576269449962495627914904616773129,
        5.7783014021321174033127313583984435707563420692168e-09,
        0.99987905894904904879894047553534619510173797607422};
    EXPECT_EQ(*seqregions_total_lh->at(10).likelihood, lh_value1_10);
    SeqRegion::LHType lh_value1_12{0.00012109700727777851427414274043670161518093664199114,
        3.7920254781824222995438774532709486075887639344728e-09,
        0.99987887895812166405562493309844285249710083007812,
        2.0242575103081267528512744062821337998059334495338e-08};
    EXPECT_EQ(*seqregions_total_lh->at(12).likelihood, lh_value1_12);
    SeqRegion::LHType lh_value1_14{7.2185298165439574223086912192976286051226963991212e-10,
        0.00012093455079658274576269449962495627914904616773129,
        5.7783014021321174033127313583984435707563420692168e-09,
        0.99987905894904904879894047553534619510173797607422};
    EXPECT_EQ(*seqregions_total_lh->at(14).likelihood, lh_value1_14);
    SeqRegion::LHType lh_value1_20{0.999878969078505264178602374158799648284912109375,
        3.79202547818242147236326490024327373618007186451e-09,
        0.00012100688689399335864343987267943703045602887868881,
        2.0242575103081260911067843638599939026789797935635e-08};
    EXPECT_EQ(*seqregions_total_lh->at(20).likelihood, lh_value1_20);
    // ----- Test 1 on a more complex seqregions -----
    
    // ----- Test 2 on a more complex seqregions -----
    seqregions2->computeTotalLhAtRoot<4>(seqregions_total_lh, tree.model, 1e-5);
    EXPECT_EQ(seqregions_total_lh->size(), 26);
    SeqRegion::LHType lh_value2_17{1.0426241904738277496250383261089389463904808508232e-06,
        0.00036249067899242780116039752691392550332238897681236,
        6.3013501138813336668181852573411561024840921163559e-06,
        0.99963016534670312562838034864398650825023651123047};
    EXPECT_EQ(*seqregions_total_lh->at(17).likelihood, lh_value2_17);
    SeqRegion::LHType lh_value2_21{0.00042526525345440798929128045635650323674781247973442,
        2.4910168754262691925165963680033343052855343557894e-06,
        0.99955743551688391868026428710436448454856872558594,
        1.4808212786277890710097057680449950112233636900783e-05};
    EXPECT_EQ(*seqregions_total_lh->at(21).likelihood, lh_value2_21);
    // ----- Test 2 on a more complex seqregions -----
    
    // ----- Test 3 on a more complex seqregions -----
    seqregions3->computeTotalLhAtRoot<4>(seqregions_total_lh, tree.model, 0.03762);
    EXPECT_EQ(seqregions_total_lh->size(), 27);
    SeqRegion::LHType lh_value3_10{0.0036579520570985094192473230378936932538636028766632,
        0.93983378575187570547200266446452587842941284179688,
        0.0036637083625175163349718676641941783600486814975739,
        0.052844553828508410153741436943164444528520107269287};
    EXPECT_EQ(*seqregions_total_lh->at(10).likelihood, lh_value3_10);
    // ----- Test 3 on a more complex seqregions -----
}

/*
 Test merge_N_O(const RealNumType lower_plength, const SeqRegion& reg_o,
 const ModelBase* model,
 const PositionType end_pos, const StateType num_states,
 SeqRegions& merged_target)
 */
TEST(SeqRegions, merge_N_O)
{
    std::unique_ptr<Params> params = ParamsBuilder().build();
    Alignment aln = loadAln5K();
    Model model(cmaple::ModelBase::GTR);
    Tree tree(&aln, &model);
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    
    // ----- Test 1 -----
    SeqRegions merged_regions;
    const PositionType end_pos = 8623;
    RealNumType lower_plength = -1;
    auto new_lh = cmaple::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value{0.27595781923293211113090706021466758102178573608398,
        8.3817209383128356075479820086471249851456377655268e-06,
        0.72401575181023847260775028189527802169322967529297,
        1.8047235891267821277644811672757896303664892911911e-05};
    (*new_lh) = new_lh_value;
    SeqRegion seqregion1(TYPE_O, 243, -1, -1, std::move(new_lh));
    merge_N_O<4>(lower_plength, seqregion1, tree.model,
              end_pos, merged_regions);
    EXPECT_EQ(merged_regions.size(), 1);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge{0.3675361745020691,0.0000068532680661245309,
        0.63243117236202639,0.000025799867838435163};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    seqregion1.plength_observation2root = 0.311; // plength_observation2root
    // does NOT affect the likelihood of new merged region
    merge_N_O<4>(lower_plength, seqregion1, tree.model,
              end_pos, merged_regions);
    EXPECT_EQ(merged_regions.size(), 2);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    seqregion1.plength_observation2node = 0; // plength_observation2node = 0
    // does NOT affect the likelihood of new merged region
    merge_N_O<4>(lower_plength, seqregion1, tree.model,
              end_pos, merged_regions);
    EXPECT_EQ(merged_regions.size(), 3);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge);
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    seqregion1.plength_observation2node = 0.1252; // plength_observation2node > 0
    // affects the likelihood of new merged region
    merge_N_O<4>(lower_plength, seqregion1, tree.model,
              end_pos, merged_regions);
    EXPECT_EQ(merged_regions.size(), 4);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge1{0.37599425541307963,0.0099624992407331657,
        0.55990605200203047,0.054137193344156891};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge1);
    // ----- Test 4 -----
    
    // ----- Test 5 -----
    lower_plength = 0; // lower_plength = 0; plength_observation2node > 0
    merge_N_O<4>(lower_plength, seqregion1, tree.model,
              end_pos, merged_regions);
    EXPECT_EQ(merged_regions.size(), 5);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge1);
    // ----- Test 5 -----
    
    /*// ----- Test 6 -----
    lower_plength = 0; // lower_plength = 0; plength_observation2node = 0
    seqregion1.plength_observation2node = 0;
    merge_N_O<4>(lower_plength, seqregion1, tree.model,
              end_pos, merged_regions);
    EXPECT_EQ(merged_regions.size(), 6);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge);
    // ----- Test 6 -----
    
    // ----- Test 7 -----
    lower_plength = 0; // lower_plength = 0; plength_observation2node = -1
    seqregion1.plength_observation2node = -1;
    merge_N_O<4>(lower_plength, seqregion1, tree.model,
              end_pos, merged_regions);
    EXPECT_EQ(merged_regions.size(), 7);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge);
    // ----- Test 7 -----
    
    // ----- Test 8 -----
    lower_plength = 0.02331; // lower_plength > 0; plength_observation2node = -1
    seqregion1.plength_observation2node = -1;
    merge_N_O<4>(lower_plength, seqregion1, tree.model,
              end_pos, merged_regions);
    EXPECT_EQ(merged_regions.size(), 8);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge2{0.36911091784101202,0.0018604164279931914,
     0.61892829252503356,0.010100373205961303};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge2);
    // ----- Test 8 -----
    
    // ----- Test 9 -----
    lower_plength = 0.02331; // lower_plength > 0; plength_observation2node = 0
    seqregion1.plength_observation2node = 0;
    merge_N_O<4>(lower_plength, seqregion1, tree.model,
              end_pos, merged_regions);
    EXPECT_EQ(merged_regions.size(), 9);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge2);
    // ----- Test 9 -----
    
    // ----- Test 10 -----
    lower_plength = 0.02331; // lower_plength > 0; plength_observation2node > 0
    seqregion1.plength_observation2node = 0.1e-6;
    merge_N_O<4>(lower_plength, seqregion1, tree.model,
              end_pos, merged_regions);
    EXPECT_EQ(merged_regions.size(), 10);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge3{0.36911092459666772,0.0018604243797870988,
     0.61892823459762114,0.010100416425924146};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge3);
    // ----- Test 10 -----*/
}

/*
 Test merge_N_RACGT(const SeqRegion& reg_racgt, const RealNumType lower_plength, const PositionType end_pos,
 const RealNumType threshold_prob, SeqRegions& merged_regions)
 */
TEST(SeqRegions, merge_N_RACGT)
{
    Alignment aln = loadAln5K();
    Model model(cmaple::ModelBase::GTR);
    Tree(&aln, &model);
    
    // dummy variables
    std::unique_ptr<Params> params = ParamsBuilder().build();
    
    // ----- Test 1 -----
    SeqRegions merged_regions;
    const PositionType end_pos = 6543;
    RealNumType lower_plength = -1;
    SeqRegion seqregion1(TYPE_R, 3442, -1, -1);
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 2);
    /*EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);*/
    EXPECT_TRUE(merged_regions.back().type == seqregion1.type);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    lower_plength = -1;
    seqregion1.plength_observation2node = -1;
    seqregion1.plength_observation2root = 0.3231;
    seqregion1.type = TYPE_R;
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
    EXPECT_TRUE(merged_regions.size() == 3);
    // EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_FALSE(merged_regions.back().plength_observation2node != -1);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    lower_plength = -1;
    seqregion1.plength_observation2node = 0;
    seqregion1.plength_observation2root = -1;
    seqregion1.type = 3;
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 10);
    EXPECT_LT(merged_regions.back().plength_observation2root, 0);
    EXPECT_LT(merged_regions.back().plength_observation2node, 0);
    EXPECT_TRUE(merged_regions.back().type == seqregion1.type);
    // ----- Test 10 -----
    
    /*// ----- Test 11 -----
    lower_plength = 0;
    seqregion1.plength_observation2node = -1;
    seqregion1.plength_observation2root = 0;
    seqregion1.type = TYPE_R;
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_N_RACGT(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 27);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0.1001211);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 27 -----*/
}

/*
 Test merge_O_N(const SeqRegion& reg_o, const RealNumType upper_plength,
 const PositionType end_pos, const ModelBase* model, const StateType num_states,
 SeqRegions& merged_regions)
 */
TEST(SeqRegions, merge_O_N)
{
    std::unique_ptr<Params> params = ParamsBuilder().build();
    Alignment aln = loadAln5K();
    Model model(cmaple::ModelBase::GTR);
    Tree tree(&aln, &model);
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    
    // ----- Test 1 -----
    SeqRegions merged_regions;
    const PositionType end_pos = 29180;
    RealNumType upper_plength = -1;
    auto new_lh = cmaple::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value{7.4346402191731947664294494204639818235591519623995e-06,
        0.66666418845326025355291221785591915249824523925781,
        7.4346402191731947664294494204639818235591519623995e-06,
        0.33332094226630137878686355179524980485439300537109};
    (*new_lh) = new_lh_value;
    SeqRegion seqregion1(TYPE_O, 243, -1, -1, std::move(new_lh));
    merge_O_N<4>(seqregion1, upper_plength, end_pos, tree.model, merged_regions);
    EXPECT_EQ(merged_regions.size(), 1);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge{0.0000074346402191731948,0.66666418845326025,
        0.0000074346402191731948,0.33332094226630138};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    seqregion1.plength_observation2root = 0.311; // plength_observation2root
    // does NOT affect the likelihood of new merged region
    merge_O_N<4>(seqregion1, upper_plength, end_pos, tree.model, merged_regions);
    EXPECT_EQ(merged_regions.size(), 2);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    seqregion1.plength_observation2node = 0; // plength_observation2node = 0
    // does NOT affect the likelihood of new merged region
    merge_O_N<4>(seqregion1, upper_plength, end_pos, tree.model, merged_regions);
    EXPECT_EQ(merged_regions.size(), 3);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge);
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    seqregion1.plength_observation2node = 0.1252; // plength_observation2node > 0
    // affects the likelihood of new merged region
    merge_O_N<4>(seqregion1, upper_plength, end_pos, tree.model, merged_regions);
    EXPECT_EQ(merged_regions.size(), 4);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge1{0.011218665480774669,0.5673625983091064,
        0.024370520108437949,0.39704821610168095};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge1);
    // ----- Test 4 -----
    
    // ----- Test 5 -----
    upper_plength = 0; // lower_plength = 0; plength_observation2node > 0
    merge_O_N<4>(seqregion1, upper_plength, end_pos, tree.model, merged_regions);
    EXPECT_EQ(merged_regions.size(), 5);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge1);
    // ----- Test 5 -----
    
    /*// ----- Test 6 -----
    upper_plength = 0; // lower_plength = 0; plength_observation2node = 0
    seqregion1.plength_observation2node = 0;
    merge_O_N<4>(seqregion1, upper_plength, end_pos, tree.model, merged_regions);
    EXPECT_EQ(merged_regions.size(), 6);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge);
    // ----- Test 6 -----
    
    // ----- Test 7 -----
    upper_plength = 0; // lower_plength = 0; plength_observation2node = -1
    seqregion1.plength_observation2node = -1;
    merge_O_N<4>(seqregion1, upper_plength, end_pos, tree.model, merged_regions);
    EXPECT_EQ(merged_regions.size(), 7);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge);
    // ----- Test 7 -----
    
    // ----- Test 8 -----
    upper_plength = 0.02331; // lower_plength > 0; plength_observation2node = -1
    seqregion1.plength_observation2node = -1;
    merge_O_N<4>(seqregion1, upper_plength, end_pos, tree.model, merged_regions);
    EXPECT_EQ(merged_regions.size(), 8);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge2{0.0020947652384088584,0.64817600901028716,
     0.00454340526533243,0.34518582048597146};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge2);
    // ----- Test 8 -----
    
    // ----- Test 9 -----
    upper_plength = 0.02331; // lower_plength > 0; plength_observation2node = 0
    seqregion1.plength_observation2node = 0;
    merge_O_N<4>(seqregion1, upper_plength, end_pos, tree.model, merged_regions);
    EXPECT_EQ(merged_regions.size(), 9);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge2);
    // ----- Test 9 -----
    
    // ----- Test 10 -----
    upper_plength = 0.02331; // lower_plength > 0; plength_observation2node > 0
    seqregion1.plength_observation2node = 0.1e-6;
    merge_O_N<4>(seqregion1, upper_plength, end_pos, tree.model, merged_regions);
    EXPECT_EQ(merged_regions.size(), 10);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge3{0.0020947741930660798,0.6481759296959182,
     0.0045434247246658715,0.34518587138634999};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge3);
    // ----- Test 10 -----*/
}

/*
 Test merge_RACGT_N(const SeqRegion& reg_n, const RealNumType upper_plength,
 const PositionType end_pos,
 const RealNumType threshold_prob, SeqRegions& merged_regions)
 */
TEST(SeqRegions, merge_RACGT_N)
{
    std::unique_ptr<Params> params = ParamsBuilder().build();
    Alignment aln = loadAln5K();
    Model model(cmaple::ModelBase::GTR);
    Tree tree(&aln, &model);
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    
    // ----- Test 1 -----
    SeqRegions merged_regions;
    const PositionType end_pos = 6543;
    RealNumType lower_plength = -1;
    SeqRegion seqregion1(TYPE_R, 3442, -1, -1);
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
    EXPECT_GE(merged_regions.size(), 1);
    EXPECT_LE(merged_regions.back().plength_observation2root, -1);
    EXPECT_LE(merged_regions.back().plength_observation2node, -1);
    EXPECT_FALSE(merged_regions.back().type != seqregion1.type);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    lower_plength = -1;
    seqregion1.plength_observation2node = -1;
    seqregion1.plength_observation2root = 0;
    seqregion1.type = 1; // C
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
    EXPECT_TRUE(merged_regions.size() == 10);
    /* EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);*/
    EXPECT_TRUE(merged_regions.back().type == seqregion1.type);
    // ----- Test 10 -----
    
    /*// ----- Test 11 -----
    lower_plength = 0;
    seqregion1.plength_observation2node = -1;
    seqregion1.plength_observation2root = 0;
    seqregion1.type = TYPE_R;
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
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
    merge_RACGT_N(seqregion1, lower_plength, end_pos, params->threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 27);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0.1001211);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0.1);
    EXPECT_EQ(merged_regions.back().type, seqregion1.type);
    // ----- Test 27 -----*/
}

/*
 Test merge_Zero_Distance(const SeqRegion& seq1_region, const SeqRegion& seq2_region,
 const RealNumType total_blength_1, const RealNumType total_blength_2,
 const PositionType end_pos, const RealNumType threshold_prob, const StateType num_states,
 SeqRegions* &merged_regions)
 */
TEST(SeqRegions, merge_Zero_Distance)
{
    std::unique_ptr<Params> params = ParamsBuilder().build();
    Alignment aln = loadAln5K();
    Model model(cmaple::ModelBase::GTR);
    Tree tree(&aln, &model);
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    
    // ----- Test 1 -----
    std::unique_ptr<SeqRegions> merged_regions_ptr = cmaple::make_unique<SeqRegions>();
    const RealNumType threshold_prob = params->threshold_prob;
    const PositionType end_pos = 6543;
    RealNumType total_blength_1 = -1;
    RealNumType total_blength_2 = -1;
    SeqRegion seqregion1(TYPE_R, 3442, -1, -1);
    SeqRegion seqregion2(0, 234, -1, 0);
    EXPECT_TRUE(merge_Zero_Distance(seqregion1, seqregion2, total_blength_1, total_blength_2,
                    end_pos, threshold_prob, num_states, merged_regions_ptr));
    EXPECT_EQ(merged_regions_ptr, nullptr);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    merged_regions_ptr = cmaple::make_unique<SeqRegions>();
    total_blength_1 = -1;
    total_blength_2 = 0;
    seqregion1.type = TYPE_R;
    seqregion2.type = 0;
    EXPECT_TRUE(merge_Zero_Distance(seqregion1, seqregion2, total_blength_1,
                total_blength_2, end_pos, threshold_prob, num_states, merged_regions_ptr));
    EXPECT_EQ(merged_regions_ptr, nullptr);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    merged_regions_ptr = cmaple::make_unique<SeqRegions>();
    total_blength_1 = -1;
    total_blength_2 = 1e-3;
    seqregion1.type = 3;
    seqregion2.type = 0;
    EXPECT_TRUE(merge_Zero_Distance(seqregion1, seqregion2, total_blength_1, total_blength_2,
                    end_pos, threshold_prob, num_states, merged_regions_ptr));
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
    EXPECT_TRUE(merge_Zero_Distance(seqregion1, seqregion2, total_blength_1, total_blength_2,
                                    end_pos, threshold_prob, num_states, merged_regions_ptr));
    EXPECT_EQ(merged_regions_ptr, nullptr);
    // ----- Test 4 -----
    
    /*// ----- Test 5 -----
    merged_regions_ptr = cmaple::make_unique<SeqRegions>();
    total_blength_1 = 0;
    total_blength_2 = 0;
    seqregion1.type = 0;
    seqregion2.type = 0;
    EXPECT_TRUE(merge_Zero_Distance(seqregion1, seqregion2, total_blength_1, total_blength_2,
     end_pos, threshold_prob, num_states, merged_regions_ptr));
    EXPECT_EQ(merged_regions_ptr, nullptr);
    // ----- Test 5 -----
    
    // ----- Test 6 -----
    merged_regions_ptr = cmaple::make_unique<SeqRegions>();
    total_blength_1 = 0;
    total_blength_2 = 1e-10;
    seqregion1.type = TYPE_R;
    seqregion2.type = TYPE_R;
    EXPECT_TRUE(merge_Zero_Distance(seqregion1, seqregion2, total_blength_1, total_blength_2,
     end_pos, threshold_prob, num_states, merged_regions_ptr));
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
    EXPECT_TRUE(merge_Zero_Distance(seqregion1, seqregion2, total_blength_1, total_blength_2,
     end_pos, threshold_prob, num_states, merged_regions_ptr));
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
    EXPECT_FALSE(merge_Zero_Distance(seqregion1, seqregion2, total_blength_1, total_blength_2,
     end_pos, threshold_prob, num_states, merged_regions_ptr));
    EXPECT_TRUE(merged_regions_ptr != nullptr);
    EXPECT_EQ(merged_regions_ptr->size(), 2);
    // ----- Test 8 -----
    
    // ----- Test 9 -----
    total_blength_1 = 0.1321;
    total_blength_2 = 1e-5;
    seqregion1.type = TYPE_O;
    seqregion2.type = TYPE_R;
    EXPECT_FALSE(merge_Zero_Distance(seqregion1, seqregion2, total_blength_1,
     total_blength_2, end_pos, threshold_prob, num_states, merged_regions_ptr));
    EXPECT_TRUE(merged_regions_ptr != nullptr);
    EXPECT_EQ(merged_regions_ptr->size(), 2);
    // ----- Test 9 -----*/
}

/*
 Test merge_O_ORACGT<4>(const SeqRegion& seq1_region, const SeqRegion& seq2_region,
 const RealNumType total_blength_1, const RealNumType total_blength_2,
 const PositionType end_pos, const RealNumType threshold_prob, const ModelBase* model,
 const Alignment& aln, SeqRegions& merged_regions)
 */
TEST(SeqRegions, merge_O_ORACGT)
{
    Alignment aln = loadAln5K();
    Model model(cmaple::ModelBase::GTR);
    Tree tree(&aln, &model);
    std::unique_ptr<Params>& params = tree.params;
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    
    // ----- Test 1 -----
    SeqRegions merged_regions;
    auto new_lh = cmaple::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value{0.27595781923293211113090706021466758102178573608398,
        8.3817209383128356075479820086471249851456377655268e-06,
        0.72401575181023847260775028189527802169322967529297,
        1.8047235891267821277644811672757896303664892911911e-05};
    (*new_lh) = new_lh_value;
    SeqRegion seqregion1(TYPE_O, 243, -1, -1, std::move(new_lh));
    new_lh = cmaple::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value1{0.27595781923293211113090706021466758102178573608398,
        8.3817209383128356075479820086471249851456377655268e-06,
        0.72401575181023847260775028189527802169322967529297,
        1.8047235891267821277644811672757896303664892911911e-05};
    (*new_lh) = new_lh_value1;
    SeqRegion seqregion2(TYPE_O, 243, -1, 1e-3, std::move(new_lh));
    RealNumType total_blength_1 = -1;
    RealNumType total_blength_2 = -1;
    const RealNumType threshold_prob = params->threshold_prob;
    const PositionType end_pos = 3213;
    merge_O_ORACGT<4>(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos,
                      threshold_prob, tree.model, tree.aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 1);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge{0.12684687976595482,1.1702018350525203E-10,
        0.87315311957450514,5.425200212298543E-10};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    total_blength_1 = -1;
    total_blength_2 = 1e-5;
    merge_O_ORACGT<4>(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos,
                      threshold_prob, tree.model, tree.aln, merged_regions);
    EXPECT_TRUE(merged_regions.size() == 2);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    /*EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);*/
    SeqRegion::LHType new_lh_value_merge1{0.12684809781766063,1.3059895886789165E-10,
        0.87315190141833243,6.3340794416577515E-10};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge1);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    total_blength_1 = 0.01;
    total_blength_2 = -1;
    merge_O_ORACGT<4>(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos,
                      threshold_prob, tree.model, tree.aln, merged_regions);
    EXPECT_FALSE(merged_regions.size() > 3);
    EXPECT_FALSE(merged_regions.back().type != TYPE_O);
    // EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_FALSE(merged_regions.back().plength_observation2node != 0);
    SeqRegion::LHType new_lh_value_merge2{0.12842468075362262,1.1709109654051516E-8,
        0.87157516057518813,1.4696207956334386E-7};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge2);
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    total_blength_1 = 0.011321;
    total_blength_2 = 1e-3;
    merge_O_ORACGT<4>(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos,
                      threshold_prob, tree.model, tree.aln, merged_regions);
    // EXPECT_EQ(merged_regions.size(), 4);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    // EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_TRUE(merged_regions.back().plength_observation2node == 0);
    SeqRegion::LHType new_lh_value_merge3{0.12875794255340739,1.6716816801441426E-7,
        0.87123893271963804,0.0000029575587865296985};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge3);
    // ----- Test 4 -----
    
    /*// ----- Test 5 -----
    total_blength_1 = -1;
    total_blength_2 = -1;
    SeqRegion seqregion3(TYPE_R, 4324, 0.121, 0);
    merge_O_ORACGT<4>(seqregion1, seqregion3, total_blength_1, total_blength_2, end_pos,
     threshold_prob, tree.model, tree.aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 5);
    EXPECT_EQ(merged_regions.back().type, seqregion3.type);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    // ----- Test 5 -----
    
    // ----- Test 6 -----
    total_blength_1 = -1;
    total_blength_2 = 0; // either total_blength_2 = -1 or = 0 leads to the same result
    seqregion3.type = 2;
    merge_O_ORACGT<4>(seqregion1, seqregion3, total_blength_1, total_blength_2, end_pos,
     threshold_prob, tree.model, tree.aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 6);
    EXPECT_EQ(merged_regions.back().type, seqregion3.type);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    // ----- Test 6 -----
    
    // ----- Test 7 -----
    total_blength_1 = -1;
    total_blength_2 = 0;
    seqregion3.type = 0;
    merge_O_ORACGT<4>(seqregion1, seqregion3, total_blength_1, total_blength_2,
     end_pos, threshold_prob, tree.model, tree.aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 7);
    EXPECT_EQ(merged_regions.back().type, seqregion3.type);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    // ----- Test 7 -----
    
    // ----- Test 8 -----
    total_blength_1 = 0.231;
    total_blength_2 = 0;
    seqregion3.type = 0;
    merge_O_ORACGT<4>(seqregion1, seqregion3, total_blength_1, total_blength_2, end_pos,
     threshold_prob, tree.model, tree.aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 8);
    EXPECT_EQ(merged_regions.back().type, seqregion3.type);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    // ----- Test 8 -----
    
    // ----- Test 9 -----
    total_blength_1 = 0.231;
    total_blength_2 = 0.3121;
    seqregion3.type = 2;
    merge_O_ORACGT<4>(seqregion1, seqregion3, total_blength_1, total_blength_2, end_pos,
     threshold_prob, tree.model, tree.aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 9);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge4{0.045597223147683601,0.0013901861006516344,
     0.92036965105659874,0.032642939695066098};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge4);
    // ----- Test 9 -----*/
}

/*
 Test merge_RACGT_O<4>(const SeqRegion& seq2_region, const RealNumType total_blength_2,
 const PositionType end_pos, SeqRegion::LHType& new_lh, const RealNumType threshold_prob,
 const ModelBase* model, const Alignment& aln, SeqRegions& merged_regions)
 */
TEST(SeqRegions, merge_RACGT_O)
{
    Alignment aln = loadAln5K();
    Model model(cmaple::ModelBase::GTR);
    Tree tree(&aln, &model);
    std::unique_ptr<Params>& params = tree.params;
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    
    // ----- Test 1 -----
    SeqRegions merged_regions;
    auto new_lh = cmaple::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value{0.49999442406129068761089229155913926661014556884766,
        5.5759387092817078802929955938516570768115343526006e-06,
        0.49999442406129068761089229155913926661014556884766,
        5.5759387092817078802929955938516570768115343526006e-06};
    (*new_lh) = new_lh_value;
    SeqRegion seqregion1(TYPE_O, 243, -1, -1, std::move(new_lh));
    RealNumType total_blength_1 = -1;
    const PositionType end_pos = 3213;
    const RealNumType threshold_prob = params->threshold_prob;
    new_lh = cmaple::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value1{0.99997490659690790870683940738672390580177307128906,
        2.5092097243268992381179730011275808010395849123597e-05,
        6.5292438309787796439726148030194621818544931102224e-10,
        6.5292438309787796439726148030194621818544931102224e-10};
    (*new_lh) = new_lh_value1;
    merge_RACGT_O<4>(seqregion1, total_blength_1, end_pos, *new_lh,
                     threshold_prob, tree.model, tree.aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 1);
    EXPECT_EQ(merged_regions.back().type, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    SeqRegion::LHType new_lh_value2{0.24998327274350143345493791002809302881360054016113,
        0.74997769674260927885711680573876947164535522460938,
        1.9515256944661845116837164959555650511902058497071e-05,
        1.9515256944661845116837164959555650511902058497071e-05};
    (*new_lh) = new_lh_value2;
    merge_RACGT_O<4>(seqregion1, total_blength_1, end_pos, *new_lh,
                     threshold_prob, tree.model, tree.aln, merged_regions);
    EXPECT_LT(merged_regions.size(), 3);
    // EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_TRUE(merged_regions.back().plength_observation2root == 0);
    // EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge{0.99988848806554697,
        0.000033453518158479732,0.000078057545796795158,8.704978900055221E-10};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    total_blength_1 = 0; // either total_blength_1 = 0 or = -1 does not affect the merged region
    (*new_lh) = new_lh_value2;
    merge_RACGT_O<4>(seqregion1, total_blength_1, end_pos, *new_lh,
                     threshold_prob, tree.model, tree.aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 3);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge);
    // ----- Test 3 -----
    
    /*// ----- Test 4 -----
    total_blength_1 = 1e-4;
    (*new_lh) = new_lh_value2;
    merge_RACGT_O<4>(seqregion1, total_blength_1, end_pos, *new_lh,
     threshold_prob, tree.model, tree.aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 4);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge1{0.99983015675026043,
     0.000091790366085780112,0.000078048395545374228,4.4881083285432782E-9};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge1);
    // ----- Test 4 -----
    
    // ----- Test 5 -----
    seqregion1.plength_observation2root = 0.231; // neither plength_observation2root
     // nor plength_observation2node affect the merged region
    total_blength_1 = 1e-4;
    (*new_lh) = new_lh_value2;
    merge_RACGT_O<4>(seqregion1, total_blength_1, end_pos, *new_lh,
     threshold_prob, tree.model, tree.aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 5);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge1);
    // ----- Test 5 -----
    
    // ----- Test 6 -----
    seqregion1.plength_observation2node = 1e-7; // neither plength_observation2root
     // nor plength_observation2node affect the merged region
    total_blength_1 = 1e-4;
    (*new_lh) = new_lh_value2;
    merge_RACGT_O<4>(seqregion1, total_blength_1, end_pos, *new_lh, threshold_prob,
     tree.model, tree.aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 6);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge1);
    // ----- Test 6 -----*/
}

/*
 Test merge_RACGT_RACGT<4>(const SeqRegion& seq2_region, const RealNumType total_blength_2,
 const PositionType end_pos, SeqRegion::LHType& new_lh, const ModelBase* model,
 const Alignment& aln, SeqRegions& merged_regions)
 */
TEST(SeqRegions, merge_RACGT_RACGT)
{
    Alignment aln = loadAln5K();
    Model model(cmaple::ModelBase::GTR);
    Tree tree(&aln, &model);
    std::unique_ptr<Params>& params = tree.params;
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    
    // ----- Test 1 -----
    SeqRegions merged_regions;
    SeqRegion seqregion1(TYPE_R, 4231);
    RealNumType total_blength_1 = -1;
    const PositionType end_pos = 3213;
    const RealNumType threshold_prob = params->threshold_prob;
    auto new_lh = cmaple::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value1{5.5759387092817078802929955938516570768115343526006e-06,
        0.49999442406129068761089229155913926661014556884766,
        5.5759387092817078802929955938516570768115343526006e-06,
        0.49999442406129068761089229155913926661014556884766};
    (*new_lh) = new_lh_value1;
    merge_RACGT_RACGT<4>(seqregion1, total_blength_1, end_pos, *new_lh,
                         tree.model, tree.aln, merged_regions);
    EXPECT_GT(merged_regions.size(), 0);
    EXPECT_TRUE(merged_regions.back().type == TYPE_O);
    EXPECT_FALSE(merged_regions.back().plength_observation2root != 0);
    // EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge{0,0,0,1};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    SeqRegion::LHType new_lh_value2{8.3640013382402129011325767060647251582850003615022e-06,
        8.3640013382402129011325767060647251582850003615022e-06,
        0.24998327199732350845096107150311581790447235107422,
        0.74999999999999988897769753748434595763683319091797};
    (*new_lh) = new_lh_value2;
    merge_RACGT_RACGT<4>(seqregion1, total_blength_1, end_pos, *new_lh,
                         tree.model, tree.aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 2);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    total_blength_1 = 0; // either total_blength_1 = 0 or = -1 does not affect the merged region
    (*new_lh) = new_lh_value2;
    merge_RACGT_RACGT<4>(seqregion1, total_blength_1, end_pos, *new_lh,
                         tree.model, tree.aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 3);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge);
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    total_blength_1 = 1e-4;
    (*new_lh) = new_lh_value2;
    merge_RACGT_RACGT<4>(seqregion1, total_blength_1, end_pos, *new_lh,
                         tree.model, tree.aln, merged_regions);
    EXPECT_GT(merged_regions.size(), 3);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_FALSE(merged_regions.back().plength_observation2root != 0);
    EXPECT_TRUE(merged_regions.back().plength_observation2node == 0);
    SeqRegion::LHType new_lh_value_merge1{8.8777603721537831E-11,
        1.5545501848967752E-9,0.000021239862924016271,0.99997875849374817};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge1);
    // ----- Test 4 -----
    
    // ----- Test 5 -----
    seqregion1.plength_observation2root = 0.231; // neither plength_observation2root
    // nor plength_observation2node affect the merged region
    total_blength_1 = 1e-4;
    (*new_lh) = new_lh_value2;
    merge_RACGT_RACGT<4>(seqregion1, total_blength_1, end_pos, *new_lh,
                         tree.model, tree.aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 5);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge1);
    // ----- Test 5 -----
    
    /*// ----- Test 6 -----
    seqregion1.plength_observation2node = 1e-7; // neither plength_observation2root
     //  nor plength_observation2node affect the merged region
    total_blength_1 = 1e-4;
    (*new_lh) = new_lh_value2;
    merge_RACGT_RACGT<4>(seqregion1, total_blength_1, end_pos, *new_lh,
     tree.model, tree.aln, merged_regions);
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
    merge_RACGT_RACGT<4>(seqregion2, total_blength_1, end_pos, *new_lh,
     tree.model, tree.aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 7);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge2{0,1,0,0};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge2);
    // ----- Test 7 -----
    
    // ----- Test 8 -----
    (*new_lh) = new_lh_value2;
    merge_RACGT_RACGT<4>(seqregion2, total_blength_1, end_pos, *new_lh,
     tree.model, tree.aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 8);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge2);
    // ----- Test 8 -----
    
    // ----- Test 9 -----
    total_blength_1 = 0; // either total_blength_1 = 0 or = -1 does not affect the merged region
    (*new_lh) = new_lh_value2;
    merge_RACGT_RACGT<4>(seqregion2, total_blength_1, end_pos, *new_lh,
     tree.model, tree.aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 9);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge2);
    // ----- Test 9 -----
    
    // ----- Test 10 -----
    total_blength_1 = 1e-4;
    (*new_lh) = new_lh_value2;
    merge_RACGT_RACGT<4>(seqregion2, total_blength_1, end_pos, *new_lh,
     tree.model, tree.aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 10);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge3{7.0898429515507956E-7,
     0.11874116559913409,0.032309155233313958,0.84894897018325677};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge3);
    // ----- Test 10 -----
    
    // ----- Test 11 -----
    seqregion1.plength_observation2root = 0.231; // neither plength_observation2root
     // nor plength_observation2node affect the merged region
    total_blength_1 = 1e-4;
    (*new_lh) = new_lh_value2;
    merge_RACGT_RACGT<4>(seqregion2, total_blength_1, end_pos, *new_lh, tree.model,
     tree.aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 11);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge3);
    // ----- Test 11 -----
    
    // ----- Test 12 -----
    seqregion1.plength_observation2node = 1e-7; // neither plength_observation2root
     // nor plength_observation2node affect the merged region
    total_blength_1 = 1e-4;
    (*new_lh) = new_lh_value2;
    merge_RACGT_RACGT<4>(seqregion2, total_blength_1, end_pos, *new_lh,
     tree.model, tree.aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 12);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge3);
    // ----- Test 12 -----*/
}

/*
 Test merge_RACGT_ORACGT<4>(const SeqRegion& seq1_region, const SeqRegion& seq2_region,
 const RealNumType total_blength_1, const RealNumType total_blength_2,
 const RealNumType upper_plength, const PositionType end_pos, const RealNumType threshold_prob,
 const ModelBase* model, const Alignment& aln, SeqRegions& merged_regions)
 */
TEST(SeqRegions, merge_RACGT_ORACGT)
{
    Alignment aln = loadAln5K();
    Model model(cmaple::ModelBase::GTR);
    Tree tree(&aln, &model);
    std::unique_ptr<Params>& params = tree.params;
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    
    // ----- Test 1 -----
    SeqRegions merged_regions;
    const RealNumType threshold_prob = params->threshold_prob;
    RealNumType upper_plength = -1;
    RealNumType total_blength_1 = -1;
    RealNumType total_blength_2 = -1;
    const PositionType end_pos = 5432;
    SeqRegion seqregion1(TYPE_R, 2131);
    auto new_lh = cmaple::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value1{5.5759387092817078802929955938516570768115343526006e-06,
        0.49999442406129068761089229155913926661014556884766,
        5.5759387092817078802929955938516570768115343526006e-06,
        0.49999442406129068761089229155913926661014556884766};
    (*new_lh) = new_lh_value1;
    SeqRegion seqregion2(TYPE_O, 543, -1, -1, std::move(new_lh));
    merge_RACGT_ORACGT<4>(seqregion1, seqregion2, total_blength_1, total_blength_2,
        upper_plength, end_pos, threshold_prob, tree.model, tree.aln, merged_regions);
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
    merge_RACGT_ORACGT<4>(seqregion1, seqregion2, total_blength_1, total_blength_2,
        upper_plength, end_pos, threshold_prob, tree.model, tree.aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 2);
    EXPECT_FALSE(merged_regions.back().type != TYPE_O);
    EXPECT_LE(merged_regions.back().plength_observation2root, 0);
    EXPECT_GE(merged_regions.back().plength_observation2root, 0);
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
    merge_RACGT_ORACGT<4>(seqregion1, seqregion2, total_blength_1, total_blength_2,
        upper_plength, end_pos, threshold_prob, tree.model, tree.aln, merged_regions);
    EXPECT_TRUE(merged_regions.size() == 3);
    EXPECT_TRUE(merged_regions.back().type == TYPE_O);
    EXPECT_TRUE(merged_regions.back().plength_observation2root == 0);
    // EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge1{0.0000043308850656905248,
        0.11650995967890213,0.067956494829030142,0.81552921460700212};
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
    merge_RACGT_ORACGT<4>(seqregion1, seqregion2, total_blength_1, total_blength_2,
        upper_plength, end_pos, threshold_prob, tree.model, tree.aln, merged_regions);
    EXPECT_GT(merged_regions.size(), 3);
    EXPECT_TRUE(merged_regions.back().type == TYPE_O);
    /*EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);*/
    SeqRegion::LHType new_lh_value_merge2{3.7896075858290732E-7,
        0.0000019907504571924796,0.00022095280891100915,0.99977667747987319};
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
    merge_RACGT_ORACGT<4>(seqregion1, seqregion2, total_blength_1, total_blength_2,
        upper_plength, end_pos, threshold_prob, tree.model, tree.aln, merged_regions);
    EXPECT_FALSE(merged_regions.size() > 5);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_LE(merged_regions.back().plength_observation2node, 0);
    EXPECT_GE(merged_regions.back().plength_observation2node, 0);
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
    merge_RACGT_ORACGT<4>(seqregion1, seqregion2, total_blength_1, total_blength_2,
        upper_plength, end_pos, threshold_prob, tree.model, tree.aln, merged_regions);
    EXPECT_FALSE(merged_regions.size() < 6);
    EXPECT_FALSE(merged_regions.back().plength_observation2root != 0);
    EXPECT_TRUE(merged_regions.back().type == TYPE_O);
    // EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge4{1.8321558474427985E-7,
        2.6859491508084263E-8,0.99999903717347061,7.5275145311095654E-7};
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
    merge_RACGT_ORACGT<4>(seqregion1, seqregion2, total_blength_1, total_blength_2,
        upper_plength, end_pos, threshold_prob, tree.model, tree.aln, merged_regions);
    // EXPECT_EQ(merged_regions.size(), 7);
    EXPECT_TRUE(merged_regions.back().type == TYPE_O);
    // EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_TRUE(merged_regions.back().plength_observation2node == 0);
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
    merge_RACGT_ORACGT<4>(seqregion1, seqregion2, total_blength_1, total_blength_2,
        upper_plength, end_pos, threshold_prob, tree.model, tree.aln, merged_regions);
    // EXPECT_EQ(merged_regions.size(), 8);
    // EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_TRUE(merged_regions.back().plength_observation2root == 0);
    EXPECT_FALSE(merged_regions.back().plength_observation2node != 0);
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
    merge_RACGT_ORACGT<4>(seqregion1, seqregion2, total_blength_1, total_blength_2,
        upper_plength, end_pos, threshold_prob, tree.model, tree.aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 9);
    EXPECT_LE(merged_regions.back().type, TYPE_O);
    EXPECT_GE(merged_regions.back().type, TYPE_O);
    /*EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);*/
    SeqRegion::LHType new_lh_value_merge7{3.8512118129099356E-7,0.50512531494378898,
        3.852643996135802E-7,0.49487391467062997};
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
    merge_RACGT_ORACGT<4>(seqregion1, seqregion2, total_blength_1, total_blength_2,
        upper_plength, end_pos, threshold_prob, tree.model, tree.aln, merged_regions);
    EXPECT_GE(merged_regions.size(), 10);
    EXPECT_EQ(merged_regions.size(), 10);
    EXPECT_LE(merged_regions.back().type, TYPE_O);
    EXPECT_GE(merged_regions.back().type, TYPE_O);
    /*EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);*/
    SeqRegion::LHType new_lh_value_merge8{0,0,0,1};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge8);
    // ----- Test 10 -----
    
    /*// ----- Test 11 -----
    upper_plength = 0.123;
    total_blength_1 = -1;
    total_blength_2 = -1;
    seqregion1.type = TYPE_R;
    seqregion1.plength_observation2root = 0.432;
    seqregion1.plength_observation2node = 3e-6;
    seqregion2.type = TYPE_O;
    seqregion2.plength_observation2root = -1;
    seqregion2.plength_observation2node = -1;
    merge_RACGT_ORACGT<4>(seqregion1, seqregion2, total_blength_1, total_blength_2,
     upper_plength, end_pos, threshold_prob, tree.model, tree.aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 11);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge9{0.0000046465940568638611,0.1249996558599968,
     0.000011794893992852322,0.87498390265195358};
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
    merge_RACGT_ORACGT<4>(seqregion1, seqregion2, total_blength_1, total_blength_2,
     upper_plength, end_pos, threshold_prob, tree.model, tree.aln, merged_regions);
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
    merge_RACGT_ORACGT<4>(seqregion1, seqregion2, total_blength_1, total_blength_2,
     upper_plength, end_pos, threshold_prob, tree.model, tree.aln, merged_regions);
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
    merge_RACGT_ORACGT<4>(seqregion1, seqregion2, total_blength_1, total_blength_2,
     upper_plength, end_pos, threshold_prob, tree.model, tree.aln, merged_regions);
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
    merge_RACGT_ORACGT<4>(seqregion1, seqregion2, total_blength_1, total_blength_2,
     upper_plength, end_pos, threshold_prob, tree.model, tree.aln, merged_regions);
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
    merge_RACGT_ORACGT<4>(seqregion1, seqregion2, total_blength_1, total_blength_2,
     upper_plength, end_pos, threshold_prob, tree.model, tree.aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 16);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge13{0.000010173897266335654,0.12178798877322355,
     0.026555202184610452,0.85164663514489958};
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
    merge_RACGT_ORACGT<4>(seqregion1, seqregion2, total_blength_1, total_blength_2,
     upper_plength, end_pos, threshold_prob, tree.model, tree.aln, merged_regions);
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
    merge_RACGT_ORACGT<4>(seqregion1, seqregion2, total_blength_1, total_blength_2,
     upper_plength, end_pos, threshold_prob, tree.model, tree.aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 18);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge15{0.28592397557801535,0.30602769829658166,
     0.000011398356098642252,0.40803692776930434};
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
    merge_RACGT_ORACGT<4>(seqregion1, seqregion2, total_blength_1, total_blength_2,
     upper_plength, end_pos, threshold_prob, tree.model, tree.aln, merged_regions);
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
    merge_RACGT_ORACGT<4>(seqregion1, seqregion2, total_blength_1, total_blength_2,
     upper_plength, end_pos, threshold_prob, tree.model, tree.aln, merged_regions);
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
    merge_RACGT_ORACGT<4>(seqregion1, seqregion2, total_blength_1, total_blength_2,
     upper_plength, end_pos, threshold_prob, tree.model, tree.aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 21);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge17{3.7529836729309572E-7,0.0000019715097300387874,
     0.99011618892809094,0.0098814642638117549};
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
    merge_RACGT_ORACGT<4>(seqregion1, seqregion2, total_blength_1,
     total_blength_2, upper_plength, end_pos, threshold_prob, tree.model,
     tree.aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 22);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge18{1.9880125157816739E-9,
     0.99902052127670382,0.00097942098640511682,5.5748878587250403E-8};
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
    merge_RACGT_ORACGT<4>(seqregion1, seqregion2, total_blength_1,
     total_blength_2, upper_plength, end_pos, threshold_prob, tree.model,
     tree.aln, merged_regions);
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
    merge_RACGT_ORACGT<4>(seqregion1, seqregion2, total_blength_1,
     total_blength_2, upper_plength, end_pos, threshold_prob, tree.model, tree.aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 24);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge20{0.0000033509369075667414,
     0.059613997536883151,0.52309211091516972,0.41729054061103954};
    EXPECT_EQ(*merged_regions.back().likelihood, new_lh_value_merge20);
    // ----- Test 24 -----*/
}

/*
    Generate testing data (seqregions1, seqregions2, test_case)
 */
void genTestData(std::unique_ptr<SeqRegions>& seqregions1,
    std::unique_ptr<SeqRegions>& seqregions2, Tree& tree,
    const RealNumType& threshold_prob, const int& test_case)
{
    // dummy variables
    const PositionType seq_length = tree.aln->ref_seq.size();
    const StateType num_states = tree.aln->num_states;
    
    // pick a few sequences
    std::unique_ptr<SeqRegions> seqregions_1 =
    tree.aln->data[test_case].getLowerLhVector(seq_length, num_states, tree.aln->getSeqType());
    std::unique_ptr<SeqRegions> seqregions_2 =
    tree.aln->data[test_case * 10].getLowerLhVector(seq_length, num_states, tree.aln->getSeqType());
    std::unique_ptr<SeqRegions> seqregions_3 =
    tree.aln->data[test_case * 100].getLowerLhVector(seq_length, num_states, tree.aln->getSeqType());
    
    // compute the output regions
    seqregions_1->mergeTwoLowers<4>(seqregions1, 1e-5,
        *seqregions_2, 123e-3, tree.aln, tree.model, tree.cumulative_rate, threshold_prob);
    seqregions_2->mergeTwoLowers<4>(seqregions2, 214e-4,
        *seqregions_3, 13e-8, tree.aln, tree.model, tree.cumulative_rate, threshold_prob);
}

/*
 Test mergeUpperLower<4>(SeqRegions* &merged_regions, RealNumType upper_plength,
 const SeqRegions& lower_regions, RealNumType lower_plength, const
 Alignment& aln, const ModelBase* model, RealNumType threshold) const
 */
TEST(SeqRegions, mergeUpperLower)
{
    Alignment aln = loadAln5K();
    Model model(cmaple::ModelBase::GTR);
    Tree tree(&aln, &model);
    std::unique_ptr<Params>& params = tree.params;
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    
    std::unique_ptr<SeqRegions> seqregions1 = nullptr;
    std::unique_ptr<SeqRegions> seqregions2 = nullptr;
    std::unique_ptr<SeqRegions> merged_regions_ptr = nullptr;
    
    // dummy variables
    const RealNumType threshold_prob = params->threshold_prob;
    
    // ----- Test 1 -----
    genTestData(seqregions1, seqregions2, tree, threshold_prob, 1);
    
    seqregions1->mergeUpperLower<4>(merged_regions_ptr,
            3.3454886086112878360986772063867533688608091324568e-05,
            *seqregions2, -1, tree.aln, tree.model, threshold_prob);
    
    EXPECT_EQ(merged_regions_ptr->size(), 17);
    SeqRegion::LHType lh_value1_3{0.25393126472061372522759370440326165407896041870117,
        0.74606825145303112822858793151681311428546905517578,
        1.527175265677545701275923378109622419174229435157e-08,
        4.6855460255815036350456482400206326133229595143348e-07};
    EXPECT_EQ(*merged_regions_ptr->at(3).likelihood, lh_value1_3);
    SeqRegion::LHType lh_value1_11{5.7966916545679062813751932197858102169263361247431e-09,
        0.65834591260301689175093997619114816188812255859375,
        0.34165191092383423443479273373668547719717025756836,
        2.1706764573211083553527789291592853260226547718048e-06};
    EXPECT_EQ(*merged_regions_ptr->at(11).likelihood, lh_value1_11);
    SeqRegion::LHType lh_value1_15{4.819386899248702337066236070712340472388390821834e-10,
        4.4331834594549051557597359836046524428354587143986e-08,
        0.78291399129324290573350708655198104679584503173828,
        0.21708596389298387419053426583559485152363777160645};
    EXPECT_EQ(*merged_regions_ptr->at(15).likelihood, lh_value1_15);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    genTestData(seqregions1, seqregions2, tree, threshold_prob, 2);
    
    seqregions1->mergeUpperLower<4>(merged_regions_ptr,
        0.00023418420260279014175064382641267002327367663383484,
        *seqregions2, 3.3454886086112878360986772063867533688608091324568e-05,
        tree.aln, tree.model, threshold_prob);
    
    EXPECT_EQ(merged_regions_ptr->size(), 19);
    SeqRegion::LHType lh_value2_5{1.473143749210116913881069098529730254654168675188e-07,
        0.11305893785257446759739963226820691488683223724365,
        1.199119918869713445909431191738203636987236677669e-06,
        0.88693971571313168222872036494663916528224945068359};
    EXPECT_EQ(*merged_regions_ptr->at(5).likelihood, lh_value2_5);
    SeqRegion::LHType lh_value2_9{1.582745833551237004852788731179558112671656999737e-07,
        0.14837879648413093702785658933862578123807907104492,
        1.2671446230282412341420789775314759140201203990728e-06,
        0.85161977809666267180688237203867174685001373291016};
    EXPECT_EQ(*merged_regions_ptr->at(9).likelihood, lh_value2_9);
    SeqRegion::LHType lh_value2_11{0.18548170160962187957842672858532750979065895080566,
        7.9596768357620026274810943675563912336201610742137e-07,
        0.81451325466009993903071517706848680973052978515625,
        4.2477625948362153376989397424168259931320790201426e-06};
    EXPECT_EQ(*merged_regions_ptr->at(11).likelihood, lh_value2_11);
    SeqRegion::LHType lh_value2_13{1.582745833551237004852788731179558112671656999737e-07,
        0.14837879648413093702785658933862578123807907104492,
        1.2671446230282412341420789775314759140201203990728e-06,
        0.85161977809666267180688237203867174685001373291016};
    EXPECT_EQ(*merged_regions_ptr->at(13).likelihood, lh_value2_13);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    genTestData(seqregions1, seqregions2, tree, threshold_prob, 3);
    
    seqregions1->mergeUpperLower<4>(merged_regions_ptr,
        0.00020072931651667725661339347631439977703848853707314,
        *seqregions2, 0.00013381954434445151344394708825547013475443236529827,
                                    tree.aln, tree.model, threshold_prob);
    
    EXPECT_EQ(merged_regions_ptr->size(), 19);
    SeqRegion::LHType lh_value3_5{0.55088461186082993137347330048214644193649291992188,
        2.1517430823900132296627729644455939705949276685715e-06,
        0.44910173316225798778589251014636829495429992675781,
        1.1503233829685540224869837178101050767509150318801e-05};
    EXPECT_EQ(*merged_regions_ptr->at(5).likelihood, lh_value3_5);
    SeqRegion::LHType lh_value3_7{4.3642717210825648466670069816619736968732468085364e-07,
        0.41287703560443073103058964079536963254213333129883,
        3.4936313425953077532650631331634372145344968885183e-06,
        0.58711903433705459054436914811958558857440948486328};
    EXPECT_EQ(*merged_regions_ptr->at(7).likelihood, lh_value3_7);
    SeqRegion::LHType lh_value3_9{4.3642717210825648466670069816619736968732468085364e-07,
        0.41287703560443073103058964079536963254213333129883,
        3.4936313425953077532650631331634372145344968885183e-06,
        0.58711903433705459054436914811958558857440948486328};
    EXPECT_EQ(*merged_regions_ptr->at(9).likelihood, lh_value3_9);
    SeqRegion::LHType lh_value3_15{1.5630903977336201439418812298435190397127847461434e-10,
        0.99999996446600269983662201411789283156394958496094,
        2.4311281533080620387737091640171886719468119508747e-10,
        3.5134575447755886574289406693374915313654582860181e-08};
    EXPECT_EQ(*merged_regions_ptr->at(15).likelihood, lh_value3_15);
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    genTestData(seqregions1, seqregions2, tree, threshold_prob, 4);
    
    seqregions1->mergeUpperLower<4>(merged_regions_ptr,
        3.3454886086112878360986772063867533688608091324568e-05,
        *seqregions2, -1, tree.aln, tree.model, threshold_prob);
    
    EXPECT_EQ(merged_regions_ptr->size(), 22);
    SeqRegion::LHType lh_value4_2{4.4313503178746672948596712557763023163093230039067e-11,
        0.18034214793852129665729933094553416594862937927246,
        5.4085351209201087649692306987166123821708652030793e-10,
        0.81965785147631176066340685792965814471244812011719};
    EXPECT_EQ(*merged_regions_ptr->at(2).likelihood, lh_value4_2);
    SeqRegion::LHType lh_value4_12{4.4313503178746672948596712557763023163093230039067e-11,
        0.18034214793852129665729933094553416594862937927246,
        5.4085351209201087649692306987166123821708652030793e-10,
        0.81965785147631176066340685792965814471244812011719};
    EXPECT_EQ(*merged_regions_ptr->at(12).likelihood, lh_value4_12);
    SeqRegion::LHType lh_value4_14{0.37269372979002129975256707439257297664880752563477,
        7.6247474514365417954404194477596203027847110433868e-10,
        0.62730626634354158532858036778634414076805114746094,
        3.1039623776691381636330143009225301931053309090203e-09};
    EXPECT_EQ(*merged_regions_ptr->at(14).likelihood, lh_value4_14);
    SeqRegion::LHType lh_value4_16{4.4313503178746672948596712557763023163093230039067e-11,
        0.18034214793852129665729933094553416594862937927246,
        5.4085351209201087649692306987166123821708652030793e-10,
        0.81965785147631176066340685792965814471244812011719};
    EXPECT_EQ(*merged_regions_ptr->at(16).likelihood, lh_value4_16);
    
    // ----- Test 4 -----
    
    // ----- Test 5 -----
    genTestData(seqregions1, seqregions2, tree, threshold_prob, 5);
    
    seqregions1->mergeUpperLower<4>(merged_regions_ptr, -1, *seqregions2,
                5.0182329129169314153348369078599944259622134268284e-05,
                                    tree.aln, tree.model, threshold_prob);
    
    EXPECT_EQ(merged_regions_ptr->size(), 15);
    SeqRegion::LHType lh_value5_6{1.2467371151780668072757457637517002069227345373292e-09,
        0.64033070070634878767634745599934831261634826660156,
        1.2180634035447382694533234245848341004148096544668e-07, 0.3596691762405736514374154921824811026453971862793};
    EXPECT_EQ(*merged_regions_ptr->at(6).likelihood, lh_value5_6);
    
    // ----- Test 5 -----
    
    /*// ----- Test 6 -----
    genTestData(seqregions1, seqregions2, tree, threshold_prob, 6);
    
    seqregions1->mergeUpperLower<4>(merged_regions_ptr,
     3.3454886086112878360986772063867533688608091324568e-05, *seqregions2,
     -1, tree.aln, tree.model, threshold_prob);
    
    EXPECT_EQ(merged_regions_ptr->size(), 21);
    SeqRegion::LHType lh_value6_10{4.4313503178746672948596712557763023163093230039067e-11,
     0.18034214793852129665729933094553416594862937927246,
     5.4085351209201087649692306987166123821708652030793e-10,
     0.81965785147631176066340685792965814471244812011719};
    EXPECT_EQ(*merged_regions_ptr->at(10).likelihood, lh_value6_10);
    SeqRegion::LHType lh_value6_12{0.37269372979002129975256707439257297664880752563477,
     7.6247474514365417954404194477596203027847110433868e-10,
     0.62730626634354158532858036778634414076805114746094,
     3.1039623776691381636330143009225301931053309090203e-09};
    EXPECT_EQ(*merged_regions_ptr->at(12).likelihood, lh_value6_12);
    SeqRegion::LHType lh_value6_14{4.4313503178746672948596712557763023163093230039067e-11,
     0.18034214793852129665729933094553416594862937927246,
     5.4085351209201087649692306987166123821708652030793e-10,
     0.81965785147631176066340685792965814471244812011719};
    EXPECT_EQ(*merged_regions_ptr->at(14).likelihood, lh_value6_14);
    
    // ----- Test 6 -----
    
    // ----- Test 7 -----
    genTestData(seqregions1, seqregions2, tree, threshold_prob, 7);
    
    seqregions1->mergeUpperLower<4>(merged_regions_ptr, -1, *seqregions2,
     -1, tree.aln, tree.model, threshold_prob);
    
    EXPECT_EQ(merged_regions_ptr->size(), 14);
    
    // ----- Test 7 -----
    
    // ----- Test 8 -----
    genTestData(seqregions1, seqregions2, tree, threshold_prob, 8);
    
    seqregions1->mergeUpperLower<4>(merged_regions_ptr,
     3.3454886086112878360986772063867533688608091324568e-05,
     *seqregions2, 3.3454886086112878360986772063867533688608091324568e-05,
     tree.aln, tree.model, threshold_prob);
    
    EXPECT_EQ(merged_regions_ptr->size(), 20);
    SeqRegion::LHType lh_value8_6{6.6049407284407675729318194790540275107559864409268e-08,
     0.35179317649038233106750794831896200776100158691406,
     5.9021711114530172845890707705729383292236889246851e-07,
     0.64820616724309920719804267719155177474021911621094};
    EXPECT_EQ(*merged_regions_ptr->at(6).likelihood, lh_value8_6);
    SeqRegion::LHType lh_value8_10{8.3716410011684010534938212531874679456223020679317e-08,
     0.54955073431965328900616896135034039616584777832031,
     6.7023268396491381295707089396640476763877813937142e-07,
     0.45044851173125272092434556725493166595697402954102};
    EXPECT_EQ(*merged_regions_ptr->at(10).likelihood, lh_value8_10);
    SeqRegion::LHType lh_value8_12{0.61451650198388751977773836188134737312793731689453,
     3.7670413001044542954150596940354756014812664943747e-07,
     0.38548111099203380414124353592342231422662734985352,
     2.0103199486559325614520952335562142820890585426241e-06};
    EXPECT_EQ(*merged_regions_ptr->at(12).likelihood, lh_value8_12);
    // ----- Test 8 -----
    
    // ----- Test 9 -----
    genTestData(seqregions1, seqregions2, tree, threshold_prob, 9);
    
    seqregions1->mergeUpperLower<4>(merged_regions_ptr,
     3.3454886086112878360986772063867533688608091324568e-05,
     *seqregions2, 6.6909772172225756721973544127735067377216182649136e-05,
     tree.aln, tree.model, threshold_prob);
    
    EXPECT_EQ(merged_regions_ptr->size(), 23);
    SeqRegion::LHType lh_value9_5{9.7723101777417999379181043091008307399647492275108e-08,
     0.52049400296263548248987262923037633299827575683594,
     8.7325305698605406513781827804177027019250090233982e-07,
     0.47950502606120570803227565193083137273788452148438};
    EXPECT_EQ(*merged_regions_ptr->at(5).likelihood, lh_value9_5);
    SeqRegion::LHType lh_value9_9{9.7723101777417999379181043091008307399647492275108e-08,
     0.52049400296263548248987262923037633299827575683594,
     8.7325305698605406513781827804177027019250090233982e-07,
     0.47950502606120570803227565193083137273788452148438};
    EXPECT_EQ(*merged_regions_ptr->at(9).likelihood, lh_value9_9);
    SeqRegion::LHType lh_value9_11{0.10386392893018993321962994968998827971518039703369,
     1.2550641926514089704187146848135547827496338868514e-07,
     0.89613510914276239827103154311771504580974578857422,
     8.3642062853197228917663877401089678187418030574918e-07};
    EXPECT_EQ(*merged_regions_ptr->at(11).likelihood, lh_value9_11);
    SeqRegion::LHType lh_value9_13{9.7723101777417999379181043091008307399647492275108e-08,
     0.52049400296263548248987262923037633299827575683594,
     8.7325305698605406513781827804177027019250090233982e-07,
     0.47950502606120570803227565193083137273788452148438};
    EXPECT_EQ(*merged_regions_ptr->at(13).likelihood, lh_value9_13);
    SeqRegion::LHType lh_value9_15{5.2219032820157223355482035195146428563361951091792e-08,
     0.72557874367873564924735774184227921068668365478516,
     4.3094360831706182765808124841833137708135836874135e-07,
     0.27442077315862339892404975216777529567480087280273};
    EXPECT_EQ(*merged_regions_ptr->at(15).likelihood, lh_value9_15);
    SeqRegion::LHType lh_value9_19{0.81584940948900563917334238794865086674690246582031,
     1.2912044068271519354942737656255502542990143410861e-05,
     1.9781469095234067875855796247996920556033728644252e-05,
     0.18411789699783087659312741379835642874240875244141};
    EXPECT_EQ(*merged_regions_ptr->at(19).likelihood, lh_value9_19);
    
    // ----- Test 9 -----
    
    // ----- Test 10 -----
    genTestData(seqregions1, seqregions2, tree, threshold_prob, 10);
    
    seqregions1->mergeUpperLower<4>(merged_regions_ptr, -1, *seqregions2, -1,
     tree.aln, tree.model, threshold_prob);
    
    EXPECT_EQ(merged_regions_ptr->size(), 11);
    // ----- Test 10 -----*/
}

/*
 Test merge_N_O_TwoLowers(const SeqRegion& seq2_region, const PositionType end_pos,
 const RealNumType plength2, SeqRegions& merged_regions)
 */
TEST(SeqRegions, merge_N_O_TwoLowers)
{
    Alignment aln = loadAln5K();
    Model model(cmaple::ModelBase::GTR);
    Tree tree(&aln, &model);
    std::unique_ptr<Params>& params = tree.params;
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    
    // ----- Test 1 -----
    SeqRegions merged_regions;
    const PositionType end_pos = 8623;
    RealNumType plength2 = -1;
    auto new_lh = cmaple::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value{0.27595781923293211113090706021466758102178573608398,
        8.3817209383128356075479820086471249851456377655268e-06,
        0.72401575181023847260775028189527802169322967529297,
        1.8047235891267821277644811672757896303664892911911e-05};
    (*new_lh) = new_lh_value;
    SeqRegion seqregion1(TYPE_O, 243, -1, -1, std::move(new_lh));
    merge_N_O_TwoLowers(seqregion1, end_pos, plength2, merged_regions);
    EXPECT_TRUE(merged_regions.size() == 1);
    EXPECT_GE(merged_regions.back().plength_observation2root, -1);
    EXPECT_LE(merged_regions.back().plength_observation2node, -1);
    EXPECT_TRUE(*merged_regions.back().likelihood == new_lh_value);
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
    
    /*// ----- Test 5 -----
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
    // ----- Test 10 -----*/
}

/*
 Test merge_N_RACGT_TwoLowers(const SeqRegion& seq2_region, const PositionType end_pos,
 const RealNumType plength2, const RealNumType threshold_prob, SeqRegions& merged_regions)
 */
TEST(SeqRegions, merge_N_RACGT_TwoLowers)
{
    Alignment aln = loadAln5K();
    Model model(cmaple::ModelBase::GTR);
    Tree tree(&aln, &model);
    std::unique_ptr<Params>& params = tree.params;
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    const RealNumType threshold_prob = params->threshold_prob;
    
    // ----- Test 1 -----
    SeqRegions merged_regions;
    const PositionType end_pos = 6543;
    RealNumType plength2 = -1;
    SeqRegion seqregion1(TYPE_R, 3442, -1, -1);
    merge_N_RACGT_TwoLowers(seqregion1, end_pos, plength2, threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 1);
    EXPECT_EQ(merged_regions[merged_regions.size() - 1].plength_observation2root, -1);
    EXPECT_EQ(merged_regions[merged_regions.size() - 1].plength_observation2node, -1);
    EXPECT_EQ(merged_regions[merged_regions.size() - 1].type, seqregion1.type);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    plength2 = -1;
    seqregion1.plength_observation2node = -1;
    seqregion1.plength_observation2root = 0;
    seqregion1.type = 1; // C
    merge_N_RACGT_TwoLowers(seqregion1, end_pos, plength2, threshold_prob, merged_regions);
    EXPECT_LT(merged_regions.size(), 3);
    EXPECT_EQ(merged_regions[merged_regions.size() - 1].plength_observation2root, -1);
    EXPECT_TRUE(merged_regions.back().plength_observation2node == -1);
    EXPECT_EQ(merged_regions[merged_regions.size() - 1].type, seqregion1.type);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    plength2 = -1;
    seqregion1.plength_observation2node = -1;
    seqregion1.plength_observation2root = 0.3231;
    seqregion1.type = TYPE_R;
    merge_N_RACGT_TwoLowers(seqregion1, end_pos, plength2, threshold_prob, merged_regions);
    EXPECT_EQ(merged_regions.size(), 3);
    EXPECT_FALSE(merged_regions[merged_regions.size() - 1].plength_observation2root != -1);
    EXPECT_EQ(merged_regions.back().plength_observation2node, -1);
    EXPECT_TRUE(merged_regions[merged_regions.size() - 1].type == seqregion1.type);
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
    EXPECT_TRUE(merged_regions.size() == 10);
    EXPECT_TRUE(merged_regions[merged_regions.size() - 1].plength_observation2root == -1);
    EXPECT_TRUE(merged_regions.back().plength_observation2node == -1);
    EXPECT_EQ(merged_regions[merged_regions.size() - 1].type, seqregion1.type);
    // ----- Test 10 -----
    
    /*// ----- Test 11 -----
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
    // ----- Test 27 -----*/
}

/*
 Test merge_identicalRACGT_TwoLowers(const SeqRegion& seq1_region, const PositionType end_pos,
 RealNumType total_blength_1, RealNumType total_blength_2, const PositionType pos,
 const RealNumType threshold_prob, const ModelBase* model, RealNumType &log_lh,
 SeqRegions& merged_regions, const bool return_log_lh)
 */
TEST(SeqRegions, merge_identicalRACGT_TwoLowers)
{
    Alignment aln = loadAln5K();
    Model model(cmaple::ModelBase::GTR);
    Tree tree(&aln, &model);
    std::unique_ptr<Params>& params = tree.params;
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    const RealNumType threshold_prob = params->threshold_prob;
    
    // ----- Test 1 -----
    SeqRegions merged_regions;
    const PositionType end_pos = 6543;
    PositionType pos = 6489;
    RealNumType total_blength_1 = -1;
    RealNumType total_blength_2 = -1;
    RealNumType log_lh = 0;
    SeqRegion seqregion1(TYPE_R, 3442, -1, -1);
    merge_identicalRACGT_TwoLowers(seqregion1, end_pos, total_blength_1,
        total_blength_2, pos, threshold_prob, tree.model,
        tree.cumulative_rate, log_lh, merged_regions, true);
    EXPECT_LT(merged_regions.size(), 2);
    EXPECT_LT(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions[merged_regions.size() - 1].plength_observation2node, -1);
    EXPECT_FALSE(merged_regions.back().type != seqregion1.type);
    EXPECT_EQ(log_lh, 0);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    total_blength_1 = 0;
    total_blength_2 = 0;
    log_lh = 0;
    seqregion1.type = 0;
    merge_identicalRACGT_TwoLowers(seqregion1, end_pos, total_blength_1, total_blength_2,
        pos, threshold_prob, tree.model, tree.cumulative_rate,
                                   log_lh, merged_regions, true);
    EXPECT_GE(merged_regions.size(), 2);
    EXPECT_GE(merged_regions[merged_regions.size() - 1].plength_observation2root, -1);
    EXPECT_GE(merged_regions[merged_regions.size() - 1].plength_observation2node, -1);
    EXPECT_LE(merged_regions[merged_regions.size() - 1].type, seqregion1.type);
    EXPECT_EQ(log_lh, 0);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    total_blength_1 = -1;
    total_blength_2 = 13e-4;
    log_lh = 0;
    seqregion1.type = 3;
    merge_identicalRACGT_TwoLowers(seqregion1, end_pos, total_blength_1,
        total_blength_2, pos, threshold_prob, tree.model,
        tree.cumulative_rate, log_lh, merged_regions, true);
    EXPECT_LE(merged_regions.size(), 3);
    EXPECT_LT(merged_regions[merged_regions.size() - 1].plength_observation2root, 0);
    EXPECT_LT(merged_regions[merged_regions.size() - 1].plength_observation2node, 0);
    EXPECT_GE(merged_regions[merged_regions.size() - 1].type, seqregion1.type);
    EXPECT_EQ(log_lh, -0.0016388829346472356);
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    total_blength_1 = 21e-6;
    total_blength_2 = -1;
    log_lh = 0;
    seqregion1.type = TYPE_R;
    merge_identicalRACGT_TwoLowers(seqregion1, end_pos, total_blength_1,
        total_blength_2, pos, threshold_prob, tree.model,
        tree.cumulative_rate, log_lh, merged_regions, true);
    EXPECT_GT(merged_regions.size(), 3);
    EXPECT_TRUE(merged_regions[merged_regions.size() - 1].plength_observation2root == -1);
    EXPECT_FALSE(merged_regions.back().plength_observation2node != -1);
    EXPECT_GE(merged_regions[merged_regions.size() - 1].type, seqregion1.type);
    EXPECT_EQ(log_lh, -0.00099527007095828778);
    // ----- Test 4 -----
    
    // ----- Test 5 -----
    total_blength_1 = 1e-5;
    total_blength_2 = 0.04231;
    log_lh = 0;
    seqregion1.type = 1;
    merge_identicalRACGT_TwoLowers(seqregion1, end_pos, total_blength_1,
        total_blength_2, pos, threshold_prob, tree.model,
        tree.cumulative_rate, log_lh, merged_regions, true);
    EXPECT_LT(merged_regions.size(), 6);
    EXPECT_EQ(merged_regions.back().plength_observation2root, -1);
    EXPECT_FALSE(merged_regions[merged_regions.size() - 1].plength_observation2node != -1);
    EXPECT_FALSE(merged_regions[merged_regions.size() - 1].type != seqregion1.type);
    EXPECT_EQ(log_lh, -0.067217084471974248);
    // ----- Test 5 -----
}

/*
 Test merge_O_O_TwoLowers<4>(const SeqRegion& seq2_region, RealNumType total_blength_2,
 const PositionType end_pos, const Alignment& aln, const ModelBase* model,
 const RealNumType threshold_prob, RealNumType &log_lh, SeqRegion::LHType& new_lh,
 SeqRegions* merged_regions, const bool return_log_lh)
 */
TEST(SeqRegions, merge_O_O_TwoLowers)
{
    Alignment aln = loadAln5K();
    Model model(cmaple::ModelBase::GTR);
    Tree tree(&aln, &model);
    std::unique_ptr<Params>& params = tree.params;
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    const RealNumType threshold_prob = params->threshold_prob;
    
    // ----- Test 1 -----
    std::unique_ptr<SeqRegions> merged_regions = cmaple::make_unique<SeqRegions>();
    RealNumType log_lh = 0;
    auto new_lh = cmaple::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value{0.12497560486006949187487435892762732692062854766846,
        9.7580559722090576469559833339140197949745925143361e-06,
        0.87500487902798607109389195102266967296600341796875,
        9.7580559722090576469559833339140197949745925143361e-06};
    (*new_lh) = new_lh_value;
    SeqRegion seqregion1(TYPE_O, 243, -1, -1, std::move(new_lh));
    auto new_lh1 = cmaple::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value1{0.59998393993535814594508792652050033211708068847656,
        2.6766774402933635499746492514283602304203668609262e-05,
        2.6766774402933635499746492514283602304203668609262e-05,
        0.39996252651583585890904259940725751221179962158203};
    (*new_lh1) = new_lh_value1;
    RealNumType total_blength_2 = -1;
    const PositionType end_pos = 3213;
    EXPECT_TRUE(merge_O_O_TwoLowers<4>(seqregion1, total_blength_2,
                end_pos, tree.aln, tree.model, threshold_prob,
                    log_lh, *new_lh1, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 1);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge{0.99963572952385693,
        3.4820599267114684E-9,0.00031223631362821728,
        0.000052030680454886103};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge);
    EXPECT_EQ(log_lh, -2.5901247759055703);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    log_lh = 0;
    total_blength_2 = 1e-5;
    (*new_lh1) = new_lh_value1;
    EXPECT_TRUE(merge_O_O_TwoLowers<4>(seqregion1, total_blength_2,
                end_pos, tree.aln, tree.model, threshold_prob,
                log_lh, *new_lh1, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 2);
    EXPECT_TRUE(merged_regions->back().type == TYPE_O);
    /*EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);*/
    SeqRegion::LHType new_lh_value_merge1{0.99961708509029623,
        3.8289366877191551E-9,0.00031222411049282615,
        0.000070686970274545192};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge1);
    EXPECT_EQ(log_lh, -2.5900955748485903);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    log_lh = 0;
    total_blength_2 = 0;
    (*new_lh1) = new_lh_value1;
    EXPECT_TRUE(merge_O_O_TwoLowers<4>(seqregion1, total_blength_2, end_pos,
                tree.aln, tree.model, threshold_prob, log_lh,
                *new_lh1, merged_regions, true));
    EXPECT_TRUE(merged_regions->size() == 3);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    // EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_FALSE(merged_regions->back().plength_observation2node != 0);
    SeqRegion::LHType new_lh_value_merge2{0.99963572952385693,
        3.4820599267114684E-9,0.00031223631362821728,0.000052030680454886103};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge2);
    EXPECT_EQ(log_lh, -2.5901247759055703);
    // ----- Test 3 -----
    
    /*// ----- Test 4 -----
    log_lh = 0;
    total_blength_2 = -1;
    SeqRegion::LHType new_lh_value2{1.1151877415789591228433182135137968771232408471406e-05,
     1.2436575683940661210420103588707597258578019250308e-10,
     0.99998884787385255989988763758447021245956420898438,
     1.2436575683940661210420103588707597258578019250308e-10};
    (*new_lh1) = new_lh_value2;
    EXPECT_TRUE(merge_O_O_TwoLowers<4>(seqregion1, total_blength_2,
     end_pos, tree.aln, tree.model, threshold_prob, log_lh, *new_lh1,
     merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 4);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge3{0.0000015928207737367849,
     1.3869404003893418E-15,0.99999840717922339,1.3869404003893418E-15};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge3);
    EXPECT_EQ(log_lh, -0.1335353759743724);
    // ----- Test 4 -----
    
    // ----- Test 5 -----
    log_lh = 0;
    total_blength_2 = 0;
    SeqRegion::LHType new_lh_value3{0.80000178432426516383912940000300295650959014892578,
     8.9216213262436184458807272856795123061601771041751e-06,
     0.19998037243308225408000566858390811830759048461914,
     8.9216213262436184458807272856795123061601771041751e-06};
    (*new_lh1) = new_lh_value3;
    EXPECT_TRUE(merge_O_O_TwoLowers<4>(seqregion1, total_blength_2,
     end_pos, tree.aln, tree.model, threshold_prob, log_lh, *new_lh1,
     merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 5);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge4{0.36361313457027911,
     3.1661424484350952E-10,0.63638686479649231,3.1661424484350952E-10};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge4);
    EXPECT_EQ(log_lh, -1.2911132491064328);
    // ----- Test 5 -----
    
    // ----- Test 6 -----
    log_lh = 0;
    total_blength_2 = 121e-5;
    SeqRegion::LHType new_lh_value4{0.99999721161091070786852696983260102570056915283203,
     2.7982277946960444848486499948455024505689081593118e-10,
     2.7878294439446026602288583595701254580490058287978e-06,
     2.7982277946960444848486499948455024505689081593118e-10};
    (*new_lh1) = new_lh_value4;
    EXPECT_TRUE(merge_O_O_TwoLowers<4>(seqregion1, total_blength_2,
     end_pos, tree.aln, tree.model, threshold_prob, log_lh,
     *new_lh1, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 6);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge5{0.99998052979430418,
     2.8492225625088802E-13,0.000019470204442126495,
     9.6862198410969555E-13};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge5);
    EXPECT_EQ(log_lh, -2.0783443388603096);
    // ----- Test 6 -----
    
    // ----- Test 7 -----
    log_lh = 0;
    total_blength_2 = 121e-5;
    SeqRegion::LHType new_lh_value5{0.99999442305078967141440671184682287275791168212891,
     5.5758298847802999098784669518291678969035274349153e-06,
     5.5966272240860645628430121287970669397004996881151e-10,
     5.5966272240860645628430121287970669397004996881151e-10};
    *seqregion1.likelihood = new_lh_value5;
    EXPECT_TRUE(merge_O_O_TwoLowers<4>(seqregion1, total_blength_2,
     end_pos, tree.aln, tree.model, threshold_prob, log_lh, *new_lh1,
     merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 7);
    EXPECT_EQ(merged_regions->back().type, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -0.00043445905798766815652);
    // ----- Test 7 -----
    
    // ----- Test 8 -----
    log_lh = 0;
    total_blength_2 = -1;
    SeqRegion::LHType new_lh_value6{0.80000044606912734668213715849560685455799102783203,
     2.2303456366633140507259747825630213924341660458595e-06,
     2.2303456366633140507259747825630213924341660458595e-06,
     0.1999950932395993530299449503218056634068489074707};
    *seqregion1.likelihood = new_lh_value6;
    EXPECT_TRUE(merge_O_O_TwoLowers<4>(seqregion1, total_blength_2, end_pos,
     tree.aln, tree.model, threshold_prob, log_lh, *new_lh1, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 8);
    EXPECT_EQ(merged_regions->back().type, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    SeqRegion::LHType new_lh_value_merge7{0.12684809781766063,
     1.3059895886789165E-10,0.87315190141833243,6.3340794416577515E-10};
    EXPECT_EQ(log_lh, -0.22314300087916599802);
    // ----- Test 8 -----*/
}

/*
 Test  merge_O_RACGT_TwoLowers<4>(const SeqRegion& seq2_region,
 RealNumType total_blength_2, const PositionType end_pos,
 const Alignment& aln, const ModelBase* model, const RealNumType threshold_prob,
 RealNumType &log_lh, SeqRegion::LHType& new_lh, RealNumType& sum_lh,
 SeqRegions* merged_regions, const bool return_log_lh)
 */
TEST(SeqRegions, merge_O_RACGT_TwoLowers)
{
    Alignment aln = loadAln5K();
    Model model(cmaple::ModelBase::GTR);
    Tree tree(&aln, &model);
    std::unique_ptr<Params>& params = tree.params;
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    const RealNumType threshold_prob = params->threshold_prob;
    
    // ----- Test 1 -----
    std::unique_ptr<SeqRegions> merged_regions = cmaple::make_unique<SeqRegions>();
    RealNumType log_lh = 0;
    RealNumType sum_lh = 0;
    SeqRegion seqregion1(TYPE_R, 243, -1, -1);
    auto new_lh1 = cmaple::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value1{0.59998393993535814594508792652050033211708068847656,
        2.6766774402933635499746492514283602304203668609262e-05,
        2.6766774402933635499746492514283602304203668609262e-05,
        0.39996252651583585890904259940725751221179962158203};
    (*new_lh1) = new_lh_value1;
    RealNumType total_blength_2 = -1;
    const PositionType end_pos = 243;
    EXPECT_TRUE(merge_O_RACGT_TwoLowers<4>(seqregion1, total_blength_2,
                end_pos, tree.aln, tree.model, threshold_prob, log_lh,
                *new_lh1, sum_lh, merged_regions, true));
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
    EXPECT_TRUE(merge_O_RACGT_TwoLowers<4>(seqregion1, total_blength_2,
                end_pos, tree.aln, tree.model, threshold_prob, log_lh,
                        *new_lh1, sum_lh, merged_regions, true));
    // EXPECT_EQ(merged_regions->size(), 2);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_FALSE(merged_regions->back().plength_observation2root != 0);
    // EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge1{0.011816164293335018,
        0.88299798641181038,8.0375755307632964E-7,0.10518504553730158};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge1);
    EXPECT_EQ(log_lh, -10.403932725076588);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    log_lh = 0;
    total_blength_2 = 0;
    seqregion1.type = 0;
    (*new_lh1) = new_lh_value1;
    EXPECT_TRUE(merge_O_RACGT_TwoLowers<4>(seqregion1, total_blength_2,
                end_pos, tree.aln, tree.model, threshold_prob,
                log_lh, *new_lh1, sum_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 3);
    EXPECT_TRUE(merged_regions->back().type == 0);
    EXPECT_EQ(log_lh, -0.510852390898630326354634689778);
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    log_lh = 0;
    total_blength_2 = -1;
    seqregion1.type = 2;
    SeqRegion::LHType new_lh_value2{1.1151877415789591228433182135137968771232408471406e-05,
        1.2436575683940661210420103588707597258578019250308e-10,
        0.99998884787385255989988763758447021245956420898438,
        1.2436575683940661210420103588707597258578019250308e-10};
    (*new_lh1) = new_lh_value2;
    EXPECT_TRUE(merge_O_RACGT_TwoLowers<4>(seqregion1, total_blength_2,
                end_pos, tree.aln, tree.model, threshold_prob, log_lh,
                *new_lh1, sum_lh, merged_regions, true));
    EXPECT_FALSE(merged_regions->size() != 4);
    EXPECT_TRUE(merged_regions->back().type == 2);
    EXPECT_TRUE(merged_regions->back().plength_observation2root == -1);
    EXPECT_TRUE(log_lh == -1.1152188332861237853011436571559755748239695094526e-05);
    // ----- Test 4 -----
    
    /*// ----- Test 5 -----
    log_lh = 0;
    total_blength_2 = 0;
    seqregion1.type = TYPE_R;
    SeqRegion::LHType new_lh_value3{0.80000178432426516383912940000300295650959014892578,
     8.9216213262436184458807272856795123061601771041751e-06,
     0.19998037243308225408000566858390811830759048461914,
     8.9216213262436184458807272856795123061601771041751e-06};
    (*new_lh1) = new_lh_value3;
    EXPECT_TRUE(merge_O_RACGT_TwoLowers<4>(seqregion1, total_blength_2,
     end_pos, tree.aln, tree.model, threshold_prob, log_lh, *new_lh1,
     sum_lh, merged_regions, true));
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
    SeqRegion::LHType new_lh_value4{0.99999721161091070786852696983260102570056915283203,
     2.7982277946960444848486499948455024505689081593118e-10,
     2.7878294439446026602288583595701254580490058287978e-06,
     2.7982277946960444848486499948455024505689081593118e-10};
    (*new_lh1) = new_lh_value4;
    EXPECT_TRUE(merge_O_RACGT_TwoLowers<4>(seqregion1, total_blength_2,
     end_pos, tree.aln, tree.model, threshold_prob, log_lh, *new_lh1,
     sum_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 6);
    EXPECT_EQ(merged_regions->back().type, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -0.0004122066213375955);
    // ----- Test 6 -----
    
    // ----- Test 7 -----
    log_lh = 0;
    total_blength_2 = 121e-5;
    seqregion1.type = 3;
    EXPECT_TRUE(merge_O_RACGT_TwoLowers<4>(seqregion1, total_blength_2,
     end_pos, tree.aln, tree.model, threshold_prob, log_lh, *new_lh1,
     sum_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 7);
    EXPECT_EQ(merged_regions->back().type, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -9.2478945382351511);
    // ----- Test 7 -----
    
    // ----- Test 8 -----
    log_lh = 0;
    total_blength_2 = -1;
    seqregion1.type = TYPE_R;
    EXPECT_TRUE(merge_O_RACGT_TwoLowers<4>(seqregion1, total_blength_2,
     end_pos, tree.aln, tree.model, threshold_prob, log_lh, *new_lh1,
     sum_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 8);
    EXPECT_EQ(merged_regions->back().type, TYPE_R);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    SeqRegion::LHType new_lh_value_merge7{0.12684809781766063,
     1.3059895886789165E-10,0.87315190141833243,6.3340794416577515E-10};
    EXPECT_EQ(log_lh, -28.181320607597758);
    // ----- Test 8 -----*/
}

/*
 Test merge_O_ORACGT_TwoLowers<4>(const SeqRegion& seq1_region,
 const SeqRegion& seq2_region, RealNumType total_blength_1,
 RealNumType total_blength_2, const PositionType end_pos, const Alignment& aln,
 const ModelBase* model, const RealNumType threshold_prob,
 RealNumType &log_lh, SeqRegions* merged_regions, const bool return_log_lh)
 */
TEST(SeqRegions, merge_O_ORACGT_TwoLowers)
{
    Alignment aln = loadAln5K();
    Model model(cmaple::ModelBase::GTR);
    Tree tree(&aln, &model);
    std::unique_ptr<Params>& params = tree.params;
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    const RealNumType threshold_prob = params->threshold_prob;
    
    // ----- Test 1 -----
    std::unique_ptr<SeqRegions> merged_regions = cmaple::make_unique<SeqRegions>();
    RealNumType log_lh = 0;
    auto new_lh = cmaple::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value{0.12497560486006949187487435892762732692062854766846,
        9.7580559722090576469559833339140197949745925143361e-06,
        0.87500487902798607109389195102266967296600341796875,
        9.7580559722090576469559833339140197949745925143361e-06};
    (*new_lh) = new_lh_value;
    SeqRegion seqregion1(TYPE_O, 412, -1, -1, std::move(new_lh));
    auto new_lh1 = cmaple::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value1{0.59998393993535814594508792652050033211708068847656,
        2.6766774402933635499746492514283602304203668609262e-05,
        2.6766774402933635499746492514283602304203668609262e-05,
        0.39996252651583585890904259940725751221179962158203};
    (*new_lh1) = new_lh_value1;
    SeqRegion seqregion2(TYPE_O, 412, -1, -1, std::move(new_lh1));
    RealNumType total_blength_1 = -1;
    RealNumType total_blength_2 = -1;
    const PositionType end_pos = 412;
    EXPECT_TRUE(merge_O_ORACGT_TwoLowers<4>(seqregion1, seqregion2,
                total_blength_1, total_blength_2, end_pos, tree.aln, tree.model,
                threshold_prob, log_lh, merged_regions, true));
    EXPECT_TRUE(merged_regions->size() == 1);
    EXPECT_FALSE(merged_regions->back().type != TYPE_O);
    EXPECT_TRUE(merged_regions->back().plength_observation2root == 0);
    // EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge{0.99963572952385693,
        3.4820599267114684E-9,0.00031223631362821728,
        0.000052030680454886103};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge);
    EXPECT_EQ(log_lh, -2.5901247759055703);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    log_lh = 0;
    total_blength_1 = -1;
    total_blength_2 = 32e-5;
    seqregion2.type = TYPE_R;
    EXPECT_TRUE(merge_O_ORACGT_TwoLowers<4>(seqregion1, seqregion2,
                total_blength_1, total_blength_2, end_pos, tree.aln,
                tree.model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_GE(merged_regions->size(), 2);
    EXPECT_FALSE(merged_regions->back().type != TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_GE(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge1{0.9993205977222771,
        2.4282817638871146E-9,0.00067939799763484592,
        1.851806273068163E-9};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge1);
    EXPECT_EQ(log_lh, -2.0790653485332990513256845588330179452896118164062);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    log_lh = 0;
    total_blength_1 = 143e-7;
    total_blength_2 = 0;
    SeqRegion::LHType new_lh_value2{2.7879382638950842936132173272012479969816922675818e-06,
        0.49999721206173608489820026079542003571987152099609,
        2.7879382638950842936132173272012479969816922675818e-06,
        0.49999721206173608489820026079542003571987152099609};
    (*seqregion1.likelihood) = new_lh_value2;
    seqregion2.type = 0;
    EXPECT_TRUE(merge_O_ORACGT_TwoLowers<4>(seqregion1, seqregion2,
                total_blength_1, total_blength_2, end_pos, tree.aln, tree.model,
                threshold_prob, log_lh, merged_regions, true));
    EXPECT_FALSE(merged_regions->size() != 3);
    EXPECT_FALSE(merged_regions->back().type != 0);
    EXPECT_EQ(log_lh, -12.484754332344007110577877028845250606536865234375);
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    log_lh = 0;
    total_blength_1 = 0;
    total_blength_2 = 12e-10;
    seqregion2.type = 3;
    EXPECT_TRUE(merge_O_ORACGT_TwoLowers<4>(seqregion1, seqregion2,
                total_blength_1, total_blength_2, end_pos, tree.aln, tree.model,
                threshold_prob, log_lh, merged_regions, true) == true);
    EXPECT_EQ(merged_regions->size(), 4);
    EXPECT_TRUE(merged_regions->back().type == 3);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_FALSE(log_lh != -0.69315275629224570863584631297271698713302612304688);
    // ----- Test 4 -----
    
    // ----- Test 5 -----
    log_lh = 0;
    total_blength_1 = 1e-5;
    total_blength_2 = 68e-4;
    SeqRegion::LHType new_lh_value3{8.3640013382402129011325767060647251582850003615022e-06,
        0.24998327199732350845096107150311581790447235107422,
        8.3640013382402129011325767060647251582850003615022e-06,
        0.74999999999999988897769753748434595763683319091797};
    (*seqregion1.likelihood) = new_lh_value3;
    seqregion2.type = TYPE_R;
    EXPECT_TRUE(merge_O_ORACGT_TwoLowers<4>(seqregion1, seqregion2,
                total_blength_1, total_blength_2, end_pos, tree.aln,
                tree.model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_LE(merged_regions->size(), 5);
    EXPECT_TRUE(merged_regions->back().type == TYPE_O);
    EXPECT_GE(merged_regions->back().plength_observation2root, 0);
    // EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge5{0.016447471276376950805,
        0.29913069287409010943,4.9917991548751474806e-05,
        0.68437191785798412447};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge5);
    EXPECT_EQ(log_lh, -7.5008336555086989);
    // ----- Test 5 -----
    
    /*// ----- Test 6 -----
    log_lh = 0;
    total_blength_1 = 0;
    total_blength_2 = 15e-2;
    seqregion2.type = TYPE_O;
    SeqRegion::LHType new_lh_value6{0.39998126426090857554740409796067979186773300170898,
     1.3382670779607485874503763900733588343427982181311e-05,
     0.59999197039753215943136410714942030608654022216797,
     1.3382670779607485874503763900733588343427982181311e-05};
    (*seqregion2.likelihood) = new_lh_value6;
    EXPECT_TRUE(merge_O_ORACGT_TwoLowers<4>(seqregion1, seqregion2,
     total_blength_1, total_blength_2, end_pos, tree.aln, tree.model,
     threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 6);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge6{0.000099915947241675934,
     0.10965211286069129,0.00013202204674923939,0.890115949145318};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge6);
    EXPECT_EQ(log_lh, -3.40271550630699692874259199015796184539794921875);
    // ----- Test 6 -----
    
    // ----- Test 7 -----
    log_lh = 0;
    total_blength_1 = 0;
    total_blength_2 = 0;
    SeqRegion::LHType new_lh_value4{7.4346402191731947664294494204639818235591519623995e-06,
     0.33332094226630137878686355179524980485439300537109,
     7.4346402191731947664294494204639818235591519623995e-06,
     0.66666418845326025355291221785591915249824523925781};
    (*seqregion1.likelihood) = new_lh_value4;
    seqregion2.type = 1;
    EXPECT_TRUE(merge_O_ORACGT_TwoLowers<4>(seqregion1, seqregion2,
     total_blength_1, total_blength_2, end_pos, tree.aln, tree.model,
     threshold_prob, log_lh, merged_regions, true));
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
    SeqRegion::LHType new_lh_value7{2.4874076016223502201974116939081453636628538106379e-10,
     7.4343638397288813108904244331132105116921593435109e-06,
     2.4874076016223502201974116939081453636628538106379e-10,
     0.99999256513867862405930964087019674479961395263672};
    (*seqregion2.likelihood) = new_lh_value7;
    EXPECT_TRUE(merge_O_ORACGT_TwoLowers<4>(seqregion1, seqregion2,
     total_blength_1, total_blength_2, end_pos, tree.aln, tree.model,
     threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 8);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge8{2.7739677142300221E-15,
     0.0000037170713771481666,2.7739677142300221E-15,
     0.9999962829286173};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge8);
    EXPECT_EQ(log_lh, -0.40547254324585230156330339923442807048559188842773);
    // ----- Test 8 -----
    
    // ----- Test 9 -----
    log_lh = 0;
    total_blength_1 = 437e-6;
    total_blength_2 = -1;
    SeqRegion::LHType new_lh_value5{0.49999442406129068761089229155913926661014556884766,
     5.5759387092817078802929955938516570768115343526006e-06,
     0.49999442406129068761089229155913926661014556884766,
     5.5759387092817078802929955938516570768115343526006e-06};
    (*seqregion1.likelihood) = new_lh_value5;
    seqregion2.type = TYPE_R;
    EXPECT_TRUE(merge_O_ORACGT_TwoLowers<4>(seqregion1, seqregion2,
     total_blength_1, total_blength_2, end_pos, tree.aln, tree.model,
     threshold_prob, log_lh, merged_regions, true));
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
    EXPECT_TRUE(merge_O_ORACGT_TwoLowers<4>(seqregion1, seqregion2,
     total_blength_1, total_blength_2, end_pos, tree.aln, tree.model,
     threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 10);
    EXPECT_EQ(merged_regions->back().type, 2);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -0.69315833249954661443581471758079715073108673095703);
    // ----- Test 10 -----*/
}

/*
 Test merge_RACGT_O_TwoLowers<4>(const SeqRegion& seq2_region,
 RealNumType total_blength_2, const PositionType end_pos, const Alignment& aln,
 const ModelBase* model, const RealNumType threshold_prob,
 SeqRegion::LHType& new_lh, RealNumType &log_lh,
 SeqRegions* merged_regions, const bool return_log_lh)
 */
TEST(SeqRegions, merge_RACGT_O_TwoLowers)
{
    Alignment aln = loadAln5K();
    Model model(cmaple::ModelBase::GTR);
    Tree tree(&aln, &model);
    std::unique_ptr<Params>& params = tree.params;
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    const RealNumType threshold_prob = params->threshold_prob;
    
    // ----- Test 1 -----
    std::unique_ptr<SeqRegions> merged_regions = cmaple::make_unique<SeqRegions>();
    RealNumType log_lh = 0;
    auto new_lh = cmaple::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value{0.12497560486006949187487435892762732692062854766846,
        9.7580559722090576469559833339140197949745925143361e-06,
        0.87500487902798607109389195102266967296600341796875,
        9.7580559722090576469559833339140197949745925143361e-06};
    (*new_lh) = new_lh_value;
    auto new_lh1 = cmaple::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value1{0.59998393993535814594508792652050033211708068847656,
        2.6766774402933635499746492514283602304203668609262e-05,
        2.6766774402933635499746492514283602304203668609262e-05,
        0.39996252651583585890904259940725751221179962158203};
    (*new_lh1) = new_lh_value1;
    SeqRegion seqregion2(TYPE_O, 412, -1, -1, std::move(new_lh1));
    RealNumType total_blength_2 = -1;
    const PositionType end_pos = 412;
    EXPECT_TRUE(merge_RACGT_O_TwoLowers<4>(seqregion2, total_blength_2,
                end_pos, tree.aln, tree.model, threshold_prob, *new_lh,
                                        log_lh, merged_regions, true));
    EXPECT_TRUE(merged_regions->size() == 1);
    EXPECT_TRUE(merged_regions->back().type == TYPE_O);
    EXPECT_TRUE(merged_regions->back().plength_observation2root == 0);
    // EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge{0.999635729523856930711644963594,
        3.48205992671146837419631824906e-09,
        0.000312236313628217279116106031012,
        5.20306804548861033901142880698e-05};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge);
    EXPECT_EQ(log_lh, -2.5901247759055703312469631782732903957366943359375);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    log_lh = 0;
    total_blength_2 = 1e-3;
    (*seqregion2.likelihood) = new_lh_value1;
    SeqRegion::LHType new_lh_value2{6.9955421312460765020272434505291649434188805400936e-11,
        6.9955421312460765020272434505291649434188805400936e-11,
        0.99999686351370775660996059741592034697532653808594,
        3.1363463814122085285697027340345854895531374495476e-06};
    (*new_lh) = new_lh_value2;
    EXPECT_TRUE(merge_RACGT_O_TwoLowers<4>(seqregion2, total_blength_2, end_pos,
                tree.aln, tree.model, threshold_prob, *new_lh, log_lh,
                merged_regions, true));
    EXPECT_GE(merged_regions->size(), 2);
    EXPECT_LE(merged_regions->size(), 2);
    EXPECT_FALSE(merged_regions->back().type != TYPE_O);
    // EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_GE(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge2{9.02597705818903812311447803357e-08,
        9.66903454131082698149291275817e-11,0.997304647644765784875175995694,
        0.00269526199877322038961358074971};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge2);
    EXPECT_EQ(log_lh, -7.6737265825075411385114421136677265167236328125);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    log_lh = 0;
    total_blength_2 = 0;
    SeqRegion::LHType new_lh_value3{4.9563776809356656518181297177427779843128519132733e-06,
        4.9563776809356656518181297177427779843128519132733e-06,
        0.11109844481259316395505010177657823078334331512451,
        0.88889164243204510373885796070680953562259674072266};
    (*seqregion2.likelihood) = new_lh_value3;
    (*new_lh) = new_lh_value;
    EXPECT_TRUE(merge_RACGT_O_TwoLowers<4>(seqregion2, total_blength_2,
                end_pos, tree.aln, tree.model, threshold_prob, *new_lh,
                log_lh, merged_regions, true));
    EXPECT_FALSE(merged_regions->size() != 3);
    EXPECT_TRUE(merged_regions->back().type == TYPE_O);
    /*EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);*/
    SeqRegion::LHType new_lh_value_merge3{6.37132369003884032904877487979e-06,
        4.9747095246365812296955045472e-10,0.999904410245811336999111063051,
        8.92179330276898546660951927478e-05};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge3);
    EXPECT_EQ(log_lh, -2.3307688028059012630421875655883923172950744628906);
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    log_lh = 0;
    total_blength_2 = 213e-7;
    (*seqregion2.likelihood) = new_lh_value1;
    SeqRegion::LHType new_lh_value4{2.6136903006998423677005767562508964374501374550164e-06,
        0.062492812351673081294745060176865081302821636199951,
        2.6136903006998423677005767562508964374501374550164e-06,
        0.93750196026772558699491355582722462713718414306641};
    (*new_lh) = new_lh_value4;
    EXPECT_TRUE(merge_RACGT_O_TwoLowers<4>(seqregion2, total_blength_2,
                end_pos, tree.aln, tree.model, threshold_prob,
                *new_lh, log_lh, merged_regions, true));
    EXPECT_LE(merged_regions->size(), 4);
    EXPECT_GE(merged_regions->size(), 4);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_TRUE(merged_regions->back().plength_observation2root == 0);
    EXPECT_FALSE(merged_regions->back().plength_observation2node != 0);
    SeqRegion::LHType new_lh_value_merge4{4.18220729036683253488498879236e-06,
        6.6470825462486196933439668022e-06,2.51442378530120592451234655152e-10,
        0.999989170458721043921457294346};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge4);
    EXPECT_EQ(log_lh, -0.98093450214178112833707245954428799450397491455078);
    // ----- Test 4 -----
    
    /*// ----- Test 5 -----
    log_lh = 0;
    total_blength_2 = -1;
    SeqRegion::LHType new_lh_value5{9.9498152920206932387357935649403392619483099679201e-10,
     8.9219993723549547142565030455330088443588465452194e-05,
     9.9498152920206932387357935649403392619483099679201e-10,
     0.99991077801631333965559633725206367671489715576172};
    (*seqregion2.likelihood) = new_lh_value5;
    (*new_lh) = new_lh_value;
    EXPECT_TRUE(merge_RACGT_O_TwoLowers<4>(seqregion2, total_blength_2,
     end_pos, tree.aln, tree.model, threshold_prob, *new_lh, log_lh,
     merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 5);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge5{1.27418556907770044148197133294e-05,
     8.92108976770603789461719368425e-05,8.92108976770603653936447807737e-05,
     0.999808836348955121131609757867};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge5);
    EXPECT_EQ(log_lh, -11.537315404593238454822312633041292428970336914062);
    // ----- Test 5 -----
    
    // ----- Test 6 -----
    log_lh = 0;
    total_blength_2 = 0;
    (*seqregion2.likelihood) = new_lh_value1;
    SeqRegion::LHType new_lh_value6{5.5759387092817078802929955938516570768115343526006e-06,
     0.49999442406129068761089229155913926661014556884766,
     5.5759387092817078802929955938516570768115343526006e-06,
     0.49999442406129068761089229155913926661014556884766};
    (*new_lh) = new_lh_value6;
    EXPECT_TRUE(merge_RACGT_O_TwoLowers<4>(seqregion2, total_blength_2,
     end_pos, tree.aln, tree.model, threshold_prob, *new_lh, log_lh,
     merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 6);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge6{1.67277228426044177343658819757e-05,
     6.6917607757947737928509723826e-05,7.46265281119078188993421229841e-10,
     0.999916353923134160197605524445};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge6);
    EXPECT_EQ(log_lh, -1.6094591028973106450195018624071963131427764892578);
    // ----- Test 6 -----
    
    // ----- Test 7 -----
    log_lh = 0;
    total_blength_2 = 3163e-7;
    SeqRegion::LHType new_lh_value7{1.2436575683940661210420103588707597258578019250308e-10,
     1.1151877415789591228433182135137968771232408471406e-05,
     1.2436575683940661210420103588707597258578019250308e-10,
     0.99998884787385255989988763758447021245956420898438};
    (*seqregion2.likelihood) = new_lh_value7;
    (*new_lh) = new_lh_value;
    EXPECT_TRUE(merge_RACGT_O_TwoLowers<4>(seqregion2, total_blength_2,
     end_pos, tree.aln, tree.model, threshold_prob, *new_lh, log_lh,
     merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 7);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge7{0.0166262913969323911089759349125,
     2.33062694270459372499967876102e-05,0.93180971485701069578766464474,
     0.0515406874766299802348434866417};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge7);
    EXPECT_EQ(log_lh, -8.5724436071309479956426002900116145610809326171875);
    // ----- Test 7 -----
    
    // ----- Test 8 -----
    log_lh = 0;
    total_blength_2 = 163e-10;
    (*seqregion2.likelihood) = new_lh_value1;
    SeqRegion::LHType new_lh_value8{0.74997769674260927885711680573876947164535522460938,
     0.24998327274350143345493791002809302881360054016113,
     1.9515256944661845116837164959555650511902058497071e-05,
     1.9515256944661845116837164959555650511902058497071e-05};
    (*new_lh) = new_lh_value8;
    EXPECT_TRUE(merge_RACGT_O_TwoLowers<4>(seqregion2, total_blength_2,
     end_pos, tree.aln, tree.model, threshold_prob, *new_lh, log_lh,
     merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 8);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge8{0.999967777775584876209791218571,
     1.48753723892170687314673652168e-05,1.16113808009755283430647568722e-09,
     1.73456908878727982935719770241e-05};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge8);
    EXPECT_EQ(log_lh, -0.79853198337463082712162076859385706484317779541016);
    // ----- Test 8 -----
    
    // ----- Test 9 -----
    log_lh = 0;
    total_blength_2 = 0;
    SeqRegion::LHType new_lh_value9{0.99997490659690790870683940738672390580177307128906,
     2.5092097243268992381179730011275808010395849123597e-05,
     6.5292438309787796439726148030194621818544931102224e-10,
     6.5292438309787796439726148030194621818544931102224e-10};
    (*seqregion2.likelihood) = new_lh_value9;
    (*new_lh) = new_lh_value;
    EXPECT_TRUE(merge_RACGT_O_TwoLowers<4>(seqregion2, total_blength_2,
     end_pos, tree.aln, tree.model, threshold_prob, *new_lh,
     log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 9);
    EXPECT_EQ(merged_regions->back().type, TYPE_R);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -2.0796618090328400363375749293481931090354919433594);
    // ----- Test 9 -----*/
}

/*
 Test merge_RACGT_RACGT_TwoLowers<4>(const SeqRegion& seq2_region,
 RealNumType total_blength_2, const PositionType end_pos,
 const Alignment& aln, const ModelBase* model, const RealNumType threshold_prob,
 SeqRegion::LHType& new_lh, RealNumType& sum_lh, RealNumType &log_lh,
 SeqRegions* merged_regions, const bool return_log_lh)
 */
TEST(SeqRegions, merge_RACGT_RACGT_TwoLowers)
{
    Alignment aln = loadAln5K();
    Model model(cmaple::ModelBase::GTR);
    Tree tree(&aln, &model);
    std::unique_ptr<Params>& params = tree.params;
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    const RealNumType threshold_prob = params->threshold_prob;
    
    // ----- Test 1 -----
    std::unique_ptr<SeqRegions> merged_regions = cmaple::make_unique<SeqRegions>();
    RealNumType log_lh = 0;
    RealNumType sum_lh = 0;
    auto new_lh = cmaple::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value{9.7580559722090576469559833339140197949745925143361e-06,
        9.7580559722090576469559833339140197949745925143361e-06,
        0.87500487902798607109389195102266967296600341796875,
        0.12497560486006949187487435892762732692062854766846};
    (*new_lh) = new_lh_value;
    SeqRegion seqregion2(TYPE_R, 412, -1, -1);
    RealNumType total_blength_2 = -1;
    const PositionType end_pos = 412;
    EXPECT_TRUE(merge_RACGT_RACGT_TwoLowers<4>(seqregion2, total_blength_2,
                    end_pos, tree.aln, tree.model, threshold_prob, *new_lh,
                    sum_lh, log_lh, merged_regions, true));
    EXPECT_TRUE(merged_regions->size() == 1);
    EXPECT_EQ(merged_regions->back().type, seqregion2.type);
    EXPECT_FALSE(merged_regions->back().plength_observation2node != -1);
    EXPECT_EQ(log_lh, -11.537417360554178102916011994238942861557006835938);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    log_lh = 0;
    sum_lh = 0;
    total_blength_2 = 3.1363463814122085285697027340345854895531374495476e-06;
    seqregion2.type = 0;
    SeqRegion::LHType new_lh_value2{6.9955421312460765020272434505291649434188805400936e-11,
        3.1363463814122085285697027340345854895531374495476e-06,
        0.99999686351370775660996059741592034697532653808594,
        6.9955421312460765020272434505291649434188805400936e-11};
    (*new_lh) = new_lh_value2;
    EXPECT_TRUE(merge_RACGT_RACGT_TwoLowers<4>(seqregion2, total_blength_2,
                    end_pos, tree.aln, tree.model, threshold_prob, *new_lh,
                    sum_lh, log_lh, merged_regions, true));
    EXPECT_FALSE(merged_regions->size() < 2);
    EXPECT_FALSE(merged_regions->back().type != TYPE_O);
    EXPECT_TRUE(merged_regions->back().plength_observation2root == 0);
    EXPECT_TRUE(merged_regions->back().plength_observation2node == 0);
    SeqRegion::LHType new_lh_value_merge2{7.35069993527212464743542108536e-05,
        1.00511327029477612984332841189e-06,
        0.999925487870280349511631357018,
        1.70965639624293546874861604245e-11};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge2);
    EXPECT_EQ(log_lh, -13.865034049156509610156717826612293720245361328125);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    log_lh = 0;
    sum_lh = 0;
    total_blength_2 = 0;
    seqregion2.type = TYPE_R;
    (*new_lh) = new_lh_value;
    EXPECT_TRUE(merge_RACGT_RACGT_TwoLowers<4>(seqregion2, total_blength_2,
                    end_pos, tree.aln, tree.model, threshold_prob, *new_lh,
                    sum_lh, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 3);
    EXPECT_FALSE(merged_regions->back().type != seqregion2.type);
    EXPECT_TRUE(merged_regions->back().plength_observation2node == -1);
    EXPECT_TRUE(log_lh == -11.537417360554178102916011994238942861557006835938);
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    log_lh = 0;
    sum_lh = 0;
    total_blength_2 = 0.062492812351673081294745060176865081302821636199951;
    seqregion2.type = 1;
    SeqRegion::LHType new_lh_value4{0.062492812351673081294745060176865081302821636199951,
        0.93750196026772558699491355582722462713718414306641,
        2.6136903006998423677005767562508964374501374550164e-06,
        2.6136903006998423677005767562508964374501374550164e-06};
    (*new_lh) = new_lh_value4;
    EXPECT_TRUE(merge_RACGT_RACGT_TwoLowers<4>(seqregion2, total_blength_2,
                end_pos, tree.aln, tree.model, threshold_prob, *new_lh,
                sum_lh, log_lh, merged_regions, true));
    EXPECT_FALSE(merged_regions->size() > 4);
    EXPECT_FALSE(merged_regions->back().type != TYPE_O);
    /*EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);*/
    SeqRegion::LHType new_lh_value_merge4{0.000276015672152719990402325311862,
        0.9997238125720525614426037464,1.76015102091620514116593592541e-08,
        1.54154284403481275470771598261e-07};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge4);
    EXPECT_EQ(log_lh, -0.16879625050215649184615074318571714684367179870605);
    // ----- Test 4 -----
    
    // ----- Test 5 -----
    log_lh = 0;
    sum_lh = 0;
    total_blength_2 = -1;
    seqregion2.type = 0;
    (*new_lh) = new_lh_value;
    EXPECT_TRUE(merge_RACGT_RACGT_TwoLowers<4>(seqregion2, total_blength_2,
                    end_pos, tree.aln, tree.model, threshold_prob,
                    *new_lh, sum_lh, log_lh, merged_regions, true));
    EXPECT_FALSE(merged_regions->size() != 5);
    EXPECT_EQ(merged_regions->back().type, seqregion2.type);
    EXPECT_TRUE(merged_regions->back().plength_observation2root == -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_FALSE(log_lh != -11.537417360554178102916011994238942861557006835938);
    // ----- Test 5 -----
    
    /*// ----- Test 6 -----
    log_lh = 0;
    sum_lh = 0;
    total_blength_2 = 0.49999442406129068761089229155913926661014556884766;
    seqregion2.type = 3;
    SeqRegion::LHType new_lh_value6{0.49999442406129068761089229155913926661014556884766,
     0.49999442406129068761089229155913926661014556884766,
     5.5759387092817078802929955938516570768115343526006e-06,
     5.5759387092817078802929955938516570768115343526006e-06};
    (*new_lh) = new_lh_value6;
    EXPECT_TRUE(merge_RACGT_RACGT_TwoLowers<4>(seqregion2, total_blength_2,
     end_pos, tree.aln, tree.model, threshold_prob, *new_lh, sum_lh,
     log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 6);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge6{0.0540225020296241503769962832848,
     0.945967079514273279094993540639,4.82257467713597815899819951091e-06,
     5.59588142544069931037494305959e-06};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge6);
    EXPECT_EQ(log_lh, -0.9987216763314137324414332397282123565673828125);
    // ----- Test 6 -----
    
    // ----- Test 7 -----
    log_lh = 0;
    sum_lh = 0;
    total_blength_2 = 1.2436575683940661210420103588707597258578019250308e-10;
    seqregion2.type = 0;
    (*new_lh) = new_lh_value;
    EXPECT_TRUE(merge_RACGT_RACGT_TwoLowers<4>(seqregion2, total_blength_2,
     end_pos, tree.aln, tree.model, threshold_prob, *new_lh, sum_lh,
     log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 7);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge7{0.99999649823344849419726187989,
     1.20936802836750174748186692148e-11,3.38363623651822093845354007258e-06,
     1.18118221219596752754965953997e-07};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge7);
    EXPECT_EQ(log_lh, -11.537413858823567736067161604296416044235229492188);
    // ----- Test 7 -----
    
    // ----- Test 8 -----
    log_lh = 0;
    sum_lh = 0;
    total_blength_2 = 1.9515256944661845116837164959555650511902058497071e-05;
    seqregion2.type = TYPE_R;
    SeqRegion::LHType new_lh_value8{0.74997769674260927885711680573876947164535522460938,
     1.9515256944661845116837164959555650511902058497071e-05,
     1.9515256944661845116837164959555650511902058497071e-05,
     0.24998327274350143345493791002809302881360054016113};
    (*new_lh) = new_lh_value8;
    EXPECT_TRUE(merge_RACGT_RACGT_TwoLowers<4>(seqregion2, total_blength_2,
     end_pos, tree.aln, tree.model, threshold_prob, *new_lh, sum_lh,
     log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 8);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge8{0.999999517409988936122999803047,
     4.93812459495446616345023798386e-11,1.54077796167036529766850756536e-10,
     4.82386551824058610632602311225e-07};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge8);
    EXPECT_EQ(log_lh, -0.28771792989182387589863765242625959217548370361328);
    // ----- Test 8 -----
    
    // ----- Test 9 -----
    log_lh = 0;
    sum_lh = 0;
    total_blength_2 = 0;
    seqregion2.type = 1;
    (*new_lh) = new_lh_value;
    EXPECT_TRUE(merge_RACGT_RACGT_TwoLowers<4>(seqregion2, total_blength_2,
     end_pos, tree.aln, tree.model, threshold_prob, *new_lh, sum_lh,
     log_lh, merged_regions, true));
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
    EXPECT_TRUE(merge_RACGT_RACGT_TwoLowers<4>(seqregion2, total_blength_2,
     end_pos, tree.aln, tree.model, threshold_prob, *new_lh,
     sum_lh, log_lh, merged_regions, true));
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
    SeqRegion::LHType new_lh_value11{2.5092097243268992381179730011275808010395849123597e-05,
     0.99997490659690790870683940738672390580177307128906,
     6.5292438309787796439726148030194621818544931102224e-10,
     6.5292438309787796439726148030194621818544931102224e-10};
    (*new_lh) = new_lh_value11;
    EXPECT_TRUE(merge_RACGT_RACGT_TwoLowers<4>(seqregion2, total_blength_2,
     end_pos, tree.aln, tree.model, threshold_prob, *new_lh,
     sum_lh, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 11);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge11{0.999997469693585716576933464239,
     2.53030640795433316882884731969e-06,
     5.15495543614262290103202928001e-15,
     1.25992117026322114735542991217e-15};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge11);
    EXPECT_EQ(log_lh, -10.592955082179518200291568064130842685699462890625);
    // ----- Test 11 -----*/
}


/*
 Test merge_RACGT_ORACGT_TwoLowers<4>(const SeqRegion& seq1_region,
 const SeqRegion& seq2_region, RealNumType total_blength_1,
 RealNumType total_blength_2, const PositionType end_pos,
 const Alignment& aln, const ModelBase* model, const RealNumType threshold_prob,
 RealNumType &log_lh, SeqRegions* merged_regions, const bool return_log_lh)
 */
TEST(SeqRegions, merge_RACGT_ORACGT_TwoLowers)
{
    Alignment aln = loadAln5K();
    Model model(cmaple::ModelBase::GTR);
    Tree tree(&aln, &model);
    std::unique_ptr<Params>& params = tree.params;
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    const RealNumType threshold_prob = params->threshold_prob;
    
    // ----- Test 1 -----
    std::unique_ptr<SeqRegions> merged_regions = cmaple::make_unique<SeqRegions>();
    RealNumType log_lh = 0;
    SeqRegion seqregion1(TYPE_R, 3223);
    auto new_lh1 = cmaple::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value1{2.6766774402933635499746492514283602304203668609262e-05,
        0.39996252651583585890904259940725751221179962158203,
        2.6766774402933635499746492514283602304203668609262e-05,
        0.59998393993535814594508792652050033211708068847656};
    (*new_lh1) = new_lh_value1;
    SeqRegion seqregion2(TYPE_O, 3223, -1, -1, std::move(new_lh1));
    RealNumType total_blength_1 = 1326e-5;
    RealNumType total_blength_2 = -1;
    const PositionType end_pos = 3223;
    EXPECT_TRUE(merge_RACGT_ORACGT_TwoLowers<4>(seqregion1, seqregion2,
                total_blength_1, total_blength_2, end_pos, tree.aln,
                tree.model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 1);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    EXPECT_EQ(log_lh, -6.7833577173188217557253665290772914886474609375);
    SeqRegion::LHType new_lh_value_merge1{0.0235298053213485250378944613203,
        0.455403998389659003809271098362,
        9.50936657505749543661116574e-05,
        0.520971102623241977269685776264};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge1);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    log_lh = 0;
    total_blength_1 = -1;
    total_blength_2 = 32e-5;
    seqregion1.type = 2;
    seqregion2.type = TYPE_R;
    EXPECT_TRUE(merge_RACGT_ORACGT_TwoLowers<4>(seqregion1, seqregion2,
                total_blength_1, total_blength_2, end_pos, tree.aln,
                tree.model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_TRUE(merged_regions->size() == 2);
    EXPECT_TRUE(merged_regions->back().type == 2);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_FALSE(log_lh != -9.23984298142250537466679816134274005889892578125);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    log_lh = 0;
    total_blength_1 = 143e-7;
    total_blength_2 = 0;
    seqregion1.type = TYPE_R;
    seqregion2.type = 1;
    EXPECT_TRUE(merge_RACGT_ORACGT_TwoLowers<4>(seqregion1, seqregion2,
                total_blength_1, total_blength_2, end_pos, tree.aln,
                tree.model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_TRUE(merged_regions->size() == 3);
    EXPECT_FALSE(merged_regions->back().type != 1);
    EXPECT_TRUE(merged_regions->back().plength_observation2node == -1);
    EXPECT_EQ(log_lh, -13.485791369259754191034517134539783000946044921875);
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    log_lh = 0;
    total_blength_1 = 0;
    total_blength_2 = 12e-10;
    seqregion1.type = 0;
    seqregion2.type = 3;
    EXPECT_TRUE(merge_RACGT_ORACGT_TwoLowers<4>(seqregion1, seqregion2,
                total_blength_1, total_blength_2, end_pos, tree.aln,
                tree.model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 4);
    EXPECT_TRUE(merged_regions->back().type == TYPE_R);
    EXPECT_FALSE(merged_regions->back().plength_observation2root != -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_TRUE(log_lh == -23.07170390644744628616535919718444347381591796875);
    // ----- Test 4 -----
    
    /*// ----- Test 5 -----
    log_lh = 0;
    total_blength_1 = 1e-5;
    total_blength_2 = 68e-4;
    seqregion1.type = 1;
    seqregion2.type = TYPE_R;
    EXPECT_TRUE(merge_RACGT_ORACGT_TwoLowers<4>(seqregion1, seqregion2,
     total_blength_1, total_blength_2, end_pos, tree.aln, tree.model,
     threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 5);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge5{0.000899932704475660738205333721851,
     0.999091155743592085336501895654,2.83758240671473494537067026877e-06,
     6.0739695257571569122004295771e-06};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge5);
    EXPECT_EQ(log_lh, -7.3204796410357877434194051602389663457870483398438);
    // ----- Test 5 -----
    
    // ----- Test 6 -----
    log_lh = 0;
    total_blength_1 = 0;
    total_blength_2 = 15e-2;
    seqregion1.type = 2;
    seqregion2.type = TYPE_O;
    SeqRegion::LHType new_lh_value6{1.3382670779607485874503763900733588343427982181311e-05,
     1.3382670779607485874503763900733588343427982181311e-05,
     0.39998126426090857554740409796067979186773300170898,
     0.59999197039753215943136410714942030608654022216797};
    (*seqregion2.likelihood) = new_lh_value6;
    EXPECT_TRUE(merge_RACGT_ORACGT_TwoLowers<4>(seqregion1, seqregion2,
     total_blength_1, total_blength_2, end_pos, tree.aln, tree.model,
     threshold_prob, log_lh, merged_regions, true));
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
    EXPECT_TRUE(merge_RACGT_ORACGT_TwoLowers<4>(seqregion1, seqregion2,
     total_blength_1, total_blength_2, end_pos, tree.aln, tree.model,
     threshold_prob, log_lh, merged_regions, true));
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
    SeqRegion::LHType new_lh_value7{0.99999256513867862405930964087019674479961395263672,
     7.4343638397288813108904244331132105116921593435109e-06,
     2.4874076016223502201974116939081453636628538106379e-10,
     2.4874076016223502201974116939081453636628538106379e-10};
    (*seqregion2.likelihood) = new_lh_value7;
    EXPECT_TRUE(merge_RACGT_ORACGT_TwoLowers<4>(seqregion1, seqregion2,
     total_blength_1, total_blength_2, end_pos, tree.aln, tree.model,
     threshold_prob, log_lh, merged_regions, true));
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
    EXPECT_TRUE(merge_RACGT_ORACGT_TwoLowers<4>(seqregion1, seqregion2,
     total_blength_1, total_blength_2, end_pos, tree.aln, tree.model,
     threshold_prob, log_lh, merged_regions, true));
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
    EXPECT_TRUE(merge_RACGT_ORACGT_TwoLowers<4>(seqregion1, seqregion2,
     total_blength_1, total_blength_2, end_pos, tree.aln, tree.model, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 10);
    EXPECT_EQ(merged_regions->back().type, 2);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -15.129806325448386772336561989504843950271606445312);
    // ----- Test 10 -----*/
}

/*
 Test merge_notN_notN_TwoLowers<4>(const SeqRegion& seq1_region,
 const SeqRegion& seq2_region, const RealNumType plength1,
 const RealNumType plength2, const PositionType end_pos,
 const PositionType pos, const Alignment& aln, const ModelBase* model,
 const RealNumType threshold_prob, RealNumType &log_lh,
 SeqRegions* merged_regions, const bool return_log_lh)
 */
TEST(SeqRegions, merge_notN_notN_TwoLowers)
{
    Alignment aln = loadAln5K();
    Model model(cmaple::ModelBase::GTR);
    Tree tree(&aln, &model);
    std::unique_ptr<Params>& params = tree.params;
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    const RealNumType threshold_prob = params->threshold_prob;
    
    // ----- Test 1 -----
    std::unique_ptr<SeqRegions> merged_regions = cmaple::make_unique<SeqRegions>();
    RealNumType log_lh = 0;
    auto new_lh = cmaple::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value{0.42855939517854246822992081433767452836036682128906,
        0.57141073458160607234646022334345616400241851806641,
        1.4935119925781397759922616841343767646321794018149e-05,
        1.4935119925781397759922616841343767646321794018149e-05};
    (*new_lh) = new_lh_value;
    SeqRegion seqregion1(TYPE_O, 1412, -1, -1, std::move(new_lh));
    auto new_lh1 = cmaple::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value1{0.61110448521316429459915298139094375073909759521484,
        1.1926616304335032169438128579752600444408017210662e-05,
        1.1926616304335032169438128579752600444408017210662e-05,
        0.38887166155422708824218602785549592226743698120117};
    (*new_lh1) = new_lh_value1;
    SeqRegion seqregion2(TYPE_O, 1412, -1, -1, std::move(new_lh1));
    RealNumType plength1 = -1;
    RealNumType plength2 = -1;
    const PositionType end_pos = 1412;
    const PositionType pos = 1355;
    EXPECT_TRUE(merge_notN_notN_TwoLowers<4>(seqregion1, seqregion2,
                plength1, plength2, end_pos, pos, tree.aln, tree.model,
                tree.cumulative_rate, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 1);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    EXPECT_EQ(log_lh, -1.3397650685348243548844493489013984799385070800781);
    SeqRegion::LHType new_lh_value_merge1{0.999951803463153265916218970233,
        2.60206546527806729897022708364e-05,
        6.80109025377646455925279108979e-10,
        2.21752020848110232805524416611e-05};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge1);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    log_lh = 0;
    plength1 = -1;
    plength2 = 1.3382670779607485874503763900733588343427982181311e-05;
    seqregion1.type = TYPE_R;
    seqregion2.type = TYPE_R;
    EXPECT_TRUE(merge_notN_notN_TwoLowers<4>(seqregion1, seqregion2,
                plength1, plength2, end_pos, pos, tree.aln, tree.model,
                tree.cumulative_rate, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 2);
    EXPECT_TRUE(merged_regions->back().type == TYPE_R);
    EXPECT_TRUE(merged_regions->back().plength_observation2root == -1);
    EXPECT_TRUE(log_lh == -0.0007633997783012927721910112488501454208744689822197);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    log_lh = 0;
    plength1 = 3.4820599267114684E-9;
    plength2 = 3.283e-9;
    seqregion1.type = 2;
    seqregion2.type = 0;
    EXPECT_FALSE(!merge_notN_notN_TwoLowers<4>(seqregion1, seqregion2,
                    plength1, plength2, end_pos, pos, tree.aln, tree.model,
                    tree.cumulative_rate, threshold_prob, log_lh,
                    merged_regions, true));
    EXPECT_FALSE(merged_regions->size() != 3);
    EXPECT_FALSE(merged_regions->back().type != TYPE_O);
    EXPECT_FALSE(merged_regions->back().plength_observation2node != 0);
    EXPECT_TRUE(log_lh == -20.199112063767703517669360735453665256500244140625);
    SeqRegion::LHType new_lh_value_merge3{0.410245869594259182644435668408,
        6.40012426019357435068582000991e-11,
        0.589754130146332378181739386491,
        1.95407241773354552292205525266e-10};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge3);
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    log_lh = 0;
    plength1 = 0;
    plength2 = 0.00031223631362821728;
    seqregion1.type = TYPE_R;
    seqregion2.type = 3;
    EXPECT_TRUE(merge_notN_notN_TwoLowers<4>(seqregion1, seqregion2,
                    plength1, plength2, end_pos, pos, tree.aln, tree.model,
                    tree.cumulative_rate, threshold_prob, log_lh,
                    merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 4);
    EXPECT_FALSE(merged_regions->back().type != TYPE_R);
    EXPECT_FALSE(log_lh != -8.5224663158967093323781227809377014636993408203125);
    // ----- Test 4 -----
    
    // ----- Test 5 -----
    log_lh = 0;
    plength1 = 1e-5;
    plength2 = 0.000052030680454886103;
    seqregion1.type = 3;
    SeqRegion::LHType new_lh_value3{4.9383548301647924830444502664050787643645890057087e-05,
        0.71422131017848822231997019116533920168876647949219,
        4.9383548301647924830444502664050787643645890057087e-05,
        0.28567992272490849714472460618708282709121704101562};
    (*seqregion2.likelihood) = new_lh_value3;
    seqregion2.type = TYPE_O;
    EXPECT_TRUE(merge_notN_notN_TwoLowers<4>(seqregion1, seqregion2,
                    plength1, plength2, end_pos, pos, tree.aln, tree.model,
                    tree.cumulative_rate, threshold_prob, log_lh,
                    merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 5);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    EXPECT_EQ(log_lh, -1.2528228994081356262313420302234590053558349609375);
    SeqRegion::LHType new_lh_value_merge5{1.47064689584769381918621458095e-10,
        3.48425531256729660231630241185e-05,
        1.38799271102831098152142401229e-09,
        0.99996515591181689419641998029};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge5);
    // ----- Test 5 -----
    
    // ----- Test 6 -----
    log_lh = 0;
    plength1 = 0;
    plength2 = 8.3640013382402129011325767060647251582850003615022e-06;
    seqregion1.type = 1;
    seqregion2.type = TYPE_O;
    SeqRegion::LHType new_lh_value6{6.6911562987199814600764238847752096717158565297723e-06,
        0.39999063238118182095348629445652477443218231201172,
        6.6911562987199814600764238847752096717158565297723e-06,
        0.59999598530622066938633452082285657525062561035156};
    (*seqregion2.likelihood) = new_lh_value6;
    EXPECT_EQ(merge_notN_notN_TwoLowers<4>(seqregion1, seqregion2,
                    plength1, plength2, end_pos, pos, tree.aln,
                    tree.model, tree.cumulative_rate, threshold_prob,
                    log_lh, merged_regions, true), true);
    EXPECT_TRUE(merged_regions->size() == 6);
    EXPECT_TRUE(merged_regions->back().type == 1);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_FALSE(merged_regions->back().plength_observation2node != -1);
    EXPECT_TRUE(log_lh == -0.91630994861674031071174795215483754873275756835938);
    // ----- Test 6 -----
    
    // ----- Test 7 -----
    log_lh = 0;
    plength1 = 0;
    plength2 = 0;
    SeqRegion::LHType new_lh_value4{4.430579981880772584663374230178360321796837695274e-10,
        0.99999581732016740165391865957644768059253692626953,
        4.430579981880772584663374230178360321796837695274e-10,
        4.181793716486326507387246559366289488934853579849e-06};
    (*seqregion1.likelihood) = new_lh_value4;
    seqregion1.type = TYPE_O;
    seqregion2.type = 1;
    EXPECT_FALSE(!merge_notN_notN_TwoLowers<4>(seqregion1, seqregion2,
                    plength1, plength2, end_pos, pos, tree.aln, tree.model,
                    tree.cumulative_rate, threshold_prob, log_lh,
                    merged_regions, true));
    EXPECT_FALSE(merged_regions->size() != 7);
    EXPECT_FALSE(merged_regions->back().type != 1);
    EXPECT_EQ(log_lh, -4.1826885800280286031383250588966404848179081454873e-06);
    // ----- Test 7 -----
    
    // ----- Test 8 -----
    log_lh = 0;
    plength1 = 0;
    plength2 = -1;
    seqregion1.type = TYPE_R;
    seqregion2.type = TYPE_O;
    SeqRegion::LHType new_lh_value7{1.1815113248226830777267349123884135342343881802663e-09,
        0.9999702585027225865133004845120012760162353515625,
        1.1815113248226830777267349123884135342343881802663e-09,
        2.9739134254674524059092188821296076639555394649506e-05};
    (*seqregion2.likelihood) = new_lh_value7;
    EXPECT_EQ(!merge_notN_notN_TwoLowers<4>(seqregion1, seqregion2,
                plength1, plength2, end_pos, pos, tree.aln, tree.model,
                tree.cumulative_rate, threshold_prob, log_lh,
                merged_regions, true), false);
    EXPECT_TRUE(merged_regions->size() == 8);
    EXPECT_EQ(merged_regions->back().type, TYPE_R);
    EXPECT_FALSE(merged_regions->back().plength_observation2root != -1);
    EXPECT_FALSE(merged_regions->back().plength_observation2node != -1);
    EXPECT_TRUE(log_lh == -20.556471434224643957122680149041116237640380859375);
    // ----- Test 8 -----
    
    /*// ----- Test 9 -----
    log_lh = 0;
    plength1 = 3.4820599267114684E-9;
    plength2 = -1;
    SeqRegion::LHType new_lh_value5{6.6911562987199814600764238847752096717158565297723e-06,
     0.39999063238118182095348629445652477443218231201172,
     6.6911562987199814600764238847752096717158565297723e-06,
     0.59999598530622066938633452082285657525062561035156};
    (*seqregion1.likelihood) = new_lh_value5;
    seqregion1.type = TYPE_O;
    seqregion2.type = TYPE_R;
    EXPECT_TRUE(merge_notN_notN_TwoLowers<4>(seqregion1, seqregion2,
     plength1, plength2, end_pos, pos, tree.aln, tree.model,
     tree.cumulative_rate, threshold_prob, log_lh, merged_regions, true));
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
    EXPECT_TRUE(merge_notN_notN_TwoLowers<4>(seqregion1, seqregion2,
     plength1, plength2, end_pos, pos, tree.aln, tree.model,
     tree.cumulative_rate, threshold_prob, log_lh, merged_regions, true));
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
    SeqRegion::LHType new_lh_value11{4.181793716486327354420193813666628557257354259491e-06,
     0.99999581732016751267622112209210172295570373535156,
     4.430579981880772584663374230178360321796837695274e-10,
     4.430579981880772584663374230178360321796837695274e-10};
    (*seqregion1.likelihood) = new_lh_value11;
    seqregion1.type = TYPE_O;
    seqregion2.type = TYPE_R;
    EXPECT_TRUE(merge_notN_notN_TwoLowers<4>(seqregion1, seqregion2,
     plength1, plength2, end_pos, pos, tree.aln, tree.model,
     tree.cumulative_rate, threshold_prob, log_lh, merged_regions, true));
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
    EXPECT_TRUE(merge_notN_notN_TwoLowers<4>(seqregion1, seqregion2,
     plength1, plength2, end_pos, pos, tree.aln, tree.model,
     tree.cumulative_rate, threshold_prob, log_lh, merged_regions, true));
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
    SeqRegion::LHType new_lh_value13{0.33332094226630137878686355179524980485439300537109,
     7.4346402191731947664294494204639818235591519623995e-06,
     0.66666418845326025355291221785591915249824523925781,
     7.4346402191731947664294494204639818235591519623995e-06};
    (*seqregion1.likelihood) = new_lh_value13;
    seqregion1.type = TYPE_O;
    seqregion2.type = 0;
    EXPECT_TRUE(merge_notN_notN_TwoLowers<4>(seqregion1,
     seqregion2, plength1, plength2, end_pos, pos, tree.aln,
     tree.model, tree.cumulative_rate, threshold_prob,
     log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 12);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    EXPECT_EQ(log_lh, -1.0986478599204548);
    SeqRegion::LHType new_lh_value_merge13{0.99999830814790946487,
     7.1779832367552409141e-12,1.6918377797180232427e-06,
     7.1327930434698958472e-12};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge13);
    // ----- Test 13 -----
    
    // ----- Test 14 -----
    log_lh = 0;
    plength1 = 0;
    plength2 = 0.31223631362821728e-10;
    seqregion1.type = 1;
    seqregion2.type = 3;
    EXPECT_TRUE(merge_notN_notN_TwoLowers<4>(seqregion1, seqregion2,
     plength1, plength2, end_pos, pos, tree.aln, tree.model,
     tree.cumulative_rate, threshold_prob, log_lh, merged_regions, true));
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
    SeqRegion::LHType new_lh_value15{9.9493714778138372440088605107603655919312757305306e-10,
     0.99996654175612176285170562550774775445461273193359,
     9.9493714778138372440088605107603655919312757305306e-10,
     3.345625400403152548828300538730218249838799238205e-05};
    (*seqregion2.likelihood) = new_lh_value15;
    seqregion2.type = TYPE_O;
    EXPECT_TRUE(merge_notN_notN_TwoLowers<4>(seqregion1, seqregion2,
     plength1, plength2, end_pos, pos, tree.aln, tree.model,
     tree.cumulative_rate, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 14);
    EXPECT_EQ(merged_regions->back().type, TYPE_O);
    EXPECT_EQ(merged_regions->back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions->back().plength_observation2node, 0);
    EXPECT_EQ(log_lh, -11.502067070054362574182960088364779949188232421875);
    SeqRegion::LHType new_lh_value_merge15{2.13259105287923615202651250744e-06,
     0.804503715022567122971963726741,0.195330717008421000935314282287,
     0.000163435377959033966777449564667};
    EXPECT_EQ(*merged_regions->back().likelihood, new_lh_value_merge15);
    // ----- Test 15 -----
    
    // ----- Test 16 -----
    log_lh = 0;
    plength1 = 0;
    plength2 = 1.09085755475e-2;
    seqregion1.type = 2;
    seqregion2.type = TYPE_O;
    SeqRegion::LHType new_lh_value16{1.5333235876614054963434918832376752106938511133194e-05,
     0.99996236262065618660699328756891191005706787109375,
     1.1152071733605872690944446623539931806590175256133e-05,
     1.1152071733605872690944446623539931806590175256133e-05};
    (*seqregion2.likelihood) = new_lh_value16;
    EXPECT_TRUE(merge_notN_notN_TwoLowers<4>(seqregion1, seqregion2,
     plength1, plength2, end_pos, pos, tree.aln, tree.model,
     tree.cumulative_rate, threshold_prob, log_lh, merged_regions, true));
    EXPECT_EQ(merged_regions->size(), 15);
    EXPECT_EQ(merged_regions->back().type, TYPE_R);
    EXPECT_EQ(merged_regions->back().plength_observation2root, -1);
    EXPECT_EQ(merged_regions->back().plength_observation2node, -1);
    EXPECT_EQ(log_lh, -6.903698068549179112096680910326540470123291015625);
    // ----- Test 16 -----*/
}

/*
 Test mergeTwoLowers<4>(SeqRegions* &merged_regions,
 RealNumType plength1, const SeqRegions* const regions2,
 RealNumType plength2, const Alignment& aln, const ModelBase* model,
 RealNumType threshold_prob, bool return_log_lh)
 */
TEST(SeqRegions, mergeTwoLowers)
{
    Alignment aln = loadAln5K();
    Model model(cmaple::ModelBase::GTR);
    Tree tree(&aln, &model);
    std::unique_ptr<Params>& params = tree.params;
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    std::unique_ptr<SeqRegions> seqregions1 = nullptr;
    std::unique_ptr<SeqRegions> seqregions2 = nullptr;
    std::unique_ptr<SeqRegions> merged_regions_ptr = nullptr;
    
    // dummy variables
    const RealNumType threshold_prob = params->threshold_prob;
    
    // ----- Test 1 -----
    genTestData(seqregions1, seqregions2, tree, threshold_prob, 1);

    seqregions1->mergeTwoLowers<4>(merged_regions_ptr,
                                   3.3454886086112878360986772063867533688608091324568e-05,
                                   *seqregions2, -1, tree.aln, tree.model,
                                   tree.cumulative_rate, threshold_prob);
    
    EXPECT_EQ(merged_regions_ptr->size(), 17);
    SeqRegion::LHType lh_value1_3{0.17283666036506234453540287177020218223333358764648,
        0.82716302666692920197988314612302929162979125976562,
        1.5848957658592456808424066543733443879204969562124e-08,
        2.9711905061841446462557632390844020164877292700112e-07};
    EXPECT_EQ(*merged_regions_ptr->at(3).likelihood, lh_value1_3);
    SeqRegion::LHType lh_value1_11{3.6381584796761232511871668248793279532016242683312e-09,
        0.67305169363907879631625519323279149830341339111328,
        0.32694703347294445938686635599879082292318344116211,
        1.2692498181937419931533585909511074873989855404943e-06};
    EXPECT_EQ(*merged_regions_ptr->at(11).likelihood, lh_value1_11);
    SeqRegion::LHType lh_value1_15{3.4523387458227984010396742980226678088051528447977e-10,
        5.1728534876242736031974253901183358195225991948973e-08,
        0.855121443561348115736109321005642414093017578125,
        0.14487850436488311500760062244808068498969078063965};
    EXPECT_EQ(*merged_regions_ptr->at(15).likelihood, lh_value1_15);
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    genTestData(seqregions1, seqregions2, tree, threshold_prob, 2);
    
    seqregions1->mergeTwoLowers<4>(merged_regions_ptr,
                                   0.00023418420260279014175064382641267002327367663383484,
                                   *seqregions2,
                                   3.3454886086112878360986772063867533688608091324568e-05,
                                   tree.aln, tree.model, tree.cumulative_rate, threshold_prob);
    
    EXPECT_EQ(merged_regions_ptr->size(), 19);
    SeqRegion::LHType lh_value2_5{1.3692647640882370534384294448460028803538079955615e-07,
        0.17101029958656552287799001987878000363707542419434,
        1.7000484440343382920856089499106467144429188920185e-06,
        0.82898786343851404989635511810774914920330047607422};
    EXPECT_EQ(*merged_regions_ptr->at(5).likelihood, lh_value2_5);
    SeqRegion::LHType lh_value2_9{1.5290708402721638054255510075912782852469717909116e-07,
        0.23349690888349838857607210229616612195968627929688,
        1.8665290161198829580602295483138242104814708000049e-06,
        0.76650107168040160221522683059447444975376129150391};
    EXPECT_EQ(*merged_regions_ptr->at(9).likelihood, lh_value2_9);
    SeqRegion::LHType lh_value2_11{0.12994397914766378510087463382660644128918647766113,
        9.0832943898191744596797318408998300753864896250889e-07,
        0.87005234005882758907546303817071020603179931640625,
        2.772464069444146518881666799161145320340438047424e-06};
    EXPECT_EQ(*merged_regions_ptr->at(11).likelihood, lh_value2_11);
    SeqRegion::LHType lh_value2_13{1.5290708402721638054255510075912782852469717909116e-07,
        0.23349690888349838857607210229616612195968627929688,
        1.8665290161198829580602295483138242104814708000049e-06,
        0.76650107168040160221522683059447444975376129150391};
    EXPECT_EQ(*merged_regions_ptr->at(13).likelihood, lh_value2_13);
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    genTestData(seqregions1, seqregions2, tree, threshold_prob, 3);
    
    seqregions1->mergeTwoLowers<4>(merged_regions_ptr,
                                   0.00020072931651667725661339347631439977703848853707314,
                                   *seqregions2,
                                   0.00013381954434445151344394708825547013475443236529827,
                                   tree.aln, tree.model,
                                   tree.cumulative_rate, threshold_prob);
    
    EXPECT_EQ(merged_regions_ptr->size(), 19);
    SeqRegion::LHType lh_value3_5{0.44582344618258057389326154407171998172998428344727,
        2.8365179941736863494035521260361321083109942264855e-06,
        0.55416504423022705516643782175378873944282531738281,
        8.6730691981839275818714704580081331641849828884006e-06};
    EXPECT_EQ(*merged_regions_ptr->at(5).likelihood, lh_value3_5);
    SeqRegion::LHType lh_value3_7{3.5786639314319516086646010867566847224452430964448e-07,
        0.55147105941417706720386604501982219517230987548828,
        4.3679568805966345593672257863193664206846733577549e-06,
        0.44852421476254927812377104601182509213685989379883};
    EXPECT_EQ(*merged_regions_ptr->at(7).likelihood, lh_value3_7);
    SeqRegion::LHType lh_value3_9{3.5786639314319516086646010867566847224452430964448e-07,
        0.55147105941417706720386604501982219517230987548828,
        4.3679568805966345593672257863193664206846733577549e-06,
        0.44852421476254927812377104601182509213685989379883};
    EXPECT_EQ(*merged_regions_ptr->at(9).likelihood, lh_value3_9);
    SeqRegion::LHType lh_value3_15{9.6067589567314981142138672126824273933554110271871e-11,
        0.99999997772928284067717186189838685095310211181641,
        2.279165015439628010410746761585116387793803482964e-10,
        2.1946733076957723211148276395987544162835547467694e-08};
    EXPECT_EQ(*merged_regions_ptr->at(15).likelihood, lh_value3_15);
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    genTestData(seqregions1, seqregions2, tree, threshold_prob, 4);
    
    seqregions1->mergeTwoLowers<4>(merged_regions_ptr,
                                   3.3454886086112878360986772063867533688608091324568e-05,
                                   *seqregions2, -1, tree.aln, tree.model,
                                   tree.cumulative_rate, threshold_prob);
    
    EXPECT_EQ(merged_regions_ptr->size(), 22);
    SeqRegion::LHType lh_value4_2{4.1908434264100559969947870306712341473276417502802e-11,
        0.27781484807075007559262758149998262524604797363281,
        7.7989680049056750533218225230345024834299749727506e-10,
        0.72218515110744463392222769471118226647377014160156};
    EXPECT_EQ(*merged_regions_ptr->at(2).likelihood, lh_value4_2);
    SeqRegion::LHType lh_value4_12{4.1908434264100559969947870306712341473276417502802e-11,
        0.27781484807075007559262758149998262524604797363281,
        7.7989680049056750533218225230345024834299749727506e-10,
        0.72218515110744463392222769471118226647377014160156};
    EXPECT_EQ(*merged_regions_ptr->at(12).likelihood, lh_value4_12);
    SeqRegion::LHType lh_value4_14{0.28039696642765232770244665516656823456287384033203,
        9.3441446340325357839424550836790808738818725487363e-10,
        0.71960303046228701884245992914657108485698699951172,
        2.1756462036521902336810207540711628593221860228368e-09};
    EXPECT_EQ(*merged_regions_ptr->at(14).likelihood, lh_value4_14);
    SeqRegion::LHType lh_value4_16{4.1908434264100559969947870306712341473276417502802e-11,
        0.27781484807075007559262758149998262524604797363281,
        7.7989680049056750533218225230345024834299749727506e-10,
        0.72218515110744463392222769471118226647377014160156};
    EXPECT_EQ(*merged_regions_ptr->at(16).likelihood, lh_value4_16);
    // ----- Test 4 -----
    
    /*// ----- Test 5 -----
    genTestData(seqregions1, seqregions2, tree, threshold_prob, 5);
    
    seqregions1->mergeTwoLowers<4>(merged_regions_ptr, -1,
     *seqregions2, 5.0182329129169314153348369078599944259622134268284e-05,
     tree.aln, tree.model, tree.cumulative_rate, threshold_prob);
    
    EXPECT_EQ(merged_regions_ptr->size(), 15);
    SeqRegion::LHType lh_value5_6{1.2467371151780668072757457637517002069227345373292e-09,
     0.64033070070634878767634745599934831261634826660156,
     1.2180634035447382694533234245848341004148096544668e-07,
     0.3596691762405736514374154921824811026453971862793};
    EXPECT_EQ(*merged_regions_ptr->at(6).likelihood, lh_value5_6);
    // ----- Test 5 -----
    
    // ----- Test 6 -----
    genTestData(seqregions1, seqregions2, tree, threshold_prob, 6);
    
    seqregions1->mergeTwoLowers<4>(merged_regions_ptr,
     3.3454886086112878360986772063867533688608091324568e-05,
     *seqregions2, -1, tree.aln, tree.model, tree.cumulative_rate, threshold_prob);
    
    EXPECT_EQ(merged_regions_ptr->size(), 21);
    SeqRegion::LHType lh_value6_10{4.1908434264100559969947870306712341473276417502802e-11,
     0.27781484807075007559262758149998262524604797363281,
     7.7989680049056750533218225230345024834299749727506e-10,
     0.72218515110744463392222769471118226647377014160156};
    EXPECT_EQ(*merged_regions_ptr->at(10).likelihood, lh_value6_10);
    SeqRegion::LHType lh_value6_12{0.28039696642765232770244665516656823456287384033203,
     9.3441446340325357839424550836790808738818725487363e-10,
     0.71960303046228701884245992914657108485698699951172,
     2.1756462036521902336810207540711628593221860228368e-09};
    EXPECT_EQ(*merged_regions_ptr->at(12).likelihood, lh_value6_12);
    SeqRegion::LHType lh_value6_14{4.1908434264100559969947870306712341473276417502802e-11,
     0.27781484807075007559262758149998262524604797363281,
     7.7989680049056750533218225230345024834299749727506e-10,
     0.72218515110744463392222769471118226647377014160156};
    EXPECT_EQ(*merged_regions_ptr->at(14).likelihood, lh_value6_14);
    // ----- Test 6 -----
    
    // ----- Test 7 -----
    genTestData(seqregions1, seqregions2, tree, threshold_prob, 7);
    
    seqregions1->mergeTwoLowers<4>(merged_regions_ptr, -1,
     *seqregions2, -1, tree.aln, tree.model, tree.cumulative_rate, threshold_prob);
    
    EXPECT_EQ(merged_regions_ptr->size(), 14);
    // ----- Test 7 -----
    
    // ----- Test 8 -----
    genTestData(seqregions1, seqregions2, tree, threshold_prob, 8);
    
    seqregions1->mergeTwoLowers<4>(merged_regions_ptr,
     3.3454886086112878360986772063867533688608091324568e-05,
     *seqregions2, 3.3454886086112878360986772063867533688608091324568e-05,
     tree.aln, tree.model, tree.cumulative_rate, threshold_prob);
    
    EXPECT_EQ(merged_regions_ptr->size(), 20);
    SeqRegion::LHType lh_value8_6{4.8049629706802666730348423023214121663215792068513e-08,
     0.41416031851624868220795860906946472823619842529297,
     6.5552353400052022819506376133391611915612884331495e-07,
     0.58583897791058758830473607304156757891178131103516};
    EXPECT_EQ(*merged_regions_ptr->at(6).likelihood, lh_value8_6);
    SeqRegion::LHType lh_value8_10{6.3671417399935346464054039278962493497715513512958e-08,
     0.68082356312456637770225142958224751055240631103516,
     7.7723376147377317832898226818150178019095619674772e-07,
     0.31917559597025468853814800240797922015190124511719};
    EXPECT_EQ(*merged_regions_ptr->at(10).likelihood, lh_value8_10);
    SeqRegion::LHType lh_value8_12{0.511129132597716306918300688266754150390625,
     5.1037610231937455663000276911978048133278207387775e-07,
     0.48886879922199172332497596471512224525213241577148,
     1.5578041896004034111807489901280199262600945075974e-06};
    EXPECT_EQ(*merged_regions_ptr->at(12).likelihood, lh_value8_12);
    // ----- Test 8 -----
    
    // ----- Test 9 -----
    genTestData(seqregions1, seqregions2, tree, threshold_prob, 9);
    
    seqregions1->mergeTwoLowers<4>(merged_regions_ptr,
     3.3454886086112878360986772063867533688608091324568e-05,
     *seqregions2, 6.6909772172225756721973544127735067377216182649136e-05,
     tree.aln, tree.model, tree.cumulative_rate, threshold_prob);
    
    EXPECT_EQ(merged_regions_ptr->size(), 23);
    SeqRegion::LHType lh_value9_5{6.7956150015771003729053104797097617506551614496857e-08,
     0.58574313490880203225685818324564024806022644042969,
     9.2710091393484677193322779364947905378357972949743e-07,
     0.41425587003413394393547264371591154485940933227539};
    EXPECT_EQ(*merged_regions_ptr->at(5).likelihood, lh_value9_5);
    SeqRegion::LHType lh_value9_9{6.7956150015771003729053104797097617506551614496857e-08,
     0.58574313490880203225685818324564024806022644042969,
     9.2710091393484677193322779364947905378357972949743e-07,
     0.41425587003413394393547264371591154485940933227539};
    EXPECT_EQ(*merged_regions_ptr->at(9).likelihood, lh_value9_9);
    SeqRegion::LHType lh_value9_11{0.10110884839269537738282167538272915408015251159668,
     1.8423977651254916097878591187131380735308994189836e-07,
     0.89889018932091035996023720144876278936862945556641,
     7.7804661786507317979504462848727719403996161418036e-07};
    EXPECT_EQ(*merged_regions_ptr->at(11).likelihood, lh_value9_11);
    SeqRegion::LHType lh_value9_13{6.7956150015771003729053104797097617506551614496857e-08,
     0.58574313490880203225685818324564024806022644042969,
     9.2710091393484677193322779364947905378357972949743e-07,
     0.41425587003413394393547264371591154485940933227539};
    EXPECT_EQ(*merged_regions_ptr->at(13).likelihood, lh_value9_13);
    SeqRegion::LHType lh_value9_15{5.0746127850612244488463607869008220596640512667364e-08,
     0.75113284338560548647478753991890698671340942382812,
     6.2643606210000642574189759492764295600863988511264e-07,
     0.24886647943220452372514728267560712993144989013672};
    EXPECT_EQ(*merged_regions_ptr->at(15).likelihood, lh_value9_15);
    SeqRegion::LHType lh_value9_19{0.8262325268812187317735151736997067928314208984375,
     2.1300026724482188061459558836574501583527307957411e-05,
     3.0545199061262617574834254963178636899101547896862e-05,
     0.17371562789299538343001927387376781553030014038086};
    EXPECT_EQ(*merged_regions_ptr->at(19).likelihood, lh_value9_19);
    // ----- Test 9 -----
    
    // ----- Test 10 -----
    genTestData(seqregions1, seqregions2, tree, threshold_prob, 10);
    
    seqregions1->mergeTwoLowers<4>(merged_regions_ptr, -1, *seqregions2, -1, tree.aln,
     tree.model, tree.cumulative_rate, threshold_prob);
    
    EXPECT_EQ(merged_regions_ptr->size(), 11);
    // ----- Test 10 -----*/
}
