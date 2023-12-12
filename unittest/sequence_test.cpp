#include "gtest/gtest.h"
#include "../alignment/sequence.h"

using namespace cmaple;

/*
 Test the default Sequence() and Sequence(std::string n_seq_name) constructors
 */
TEST(Sequence, constructor_1)
{
    Sequence sequence1;
    Sequence sequence2("");
    Sequence sequence3("sequence name");
    
    EXPECT_EQ(sequence1.seq_name, "");
    EXPECT_EQ(sequence1.size(), 0);
    
    EXPECT_EQ(sequence2.seq_name, "");
    EXPECT_EQ(sequence2.size(), 0);
    
    EXPECT_EQ(sequence3.seq_name, "sequence name");
    EXPECT_EQ(sequence3.size(), 0);
}

/*
 Test Sequence(std::string n_seq_name, vector<Mutation> n_mutations) constructor
*/
TEST(Sequence, constructor_2)
{
    // detect the path to the example directory
    std::string example_dir = "../../example/";
    if (!fileExists(example_dir + "example.maple"))
        example_dir = "../example/";
    
    std::vector<Mutation> mutations1;
    Sequence sequence1("sequence name 1", std::move(mutations1));
    EXPECT_EQ(sequence1.seq_name, "sequence name 1");
    EXPECT_EQ(sequence1.size(), 0);
    
    Alignment aln(example_dir + "test_5K.maple");
    Sequence sequence2("sequence name 2", std::move(aln.data[0]));
    EXPECT_EQ(sequence2.seq_name, "sequence name 2");
    EXPECT_EQ(sequence2.size(), 5);
    EXPECT_EQ(sequence2[0].type, 1);
    EXPECT_EQ(sequence2[2].position, 14407);
    EXPECT_EQ(sequence2[4].getLength(), 1);
}

/*
 Test operators = and Move Ctor
*/
TEST(Sequence, operators)
{
    // detect the path to the example directory
    std::string example_dir = "../../example/";
    if (!fileExists(example_dir + "example.maple"))
        example_dir = "../example/";
    
    Alignment aln(example_dir + "test_5K.maple");
    Sequence sequence1("sequence name 1", std::move(aln.data[10]));
    
    Sequence sequence2(move(sequence1));
    EXPECT_EQ(sequence1.seq_name, "");
    EXPECT_EQ(sequence1.size(), 0);
    EXPECT_EQ(sequence2.seq_name, "sequence name 1");
    EXPECT_EQ(sequence2.size(), 6);
    EXPECT_EQ(sequence2[0].type, 1);
    EXPECT_EQ(sequence2[2].position, 3036);
    EXPECT_EQ(sequence2[4].getLength(), 1);
    
    Sequence sequence3 = std::move(sequence2);
    EXPECT_EQ(sequence2.seq_name, "");
    EXPECT_EQ(sequence2.size(), 0);
    EXPECT_EQ(sequence3.seq_name, "sequence name 1");
    EXPECT_EQ(sequence3.size(), 6);
    EXPECT_EQ(sequence3[0].type, 1);
    EXPECT_EQ(sequence3[2].position, 3036);
    EXPECT_EQ(sequence3[4].getLength(), 1);
}

/*
 Test SeqRegions* getLowerLhVector(PositionType sequence_length, StateType num_states, SeqType seq_type)
*/
TEST(Sequence, getLowerLhVector)
{
    // detect the path to the example directory
    std::string example_dir = "../../example/";
    if (!fileExists(example_dir + "example.maple"))
        example_dir = "../example/";
    
    Sequence sequence1;
    std::unique_ptr<SeqRegions> seqregions1 = sequence1.getLowerLhVector(30000, 4, cmaple::SeqRegion::SEQ_DNA);
    EXPECT_EQ(seqregions1->size(), 1);
    SeqRegion& seqregion0 = seqregions1->data()[0];
    EXPECT_EQ(seqregion0.type, TYPE_R);
    EXPECT_EQ(seqregion0.position, 30000 - 1);
    EXPECT_EQ(seqregion0.plength_observation2root, -1);
    EXPECT_EQ(seqregion0.plength_observation2node, -1);
    EXPECT_EQ(seqregion0.likelihood, nullptr);
    
    Alignment aln(example_dir + "test_5K.maple");
    std::unique_ptr<SeqRegions> seqregions2 = aln.data[2].getLowerLhVector(aln.ref_seq.size(), aln.num_states, aln.getSeqType());
    EXPECT_EQ(seqregions2->size(), 11);
    EXPECT_EQ(seqregions2->data()[0].type, TYPE_R);
    EXPECT_EQ(seqregions2->data()[1].position, 240);
    EXPECT_EQ(seqregions2->data()[2].plength_observation2root, -1);
    EXPECT_EQ(seqregions2->data()[3].plength_observation2node, -1);
    EXPECT_EQ(seqregions2->data()[4].likelihood, nullptr);
    EXPECT_EQ(seqregions2->data()[5].likelihood, nullptr);
    EXPECT_EQ(seqregions2->data()[7].type, 0);
    EXPECT_EQ(seqregions2->data()[10].position, 29890);
}
