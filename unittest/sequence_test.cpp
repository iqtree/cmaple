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
    std::string name = "sequence name";
    Sequence sequence3(name);
    
    EXPECT_EQ(sequence1.seq_name, "");
    EXPECT_EQ(sequence1.size(), 0);
    
    EXPECT_EQ(sequence2.seq_name, "");
    EXPECT_EQ(sequence2.size(), 0);
    
    EXPECT_EQ(sequence3.seq_name, name);
    EXPECT_EQ(sequence3.size(), 0);
}

/*
 Test Sequence(std::string n_seq_name, vector<Mutation> n_mutations) constructor
*/
TEST(Sequence, constructor_2)
{
    std::vector<Mutation> mutations1;
    std::string name1 = "sequence name 1";
    Sequence sequence1(name1, mutations1);
    EXPECT_EQ(sequence1.seq_name, name1);
    EXPECT_EQ(sequence1.size(), 0);
    
    mutations1.emplace_back(1, 132, 1);
    mutations1.emplace_back(TYPE_N, 5434, 943);
    mutations1.emplace_back(TYPE_O, 9563, 1);
    mutations1.emplace_back(TYPE_R, 38432, 45);
    mutations1.emplace_back(TYPE_DEL, 153, 3);
    std::string name2 = "sequence name 2";
    Sequence sequence2(name2, mutations1);
    EXPECT_EQ(sequence2.seq_name, name2);
    EXPECT_EQ(sequence2.size(), 5);
    EXPECT_EQ(sequence2[0].type, 1);
    EXPECT_EQ(sequence2[2].position, 9563);
    EXPECT_EQ(sequence2[4].getLength(), 3);
}

/*
 Test operators = and Move Ctor
*/
TEST(Sequence, operators)
{
    std::vector<Mutation> mutations1;
    mutations1.emplace_back(1, 132, 1);
    mutations1.emplace_back(TYPE_N, 5434, 943);
    mutations1.emplace_back(TYPE_O, 9563, 1);
    mutations1.emplace_back(TYPE_R, 38432, 45);
    mutations1.emplace_back(TYPE_DEL, 153, 3);
    std::string name1 = "sequence name 1";
    Sequence sequence1(name1, mutations1);
    
    Sequence sequence2(move(sequence1));
    EXPECT_EQ(sequence1.seq_name, "");
    EXPECT_EQ(sequence1.size(), 0);
    EXPECT_EQ(sequence2.seq_name, name1);
    EXPECT_EQ(sequence2.size(), 5);
    EXPECT_EQ(sequence2[0].type, 1);
    EXPECT_EQ(sequence2[2].position, 9563);
    EXPECT_EQ(sequence2[4].getLength(), 3);
    
    Sequence sequence3 = std::move(sequence2);
    EXPECT_EQ(sequence2.seq_name, "");
    EXPECT_EQ(sequence2.size(), 0);
    EXPECT_EQ(sequence3.seq_name, name1);
    EXPECT_EQ(sequence3.size(), 5);
    EXPECT_EQ(sequence3[0].type, 1);
    EXPECT_EQ(sequence3[2].position, 9563);
    EXPECT_EQ(sequence3[4].getLength(), 3);
}

/*
 Test SeqRegions* getLowerLhVector(PositionType sequence_length, StateType num_states, SeqType seq_type)
*/
TEST(Sequence, getLowerLhVector)
{
    Sequence sequence1;
    std::unique_ptr<SeqRegions> seqregions1 = sequence1.getLowerLhVector(30000, 4, cmaple::SeqRegion::SEQ_DNA);
    EXPECT_EQ(seqregions1->size(), 1);
    SeqRegion& seqregion0 = seqregions1->data()[0];
    EXPECT_EQ(seqregion0.type, TYPE_R);
    EXPECT_EQ(seqregion0.position, 30000 - 1);
    EXPECT_EQ(seqregion0.plength_observation2root, -1);
    EXPECT_EQ(seqregion0.plength_observation2node, -1);
    EXPECT_EQ(seqregion0.likelihood, nullptr);
    
    Sequence sequence2;
    sequence2.emplace_back(1, 132, 1);
    sequence2.emplace_back(TYPE_DEL, 153, 3);
    sequence2.emplace_back(1+4+3, 154, 1);
    sequence2.emplace_back(TYPE_N, 5434, 943);
    sequence2.emplace_back(1+2+8+3, 6390, 1);
    sequence2.emplace_back(0, 9563, 1);
    sequence2.emplace_back(TYPE_N, 15209, 8);
    sequence2.emplace_back(3, 28967, 1);
    sequence2.emplace_back(1+2+4+3, 28968, 1);
    
    std::unique_ptr<SeqRegions> seqregions2 = sequence2.getLowerLhVector(30000, 4, cmaple::SeqRegion::SEQ_DNA);
    EXPECT_EQ(seqregions2->size(), 17);
    EXPECT_EQ(seqregions2->data()[0].type, TYPE_R);
    EXPECT_EQ(seqregions2->data()[1].position, 132);
    EXPECT_EQ(seqregions2->data()[2].plength_observation2root, -1);
    EXPECT_EQ(seqregions2->data()[3].plength_observation2node, -1);
    EXPECT_FALSE(seqregions2->data()[4].likelihood == nullptr);
    EXPECT_EQ(seqregions2->data()[5].likelihood, nullptr);
    EXPECT_EQ(seqregions2->data()[8].type, TYPE_O);
    EXPECT_EQ(seqregions2->data()[12].position, 15209 + 8 - 1);
    EXPECT_EQ(seqregions2->data()[15].getLH(0), seqregions2->data()[15].getLH(2)); // Likelihood (1.0 / 3, 1.0 / 3, 1.0 / 3, 0)
    EXPECT_EQ(seqregions2->data()[15].getLH(3), 0);
}
