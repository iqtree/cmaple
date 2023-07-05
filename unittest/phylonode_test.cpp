#include "gtest/gtest.h"
#include "../tree/phylonode.h"
#include "../model/model_dna.h"

using namespace cmaple;

/*
    Test constructors
    Also test get/setSeqNameIndex()
 */
TEST(PhyloNode, TestConstructors)
{
    // make sure it fits on a cacheline
    EXPECT_LE(sizeof(PhyloNode), 64);
    
    // default constructor PhyloNode()
    PhyloNode node1((InternalNode()));
    EXPECT_TRUE(node1.isInternal());
    // invalid access SeqNameIndex (from an internal node)
    EXPECT_DEATH(node1.getSeqNameIndex(), ".*");
    EXPECT_DEATH(node1.setSeqNameIndex(3), ".*");
    
    // PhyloNode(LeafNode&& leaf)
    LeafNode leaf1(100);
    SeqRegions seqregions1;
    seqregions1.emplace_back(TYPE_R, 382, -1, 0.321);
    seqregions1.emplace_back(TYPE_N, 654);
    seqregions1.emplace_back(0, 655, 0);
    leaf1.partial_lh_ = std::make_unique<SeqRegions>(std::move(seqregions1));
    
    PhyloNode node2(std::move(leaf1));
    EXPECT_EQ(leaf1.partial_lh_, nullptr);
    EXPECT_EQ(leaf1.seq_name_index_, 100);
    EXPECT_FALSE(node2.isInternal());
    EXPECT_EQ(node2.getSeqNameIndex(), 100);
    node2.setSeqNameIndex(200);
    EXPECT_EQ(node2.getSeqNameIndex(), 200);
    EXPECT_EQ(node2.getPartialLh(TOP)->size(), 3);
    
    PhyloNode node3(std::move(leaf1));
    EXPECT_FALSE(node3.isInternal());
    EXPECT_EQ(node3.getSeqNameIndex(), 100);
    EXPECT_EQ(node3.getPartialLh(TOP), nullptr);
}

/*
    Test get/setTotalLh()
 */
TEST(PhyloNode, TestSetGetTotalLh) {
    PhyloNode node((InternalNode()));
    std::unique_ptr<SeqRegions> total_lh = std::make_unique<SeqRegions>();
    node.setTotalLh(std::move(total_lh));
    EXPECT_EQ(total_lh, nullptr);
    EXPECT_EQ(node.getTotalLh()->size(), 0);
    node.getTotalLh()->emplace_back(TYPE_R, 382, -1, 0.321);
    node.getTotalLh()->emplace_back(TYPE_N, 654);
    EXPECT_EQ(node.getTotalLh()->size(), 2);
}

/*
    Test get/setMidBranchLh()
 */
TEST(PhyloNode, TestSetGetMidBranchLh) {
    PhyloNode node((InternalNode()));
    std::unique_ptr<SeqRegions> mid_branch_lh = std::make_unique<SeqRegions>();
    node.setMidBranchLh(std::move(mid_branch_lh));
    EXPECT_EQ(mid_branch_lh, nullptr);
    EXPECT_EQ(node.getMidBranchLh()->size(), 0);
    node.getMidBranchLh()->emplace_back(TYPE_R, 382, -1, 0.321);
    EXPECT_EQ(node.getMidBranchLh()->size(), 1);
}

/*
    Test get/setIsOutdated()
 */
TEST(PhyloNode, TestSetGetOutdated) {
    PhyloNode node((InternalNode()));
    node.setOutdated(true);
    EXPECT_TRUE(node.isOutdated());
    
    node.setOutdated(false);
    EXPECT_FALSE(node.isOutdated());
}

/*
    Test get/setSPRApplied()
 */
TEST(PhyloNode, TestSetGetSPRApplied) {
    PhyloNode node((InternalNode()));
    node.setSPRApplied(true);
    EXPECT_TRUE(node.isSPRApplied());
    
    node.setSPRApplied(false);
    EXPECT_FALSE(node.isSPRApplied());
}

/*
    Test get/setUpperLength()
 */
TEST(PhyloNode, TestSetGetUpperLength) {
    PhyloNode node((InternalNode()));
    float blength = 1.23;
    node.setUpperLength(blength);
    EXPECT_EQ(float(node.getUpperLength()), blength);
    
    // set upper length
    node.setUpperLength(0);
    EXPECT_EQ(node.getUpperLength(), 0);
    
    // set upper length
    node.setUpperLength(-0.5);
    EXPECT_EQ(node.getUpperLength(), -0.5);
}

/*
    Test get/setCorrespondingLength()
 */
TEST(PhyloNode, TestSetGetCorrespondingLength) {
    const int NUM_INTERNALS = 10;
    std::vector<PhyloNode> nodes;
    nodes.reserve(NUM_INTERNALS);
    for (int i = 0; i < NUM_INTERNALS; ++i)
        nodes.emplace_back(InternalNode());
    
    // test on a leaf
    LeafNode leaf1(100);
    PhyloNode node1(std::move(leaf1));
    EXPECT_EQ(node1.getCorrespondingLength(TOP, nodes), 0); //default value
    const float blength = 0.5;
    node1.setCorrespondingLength(TOP, nodes, blength);
    EXPECT_EQ(float(node1.getUpperLength()), blength);
    EXPECT_EQ(node1.getCorrespondingLength(TOP, nodes), node1.getUpperLength());
    EXPECT_EQ(node1.getCorrespondingLength(RIGHT, nodes), node1.getUpperLength());
    EXPECT_EQ(node1.getCorrespondingLength(LEFT, nodes), node1.getUpperLength());
    
    // set another blength
    const float blength1 = 0.75;
    node1.setCorrespondingLength(LEFT, nodes, blength1);
    EXPECT_EQ(float(node1.getUpperLength()), blength1);
    EXPECT_EQ(node1.getCorrespondingLength(TOP, nodes), node1.getUpperLength());
    EXPECT_EQ(node1.getCorrespondingLength(RIGHT, nodes), node1.getUpperLength());
    EXPECT_EQ(node1.getCorrespondingLength(LEFT, nodes), node1.getUpperLength());
    
    // test on an internal node
    PhyloNode node2((InternalNode()));
    node2.setNeighborIndex(TOP, Index(0, LEFT));
    node2.setNeighborIndex(RIGHT, Index(1, TOP));
    node2.setNeighborIndex(LEFT, Index(2, TOP));
    EXPECT_EQ(node2.getUpperLength(), 0);
    EXPECT_EQ(node2.getCorrespondingLength(TOP, nodes), node2.getUpperLength());
    EXPECT_EQ(node2.getCorrespondingLength(RIGHT, nodes), 0);
    EXPECT_EQ(node2.getCorrespondingLength(LEFT, nodes), 0);
    node2.setCorrespondingLength(RIGHT, nodes, blength);
    EXPECT_EQ(node2.getUpperLength(), 0);
    EXPECT_EQ(node2.getCorrespondingLength(TOP, nodes), node2.getUpperLength());
    EXPECT_EQ(node2.getCorrespondingLength(RIGHT, nodes), blength);
    EXPECT_EQ(node2.getCorrespondingLength(LEFT, nodes), 0);
    node2.setCorrespondingLength(TOP, nodes, blength1);
    EXPECT_EQ(node2.getUpperLength(), blength1);
    EXPECT_EQ(node2.getCorrespondingLength(TOP, nodes), node2.getUpperLength());
    EXPECT_EQ(node2.getCorrespondingLength(RIGHT, nodes), blength);
    EXPECT_EQ(node2.getCorrespondingLength(LEFT, nodes), 0);
}

/*
    Test get/setNode(Leaf&&);
    Also test setNeighborIndex(), getNeighborIndex();
 */
TEST(PhyloNode, TestGetSetNodeWithRvalueLeaf)
{
    // Create a leaf node
    LeafNode leaf(42);

    // Create a phylonode to test
    PhyloNode node(std::move(leaf));
    EXPECT_EQ(node.getSeqNameIndex(), 42);
    EXPECT_EQ(node.getNeighborIndex(TOP).getVectorIndex(), 0);
    EXPECT_EQ(node.getNeighborIndex(TOP).getMiniIndex(), UNDEFINED);
    node.setNeighborIndex(TOP, Index(10, LEFT));
    EXPECT_EQ(node.getNeighborIndex(TOP).getVectorIndex(), 10);
    EXPECT_EQ(node.getNeighborIndex(TOP).getMiniIndex(), LEFT);
    node.setNeighborIndex(RIGHT, Index(20, RIGHT));
    EXPECT_EQ(node.getNeighborIndex(RIGHT).getVectorIndex(), 20);
    EXPECT_EQ(node.getNeighborIndex(RIGHT).getMiniIndex(), RIGHT);
    node.setNeighborIndex(LEFT, Index(30, TOP));
    EXPECT_EQ(node.getNeighborIndex(LEFT).getVectorIndex(), 30);
    EXPECT_EQ(node.getNeighborIndex(LEFT).getMiniIndex(), TOP);

    // replace the current node by another leaf
    LeafNode leaf2(100);
    node.setNode(std::move(leaf2));
    EXPECT_EQ(node.getSeqNameIndex(), 100);
    EXPECT_EQ(node.getNode().leaf_.seq_name_index_, 100);
    
    // replace the current node by another internal (invalid)
    InternalNode internal;
    EXPECT_DEATH(node.setNode(std::move(internal)), ".*");
}

/*
    Test get/setNode(InternalNode&&);
    Also test setNeighborIndex(), getNeighborIndex();
 */
TEST(PhyloNode, TestGetSetNodeWithRvalueInternal)
{
    // Create a phylonode to test
    PhyloNode node((InternalNode()));
    node.setNeighborIndex(RIGHT, Index(10, TOP));
    EXPECT_TRUE(node.isInternal());
    EXPECT_EQ(node.getNeighborIndex(TOP).getVectorIndex(), 0);
    EXPECT_EQ(node.getNeighborIndex(TOP).getMiniIndex(), UNDEFINED);
    EXPECT_EQ(node.getNeighborIndex(RIGHT).getVectorIndex(), 10);
    EXPECT_EQ(node.getNeighborIndex(RIGHT).getMiniIndex(), TOP);
    EXPECT_EQ(node.getNeighborIndex(LEFT).getVectorIndex(), 0);
    EXPECT_EQ(node.getNeighborIndex(LEFT).getMiniIndex(), UNDEFINED);

    // replace the current node by another internal
    InternalNode internal;
    node.setNode(std::move(internal));
    EXPECT_TRUE(node.isInternal());
    node.setNeighborIndex(LEFT, Index(55, TOP));
    EXPECT_EQ(node.getNeighborIndex(TOP).getVectorIndex(), 0);
    EXPECT_EQ(node.getNeighborIndex(TOP).getMiniIndex(), UNDEFINED);
    EXPECT_EQ(node.getNeighborIndex(RIGHT).getVectorIndex(), 0);
    EXPECT_EQ(node.getNeighborIndex(RIGHT).getMiniIndex(), UNDEFINED);
    EXPECT_EQ(node.getNeighborIndex(LEFT).getVectorIndex(), 55);
    EXPECT_EQ(node.getNeighborIndex(LEFT).getMiniIndex(), TOP);
    EXPECT_EQ(node.getNode().internal_.neighbor_index3_[0].getVectorIndex(), 0);
    
    // replace the current node by another leaf (invalid)
    LeafNode leaf(0);
    EXPECT_DEATH(node.setNode(std::move(leaf)), ".*");
}

/*
    Test getLessInfoSeqs() and addLessInfoSeqs() functions
 */
TEST(PhyloNode, TestAddGetLessInfoSeqs) {
    LeafNode leaf1(0);
    PhyloNode node1(std::move(leaf1));
    EXPECT_EQ(node1.getLessInfoSeqs().size(), 0);
    node1.addLessInfoSeqs(3);
    node1.addLessInfoSeqs(1);
    node1.addLessInfoSeqs(2);
    node1.addLessInfoSeqs(4);
    EXPECT_EQ(node1.getLessInfoSeqs().size(), 4);
    
    // invalid access lessinfoseqs (from an internal node)
    PhyloNode node2((InternalNode()));
    EXPECT_DEATH(node2.getLessInfoSeqs(), ".*");
    EXPECT_DEATH(node2.addLessInfoSeqs(3), ".*");
}

/*
    Test set/getPartialLh
 */
TEST(PhyloNode, TestGetSetPartialLh)
{
    // Create a leaf node
    LeafNode leaf(42);

    // Create a phylonode to test
    PhyloNode node(std::move(leaf));
    EXPECT_EQ(node.getPartialLh(TOP), nullptr);
    EXPECT_EQ(node.getPartialLh(LEFT), nullptr);
    EXPECT_EQ(node.getPartialLh(RIGHT), nullptr);
    
    // set partial lh for a leaf
    SeqRegions seqregions1;
    seqregions1.emplace_back(TYPE_R, 382, -1, 0.321);
    seqregions1.emplace_back(TYPE_N, 654);
    seqregions1.emplace_back(0, 655, 0);
    node.setPartialLh(TOP, std::make_unique<SeqRegions>(std::move(seqregions1)));
    EXPECT_EQ(node.getPartialLh(TOP)->size(), 3);
    EXPECT_EQ(node.getPartialLh(LEFT), node.getPartialLh(TOP));
    EXPECT_EQ(node.getPartialLh(RIGHT), node.getPartialLh(TOP));
    
    // create another (internal) phylonode
    PhyloNode node2((InternalNode()));
    EXPECT_EQ(node2.getPartialLh(TOP), nullptr);
    EXPECT_EQ(node2.getPartialLh(LEFT), nullptr);
    EXPECT_EQ(node2.getPartialLh(RIGHT), nullptr);
    
    // set partial lh for an internal
    SeqRegions seqregions2;
    seqregions2.emplace_back(TYPE_R, 382, -1, 0.321);
    node2.setPartialLh(TOP, std::make_unique<SeqRegions>(std::move(seqregions2)));
    EXPECT_EQ(node2.getPartialLh(TOP)->size(), 1);
    EXPECT_EQ(node2.getPartialLh(LEFT), nullptr);
    EXPECT_EQ(node2.getPartialLh(RIGHT), nullptr);
    auto seqregions3_ptr = std::make_unique<SeqRegions>();
    seqregions3_ptr->emplace_back(TYPE_N, 654);
    seqregions3_ptr->emplace_back(0, 655, 0);
    node2.setPartialLh(RIGHT, std::move(seqregions3_ptr));
    EXPECT_EQ(node2.getPartialLh(TOP)->size(), 1);
    EXPECT_EQ(node2.getPartialLh(LEFT), nullptr);
    EXPECT_EQ(node2.getPartialLh(RIGHT)->size(), 2);
}

/*
    Test exportString()
 */
TEST(PhyloNode, TestExportString)
{
    const int NUM_SEQS = 10;
    std::unique_ptr<Alignment> aln = std::make_unique<Alignment>();
    // init NUM_SEQS
    aln->data.reserve(NUM_SEQS);
    for (int i =0 ; i < NUM_SEQS; ++i)
    {
        aln->data.emplace_back("sequence " + convertIntToString(i));
    }
    
    // test on an internal node
    PhyloNode node1((InternalNode()));
    EXPECT_EQ(node1.exportString(true, aln, false), ""); // internal node returns ""
    EXPECT_EQ(node1.exportString(false, aln, false), ""); // internal node returns ""
    
    // test on a leaf
    PhyloNode node2(LeafNode(1));
    EXPECT_EQ(node2.exportString(true, aln, false), "sequence 1:0");
    EXPECT_EQ(node2.exportString(false, aln, false), "sequence 1:0");
    
    // add a lessinfoseq
    node2.addLessInfoSeqs(3);
    node2.setUpperLength(-1);
    EXPECT_EQ(node2.exportString(true, aln, false), "(sequence 1:0,sequence 3:0):0");
    EXPECT_EQ(node2.exportString(false, aln, false), "(sequence 1:0,sequence 3:0):0");
    
    // add two more lessinfoseqs
    node2.addLessInfoSeqs(6);
    node2.addLessInfoSeqs(8);
    node2.setUpperLength(0.5);
    EXPECT_EQ(node2.exportString(true, aln, false), "(sequence 1:0,(sequence 3:0,(sequence 6:0,sequence 8:0):0):0):0.5");
    EXPECT_EQ(node2.exportString(false, aln, false), "(sequence 1:0,sequence 3:0,sequence 6:0,sequence 8:0):0.5");
}

/*
    Initialize Alignment, Model, and Parameters
 */
void initTestData(Params& params, std::unique_ptr<Alignment>& aln, std::unique_ptr<ModelBase>& model, const std::string model_name = "GTR")
{
    // Init params, aln, and model
    params.model_name = model_name;
    std::string diff_file_path("../../example/test_5K.maple");
    aln->readMapleFile(diff_file_path, "");
    model = std::make_unique<ModelDNA>(ModelDNA(params.model_name));
    // extract related info (freqs, log_freqs) of the ref sequence
    model->extractRefInfo(aln);
    // init the mutation matrix from a model name
    model->initMutationMat();
    // compute cumulative rates of the ref sequence
    model->computeCumulativeRate(aln);
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
    SeqRegion::LHType new_lh_value{7.4346402191731947664294494204639818235591519623995e-06,0.66666418845326025355291221785591915249824523925781,7.4346402191731947664294494204639818235591519623995e-06,0.33332094226630137878686355179524980485439300537109};
    (*new_lh) = new_lh_value;
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
    Test computeTotalLhAtNode()
 */
TEST(PhyloNode, TestComputeTotalLhAtNode)
{
    std::unique_ptr<Alignment> aln = std::make_unique<Alignment>();
    std::unique_ptr<ModelBase> model = nullptr;
    Params params = Params::getInstance();
    
    // Init params, aln, and model
    initTestData(params, aln, model);
    
    SeqRegions seqregions1, seqregions2;
    genTestData(seqregions1, seqregions2);
    
    // test on a root
    PhyloNode neighbor((InternalNode()));
    PhyloNode node1((InternalNode()));
    std::unique_ptr<SeqRegions> total_lh = nullptr;
    node1.setPartialLh(TOP, std::make_unique<SeqRegions>(std::move(seqregions1)));
    node1.computeTotalLhAtNode<4>(total_lh, neighbor, aln, model, params.threshold_prob, true); // deafault blength = -1
    EXPECT_EQ(total_lh->size(), 9);
    EXPECT_EQ(total_lh->at(0).type, TYPE_R);
    EXPECT_EQ(total_lh->at(1).position, 654);
    EXPECT_EQ(total_lh->at(2).plength_observation2root, 0);
    EXPECT_EQ(total_lh->at(4).plength_observation2node, 0.432);
    EXPECT_EQ(total_lh->at(4).likelihood, nullptr);
    SeqRegion::LHType lh_value{0.0000096923454377779065,0.53355990895437144,0.0000063567736994888323,0.46642404192649117};
    EXPECT_EQ(*total_lh->at(6).likelihood, lh_value);
    
    // test on a non-root node
    aln->ref_seq.resize(3500);
    const MiniIndex parent_mini = RIGHT;
    neighbor.setPartialLh(parent_mini, std::make_unique<SeqRegions>(std::move(seqregions2)));
    node1.setNeighborIndex(TOP, Index(0, parent_mini));
    node1.setUpperLength(0.123);
    node1.computeTotalLhAtNode<4>(total_lh, neighbor, aln, model, params.threshold_prob, false, 141e-5);
    EXPECT_EQ(total_lh->size(), 15);
    EXPECT_EQ(total_lh->at(0).type, TYPE_R);
    EXPECT_EQ(total_lh->at(1).position, 381);
    EXPECT_EQ(total_lh->at(2).plength_observation2root, 0);
    EXPECT_EQ(total_lh->at(3).plength_observation2node, -1);
    EXPECT_EQ(total_lh->at(5).likelihood, nullptr);
    SeqRegion::LHType lh_value1{0.20522485747557062,0.20388879902610788,0.20107341695416275,0.38981292654415878};
    EXPECT_EQ(*total_lh->at(6).likelihood, lh_value1);
    
    // default blength
    node1.setUpperLength(0.012);
    node1.computeTotalLhAtNode<4>(total_lh, neighbor, aln, model, params.threshold_prob, false);
    EXPECT_EQ(total_lh->size(), 13);
    EXPECT_EQ(total_lh->at(0).type, TYPE_R);
    EXPECT_EQ(total_lh->at(1).position, 654);
    EXPECT_EQ(total_lh->at(2).plength_observation2root, -1);
    EXPECT_EQ(total_lh->at(3).plength_observation2node, -1);
    EXPECT_EQ(total_lh->at(5).likelihood, nullptr);
    SeqRegion::LHType lh_value2{0.20054776730537416,0.20040769666419275,0.2001125356462399,0.39893200038419324};
    EXPECT_EQ(*total_lh->at(4).likelihood, lh_value2);
}


