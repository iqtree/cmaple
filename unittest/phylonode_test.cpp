#include "gtest/gtest.h"
#include "../tree/phylonode.h"
#include "../tree/tree.h"
#include "../model/model.h"

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
#ifdef DEBUG
    EXPECT_DEATH(node1.getSeqNameIndex(), ".*");
    EXPECT_DEATH(node1.setSeqNameIndex(3), ".*");
 #endif

    
    // PhyloNode(LeafNode&& leaf)
    LeafNode leaf1(100);
    SeqRegions seqregions1;
    seqregions1.resize(3);
    leaf1.partial_lh_ = cmaple::make_unique<SeqRegions>(std::move(seqregions1));
    
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
    std::unique_ptr<SeqRegions> total_lh = cmaple::make_unique<SeqRegions>();
    node.setTotalLh(std::move(total_lh));
    EXPECT_EQ(total_lh, nullptr);
    EXPECT_EQ(node.getTotalLh()->size(), 0);
    node.getTotalLh()->emplace_back(TYPE_N, 654);
    EXPECT_EQ(node.getTotalLh()->size(), 1);
}

/*
    Test get/setMidBranchLh()
 */
TEST(PhyloNode, TestSetGetMidBranchLh) {
    PhyloNode node((InternalNode()));
    std::unique_ptr<SeqRegions> mid_branch_lh = cmaple::make_unique<SeqRegions>();
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
    node.setSPRCount(10);
    EXPECT_EQ(node.getSPRCount(), 10);
    
    node.setSPRCount(0);
    EXPECT_EQ(node.getSPRCount(), 0);
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
#ifdef DEBUG
    InternalNode internal;
    EXPECT_DEATH(node.setNode(std::move(internal)), ".*");
#endif
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
#ifdef DEBUG
    LeafNode leaf(0);
    EXPECT_DEATH(node.setNode(std::move(leaf)), ".*");
#endif
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
#ifdef DEBUG
    PhyloNode node2((InternalNode()));
    EXPECT_DEATH(node2.getLessInfoSeqs(), ".*");
    EXPECT_DEATH(node2.addLessInfoSeqs(3), ".*");
#endif
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
    seqregions1.resize(3);
    node.setPartialLh(TOP, cmaple::make_unique<SeqRegions>(std::move(seqregions1)));
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
    node2.setPartialLh(TOP, cmaple::make_unique<SeqRegions>(std::move(seqregions2)));
    EXPECT_EQ(node2.getPartialLh(TOP)->size(), 1);
    EXPECT_EQ(node2.getPartialLh(LEFT), nullptr);
    EXPECT_EQ(node2.getPartialLh(RIGHT), nullptr);
    auto seqregions3_ptr = cmaple::make_unique<SeqRegions>();
    seqregions3_ptr->resize(2);
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
    std::vector<std::string> seq_names;
    // init NUM_SEQS
    seq_names.reserve(NUM_SEQS);
    for (int i =0; i < NUM_SEQS; ++i)
    {
        seq_names.emplace_back("sequence " + convertIntToString(i));
    }
    
    // test on an internal node
    PhyloNode node1((InternalNode()));
    EXPECT_EQ(node1.exportString(true, seq_names, false, false, false, ""), ""); // internal node returns ""
    EXPECT_EQ(node1.exportString(false, seq_names, false, false, false, ""), ""); // internal node returns ""
    
    // test on a leaf
    PhyloNode node2(LeafNode(1));
    EXPECT_EQ(node2.exportString(true, seq_names, false, false, false, ""), "sequence 1:0");
    EXPECT_EQ(node2.exportString(false, seq_names, false, false, false, ""), "sequence 1:0");
    
    // add a lessinfoseq
    node2.addLessInfoSeqs(3);
    node2.setUpperLength(-1);
    EXPECT_EQ(node2.exportString(true, seq_names, false, false, false, ""), "(sequence 1:0,sequence 3:0):0");
    EXPECT_EQ(node2.exportString(false, seq_names, false, false, false, ""), "(sequence 1:0,sequence 3:0):0");
    
    // add two more lessinfoseqs
    node2.addLessInfoSeqs(6);
    node2.addLessInfoSeqs(8);
    node2.setUpperLength(0.5);
    EXPECT_EQ(node2.exportString(true, seq_names, false, false, false, ""),
        "(((sequence 1:0,sequence 3:0):0,sequence 6:0):0,sequence 8:0):0.5");
    EXPECT_EQ(node2.exportString(false, seq_names, false, false, false, ""),
        "(sequence 1:0,sequence 3:0,sequence 6:0,sequence 8:0):0.5");
}

/*
    Test computeTotalLhAtNode()
 */
TEST(PhyloNode, TestComputeTotalLhAtNode)
{
    // detect the path to the example directory
    std::string example_dir = "../../example/";
    if (!fileExists(example_dir + "example.maple"))
        example_dir = "../example/";
    
    Alignment aln(example_dir + "test_5K.maple");
    Model model(aln.ref_seq.size(), false, false, cmaple::ModelBase::GTR);
    Tree tree(&aln, &model);
    std::unique_ptr<Params> params = ParamsBuilder().build();
    
    std::unique_ptr<SeqRegions> seqregions1 = aln.data[0]
        .getLowerLhVector(aln.ref_seq.size(), aln.num_states, aln.getSeqType());
    std::unique_ptr<SeqRegions> seqregions2 = aln.data[10]
        .getLowerLhVector(aln.ref_seq.size(), aln.num_states, aln.getSeqType());
    std::unique_ptr<SeqRegions> seqregions3 = aln.data[100]
        .getLowerLhVector(aln.ref_seq.size(), aln.num_states, aln.getSeqType());
    
    // test on a root
    std::unique_ptr<SeqRegions> merge_regions1 = nullptr;
    seqregions1->mergeTwoLowers<4>(merge_regions1, 1e-5, *seqregions2, 123e-3, tree.aln,
            tree.model, tree.cumulative_rate, params->threshold_prob);
    PhyloNode neighbor((InternalNode()));
    PhyloNode node1((InternalNode()));
    std::unique_ptr<SeqRegions> total_lh = nullptr;
    node1.setPartialLh(TOP, cmaple::make_unique<SeqRegions>(std::move(merge_regions1)));
    node1.computeTotalLhAtNode<4>(total_lh, neighbor, tree.aln,
        tree.model, params->threshold_prob, true); // deafault blength = -1
    EXPECT_EQ(total_lh->size(), 13);
    EXPECT_EQ(total_lh->at(0).type, TYPE_R);
    EXPECT_EQ(total_lh->at(1).position, 240);
    EXPECT_EQ(total_lh->at(2).plength_observation2root, -1);
    EXPECT_EQ(total_lh->at(4).plength_observation2node, -1);
    EXPECT_EQ(total_lh->at(4).likelihood, nullptr);
    SeqRegion::LHType lh_value{0.00007791197759195882027640628342268769301881548017,
        0.99990842564971482708813255158020183444023132324219,
        0.00000303392141268825073045778406566341800498776138,
        0.00001062845128063864023817020748596817725228902418};
    EXPECT_EQ(*total_lh->at(3).likelihood, lh_value);
    
    // test on a non-root node
    std::unique_ptr<SeqRegions> merge_regions2 = nullptr;
    seqregions1->mergeTwoLowers<4>(merge_regions2, 14e-6, *seqregions3, 22e-5,
            tree.aln, tree.model, tree.cumulative_rate, params->threshold_prob);
    const MiniIndex parent_mini = RIGHT;
    neighbor.setPartialLh(parent_mini, cmaple::make_unique<SeqRegions>(std::move(merge_regions2)));
    node1.setNeighborIndex(TOP, Index(0, parent_mini));
    node1.setUpperLength(0.013);
    node1.computeTotalLhAtNode<4>(total_lh, neighbor, tree.aln, tree.model, params->threshold_prob, false, 141e-5);
    EXPECT_EQ(total_lh->size(), 19);
    EXPECT_EQ(total_lh->at(0).type, TYPE_R);
    EXPECT_EQ(total_lh->at(1).position, 240);
    EXPECT_EQ(total_lh->at(2).plength_observation2root, -1);
    EXPECT_EQ(total_lh->at(3).plength_observation2node, 0);
    EXPECT_EQ(total_lh->at(4).likelihood, nullptr);
    SeqRegion::LHType lh_value1{0.00000017075294482921272644353357077900978922002651,
        0.99997870393335785976773877337109297513961791992188,
        0.00000016973979605241345684048996731579928010091862,
        0.00002095557390147909605977206981552996012396761216};
    EXPECT_EQ(*total_lh->at(5).likelihood, lh_value1);
    
    // default blength
    node1.setUpperLength(0.0112);
    node1.computeTotalLhAtNode<4>(total_lh, neighbor, tree.aln, tree.model, params->threshold_prob, false);
    EXPECT_EQ(total_lh->size(), 13);
    EXPECT_EQ(total_lh->at(0).type, TYPE_R);
    EXPECT_EQ(total_lh->at(1).position, 240);
    EXPECT_EQ(total_lh->at(2).plength_observation2root, -1);
    EXPECT_EQ(total_lh->at(3).plength_observation2node, 0);
    EXPECT_EQ(total_lh->at(5).likelihood, nullptr);
    SeqRegion::LHType lh_value2{0.00000005304243814999847067991444317207327951990692,
        0.99999984718379919534925193147500976920127868652344,
        0.00000000314931216745035814445441008362514684337796,
        0.00000009662445065181004761007619837179238864166564};
    EXPECT_EQ(*total_lh->at(3).likelihood, lh_value2);
}


