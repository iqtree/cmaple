#include "gtest/gtest.h"
#include "../maple/cmaple.h"

using namespace cmaple;

TEST(CMapleTest, SetAlignmentTest) {
    CMaple cmaple;
    
    // Set up test data
    std::string aln_filename = "../../example/test_100.fa";
    std::string format = "FASTA";
    std::string seqtype = "DNA";

    // Call the method under test
    int status = cmaple.setAlignment(aln_filename, format, seqtype);

    // Check the result
    ASSERT_EQ(status, 0);
    // Add additional assertions if needed
}

TEST(CMapleTest, SetModelTest) {
    CMaple cmaple;

    // Set up test data and configuration
    std::string model_name = "GTR";

    // Call the method under test
    int status = cmaple.setModel(model_name);

    // Check the result
    ASSERT_EQ(status, 0);
    // Add additional assertions if needed
}

TEST(CMapleTest, RunInferenceTest) {
    CMaple cmaple;

    // Set up test data and configuration
    cmaple.setAlignment("../../example/test_100.maple");
    cmaple.setModel("GTR");
    bool force_rerun = false;
    std::string tree_type = "BIN";

    // Call the method under test
    cmaple.overwriteOutputs(true);
    int status = cmaple.runInference(force_rerun, tree_type);

    // Check the result
    ASSERT_EQ(status, 0);
    // Add additional assertions if needed
}

TEST(CMapleTest, ComputeBranchSupportsTest) {
    CMaple cmaple;

    // Set up test data and configuration
    cmaple.setAlignment("../../example/test_100.maple");
    cmaple.setModel("GTR");
    bool force_rerun = false;
    int num_threads = 4;
    int num_replicates = 1000;
    double epsilon = 0.05;

    // Call the method under test
    cmaple.overwriteOutputs(true);
    int status = cmaple.computeBranchSupports(force_rerun, num_threads, num_replicates, epsilon);

    // Check the result
    ASSERT_EQ(status, 0);
    // Add additional assertions if needed
}

TEST(CMapleTest, SetInputTreeTest) {
    CMaple cmaple;

    // Set up test data
    std::string tree_filename = "../../example/test_100.treefile";

    // Call the method under test
    int status = cmaple.setInputTree(tree_filename);

    // Check the result
    ASSERT_EQ(status, 0);
    // Add additional assertions if needed
}

TEST(CMapleTest, ExtractMapleTest) {
    CMaple cmaple;

    // Set up test data
    std::string aln_filename = "../../example/test_100.fa";
    std::string output_filename = "output.maple";

    // Call the method under test
    cmaple.overwriteOutputs(true);
    int status = cmaple.extractMaple(aln_filename, output_filename);

    // Check the result
    ASSERT_EQ(status, 0);
    // Add additional assertions if needed
}

TEST(CMapleTest, ExtractFASTATest) {
    CMaple cmaple;

    // Set up test data
    std::string aln_filename = "../../example/test_100.maple";
    std::string output_filename = "output.fa";

    // Call the method under test
    int status = cmaple.extractFASTA(aln_filename, output_filename);

    // Check the result
    ASSERT_EQ(status, 0);
    // Add additional assertions if needed
}

TEST(CMapleTest, GetTreeStringTest) {
    CMaple cmaple;

    // Set up test data
    cmaple.setAlignment("../../example/test_100.maple");
    cmaple.overwriteOutputs(true);
    cmaple.runInference();
    std::string tree_type = "BIN";
    bool show_branch_supports = true;

    // Call the method under test
    std::string treeString = cmaple.getTreeString(tree_type, show_branch_supports);

    // Check the result
    // Add assertions to validate the treeString if needed
}

TEST(CMapleTest, GetModelStringTest) {
    CMaple cmaple;

    // Set up test data
    cmaple.setModel("GTR");

    // Call the method under test
    std::string modelString = cmaple.getModelString();

    // Check the result
    // Add assertions to validate the modelString if needed
}

// Add more tests for the remaining methods of the CMaple class



TEST(CMapleTest, GetVersionTest) {
    CMaple cmaple;

    // Call the method under test
    std::string version = cmaple.getVersion();

    // Check the result
    // Add assertions to validate the version if needed
}

TEST(CMapleTest, GetCitationsTest) {
    CMaple cmaple;

    // Call the method under test
    std::string citations = cmaple.getCitations();

    // Check the result
    // Add assertions to validate the citations if needed
}

TEST(CMapleTest, SetRefTest) {
    CMaple cmaple;

    // Set up test data
    std::string ref_filename = "../../example/ref.fa";

    // Call the method under test
    int status = cmaple.setRef(ref_filename);

    // Check the result
    ASSERT_EQ(status, 0);
    // Add additional assertions if needed
}

TEST(CMapleTest, SetTreeSearchTypeTest) {
    CMaple cmaple;

    // Set up test data
    std::string tree_search_type = "FAST";

    // Call the method under test
    int status = cmaple.setTreeSearchType(tree_search_type);

    // Check the result
    ASSERT_EQ(status, 0);
    // Add additional assertions if needed
}

TEST(CMapleTest, SetShortRangeTreeSearchTest) {
    CMaple cmaple;

    // Set up test data
    bool enable = true;

    // Call the method under test
    int status = cmaple.setShortRangeTreeSearch(enable);

    // Check the result
    ASSERT_EQ(status, 0);
    // Add additional assertions if needed
}

TEST(CMapleTest, SetPrefixTest) {
    CMaple cmaple;

    // Set up test data
    std::string prefix = "output_";

    // Call the method under test
    int status = cmaple.setPrefix(prefix);

    // Check the result
    ASSERT_EQ(status, 0);
    // Add additional assertions if needed
}

TEST(CMapleTest, EnableOptimizingBlengthsTest) {
    CMaple cmaple;

    // Set up test data
    bool enable = true;

    // Call the method under test
    int status = cmaple.enableOptimizingBlengths(enable);

    // Check the result
    ASSERT_EQ(status, 0);
    // Add additional assertions if needed
}

TEST(CMapleTest, SetMinBlengthTest) {
    CMaple cmaple;

    // Set up test data
    double min_blength = 1e-9;

    // Call the method under test
    int status = cmaple.setMinBlength(min_blength);

    // Check the result
    ASSERT_EQ(status, 0);
    // Add additional assertions if needed
}

TEST(CMapleTest, SetThreshProbTest) {
    CMaple cmaple;

    // Set up test data
    double thresh_prob = 1e-6;

    // Call the method under test
    int status = cmaple.setThreshProb(thresh_prob);

    // Check the result
    ASSERT_EQ(status, 0);
    // Add additional assertions if needed
}

TEST(CMapleTest, OverwriteOutputsTest) {
    CMaple cmaple;

    // Set up test data
    bool enable = true;

    // Call the method under test
    int status = cmaple.overwriteOutputs(enable);

    // Check the result
    ASSERT_EQ(status, 0);
    // Add additional assertions if needed
}

// Add more tests for the remaining methods of the CMaple class



TEST(CMapleTest, SetRandomSeedTest) {
    CMaple cmaple;

    // Set up test data
    int seed = 1234;

    // Call the method under test
    int status = cmaple.setRandomSeed(seed);

    // Check the result
    ASSERT_EQ(status, 0);
    // Add additional assertions if needed
}

TEST(CMapleTest, GetSettingsTest) {
    CMaple cmaple;

    // Call the method under test
    cmaple::Params& settings = cmaple.getSettings();

    // Check the result
    // Add assertions to validate the settings if needed
}
