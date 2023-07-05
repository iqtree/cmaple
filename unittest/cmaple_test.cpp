#include "gtest/gtest.h"
#include "../maple/cmaple.h"

using namespace cmaple;

TEST(CMapleTest, DefaultConstructor) {
    CMaple cmaple;
    cmaple::Params& params = cmaple.getSettings();
    
    EXPECT_EQ(params.aln_path, "");
    EXPECT_EQ(params.maple_path, "");
    EXPECT_EQ(params.ref_path, "");
    EXPECT_EQ(params.aln_format, IN_UNKNOWN);
    EXPECT_FALSE(params.only_extract_maple);
    EXPECT_EQ(params.hamming_weight, 1000);
    EXPECT_EQ(params.model_name, "GTR");
    EXPECT_TRUE(params.optimize_blength);
    EXPECT_FALSE(params.overwrite_output);
    EXPECT_FALSE(params.strict_stop_seeking_placement_sample);
    EXPECT_FALSE(params.strict_stop_seeking_placement_subtree);
    EXPECT_TRUE(params.strict_stop_seeking_placement_subtree_short_search);
    EXPECT_EQ(params.output_aln, nullptr);
    EXPECT_TRUE(params.export_binary_tree);
    EXPECT_FALSE(params.short_range_topo_search);
    EXPECT_EQ(params.output_testing, nullptr);
    EXPECT_FALSE(params.compute_aLRT_SH);
    EXPECT_EQ(params.aLRT_SH_replicates, 1000);
    EXPECT_EQ(params.aLRT_SH_epsilon, 0.05);
    EXPECT_EQ(params.num_threads, 1);
    EXPECT_EQ(params.input_treefile, "");
    EXPECT_EQ(params.output_prefix, "");
    EXPECT_FALSE(params.allow_replace_input_tree);
    EXPECT_EQ(params.seq_type, SEQ_UNKNOWN);
    EXPECT_EQ(params.tree_search_type, PARTIAL_TREE_SEARCH);
}

TEST(CMapleTest, DefaultParamsConstructor) {
    // Test case 1
    cmaple::Params default_params;
    cmaple::Params params1;
    CMaple cmaple1;
    cmaple::Params& params = cmaple1.getSettings();
    params = std::move(params1);
    cmaple::Params& params1_ref = cmaple1.getSettings();
    EXPECT_EQ(params1_ref.aln_path, default_params.aln_path);
    EXPECT_EQ(params1_ref.maple_path, default_params.aln_path);
    EXPECT_EQ(params1_ref.ref_path, default_params.aln_path);
    
    // Test case 2
    cmaple::Params params2;
    params2.aln_path = "../../example/test_100.fa";
    params2.model_name = "GTR";
    params2.num_threads = 4;
    CMaple cmaple2(std::move(params2));
    cmaple::Params& params2_ref = cmaple2.getSettings();
    EXPECT_EQ(params2_ref.aln_path, "../../example/test_100.fa");
    EXPECT_EQ(params2_ref.model_name, "GTR");
    EXPECT_EQ(params2_ref.num_threads, 4);

    // Test case 3
    cmaple::Params params3;
    params3.aln_path = "../../example/test_100.fa";
    params3.model_name = "JC";
    params3.num_threads = 2;
    
    CMaple cmaple3(std::move(params3));
    // Test that the original params object has been moved
    EXPECT_EQ(params3.aln_path, "");
    EXPECT_EQ(params3.model_name, "");
    EXPECT_EQ(params3.num_threads, 2);
}

TEST(CMapleTest, ValidAlignmentConstructor) {
    // Test case 1
    std::string aln_filename = "alignment.fasta";
    std::string format = "FASTA";
    std::string seqtype = "AA";
    CMaple cmaple1(aln_filename, format, seqtype);
    cmaple::Params& params1 = cmaple1.getSettings();
    EXPECT_EQ(params1.aln_path, "alignment.fasta");
    EXPECT_EQ(params1.aln_format, IN_FASTA);
    EXPECT_EQ(params1.seq_type, SEQ_PROTEIN);
    
    // Test case 2
    CMaple cmaple2("alignment.fasta", "", "DNA");
    cmaple::Params& params2 = cmaple2.getSettings();
    EXPECT_EQ(params2.aln_path, "alignment.fasta");
    EXPECT_EQ(params2.aln_format, IN_UNKNOWN);
    EXPECT_EQ(params2.seq_type, SEQ_DNA);
    
    // Test case 3
    CMaple cmaple3("alignment.fasta", "invalid_format", "AA");
    cmaple::Params& params3 = cmaple3.getSettings();
    EXPECT_EQ(params3.aln_path, "alignment.fasta");
    
    // Test case 4
    CMaple cmaple4("alignment.fasta", "FASTA", "invalid_seqtype");
    cmaple::Params& params4 = cmaple4.getSettings();
    EXPECT_EQ(params4.aln_path, "alignment.fasta");
    EXPECT_EQ(params4.aln_format, IN_FASTA);
    
    // Test case 5
    CMaple cmaple5("", "FASTA", "");
    cmaple::Params& params5 = cmaple5.getSettings();
    EXPECT_EQ(params5.aln_path, "");
}


TEST(CMapleTest, SetAlignmentTest) {
    // Test case 1
    CMaple cmaple1;
    std::string aln_filename = "alignment.fasta";
    std::string format = "PHYLIP";
    std::string seqtype = "AA";
    int result1 = cmaple1.setAlignment(aln_filename, format, seqtype);
    cmaple::Params& params = cmaple1.getSettings();
    EXPECT_EQ(params.aln_path, "alignment.fasta");
    EXPECT_EQ(params.aln_format, IN_PHYLIP);
    EXPECT_EQ(params.seq_type, SEQ_PROTEIN);
    EXPECT_EQ(result1, CODE_SUCCESS);
    
    // Test case 2
    CMaple cmaple2;
    int result2 = cmaple2.setAlignment("", "FASTA", "");
    EXPECT_NE(result2, CODE_SUCCESS);
    
    // Test case 3
    CMaple cmaple3;
    int result3 = cmaple3.setAlignment("alignment.fasta", "invalid_format", "");
    EXPECT_NE(result3, CODE_SUCCESS);
    
    // Test case 4
    CMaple cmaple4;
    int result4 = cmaple4.setAlignment("alignment.fasta", "PHYLIP", "invalid_seqtype");
    EXPECT_NE(result4, CODE_SUCCESS);
}

TEST(CMapleTest, SetModelTest) {
    // Test case 1
    CMaple cmaple1;
    // Set up test data and configuration
    std::string model_name = "GTR";
    // Call the method under test
    int result1 = cmaple1.setModel(model_name);
    cmaple::Params& params = cmaple1.getSettings();
    // Check the result
    EXPECT_EQ(result1, CODE_SUCCESS);
    EXPECT_EQ(params.model_name, "GTR");
    
    // Test case 2
    int result2 = cmaple1.setModel("");
    EXPECT_NE(result2, CODE_SUCCESS);
    
    // Test case 3
    int result3 = cmaple1.setModel("JC");
    EXPECT_EQ(result3, CODE_SUCCESS);
    int result4 = cmaple1.setModel("UNREST");
    EXPECT_EQ(result4, CODE_SUCCESS);
    EXPECT_EQ(params.model_name, "UNREST");
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
    //int status = cmaple.runInference(force_rerun, tree_type);

    // Check the result
    //EXPECT_EQ(status, 0);
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
    //int status = cmaple.computeBranchSupports(force_rerun, num_threads, num_replicates, epsilon);

    // Check the result
    //EXPECT_EQ(status, 0);
    // Add additional assertions if needed
}

TEST(CMapleTest, SetInputTreeTest) {
    // Test case 1
    CMaple cmaple1;
    // Set up test data and configuration
    std::string tree_filename = "../../example/test_100.treefile";
    // Call the method under test
    int result1 = cmaple1.setInputTree(tree_filename);
    cmaple::Params& params = cmaple1.getSettings();
    // Check the result
    EXPECT_EQ(result1, CODE_SUCCESS);
    EXPECT_EQ(params.input_treefile, "../../example/test_100.treefile");
    
    // Test case 2
    int result2 = cmaple1.setInputTree("");
    EXPECT_NE(result2, CODE_SUCCESS);
    
    // Test case 3
    int result3 = cmaple1.setInputTree("tree1.nwk");
    EXPECT_EQ(result3, CODE_SUCCESS);
    int result4 = cmaple1.setInputTree("tree2.nwk");
    EXPECT_EQ(result4, CODE_SUCCESS);
    EXPECT_EQ(params.input_treefile, "tree2.nwk");
}

TEST(CMapleTest, ExtractMapleTest) {
    // Test case 1
    CMaple cmaple;
    // Set up test data
    std::string aln_filename = "../../example/test_100.fa";
    std::string output_filename = "output.maple";
    // Call the method under test
    cmaple.overwriteOutputs(true);
    int result1 = cmaple.extractMaple(aln_filename, output_filename);
    // Check the result
    EXPECT_EQ(result1, 0);
    
    // Test case 2
    int result2 = cmaple.extractMaple("", output_filename);
    EXPECT_NE(result2, CODE_SUCCESS);
    
    // Test case 3
    int result3 = cmaple.extractMaple(aln_filename, "");
    EXPECT_NE(result3, CODE_SUCCESS);
    
    // Test case 4
    CMaple cmaple2;
    cmaple2.overwriteOutputs(true);
    EXPECT_DEATH(cmaple2.extractMaple("notfound", output_filename), ".*");
    
    // Test case 5
    CMaple cmaple3;
    // Output file is existing
    EXPECT_DEATH(cmaple3.extractMaple(aln_filename, output_filename), ".*");
}

TEST(CMapleTest, ExtractFASTATest) {
    // Test case 1
    CMaple cmaple;
    // Set up test data
    std::string aln_filename = "../../example/test_100.maple";
    std::string output_filename = "output.fa";
    // Call the method under test
    cmaple.overwriteOutputs(true);
    int result1 = cmaple.extractFASTA(aln_filename, output_filename);
    // Check the result
    EXPECT_EQ(result1, 0);
    
    // Test case 2
    int result2 = cmaple.extractFASTA("", output_filename);
    EXPECT_NE(result2, CODE_SUCCESS);
    
    // Test case 3
    int result3 = cmaple.extractFASTA(aln_filename, "");
    EXPECT_NE(result3, CODE_SUCCESS);
    
    // Test case 4
    CMaple cmaple2;
    cmaple2.overwriteOutputs(true);
    EXPECT_DEATH(cmaple2.extractFASTA("notfound", output_filename), ".*");
    
    // Test case 5
    CMaple cmaple3;
    // Output file is existing
    EXPECT_DEATH(cmaple3.extractFASTA(aln_filename, output_filename), ".*");
}

TEST(CMapleTest, GetTreeAndModelStringTest) {
    // Test case 1
    CMaple cmaple1;
    std::string tree_str = cmaple1.getTreeString();
    std::string model_str = cmaple1.getModelString();
    // Example assertion:
    EXPECT_TRUE(tree_str.empty());
    EXPECT_EQ(model_str, "");
    
    // Test case 2
    CMaple cmaple2;
    const std::string tree_bin = "(((((((30:0,(31:0,(26:0,6:0):0):0):0,((25:0,(7:0,(33:0,(3:0,(21:0,(13:0,(15:0,(4:0,(1:0,(2:0,(12:0,(8:0,(9:0,(19:0,(35:0,(36:0,(27:0,(20:0,(16:0,(32:0,(18:0,(17:0,(34:0,(37:0,(23:0,(38:0,(10:0,(11:0,22:0):0):0):0):0):0):0):0):0):0):0):0):0):0):0):0):0):0):0):0):0):0):0):0):0):0):0):0):2.509116529836319387e-05,14:3.3454885851824656129e-05):8.3637214629561640322e-06):0,28:3.3454885851824656129e-05):0,5:3.3454885851824656129e-05):0,(24:1.6727442925912328064e-05,39:0.00013381954340729862452):1.6727442925912328064e-05):8.4585071817855350673e-06,(0:0,29:0.00052877567941322922707):3.2175408705370500684e-05):1.3757970009464770555e-05,(((42:0,((44:0,(46:0,(45:0,40:0):0):0):0,(43:0,47:0):7.3747803980950266123e-05):7.2997572715394198895e-05):0,41:3.3454885851824656129e-05):7.364175689872354269e-05,(48:5.7780795259532169439e-08,49:3.3454885851824656129e-05):3.3454885851824656129e-05):0.00013381954340729862452):0;";
    const std::string true_model_str2 = "\nMODEL: GTR\n\nROOT FREQUENCIES\nA\t\t\tC\t\t\tG\t\t\tT\t\t\t\n0.29912\t0.183634\t0.196179\t0.321067\t\n\nMUTATION MATRIX\n\tA\t\t\tC\t\t\tG\t\t\tT\t\t\t\nA\t-0.552468\t0.127493\t0.311648\t0.113327\t\nC\t0.207672\t-1.4537\t0.161523\t1.08451\t\nG\t0.475179\t0.151193\t-1.16635\t0.539976\t\nT\t0.10558\t0.620284\t0.329938\t-1.0558\t\n";
    // Set up test data
    cmaple2.setAlignment("../../example/test_50.maple");
    cmaple2.overwriteOutputs(true);
    cmaple2.runInference();
    // Call the method under test
    std::string treeString2 = cmaple2.getTreeString();
    std::string model_str2 = cmaple2.getModelString();
    EXPECT_EQ(treeString2, tree_bin);
    EXPECT_EQ(model_str2, true_model_str2);
    
    // Test case 3
    const std::string tree_mul = "(((((((30:0,31:0,26:0,6:0):0,((25:0,7:0,33:0,3:0,21:0,13:0,15:0,4:0,1:0,2:0,12:0,8:0,9:0,19:0,35:0,36:0,27:0,20:0,16:0,32:0,18:0,17:0,34:0,37:0,23:0,38:0,10:0,11:0,22:0):2.509116529836319387e-05,14:3.3454885851824656129e-05):8.3637214629561640322e-06):0,28:3.3454885851824656129e-05):0,5:3.3454885851824656129e-05):0,(24:1.6727442925912328064e-05,39:0.00013381954340729862452):1.6727442925912328064e-05):8.4585071817855350673e-06,(0:0,29:0.00052877567941322922707):3.2175408705370500684e-05):1.3757970009464770555e-05,(((42:0,((44:0,46:0,45:0,40:0):0,(43:0,47:0):7.3747803980950266123e-05):7.2997572715394198895e-05):0,41:3.3454885851824656129e-05):7.364175689872354269e-05,(48:5.7780795259532169439e-08,49:3.3454885851824656129e-05):3.3454885851824656129e-05):0.00013381954340729862452):0;";
    std::string treeString3 = cmaple2.getTreeString("MUL");
    EXPECT_EQ(treeString3, tree_mul);
    
    // Test case 4
    std::string treeString4 = cmaple2.getTreeString("MUL", true);
    EXPECT_EQ(treeString4, "");
    
    // Test case 5
    CMaple cmaple3;
    EXPECT_DEATH(cmaple3.getTreeString("INVALID"), ".*");
    
    // Test case 6
    const std::string true_model_str6 = "\nMODEL: GTR\n\nROOT FREQUENCIES\nA\t\t\tC\t\t\tG\t\t\tT\t\t\t\n0.299622\t0.183634\t0.196046\t0.320699\t\n\nMUTATION MATRIX\n\tA\t\t\tC\t\t\tG\t\t\tT\t\t\t\nA\t-0.541994\t0.128367\t0.299523\t0.114104\t\nC\t0.209447\t-1.46613\t0.139631\t1.11705\t\nG\t0.457769\t0.130791\t-1.13352\t0.544963\t\nT\t0.106605\t0.639629\t0.33314\t-1.07937\t\n";
    CMaple cmaple4;
    cmaple4.setAlignment("../../example/test_50.fa");
    cmaple4.overwriteOutputs(true);
    cmaple4.computeBranchSupports(false, 4);
    // Call the method under test
    std::string treeString6 = cmaple4.getTreeString("BIN", true);
    EXPECT_TRUE(treeString6.length() > 0);
    std::string model_str6 = cmaple4.getModelString();
    EXPECT_EQ(model_str6, true_model_str6);
}

TEST(CMapleTest, GetVersionTest) {
    CMaple cmaple;

    // Call the method under test
    std::string version = cmaple.getVersion();

    // Check the result
    EXPECT_EQ(version, "CMAPLE " + convertIntToString(cmaple_VERSION_MAJOR) + "." + convertIntToString(cmaple_VERSION_MINOR) + cmaple_VERSION_PATCH);
}

TEST(CMapleTest, GetCitationsTest) {
    CMaple cmaple;

    // Call the method under test
    std::string citations = cmaple.getCitations();

    // Check the result
    EXPECT_TRUE(citations.length() > 0);
}

TEST(CMapleTest, SetRefTest) {
    // Test case 1
    CMaple cmaple1;
    // Set up test data and configuration
    std::string ref_filename = "../../example/ref.fa";
    // Call the method under test
    int result1 = cmaple1.setRef(ref_filename);
    cmaple::Params& params = cmaple1.getSettings();
    // Check the result
    EXPECT_EQ(result1, CODE_SUCCESS);
    EXPECT_EQ(params.ref_path, "../../example/ref.fa");
    
    // Test case 2
    int result2 = cmaple1.setRef("");
    EXPECT_NE(result2, CODE_SUCCESS);
    
    // Test case 3
    int result3 = cmaple1.setRef("ref1.fa");
    EXPECT_EQ(result3, CODE_SUCCESS);
    int result4 = cmaple1.setRef("ref2.fa");
    EXPECT_EQ(result4, CODE_SUCCESS);
    EXPECT_EQ(params.ref_path, "ref2.fa");
}

TEST(CMapleTest, SetTreeSearchTypeTest) {
    // Test case 1
    CMaple cmaple;
    cmaple::Params& params = cmaple.getSettings();
    EXPECT_EQ(params.tree_search_type, PARTIAL_TREE_SEARCH);
    // Set up test data
    std::string tree_search_type = "FAST";
    // Call the method under test
    int result1 = cmaple.setTreeSearchType(tree_search_type);
    // Check the result
    EXPECT_EQ(result1, 0);
    EXPECT_EQ(params.tree_search_type, NO_TREE_SEARCH);
    
    // Test case 2
    int result2 = cmaple.setTreeSearchType("NORMAL");
    // Check the result
    EXPECT_EQ(result2, 0);
    EXPECT_EQ(params.tree_search_type, PARTIAL_TREE_SEARCH);
    
    // Test case 3
    int result3 = cmaple.setTreeSearchType("slow");
    // Check the result
    EXPECT_EQ(result3, 0);
    EXPECT_EQ(params.tree_search_type, COMPLETE_TREE_SEARCH);
    
    // Test case 4
    int result4 = cmaple.setTreeSearchType("");
    // Check the result
    EXPECT_NE(result4, 0);
    
    // Test case 5
    int result5 = cmaple.setTreeSearchType("invalid");
    // Check the result
    EXPECT_NE(result5, 0);
}

TEST(CMapleTest, SetShortRangeTreeSearchTest) {
    CMaple cmaple;

    // Set up test data
    bool enable = true;

    // Call the method under test
    int status = cmaple.setShortRangeTreeSearch(enable);

    // Check the result
    EXPECT_EQ(status, 0);
    // Add additional assertions if needed
}

TEST(CMapleTest, SetPrefixTest) {
    CMaple cmaple;

    // Set up test data
    std::string prefix = "output_";

    // Call the method under test
    int status = cmaple.setPrefix(prefix);

    // Check the result
    EXPECT_EQ(status, 0);
    // Add additional assertions if needed
}

TEST(CMapleTest, EnableOptimizingBlengthsTest) {
    CMaple cmaple;

    // Set up test data
    bool enable = true;

    // Call the method under test
    int status = cmaple.enableOptimizingBlengths(enable);

    // Check the result
    EXPECT_EQ(status, 0);
    // Add additional assertions if needed
}

TEST(CMapleTest, SetMinBlengthTest) {
    CMaple cmaple;

    // Set up test data
    double min_blength = 1e-9;

    // Call the method under test
    int status = cmaple.setMinBlength(min_blength);

    // Check the result
    EXPECT_EQ(status, 0);
    // Add additional assertions if needed
}

TEST(CMapleTest, SetThreshProbTest) {
    CMaple cmaple;

    // Set up test data
    double thresh_prob = 1e-6;

    // Call the method under test
    int status = cmaple.setThreshProb(thresh_prob);

    // Check the result
    EXPECT_EQ(status, 0);
    // Add additional assertions if needed
}

TEST(CMapleTest, OverwriteOutputsTest) {
    CMaple cmaple;

    // Set up test data
    bool enable = true;

    // Call the method under test
    int status = cmaple.overwriteOutputs(enable);

    // Check the result
    EXPECT_EQ(status, 0);
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
    EXPECT_EQ(status, 0);
    // Add additional assertions if needed
}

TEST(CMapleTest, GetSettingsTest) {
    CMaple cmaple;

    // Call the method under test
    cmaple::Params& settings = cmaple.getSettings();

    // Check the result
    // Add assertions to validate the settings if needed
}
