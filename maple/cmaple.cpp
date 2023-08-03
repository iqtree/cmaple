#include "cmaple.h"
using namespace std;
using namespace cmaple;

void cmaple::runCMaple(cmaple::Params &params)
{
    try
    {
        // record the start time
        auto start = getRealTime();
        
        // Initialize output filename -> use aln_ as the output prefix if users didn't specify it
        const std::string prefix = (params.output_prefix.length() ? params.output_prefix :  params.aln_path);
        const std::string output_treefile = prefix + ".treefile";
        // check whether output file is already exists
        if (!params.overwrite_output && fileExists(output_treefile))
            outError("File " + output_treefile + " already exists. Use `--overwrite` option if you really want to overwrite it.\n");
        
        // Initialize a Model
        Model model(params.sub_model, params.seq_type);
        
        // Initializa an Alignment
        // Retrieve the reference genome (if specified) from an alignment -> this feature has not yet exposed to APIs -> should be refactoring later
        const std::string ref_seq = "";
        if (params.ref_path.length() && params.ref_seqname.length())
        {
            Alignment aln_tmp;
            aln_tmp.readRefSeq(params.ref_path, params.ref_seqname);
        }
        Alignment aln(params.aln_path, ref_seq, params.aln_format, params.seq_type);
        
        // If users only want to convert the alignment to another format -> convert it and terminate
        if (params.output_aln.length() && params.output_aln_format != IN_UNKNOWN)
        {
            if (cmaple::verbose_mode > cmaple::VB_QUIET)
                std::cout << "Write the alignment to " + params.output_aln << std::endl;
            aln.write(params.output_aln, params.output_aln_format, params.overwrite_output);
            return;
        }
        
        // Initialize a Tree
        Tree tree(&aln, &model, params.input_treefile, params.fixed_blengths);
        // clone all settings
        cmaple::Params& tree_params = tree.getParams();
        tree_params = params;
        
        // Infer a phylogenetic tree
        std::string redirected_msgs = tree.doInference(params.tree_search_type, params.shallow_tree_search);
        if (cmaple::verbose_mode >= cmaple::VB_MED)
            std::cout << redirected_msgs << std::endl;
        
        // Compute branch supports (if users want to do so)
        if (params.compute_aLRT_SH)
        {
            // if users don't input a tree file, always allow CMaple to replace the ML tree by a higher-likelihood tree (if found)
            bool allow_replacing_ML_tree = true;
            // if users input a tree -> depend on the setting in params (false ~ don't allow replacing (by default)
            if (params.input_treefile.length())
                allow_replacing_ML_tree = params.allow_replace_input_tree;
            
            redirected_msgs = tree.computeBranchSupport(params.num_threads, params.aLRT_SH_replicates, params.aLRT_SH_half_epsilon + params.aLRT_SH_half_epsilon, allow_replacing_ML_tree);
            if (cmaple::verbose_mode >= cmaple::VB_MED)
                std::cout << redirected_msgs << std::endl;
            
            // write the tree file with branch supports
            ofstream out_tree_branch_supports = ofstream(prefix + ".aLRT_SH.treefile");
            out_tree_branch_supports << tree.exportNewick(params.tree_format, true);
            out_tree_branch_supports.close();
        }
        
        // output log-likelihood of the tree
        if (cmaple::verbose_mode > cmaple::VB_QUIET)
            std::cout << std::setprecision(10) << "Tree log likelihood: " << tree.computeLh() << std::endl;
        
        // Write the normal tree file
        ofstream out = ofstream(output_treefile);
        out << tree.exportNewick(params.tree_format);
        out.close();
            
        // Show model parameters
        if (cmaple::verbose_mode > cmaple::VB_QUIET)
        {
            cmaple::ModelParams model_params = model.getParams();
            std::cout << "\nMODEL: " + model_params.model_name + "\n";
            std::cout << "\nROOT FREQUENCIES\n";
            std::cout << model_params.state_freqs;
            std::cout << "\nMUTATION MATRIX\n";
            std::cout << model_params.mut_rates << std::endl;
        }
            
        // Show information about output files
        std::cout << "Analysis results written to:" << std::endl;
        std::cout << "Maximum-likelihood tree:       " << output_treefile << std::endl;
        if (params.compute_aLRT_SH)
            std::cout << "Tree with aLRT-SH values:      " << prefix + ".aLRT_SH.treefile" << std::endl;
        std::cout << "Screen log file:               " << prefix + ".log" << std::endl << std::endl;
        
        // show runtime
        auto end = getRealTime();
        if (cmaple::verbose_mode > cmaple::VB_QUIET)
            cout << "Runtime: " << end - start << "s" << endl;
    }
    catch (std::invalid_argument e)
    {
        outError(e.what());
    }
    catch (std::logic_error e)
    {
        outError(e.what());
    }
    catch (ios::failure e)
    {
        outError(e.what());
    }
}

std::string cmaple::getVersion()
{
    return "CMAPLE " + convertIntToString(cmaple_VERSION_MAJOR) + "." + convertIntToString(cmaple_VERSION_MINOR) + cmaple_VERSION_PATCH;
}

std::string cmaple::getCitations()
{
    return "[Citations]";
}

void cmaple::testing(cmaple::Params& params)
{
    cmaple::verbose_mode = VB_DEBUG;
    
    // -------- Test write and read tree -----------
    Alignment aln10("test_100.maple");
    Model model10(GTR);
    Tree tree10(&aln10, &model10);
    std::cout << tree10.doInference() << std::endl;
    std::stringstream tree10_stream(tree10.exportNewick(MUL_TREE));
    //tree10_stream << tree10;
    
    Tree tree10_1(&aln10, &model10, tree10_stream);
    std::stringstream tree10_1_stream(tree10.exportNewick(MUL_TREE));
    //tree10_1_stream << tree10_1;
    
    std::cout << tree10_stream.str() << std::endl;
    std::cout << tree10_1_stream.str() << std::endl;
    string str1 = tree10_stream.str();
    string str2 = tree10_1_stream.str();
    for (auto i = 0; i < str1.length(); ++i)
        if (str1[i] != str2[i])
            std::cout << i << std::endl;
    if (tree10_stream.str() != tree10_1_stream.str())
        std::cout << "Mismatch tree strings" << std:: endl;
    // -------- Test write and read tree -----------
    
    // Initialize a Model
    Model model(params.sub_model, params.seq_type);
    
    Alignment aln_empty;
    // Initialize a Tree
    //Tree tree_empty(&aln_empty, &model);
    
    // Initialize an Alignment
    Alignment aln(params.aln_path, "", params.aln_format, params.seq_type);
    
    // Initialize a Tree
    Tree tree(&aln, &model);
    
    // Infer a phylogenetic tree
    std::string redirected_msgs = tree.doInference(params.tree_search_type, params.shallow_tree_search);
    if (cmaple::verbose_mode >= cmaple::VB_MED)
        std::cout << redirected_msgs << std::endl;
    
    // Compute branch supports (if users want to do so)
    if (params.compute_aLRT_SH)
    {
        redirected_msgs = tree.computeBranchSupport(params.num_threads, params.aLRT_SH_replicates, params.aLRT_SH_half_epsilon + params.aLRT_SH_half_epsilon, true);
        if (cmaple::verbose_mode >= cmaple::VB_MED)
            std::cout << redirected_msgs << std::endl;
    }
    
    std::cout << tree << std::endl;
    
    // -------- Test changeAln -----------
    // tree.changeAln(&aln_empty); // PASSED
    
    /*Alignment aln_100("test_100.maple");
    // tree.changeAln(&aln_100); // PASSED some existing leaves are not found in the new (smaller) alignment
    
    std::stringstream tree_str_1;
    tree_str_1 << tree;
    
    Alignment aln_500("test_500.maple");
    tree.changeAln(&aln_500);
    std::stringstream tree_str_2;
    tree_str_2 << tree;
    if (tree_str_1.str() != tree_str_2.str())
        std::cout << "Mismatch tree strings" << std:: endl;
    Model model2(JC);
    tree.changeModel(&model2);
    tree.doInference();
    tree.computeLh();
    tree.computeBranchSupport(4);*/
    
    // -------- Test changeAln -----------
    
    // -------- Test re-read alignment -----------
    std::stringstream tree_str_1;
    tree_str_1 << tree.exportNewick(MUL_TREE);
    
    Model model3(GTR);
    Tree tree2(&aln,&model3, tree_str_1, true);
    std::stringstream tree_str_2;
    tree_str_2 << tree2.exportNewick(MUL_TREE);
    if (tree_str_1.str() != tree_str_2.str())
        std::cout << "Mismatch tree strings" << std:: endl;
    
    str1 = tree_str_1.str();
    str2 = tree_str_2.str();
    std::cout << str1 << std::endl;
    std::cout << str2 << std::endl;
    for (auto i = 0; i < str1.length(); ++i)
        if (str1[i] != str2[i])
            std::cout << i << std::endl;
    
    tree_str_1 >> tree2;
    std::stringstream tree_str_2_6;
    tree_str_2_6 << tree2.exportNewick(MUL_TREE);
    if (tree_str_2_6.str() != tree_str_2.str())
        std::cout << "Mismatch tree strings" << std:: endl;
    
    aln.read("test_200.maple");
    std::cout << tree2.doInference(FAST_TREE_SEARCH) << std::endl;
    std::stringstream tree_str_2_5;
    tree_str_2_5 << tree2.exportNewick(MUL_TREE);
    if (tree_str_2_5.str() != tree_str_2.str())
        std::cout << "Mismatch tree strings" << std:: endl;
    
    // aln.read("test_100.maple"); // PASSED - unable to change to a smaller alignment
    aln.read("test_500.maple");
    Model model2(JC);
    /*
    std::stringstream tree_str_2;
    tree.changeModel(&model2);
    tree_str_2 << tree;
    if (tree_str_1.str() != tree_str_2.str())
        std::cout << "Mismatch tree strings" << std:: endl;*/
    
    std::cout << tree.computeLh() << std::endl;
    std::stringstream tree_str_3;
    tree_str_3 << tree;
    if (tree_str_1.str() != tree_str_3.str())
        std::cout << "Mismatch tree strings" << std:: endl;
    
    tree_str_2_5 >> tree;
    std::cout << tree.computeLh() << std::endl;
    
    std::stringstream tree_str_4;
    std::cout << tree.doInference() << std::endl;
    tree_str_4 << tree;
    if (tree_str_1.str() != tree_str_4.str())
        std::cout << "Mismatch tree strings" << std:: endl;
    tree.computeLh();
    
    std::cout << tree.computeLh() << std::endl;
    tree_str_2_5 >> tree;
    std::cout << tree.computeLh() << std::endl;
    
    /*std::stringstream tree_str_2_2;
    tree_str_2_2 << tree2;
    if (tree_str_1.str() != tree_str_2_2.str())
        std::cout << "Mismatch tree strings" << std:: endl;
    std::cout << tree2.doInference(FAST_TREE_SEARCH, false) << std::endl;
    std::stringstream tree_str_2_3;
    tree_str_2_3 << tree2;
    if (tree_str_1.str() != tree_str_2_3.str())
        std::cout << "Mismatch tree strings" << std:: endl;*/
    
    
    tree.computeBranchSupport(4);
    
    // -------- Test re-read alignment -----------
    
    std::cout << "DONE" << std::endl;
    
    /*CMaple cmaple;
    cmaple::Params& cmaple_params = cmaple.getSettings();
    cmaple_params = params;
    cmaple.inferTree();*/
    /*Model model(params.sub_model);
    // with an alignment file
    Alignment aln1(params.aln_path);
    // with an alignment file
    // Alignment aln2(""); // tested PASS
    // Alignment aln3("notfound"); // tested PASS
    // with a stream
    const std::string aln_filename = params.aln_path;
    std::ifstream aln_stream;
    try {
        aln_stream.exceptions(ios::failbit | ios::badbit);
        aln_stream.open(aln_filename);
    } catch (ios::failure) {
        outError(ERR_READ_INPUT, aln_filename);
    }
    Alignment aln(aln_stream);
    aln_stream.close();
    // Write alignment to file in MAPLE format
    aln.write("output.maple", "MAPLE", true);
    // write aln to a file in FASTA format
    aln.write("output.fa", "FASTA", true);
    aln.write("output.phy", "PHYLIP", true);
    // Read ref_seq from an alignment file (not yet exposed to APIs)
    ASSERT(params.ref_path.length() && params.ref_seqname.length());
     Alignment aln_base;
    std::string ref_seq = aln_base.readRefSeq(params.ref_path, params.ref_seqname);
    Alignment aln2("output.fa", ref_seq);
    aln2.write("output1.maple", "MAPLE", true);
    Alignment aln3("output.fa");
    aln3.write("output2.maple", "MAPLE", true);
    Alignment aln4("output.phy", ref_seq);
    aln4.write("output_phy1.maple", "MAPLE", true);
    Alignment aln5("output.phy");
    aln5.write("output_phy2.maple", "MAPLE", true);
    Alignment aln6("output2.maple");
    aln4.write("output3.fa", "FASTA", true);
    // aln5.write("output_phy3.maple", "INVALID", true);
    // aln.write("output.maple");
    // without tree file
    Tree tree1(aln, model);
    // with a tree file
    const std::string tree_filename = params.input_treefile;
    Tree tree2(aln, model, tree_filename, true);
    // Test keeping blengths fixed when inputting an empty tree
    Tree tree3(aln, model, "", true);
    // Test keeping blengths fixed when inputting a tree topology without branch lengths
    //Tree tree4(aln, model, "topo.treefile", true);
    // Test keeing blengths fixed (successfully)
    Tree tree5(aln, model, "test_200_5.diff.treefile", true);
    std::cout << tree5.doInference("Normal", true) << std::endl;
    
    // with a tree stream
    std::ifstream tree_stream;
    try {
        tree_stream.exceptions(ios::failbit | ios::badbit);
        tree_stream.open(tree_filename);
    } catch (ios::failure) {
        outError(ERR_READ_INPUT, tree_filename);
    }
    Tree tree(aln, model, tree_stream);
    tree_stream.close();
    std::cout << tree.doInference("FAST", true) << std::endl;
    
    std::cout << tree.doInference("MORE_accurate") << std::endl;

    // std::cout << tree.exportNewick("BIN", true) << std::endl;
    //Tree tree(aln, model, "");
    std::cout << "Tree likelihood: " << tree.computeLh() << std::endl;
    std::cout << tree.computeBranchSupport(8, 100, 0.1, false) << std::endl;
    //std::cout << tree.computeBranchSupport(8, 100, 0.1, true) << std::endl;
    std::cout << tree.exportNewick("BIN", true) << std::endl;
    tree.doInference();
    std::cout << tree.exportNewick() << std::endl;
    tree.computeBranchSupport(8, 100);
    std::cout << tree.exportNewick("BIN", true) << std::endl;
    std::cout << "Tree likelihood: " << tree.computeLh() << std::endl;*/

    // Create an alignment from a file
    /*cmaple::Alignment aln(params.aln_path);

    // Create a model
    cmaple::Model model("GTR");

    // Create a tree, attach the alignment and model to the tree
    cmaple::Tree tree(aln, model);

    // Infer a phylogenetic tree from the alignment and the model
    cout << tree.doInference() << endl;

    // Compute the branch supports for the inferred tree
    cout << tree.computeBranchSupport() << endl;

    // Compute the likelihood of the tree
    cout << "- Tree log likelihood: " << tree.computeLh() << endl;

    // Export the tree (with branch supports) in NEWICK format
    cout << "- Tree: " << tree.exportNewick("BIN", true) << endl;*/
    
    // Create an alignment from a file
    /*cmaple::Alignment aln(params.aln_path);

    // Create a model
    cmaple::Model model(GTR);

    // Create a tree, attach the alignment and model to the tree
    cmaple::Tree tree(&aln, &model);
    
    // change aln
    cmaple::Alignment aln1("test_100.maple");
    tree.changeAln(&aln1);
    
    // infer
    cout << tree.doInference()<< endl;
    
    cout << tree.exportNewick() << endl;
    
    // change model
    cmaple::Model model1(JC);
    tree.changeModel(&model1);
    
    // change model again
    cmaple::Model model2(GTR);
    tree.changeModel(&model2);
    
    // change aln
    cmaple::Alignment aln2("test_200.maple");
    tree.changeAln(&aln2);
    
    // infer again
    cout << tree.doInference()<< endl;
    
    cout << tree.exportNewick() << endl;
    
    // change aln
    // cmaple::Alignment aln3("test_100.maple");
    // tree.changeAln(aln3);
    
    // change aln with different seqtypes
    Alignment aln4("test_aa/test_aa_100K.maple");
    tree.changeAln(&aln4);*/
}
