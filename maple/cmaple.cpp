#include "cmaple.h"
using namespace std;
using namespace cmaple;

void cmaple::runCMaple(cmaple::Params &params)
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
    Model model(params.model_name);
    
    // Initializa an Alignment
    // Retrieve the reference genome (if specified) from an alignment -> this feature has not yet exposed to APIs -> should be refactoring later
    const std::string ref_seq = "";
    if (params.ref_path.length() && params.ref_seqname.length())
    {
        AlignmentBase aln_tmp;
        aln_tmp.readRefSeq(params.ref_path, params.ref_seqname);
    }
    Alignment aln(params.aln_path, ref_seq, params.aln_format_str, params.seq_type_str);
    
    // If users only want to convert the alignment to another format -> convert it and terminate
    if (params.output_aln.length() && params.output_aln_format.length())
    {
        std::cout << "Write the alignment to " + params.output_aln + " in " + params.output_aln_format + " format." << std::endl;
        aln.write(params.output_aln, params.output_aln_format, params.overwrite_output);
        return;
    }
    
    // Initialize a Tree
    Tree tree(aln, model, params.input_treefile, params.fixed_blengths);
    // clone all settings
    cmaple::Params& tree_params = tree.getParams();
    tree_params = params;
    
    // Infer a phylogenetic tree
    std::cout << tree.infer(params.tree_search_type_str, params.shallow_tree_search) << std::endl;
    
    // Compute branch supports (if users want to do so)
    if (params.compute_aLRT_SH)
    {
        // if users don't input a tree file, always allow CMaple to replace the ML tree by a higher-likelihood tree (if found)
        bool allow_replacing_ML_tree = true;
        // if users input a tree -> depend on the setting in params (false ~ don't allow replacing (by default)
        if (params.input_treefile.length())
            allow_replacing_ML_tree = params.allow_replace_input_tree;
        std::cout << tree.computeBranchSupports(params.num_threads, params.aLRT_SH_replicates, params.aLRT_SH_half_epsilon + params.aLRT_SH_half_epsilon, allow_replacing_ML_tree) << std::endl;
        
        // write the tree file with branch supports
        ofstream out_tree_branch_supports = ofstream(prefix + ".aLRT_SH.treefile");
        out_tree_branch_supports << tree.exportString(params.tree_format, true);
        out_tree_branch_supports.close();
    }
    
    // Write the normal tree file
    ofstream out = ofstream(output_treefile);
    out << tree.exportString(params.tree_format);
    out.close();
        
    // Show model parameters
    std::map<std::string, std::string> model_params = model.getParams();
    std::cout << "\nMODEL: " + model_params[cmaple::MODEL_NAME] + "\n";
    std::cout << "\nROOT FREQUENCIES\n";
    std::cout << model_params[cmaple::MODEL_FREQS];
    std::cout << "\nMUTATION MATRIX\n";
    std::cout << model_params[cmaple::MODEL_RATES] << std::endl;
        
    // Show information about output files
    std::cout << "Analysis results written to:" << std::endl;
    std::cout << "Maximum-likelihood tree:       " << output_treefile << std::endl;
    if (params.compute_aLRT_SH)
        std::cout << "Tree with aLRT-SH values:      " << prefix + ".aLRT_SH.treefile" << std::endl;
    std::cout << "Screen log file:               " << prefix + ".log" << std::endl << std::endl;
    
    // show runtime
    auto end = getRealTime();
    cout << "Runtime: " << end - start << "s" << endl;
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
    /*CMaple cmaple;
    cmaple::Params& cmaple_params = cmaple.getSettings();
    cmaple_params = params;
    cmaple.inferTree();*/
    Model model(params.model_name);
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
    AlignmentBase aln_base;
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
    std::cout << tree5.infer("Normal", true) << std::endl;
    
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
    std::cout << tree.infer("FAST", true) << std::endl;
    
    std::cout << tree.infer("MORE_accurate") << std::endl;

    // std::cout << tree.exportString("BIN", true) << std::endl;
    //Tree tree(aln, model, "");
    std::cout << "Tree likelihood: " << tree.computeLh() << std::endl;
    std::cout << tree.computeBranchSupports(8, 100, 0.1, false) << std::endl;
    //std::cout << tree.computeBranchSupports(8, 100, 0.1, true) << std::endl;
    std::cout << tree.exportString("BIN", true) << std::endl;
    tree.infer();
    std::cout << tree.exportString() << std::endl;
    tree.computeBranchSupports(8, 100);
    std::cout << tree.exportString("BIN", true) << std::endl;
    std::cout << "Tree likelihood: " << tree.computeLh() << std::endl;
}
