#include "cmaple.h"
using namespace std;
using namespace cmaple;

CMaple::CMaple():tree(cmaple::Params()),status(NEW_INSTANCE) {};

CMaple::CMaple(cmaple::Params&& params):tree(std::move(params)),status(NEW_INSTANCE) {};

CMaple::CMaple(const std::string& aln_filename, const std::string& format, const std::string& seqtype):tree(cmaple::Params()),status(NEW_INSTANCE)
{
    setAlignment(aln_filename, format, seqtype);
}

void CMaple::resetStatus()
{
    // reset status
    status = NEW_INSTANCE;
    
    // reset tree
    tree.resetTree();
}

int CMaple::setAlignment(const std::string& aln_filename, const std::string& format, const std::string& seqtype)
{
    // set aln_filename
    if (aln_filename.length())
        tree.params->aln_path = aln_filename;
    else
        return CODE_ERROR_1;
    
    // set aln format (if specified)
    if (format.length())
    {
        tree.params->aln_format = tree.aln->getAlignmentFormat(format);
        if (tree.params->aln_format == IN_UNKNOWN)
        {
            outError("Unsupported alignment format " + format + ". Please use MAPLE, FASTA, or PHYLIP");
            return CODE_ERROR_1;
        }
    }
    
    // set sequence type
    if (seqtype.length())
    {
        tree.params->seq_type = tree.aln->getSeqType(seqtype);
        if (tree.params->seq_type == SEQ_UNKNOWN)
        {
            outError("Unknown sequence type " + seqtype + ", please use DNA or AA");
            return CODE_ERROR_1;
        }
    }
    
    // success
    return CODE_SUCCESS;
}

int CMaple::setModel(const std::string& model_name)
{
    // set model (if specified)
    if (model_name.length())
    {
        tree.params->model_name = model_name;
        return CODE_SUCCESS;
    }
    
    return CODE_ERROR_1;
}

int CMaple::runInference(const bool force_rerun, const std::string& tree_type)
{
    // If this instance is not NEW -> already run inference/computing branch supports
    if (status != NEW_INSTANCE)
    {
        if (force_rerun)
        {
            std::cout << "Rerun the inference" << std::endl;
            resetStatus();
        }
        else
        {
            std::cout << "Inference has been done. Use runInference(TRUE, <tree_type>) if you want to rerun the inference!" << std::endl;
            return CODE_ERROR_1;
        }
    }
    
    // validate inputs
    if ((!tree.params->maple_path.length()) && (!tree.params->aln_path.length()))
        outError("Please specify an alignment file via setAlignment(...)");
    
    // update status
    status = INFERENCE_DONE;
    
    // Show the random seed number
    cout << "Seed:    " << tree.params->ran_seed <<  " " << std::endl;
    
    // Show the number of threads
    setNumThreads(tree.params->num_threads);

    // Run the inference
    auto start = getRealTime();
    
    // load input data
    loadInput();
    
    // terminate if the user only wants to export a MAPLE file from an alignment
    // or only want to reconstruct an aln from a MAPLE file
    if (tree.params->only_extract_maple || tree.params->output_aln)
        return CODE_SUCCESS;
    
    // prepare for the inference
    preInference();
    
    // infer trees and model params
    doInference();
    
    // complete remaining stuff after the inference
    postInference();
    
    // show runtime
    auto end = getRealTime();
    cout << "Runtime: " << end - start << "s" << endl;
    
    return CODE_SUCCESS;
}

int CMaple::computeBranchSupports(const bool force_rerun, const int num_threads, const int num_replicates, const double epsilon)
{
    // validate inputs
    if (num_threads < 0)
    {
        std::cout << "Number of threads cannot be negative!" << std::endl;
        return CODE_ERROR_1;
    }
    if (num_replicates <= 0)
    {
        std::cout << "Number of replicates must be positive!" << std::endl;
        return CODE_ERROR_1;
    }
    if (epsilon < 0)
    {
        std::cout << "Epsilon cannot be negative!" << std::endl;
        return CODE_ERROR_1;
    }
    
    // If the branch supports have already been computed -> terminate with an error or recompute them (if users want to do so)
    if (status == BRANCH_SUPPORT_DONE)
    {
        // Terminate with an error
        if (!force_rerun)
        {
            std::cout << "Branch supports have already been computed. Use computeBranchSupports(TRUE, <...>) if you want to recompute them!" << std::endl;
            return CODE_ERROR_1;
        }
        // Recompute branch supports (if users want to do so) -> only need to reset the status (delete all objects (e.g., tree), except params)
        else
        {
            std::cout << "Recompute branch supports" << std::endl;
            resetStatus();
        }
    }
    
    // We now need to compute branch supports
    // backup the current option of compute_aLRT_SH
    const bool compute_aLRT_SH = tree.params->compute_aLRT_SH;
    
    // turn on the flag to compute branch supports
    tree.params->compute_aLRT_SH = true;
    
    // Run the entier inference process (if the inference has yet run)
    if (status == NEW_INSTANCE)
        runInference();
    // Otherwise, only need to compute the branch supports (if the inference has already been run)
    else if (status == INFERENCE_DONE)
    {
        // Show the random seed number
        cout << "Seed:    " << tree.params->ran_seed <<  " " << std::endl;
        
        // Show the number of threads
        const int act_num_threads = num_threads == 1 ? tree.params->num_threads : num_threads; // if users specify num_threads when calling this function -> use that value. Otherwise, use the setting in params instance.
        setNumThreads(act_num_threads);
        
        // Only compute the branch supports
        postInference();
    }
        
    // restore compute_aLRT_SH
    tree.params->compute_aLRT_SH = compute_aLRT_SH;
    
    return CODE_SUCCESS;
}

int CMaple::setInputTree(const std::string& tree_filename)
{
    // set input tree file (if specified)
    if (tree_filename.length())
    {
        tree.params->input_treefile = tree_filename;
        return CODE_SUCCESS;
    }
    
    return CODE_ERROR_1;
}

int CMaple::extractMaple(const std::string& aln_filename, const std::string& output_filename)
{
    if (aln_filename.length() && output_filename.length())
    {
        tree.aln->extractMapleFile(aln_filename, output_filename, *tree.params);
        return CODE_SUCCESS;
    }
    
    return CODE_ERROR_1;
}

int CMaple::extractFASTA(const std::string& aln_filename, const std::string& output_filename)
{
    if (aln_filename.length() && output_filename.length())
    {
        tree.aln->reconstructAln(aln_filename, output_filename);
        return CODE_SUCCESS;
    }
    
    return CODE_ERROR_1;
}

std::string CMaple::getTreeString(const std::string& tree_type, const bool show_branch_supports)
{
    return tree.exportTreeString(tree_type, show_branch_supports);
}

std::string CMaple::getModelString()
{
    return tree.exportModelString();
}

std::string CMaple::getVersion()
{
    return "CMAPLE " + convertIntToString(cmaple_VERSION_MAJOR) + "." + convertIntToString(cmaple_VERSION_MINOR) + cmaple_VERSION_PATCH;
}

std::string CMaple::getCitations()
{
    return "[Citations]";
}

int CMaple::setRef(const std::string& ref_filename)
{
    // set ref_filename(if specified)
    if (ref_filename.length())
    {
        tree.params->ref_path = ref_filename;
        return CODE_SUCCESS;
    }
    
    return CODE_ERROR_1;
}

int CMaple::setTreeSearchType(const std::string& tree_search_type)
{
    // set tree_search_type (if specified)
    if (tree_search_type.length())
    {
        tree.params->tree_search_type = parseTreeSearchType(tree_search_type);
        if (tree.params->tree_search_type == UNKNOWN_TREE_SEARCH)
        {
            outError("Unknown tree search type " + tree_search_type + ". Please use <FAST|NORMAL|SLOW>");
            return CODE_ERROR_1;
        }
        return CODE_SUCCESS;
    }
    
    return CODE_ERROR_1;
}

int CMaple::setShortRangeTreeSearch(const bool enable)
{
    tree.params->short_range_topo_search = enable;
    return CODE_SUCCESS;
}

int CMaple::setPrefix(const std::string& prefix)
{
    // set prefix (if specified)
    if (prefix.length())
    {
        tree.params->output_prefix = prefix;
        return CODE_SUCCESS;
    }
    
    return CODE_ERROR_1;
}

int CMaple::enableOptimizingBlengths(const bool enable)
{
    tree.params->optimize_blength = enable;
    return CODE_SUCCESS;
}

int CMaple::setMinBlength(const double min_blength)
{
    // set min_blength (if specified)
    if (min_blength > 0)
    {
        tree.params->fixed_min_blength = min_blength;
        return CODE_SUCCESS;
    }
    else
        outError("min_blength must be positive!");
    
    return CODE_ERROR_1;
}

int CMaple::setThreshProb(const double thresh_prob)
{
    // set thresh_prob (if specified)
    if (thresh_prob > 0)
    {
        tree.params->threshold_prob = thresh_prob;
        return CODE_SUCCESS;
    }
    else
        outError("thresh_prob must be positive!");
    
    return CODE_ERROR_1;
}

int CMaple::overwriteOutputs(const bool enable)
{
    tree.params->overwrite_output = enable;
    return CODE_SUCCESS;
}

int CMaple::setRandomSeed(const int seed)
{
    tree.params->ran_seed = seed;
    return CODE_SUCCESS;
}

cmaple::Params& CMaple::getSettings()
{
    return *tree.params;
}


void CMaple::loadInput()
{
    cmaple::Params& params = *tree.params;
    ASSERT(params.aln_path.length());
    
    // Synchronize seq_type
    tree.aln->setSeqType(params.seq_type);
    
    // detect alignment format (if not specified)
    if (params.aln_format == IN_UNKNOWN)
        params.aln_format = detectInputFile(params.aln_path.c_str());
    
    // if alignment is in PHYLIP or FASTA format -> convert to MAPLE format
    if (params.aln_format != IN_MAPLE)
    {
        // record the starting time
        auto start = getRealTime();
        
        // prepare output (MAPLE) file
        // init maple_path if it's blank
        if (!params.maple_path.length())
            params.maple_path = params.aln_path + ".maple";
        
        tree.aln->extractMapleFile(params.aln_path, params.maple_path, params, params.only_extract_maple);
        
        // record the end time and show the runtime
        auto end = getRealTime();
        cout << "The input alignment is converted into MAPLE format at " << params.maple_path << endl;
        cout << " - Converting time: " << end-start << endl;
    }
    // otherwise, alignment is in MAPLE format -> read it
    else
    {
        // update input file path
        params.maple_path = params.aln_path;
        params.aln_path = "";
        
        if (params.only_extract_maple)
            outError("To export a MAPLE file, please supply an alignment (in FASTA or PHYLIP format) via -aln <ALIGNMENT>");
        
        // only want to reconstruct the aln file from the MAPLE file
        if (params.output_aln)
            tree.aln->reconstructAln(params.maple_path, params.output_aln);
        // otherwise, read the MAPLE file
        else
            tree.aln->readMapleFile(params.maple_path, params.ref_path);
    }
}

template <const StateType num_states>
void CMaple::loadInputTree()
{
    // read tree from the input treefile
    tree.readTree(tree.params->input_treefile);
    
    // calculate all lower, upper left/right likelihoods
    tree.refreshAllLhs<num_states>(true);
    
    // update model params
    /*cout << " - Model params before updating: " << endl;
     tree.showModelParams();*/
    tree.updateModelParams<num_states>();
    /*cout << " - Model params after updating: " << endl;
     tree.showModelParams();*/
    
    // refresh all lower after updating model params
    tree.performDFS<&Tree::updateLowerLh<num_states>>();
    
    // set outdated = false at all nodes to avoid considering SPR moves at those nodes
    if (tree.params->tree_search_type == PARTIAL_TREE_SEARCH)
        tree.resetSPRFlags(false, true, false);
}

void CMaple::preInference()
{
    // validate input
    ASSERT(tree.aln->ref_seq.size() > 0 && "Reference sequence is not found!");
    ASSERT(tree.aln->data.size() >= 3 && "The number of input sequences must be at least 3! Please check and try again!");
    
    // use maple_path as the output prefix if users didn't specify it
    if (!tree.params->output_prefix.length())
        tree.params->output_prefix = tree.params->maple_path;
    
    // check whether output file is already exists
    string output_file(tree.params->output_prefix);
    output_file += ".treefile";
    if (!tree.params->overwrite_output && fileExists(output_file))
        outError("File " + output_file + " already exists. Use `--overwrite` option if you really want to redo the analysis and overwrite all output files.\n");
    
    
    // setup function pointers in tree
    tree.setup();
    
    // sort sequences by their distances to the reference sequence
    tree.aln->sortSeqsByDistances(tree.params->hamming_weight);
    
    // extract related info (freqs, log_freqs) of the ref sequence
    tree.model->extractRefInfo(tree.aln);
    
    // init the mutation matrix from a model name
    tree.model->initMutationMat();
    
    // compute cumulative rates of the ref sequence
    tree.model->computeCumulativeRate(tree.aln);
    
    // compute thresholds for approximations
    tree.params->threshold_prob2 = tree.params->threshold_prob * tree.params->threshold_prob;
    
    // setup function pointers for CMaple
    setupFuncPtrs(tree.aln->num_states);
}

template <const StateType num_states>
void CMaple::buildInitialTree()
{
    // record the start time
    auto start = getRealTime();
    
    // dummy variables
    std::unique_ptr<Alignment>& aln = tree.aln;
    std::unique_ptr<Model>& model = tree.model;
    const PositionType seq_length = aln->ref_seq.size();
    const PositionType num_seqs = aln->data.size();
    const bool with_input_tree = tree.params->input_treefile.length();
    PositionType num_new_sequences = num_seqs;
    Sequence* sequence = &tree.aln->data.front();
    tree.nodes.reserve(num_seqs + num_seqs);
    NumSeqsType i = 0;
    
    // if users don't input a tree -> create the root from the first sequence
    if (!with_input_tree)
    {
        // place the root node
        tree.root_vector_index = 0;
        tree.nodes.emplace_back(LeafNode(0));
        PhyloNode& root = tree.nodes[0];
        root.setPartialLh(TOP, std::move(sequence->getLowerLhVector(seq_length, num_states, aln->getSeqType())));
        root.getPartialLh(TOP)->computeTotalLhAtRoot<num_states>(root.getTotalLh(), model);
        root.setUpperLength(0);
        
        // move to the next sequence in the alignment
        ++sequence;
        ++i;
    }
    
    // iteratively place other samples (sequences)
    for (; i < num_seqs; ++i, ++sequence)
    {
        // don't add sequence that was already added in the input tree
        if (with_input_tree && sequence->is_added)
        {
            --num_new_sequences;
            continue;
        }
        
        // get the lower likelihood vector of the current sequence
        std::unique_ptr<SeqRegions> lower_regions = sequence->getLowerLhVector(seq_length, num_states, aln->getSeqType());
        
        // update the mutation matrix from empirical number of mutations observed from the recent sequences
        if (i % tree.params->mutation_update_period == 0)
            tree.model->updateMutationMatEmpirical(aln);
        
        // NHANLT: debug
        //if ((*sequence)->seq_name == "39")
        //    cout << "debug" <<endl;
        
        // seek a position for new sample placement
        Index selected_node_index;
        RealNumType best_lh_diff = MIN_NEGATIVE;
        bool is_mid_branch = false;
        RealNumType best_up_lh_diff = MIN_NEGATIVE;
        RealNumType best_down_lh_diff = MIN_NEGATIVE;
        Index best_child_index;
        tree.seekSamplePlacement<num_states>(Index(tree.root_vector_index, TOP), i, lower_regions, selected_node_index, best_lh_diff, is_mid_branch, best_up_lh_diff, best_down_lh_diff, best_child_index);
        
        // if new sample is not less informative than existing nodes (~selected_node != NULL) -> place the new sample in the existing tree
        if (selected_node_index.getMiniIndex() != UNDEFINED)
        {
            // place new sample as a descendant of a mid-branch point
            if (is_mid_branch)
                tree.placeNewSampleMidBranch<num_states>(selected_node_index, lower_regions, i, best_lh_diff);
            // otherwise, best lk so far is for appending directly to existing node
            else
                tree.placeNewSampleAtNode<num_states>(selected_node_index, lower_regions, i, best_lh_diff, best_up_lh_diff, best_down_lh_diff, best_child_index);
        }
        
        // NHANLT: debug
        //cout << "Added node " << (*sequence)->seq_name << endl;
        //cout << (*sequence)->seq_name << endl;
        //cout << tree.exportTreeString() << endl;
        
        //if ((*sequence)->seq_name == "2219")
        //{
            //cout << tree.exportTreeString() << endl;
            //string output_file(tree.params->output_prefix);
            //exportOutput(output_file + "_init.treefile");
            //exit(0);
        //}
    }
    
    // flag denotes whether there is any new nodes added
    // show the number of new sequences added to the tree
    if (num_new_sequences > 0)
    {
        std::cout << num_new_sequences << " sequences have been added to the tree." << std::endl;
        
        // traverse the intial tree from root to re-calculate all likelihoods regarding the latest/final estimated model parameters
        tree.refreshAllLhs<num_states>();
    }
    else
        std::cout << "All sequences were presented in the input tree. No new sequence has been added!" << std::endl;
    
    // show the runtime for building an initial tree
    auto end = getRealTime();
    cout << " - Time spent on building an initial tree: " << std::setprecision(3) << end - start << endl;
}

template <const StateType num_states>
void CMaple::optimizeTree()
{
    // tree.params->debug = true;
    string output_file(tree.params->output_prefix);
    exportOutput(output_file + "_init.treefile");
    
    // run a short range search for tree topology improvement (if neccessary)
    // NOTES: don't apply short range search when users input a tree because only a few new sequences were added -> we only apply a deep SPR search
    if (tree.params->short_range_topo_search && tree.params->tree_search_type != NO_TREE_SEARCH && !(tree.params->input_treefile.length() && tree.params->tree_search_type == PARTIAL_TREE_SEARCH))
    {
        // apply short-range SPR search
        optimizeTreeTopology<num_states>(true);
        exportOutput(output_file + "_short_search.treefile");
        
        // reset the SPR flags so that we can start a deeper SPR search later
        tree.resetSPRFlags(false, true, true);
    }
    
    // output log-likelihood of the tree
    std::cout << std::setprecision(10) << "Tree log likelihood (before topo-opt): " << tree.calculateTreeLh<num_states>() << std::endl;
    
    // run a normal search for tree topology improvement
    if (tree.params->tree_search_type != NO_TREE_SEARCH)
    {
        optimizeTreeTopology<num_states>();
        exportOutput(output_file + "_topo.treefile");
    }
    
    // traverse the tree from root to re-calculate all likelihoods after optimizing the tree topology
    tree.refreshAllLhs<num_states>();
    
    // output log-likelihood of the tree
    std::cout << std::setprecision(10) << "Tree log likelihood (before optimizing branch lengths): " << tree.calculateTreeLh<num_states>() << std::endl;
    
    // do further optimization on branch lengths (if needed)
    if (tree.params->optimize_blength)
        optimizeBranchLengthsOfTree<num_states>();
    
    // NhanLT: update the model params
    if (tree.aln->getSeqType() == SEQ_DNA)
    {
        tree.model->initMutationMat();
        tree.updateModelParams<num_states>();
    }
    
    // traverse the tree from root to re-calculate all lower likelihoods after optimizing branch lengths
    tree.performDFS<&Tree::updateLowerLh<num_states>>();
    
    // output log-likelihood of the tree
    std::cout << std::setprecision(10) << "Tree log likelihood (after updating model): " << tree.calculateTreeLh<num_states>() << std::endl;
}

template <const StateType num_states>
void CMaple::optimizeTreeTopology(bool short_range_search)
{
    // record the start time
    auto start = getRealTime();
    int num_tree_improvement = short_range_search ? 1 : tree.params->num_tree_improvement;
    
    for (int i = 0; i < num_tree_improvement; ++i)
    {
        // first, set all nodes outdated
        // no need to do so anymore as new nodes were already marked as outdated
        // tree.resetSPRFlags(false, true, true);
        
        // traverse the tree from root to try improvements on the entire tree
        RealNumType improvement = tree.improveEntireTree<num_states>(short_range_search);
        
        // stop trying if the improvement is so small
        if (improvement < tree.params->thresh_entire_tree_improvement)
        {
            cout << "Small improvement, stopping topological search." << endl;
            break;
        }
        
        // run improvements only on the nodes that have been affected by some changes in the last round, and so on
        for (int j = 0; j < 20; ++j)
        {
            // forget SPR_applied flag to allow new SPR moves
            tree.resetSPRFlags(false, false, true);
            
            improvement = tree.improveEntireTree<num_states>(short_range_search);
            cout << "Tree was improved by " + convertDoubleToString(improvement) + " at subround " + convertIntToString(j + 1) << endl;
            
            // stop trying if the improvement is so small
            if (improvement < tree.params->thresh_entire_tree_improvement)
                break;
        }
            
    }
    
    // show the runtime for optimize the tree
    auto end = getRealTime();
    cout << " - Time spent on";
    cout << (short_range_search ? " a short range search for" : "");
    cout << " optimizing the tree topology: " << std::setprecision(3) << end - start << endl;
}

template <const StateType num_states>
void CMaple::optimizeBranchLengthsOfTree()
{
    // record the start time
    auto start = getRealTime();
    
    cout << "Start optimizing branch lengths" << endl;
    
    // first, set all nodes outdated
    tree.resetSPRFlags(false, true, true);
   
    // traverse the tree from root to optimize branch lengths
    PositionType num_improvement = tree.optimizeBranchLengths<num_states>();
   
    // run improvements only on the nodes that have been affected by some changes in the last round, and so on
    for (int j = 0; j < 20; ++j)
    {
        // stop trying if the improvement is so small
        if (num_improvement < tree.params->thresh_entire_tree_improvement)
        //if (num_improvement == 0)
            break;
        
        // traverse the tree from root to optimize branch lengths
        num_improvement = tree.optimizeBranchLengths<num_states>();
    }

    // show the runtime for optimize the branch lengths
    auto end = getRealTime();
    cout << " - Time spent on optimizing the branch lengths: " << std::setprecision(3) << end - start << endl;
}

void CMaple::doInference()
{
    (this->*doInferencePtr)();
}

template <const StateType num_states>
void CMaple::doInferenceTemplate()
{
    // 0. Load an input tree if users supply a treefile
    if (tree.params->input_treefile.length())
        loadInputTree<num_states>();
        
    // 1. Build an initial tree
    buildInitialTree<num_states>();
    
    // 2. Optimize the tree with SPR if there is any new nodes added to the tree
    optimizeTree<num_states>();
}

void CMaple::postInference()
{
    (this->*postInferencePtr)();
}

template <const StateType num_states>
void CMaple::postInferenceTemplate()
{
    // compute branch support if requested
    if (tree.params->compute_aLRT_SH)
        calculateBranchSupports<num_states>();
    
    // output log-likelihood of the tree
    std::cout << std::setprecision(10) << "Tree log likelihood: " << tree.calculateTreeLh<num_states>() << std::endl;
    
    // output treefile
    string output_file(tree.params->output_prefix);
    exportOutput(output_file + ".treefile");
    
    // output treefile
    if (tree.params->compute_aLRT_SH)
        exportOutput(output_file + "_aLRT_SH.treefile", true);
    
    // output model params
    std::cout << tree.exportModelString() << std::endl;
    
    // list output files
    std::cout << "Analysis results written to:" << std::endl;
    std::cout << "Maximum-likelihood tree:       " << output_file + ".treefile" << std::endl;
    if (tree.params->compute_aLRT_SH)
        std::cout << "Tree with aLRT-SH values:      " << output_file + "_aLRT_SH.treefile" << std::endl;
    std::cout << "Screen log file:               " << output_file + ".log" << std::endl << std::endl;
}

template <const StateType num_states>
void CMaple::calculateBranchSupports()
{
    // show current lh
    // std::cout << std::setprecision(10) << "Tree log likelihood (before calculating branch supports): " << tree.calculateTreeLh() << std::endl;
    
    // record the start time
    auto start = getRealTime();
    cout << "Start calculating branch supports" << endl;
    
    // update CMaple status
    status = BRANCH_SUPPORT_DONE;
    
    // calculate branch supports
    tree.calculateBranchSupports<num_states>();
    
    // show the runtime for calculating branch supports
    auto end = getRealTime();
    cout << " - Time spent on calculating branch supports: " << std::setprecision(3) << end - start << endl;
}

void CMaple::exportOutput(const string &filename, const bool show_branch_support)
{
    // open the tree file
    ofstream out = ofstream(filename);
    
    // write tree string into the tree file
    out << tree.exportTreeString(tree.params->export_binary_tree, show_branch_support) << endl;
    
    // close the output file
    out.close();
}

void CMaple::setupFuncPtrs(const StateType num_states)
{
    switch (num_states) {
        case 2:
            doInferencePtr = &CMaple::doInferenceTemplate<2>;
            postInferencePtr = &CMaple::postInferenceTemplate<2>;
            break;
        case 4:
            doInferencePtr = &CMaple::doInferenceTemplate<4>;
            postInferencePtr = &CMaple::postInferenceTemplate<4>;
            break;
        case 20:
            doInferencePtr = &CMaple::doInferenceTemplate<20>;
            postInferencePtr = &CMaple::postInferenceTemplate<20>;
            break;
            
        default:
            outError("Sorry! currently we only support DNA and Protein data!");
            break;
    }
}
