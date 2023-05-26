#include "cmaple.h"
using namespace std;
using namespace cmaple;

void CMaple::loadInput()
{
    ASSERT(tree.params->input_path);
    const InputType input_type = detectInputFile(tree.params->input_path);
    
    // extract sequences (in vectors of mutations) from an input alignment (in PHYLIP or FASTA format)
    if (input_type != IN_MAPLE)
    {
        // record the starting time
        auto start = getRealTime();
        
        tree.aln.extractDiffFile(*tree.params);
        
        // record the end time and show the runtime
        auto end = getRealTime();
        cout << "The input alignment is converted into DIFF format at " << tree.params->diff_path << endl;
        cout << " - Converting time: " << end-start << endl;
    }
    // or read sequences (in vectors of mutations) from a DIFF file
    else
    {
        // update input file path
        tree.params->diff_path = tree.params->input_path;
        tree.params->input_path = NULL;
        
        if (tree.params->only_extract_diff)
            outError("To export a Diff file, please supple an alignment via --input <ALIGNMENT>");
        
        // only want to reconstruc the aln file from the Diff file
        if (tree.params->output_aln)
            tree.aln.reconstructAln(tree.params->diff_path, tree.params->output_aln);
        // otherwise, read the Diff file
        else
            tree.aln.readDiff(tree.params->diff_path, tree.params->ref_path);
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
    
    // set outdated = true at all nodes to avoid considering SPR moves at those nodes
    if (!tree.params->redo_inference)
        tree.resetSPRFlags(false, true, false);
}

void CMaple::preInference()
{
    // validate input
    ASSERT(tree.aln.ref_seq.size() > 0 && "Reference sequence is not found!");
    ASSERT(tree.aln.data.size() >= 3 && "The number of input sequences must be at least 3! Please check and try again!");
    
    // use diff_path as the output prefix if users didn't specify it
    if (!tree.params->output_prefix)
    {
        string diff_path_str(tree.params->diff_path);
        tree.params->output_prefix = new char[diff_path_str.length() + 1];
        strcpy(tree.params->output_prefix, diff_path_str.c_str());
    }
    
    // check whether output file is already exists
    string output_file(tree.params->output_prefix);
    output_file += ".treefile";
    if (!tree.params->overwrite_output && fileExists(output_file))
        outError("File " + output_file + " already exists. Use `--overwrite` option if you really want to redo the analysis and overwrite all output files.\n");
    
    
    // setup function pointers in tree
    tree.setup();
    
    // sort sequences by their distances to the reference sequence
    tree.aln.sortSeqsByDistances(tree.params->hamming_weight);
    
    // extract related info (freqs, log_freqs) of the ref sequence
    tree.model->extractRefInfo(tree.aln);
    
    // init the mutation matrix from a model name
    tree.model->initMutationMat();
    
    // compute cumulative rates of the ref sequence
    tree.model->computeCumulativeRate(tree.aln);
    
    // compute thresholds for approximations
    tree.params->threshold_prob2 = tree.params->threshold_prob * tree.params->threshold_prob;
    
    // setup function pointers for CMaple
    setupFuncPtrs(tree.aln.num_states);
}

template <const StateType num_states>
bool CMaple::buildInitialTree()
{
    // record the start time
    auto start = getRealTime();
    
    // dummy variables
    Alignment& aln = tree.aln;
    std::unique_ptr<Model>& model = tree.model;
    const PositionType seq_length = aln.ref_seq.size();
    const PositionType num_seqs = aln.data.size();
    const bool with_input_tree = tree.params->input_treefile != NULL;
    PositionType num_new_sequences = num_seqs;
    Sequence* sequence = &tree.aln.data.front();
    tree.nodes.reserve(num_seqs + num_seqs);
    NumSeqsType i = 0;
    
    // if users don't input a tree -> create the root from the first sequence
    if (!with_input_tree)
    {
        // place the root node
        tree.root_vector_index = 0;
        tree.nodes.emplace_back(LeafNode(0));
        PhyloNode& root = tree.nodes[0];
        root.setPartialLh(TOP, std::move(sequence->getLowerLhVector(seq_length, num_states, aln.getSeqType())));
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
        std::unique_ptr<SeqRegions> lower_regions = sequence->getLowerLhVector(seq_length, num_states, aln.getSeqType());
        
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
        //cout << tree.exportTreeString() << ";" << endl;
        
        //if ((*sequence)->seq_name == "2219")
        //{
            //cout << tree.exportTreeString() << ";" << endl;
            //string output_file(tree.params->output_prefix);
            //exportOutput(output_file + "_init.treefile");
            //exit(0);
        //}
    }
    
    // flag denotes whether there is any new nodes added
    bool new_node_added = true;
    // show the number of new sequences added to the tree
    if (num_new_sequences > 0)
    {
        std::cout << num_new_sequences << " sequences have been added to the tree." << std::endl;
        
        // traverse the intial tree from root to re-calculate all likelihoods regarding the latest/final estimated model parameters
        tree.refreshAllLhs<num_states>();
    }
    else
    {
        new_node_added = false;
        std::cout << "All sequences were presented in the input tree. No new sequence has been added!" << std::endl;
    }
    
    // show the runtime for building an initial tree
    auto end = getRealTime();
    cout << " - Time spent on building an initial tree: " << std::setprecision(3) << end - start << endl;
    
    // return a flag denotes whether we should optimize the tree or not
    return new_node_added;
}

template <const StateType num_states>
void CMaple::optimizeTree(const bool new_sequences_added)
{
    // tree.params->debug = true;
    string output_file(tree.params->output_prefix);
    exportOutput(output_file + "_init.treefile");
    
    // run a short range search for tree topology improvement (if neccessary)
    // NOTES: don't apply short range search when users input a tree because only a few new sequences were added -> we only apply a deep SPR search
    if (tree.params->short_range_topo_search && !tree.params->input_treefile)
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
    if (new_sequences_added || tree.params->redo_inference)
    {
        optimizeTreeTopology<num_states>();
        exportOutput(output_file + "_topo.treefile");
    }
    
    // traverse the tree from root to re-calculate all likelihoods after optimizing the tree topology
    // tree.refreshAllLhs();
    
    // output log-likelihood of the tree
    std::cout << std::setprecision(10) << "Tree log likelihood (before optimizing branch lengths): " << tree.calculateTreeLh<num_states>() << std::endl;
    
    // do further optimization on branch lengths (if needed)
    if (tree.params->optimize_branch_length &&
        (new_sequences_added || tree.params->redo_blength))
        optimizeBranchLengthsOfTree<num_states>();
    
    // NhanLT: update the model params
    if (tree.aln.getSeqType() == SEQ_DNA)
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
    if (tree.params->input_treefile)
        loadInputTree<num_states>();
        
    // 1. Build an initial tree
    const bool new_sequences_added = buildInitialTree<num_states>();
    
    // 2. Optimize the tree with SPR if there is any new nodes added to the tree
    optimizeTree<num_states>(new_sequences_added);
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
    tree.showModelParams();
    
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
    out << tree.exportTreeString(tree.params->export_binary_tree, tree.root_vector_index, show_branch_support) << ";" << endl;
    
    // close the output file
    out.close();
}

void test()
{
    std::vector<PhyloNode> nodes;
    for (int i = 0; i < 100; i++)
        nodes.emplace_back(InternalNode());
    
    std::cout << "Fdsfsd" <<std::endl;
    // test phylonode is a leaf
    /*MiniIndex mini_index{RIGHT};
    Index tmp_index{0,mini_index};
    LeafNode leaf;
    leaf.less_info_seqs_.push_back(0);
    leaf.seq_name_index_ = 100;
    leaf.neighbor_index_ = Index(200, LEFT);
    PhyloNode phylonode1(std::move(leaf));
    std::cout << "Leaf node: " << std::endl;

    std::cout << "- seq_name_index: " << phylonode1.getSeqNameIndex() << std::endl;
    std::cout << "- neighbor_index: " << phylonode1.getNeighborIndex(mini_index) << std::endl;
    std::cout << "- (size of) less_info_seqs: " << phylonode1.getLessInfoSeqs().size() << std::endl;
    
    
    // change total_lh
    std::cout << " Update total_lh " << std::endl;
    std::unique_ptr<SeqRegions> new_total_lh = std::make_unique<SeqRegions>();
    new_total_lh->reserve(10);
    new_total_lh->emplace_back(0, 100);
    phylonode1.setTotalLh(std::move(new_total_lh));
    std::cout << "- (size of) total_lh (after updating): " << phylonode1.getTotalLh()->size() << std::endl;
    
    std::cout << " Update partial_lh " << std::endl;
    std::unique_ptr<SeqRegions> new_partial_lh = std::make_unique<SeqRegions>();
    new_partial_lh->emplace_back(0, 100);
    phylonode1.setPartialLh(mini_index, std::move(new_partial_lh));
    std::cout << "- (size of) partial_lh: " << phylonode1.getPartialLh(mini_index)->size() << std::endl;
    
    std::cout << " Update partial_lh 2nd time" << std::endl;
    std::unique_ptr<SeqRegions> new_partial_lh2 = std::make_unique<SeqRegions>();
    new_partial_lh2->emplace_back(1, 200);
    new_partial_lh2->emplace_back(2, 300);
    phylonode1.setPartialLh(mini_index, std::move(new_partial_lh2));
    std::cout << "- (size of) partial_lh: " << phylonode1.getPartialLh(mini_index)->size() << std::endl;
    
    // test phylonode is an internal node
    MiniIndex mini_index1{TOP};
    MiniIndex mini_index2{LEFT};
    MiniIndex mini_index3{RIGHT};
    InternalNode internal;
    internal.neighbor_index3_ = std::array{Index(400, LEFT),Index(500, TOP),Index(600, TOP)};
    PhyloNode phylonode2(std::move(internal));
    std::cout << "\n\nInternal node: " << std::endl;
    std::cout << "- neighbor_index (top, left, right): ";
    std::cout << phylonode2.getNeighborIndex(mini_index1) << " ";
    std::cout << phylonode2.getNeighborIndex(mini_index2) << " ";
    std::cout << phylonode2.getNeighborIndex(mini_index3) << std::endl;
    
    std::cout << " Update partial_lh at the left mininode" << std::endl;
    std::unique_ptr<SeqRegions> new_partial_lh1 = std::make_unique<SeqRegions>();
    new_partial_lh1->emplace_back(0, 100);
    new_partial_lh1->emplace_back(1, 200);
    phylonode2.setPartialLh(mini_index2, std::move(new_partial_lh1));
    std::cout << "- (size of) partial_lh at the left mininode: " << phylonode2.getPartialLh(mini_index2)->size() << std::endl;
    
    // create a phylonode as an internal node
    PhyloNode phylonode3;
    // test delete a phylonode created by a default constructor
    phylonode3.~PhyloNode();
    
    std::cout << "\n\n\nsize of a single MiniNode (old):  " << sizeof(Node) << '\n';
    // however, the actual memory allocated by the call to `new` will be larger, since the allocator used predefined block sizes:
    Node* p = new Node();
    Node* p2 = new Node();
    auto diff = (char*)p2 - (char*)p;
    std::cout << "size of a single MiniNode (old, when allocated by new): " << diff << '\n';

    std::cout << "size of a full PhyloNode (new): " << sizeof(PhyloNode) << '\n';
    std::cout << " + size of a InternalNode (new): " << sizeof(InternalNode) << '\n';
    std::cout << " + size of a LeafNode (new): " << sizeof(LeafNode) << '\n';
    // std::cout << " + size of a std::variant (new): " << sizeof(std::variant<InternalNode, LeafNode>) << '\n';
    
    const int nr_seqs = 5e5;
    const int old_nodes = nr_seqs * 3 + nr_seqs; // estimate: as many internal phylonodes (with 3 Node's each) as leaf nodes
    const int pylonodes = nr_seqs * 2; // number of internal nodes equals leaf nodes
    std::cout << "\n\nMemory usage:\n"
    "  for " << nr_seqs/1000 << "k Seqs:\n"
    //"    old nodes (" << old_nodes/1000 << "k) allocated with 'new'  : " << diff * old_nodes / 1024/ 1024 << " MB\n"
    "    new phylonodes (" << pylonodes/1000 << "k) in std::vector : " << sizeof(PhyloNode) * pylonodes / 1024/ 1024 << " MB\n";*/
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

void runCMaple(Params &params)
{
    // NHANLT: test new funtions
    // test();
    
    auto start = getRealTime();
    CMaple cmaple(params);
    
    // load input data
    cmaple.loadInput();
    
    // terminate if the user only wants to export a Diff file from an alignment
    // or only want to reconstruct an aln from a Diff file
    if (params.only_extract_diff || params.output_aln)
        return;
    
    // prepare for the inference
    cmaple.preInference();
    
    // debug
    // cmaple.tree.showModelParams();
    
    // infer trees and model params
    cmaple.doInference();
    
    // complete remaining stuff after the inference
    cmaple.postInference();
    
    // show runtime
    auto end = getRealTime();
    cout << "Runtime: " << end - start << "s" << endl;
}
