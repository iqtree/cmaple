/*
 *  cmaple.h
 *  Created on: Jan 25, 2022
 *      Author: Nhan Ly-Trong
 */

#include "cmaple.h"
using namespace std;

void CMaple::loadInput()
{
    // extract sequences (in vectors of mutations) from an input alignment (in PHYLIP or FASTA format)
    if (tree.params->aln_path)
    {
        // record the starting time
        auto start = getRealTime();
        
        tree.aln.extractDiffFile(tree.params.value());
        
        // record the end time and show the runtime
        auto end = getRealTime();
        cout << "The input alignment is converted into DIFF format at " << tree.params->diff_path << endl;
        cout << " - Converting time: " << end-start << endl;
    }
    // or read sequences (in vectors of mutations) from a DIFF file
    else
    {
        if (tree.params->only_extract_diff)
            outError("To export a Diff file, please supple an alignment via --aln <ALIGNMENT>");
        
        // only want to reconstruc the aln file from the Diff file
        if (tree.params->output_aln)
            tree.aln.reconstructAln(tree.params->diff_path, tree.params->output_aln);
        // otherwise, read the Diff file
        else
            tree.aln.readDiff(tree.params->diff_path, tree.params->ref_path);
    }
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
    if (!tree.params->redo_inference && fileExists(output_file))
        outError("File " + output_file + " already exists. Use `-redo` option if you really want to redo the analysis and overwrite all output files.\n");
    
    // sort sequences by their distances to the reference sequence
    tree.aln.sortSeqsByDistances(tree.params->hamming_weight);
    
    // extract related info (freqs, log_freqs) of the ref sequence
    tree.model.extractRefInfo(tree.aln.ref_seq, tree.aln.num_states);
    
    // init the mutation matrix from a model name
    tree.model.initMutationMat(tree.params->model_name, tree.aln.num_states);
    
    // compute cumulative rates of the ref sequence
    tree.model.computeCumulativeRate(tree.aln);
    
    // setup function pointers in tree
    tree.setup();
    
    // compute thresholds for approximations
    tree.params->threshold_prob2 = tree.params->threshold_prob * tree.params->threshold_prob;
}

void CMaple::buildInitialTree()
{
    // record the start time
    auto start = getRealTime();
    
    // dummy variables
    Alignment& aln = tree.aln;
    Model& model = tree.model;
    const StateType num_states = aln.num_states;
    const PositionType seq_length = aln.ref_seq.size();
    const PositionType num_seqs = aln.data.size();
    
    // place the root node
    tree.nodes.reserve(num_seqs + num_seqs);
    tree.root_vector_index = 0;
    tree.nodes.emplace_back(LeafNode(0));
    PhyloNode& root = tree.nodes[0];
    Sequence* sequence = &tree.aln.data.front();
    root.setPartialLh(TOP, std::move(sequence->getLowerLhVector(seq_length, num_states, aln.seq_type)));
    root.getPartialLh(TOP)->computeTotalLhAtRoot(root.getTotalLh(), aln.num_states, model);
    root.setUpperLength(0);
    
    // move to the next sequence in the alignment
    ++sequence;
    
    // iteratively place other samples (sequences)
    for (NumSeqsType i = 1; i < num_seqs; ++i, ++sequence)
    {
        // get the lower likelihood vector of the current sequence
        std::unique_ptr<SeqRegions> lower_regions = sequence->getLowerLhVector(seq_length, num_states, aln.seq_type);
        
        // update the mutation matrix from empirical number of mutations observed from the recent sequences
        if (i % tree.params->mutation_update_period == 0)
            tree.model.updateMutationMatEmpirical(aln);
        
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
        tree.seekSamplePlacement(Index(tree.root_vector_index, TOP), i, lower_regions, selected_node_index, best_lh_diff, is_mid_branch, best_up_lh_diff, best_down_lh_diff, best_child_index);
        
        // if new sample is not less informative than existing nodes (~selected_node != NULL) -> place the new sample in the existing tree
        if (selected_node_index.getMiniIndex() != UNDEFINED)
        {
            // place new sample as a descendant of a mid-branch point
            if (is_mid_branch)
                tree.placeNewSampleMidBranch(selected_node_index, lower_regions, i, best_lh_diff);
            // otherwise, best lk so far is for appending directly to existing node
            else
                tree.placeNewSampleAtNode(selected_node_index, lower_regions, i, best_lh_diff, best_up_lh_diff, best_down_lh_diff, best_child_index);
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
    
    // show the runtime for building an initial tree
    auto end = getRealTime();
    cout << " - Time spent on building an initial tree: " << end - start << endl;
    
    // traverse the intial tree from root to re-calculate all likelihoods regarding the latest/final estimated model parameters
    tree.refreshAllLhs();
}

void CMaple::optimizeTree()
{
    // tree.params->debug = true;
    string output_file(tree.params->output_prefix);
    exportOutput(output_file + "_init.treefile");
    
    // run a short range search for tree topology improvement (if neccessary)
    if (tree.params->short_range_topo_search)
    {
        optimizeTreeTopology(true);
        exportOutput(output_file + "_short_search.treefile");
    }
    
    // run a normal search for tree topology improvement
    optimizeTreeTopology();
    exportOutput(output_file + "_topo.treefile");
    
    // traverse the tree from root to re-calculate all likelihoods after optimizing the tree topology
    // tree.refreshAllLhs();
    
    // output log-likelihood of the tree
    // std::cout << "Tree log likelihood (before optimizing branch lengths): " << tree.calculateTreeLh() << std::endl;
    
    // do further optimization on branch lengths (if needed)
    if (tree.params->optimize_branch_length)
        optimizeBranchLengthsOfTree();
    
    // traverse the tree from root to re-calculate all lower likelihoods after optimizing branch lengths
    tree.performDFS<&Tree::updateLowerLh>();
}

void CMaple::optimizeTreeTopology(bool short_range_search)
{
    // record the start time
    auto start = getRealTime();
    int num_tree_improvement = short_range_search ? 1 : tree.params->num_tree_improvement;
    
    for (int i = 0; i < num_tree_improvement; ++i)
    {
        // first, set all nodes outdated
        tree.resetSPRFlags(true);
        
        // traverse the tree from root to try improvements on the entire tree
        RealNumType improvement = tree.improveEntireTree(short_range_search);
        
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
            tree.resetSPRFlags(false);
            
            improvement = tree.improveEntireTree(short_range_search);
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
    cout << " optimizing the tree topology: " << end - start << endl;
}

void CMaple::optimizeBranchLengthsOfTree()
{
    // record the start time
    auto start = getRealTime();
    
    cout << "Start optimizing branch lengths" << endl;
    
    // first, set all nodes outdated
    tree.resetSPRFlags(true);
   
    // traverse the tree from root to optimize branch lengths
    PositionType num_improvement = tree.optimizeBranchLengths();
   
    // run improvements only on the nodes that have been affected by some changes in the last round, and so on
    for (int j = 0; j < 20; ++j)
    {
        // stop trying if the improvement is so small
        if (num_improvement < tree.params->thresh_entire_tree_improvement)
        //if (num_improvement == 0)
            break;
        
        // traverse the tree from root to optimize branch lengths
        num_improvement = tree.optimizeBranchLengths();
    }

    // show the runtime for optimize the branch lengths
    auto end = getRealTime();
    cout << " - Time spent on optimizing the branch lengths: " << end - start << endl;
}

void CMaple::doInference()
{
    // if users input a tree, an alignment and only need to compute aLRT-SH
    if (tree.params->compute_aLRT_SH && tree.params->input_treefile)
    {
        tree.readTree(tree.params->input_treefile);
        
        // calculate all lower, upper left/right likelihoods
        tree.refreshAllLhs();
        
        // update model params
        /*cout << " - Model params before updating: " << endl;
        tree.showModelParams();*/
        tree.updateModelParams();
        /*cout << " - Model params after updating: " << endl;
        tree.showModelParams();*/
        
        // refresh all lower, upper left/right likelihoods after updating model params
        tree.refreshAllLhs();
    }
    // otherwise, infer a phylogenetic tree from the alignment
    else
    {
        // 1. Build an initial tree
        buildInitialTree();
        
        // 2. Optimize the tree with SPR
        optimizeTree();
    }
}

void CMaple::postInference()
{
    // output log-likelihood of the tree
    std::cout << std::setprecision(20) << "Tree log likelihood (before calculating branch supports): " << tree.calculateTreeLh() << std::endl;
    
    // compute branch support if requested
    if (tree.params->compute_aLRT_SH)
        calculateBranchSupports();
    
    // output log-likelihood of the tree
    std::cout << std::setprecision(20) << "Tree log likelihood: " << tree.calculateTreeLh() << std::endl;
    
    // output treefile
    string output_file(tree.params->output_prefix);
    exportOutput(output_file + ".treefile");
    
    // output treefile
    if (tree.params->compute_aLRT_SH)
        exportOutput(output_file + "_aLRT_SH.treefile", true);
    
    // output model params
    tree.showModelParams();
}

void CMaple::calculateBranchSupports()
{
    // record the start time
    auto start = getRealTime();
    cout << "Start calculating branch supports" << endl;
    
    // calculate branch supports
    tree.calculateBranchSupports();
    
    // show the runtime for calculating branch supports
    auto end = getRealTime();
    cout << " - Time spent on calculating branch supports: " << end - start << endl;
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
    
    // infer trees and model params
    cmaple.doInference();
    
    // complete remaining stuff after the inference
    cmaple.postInference();
    
    // show runtime
    auto end = getRealTime();
    cout << "Runtime: " << end - start << "s" << endl;
}
