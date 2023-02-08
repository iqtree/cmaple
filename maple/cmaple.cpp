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
    
    // check whether output file is already exists
    string output_file(tree.params->diff_path);
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
    
    // place the root node
    Sequence* sequence = &tree.aln.data.front();
    tree.root = new Node(sequence->seq_name);
    tree.root->partial_lh = sequence->getLowerLhVector(seq_length, num_states, aln.seq_type);
    tree.root->computeTotalLhAtNode(aln, model, tree.params->threshold_prob, true);
    
    // move to the next sequence in the alignment
    ++sequence;
    
    // iteratively place other samples (sequences)
    for (PositionType i = 1; i < (PositionType) aln.data.size(); ++i, ++sequence)
    {
        // get the lower likelihood vector of the current sequence
        SeqRegions* lower_regions = sequence->getLowerLhVector(seq_length, num_states, aln.seq_type);
        
        // update the mutation matrix from empirical number of mutations observed from the recent sequences
        if (i % tree.params->mutation_update_period == 0)
            tree.model.updateMutationMatEmpirical(aln);
        
        // NHANLT: debug
        /*if ((*sequence)->seq_name == "39")
            cout << "debug" <<endl;*/
        
        // seek a position for new sample placement
        Node* selected_node = NULL;
        RealNumType best_lh_diff = MIN_NEGATIVE;
        bool is_mid_branch = false;
        RealNumType best_up_lh_diff = MIN_NEGATIVE;
        RealNumType best_down_lh_diff = MIN_NEGATIVE;
        Node* best_child = NULL;
        tree.seekSamplePlacement(tree.root, sequence->seq_name, lower_regions, selected_node, best_lh_diff, is_mid_branch, best_up_lh_diff, best_down_lh_diff, best_child);
        
        // if new sample is not less informative than existing nodes (~selected_node != NULL) -> place the new sample in the existing tree
        if (selected_node)
        {
            // place new sample as a descendant of a mid-branch point
            if (is_mid_branch)
                tree.placeNewSampleMidBranch(selected_node, lower_regions, sequence->seq_name, best_lh_diff);
            // otherwise, best lk so far is for appending directly to existing node
            else
                tree.placeNewSampleAtNode(selected_node, lower_regions, sequence->seq_name, best_lh_diff, best_up_lh_diff, best_down_lh_diff, best_child);
        }
        else
            delete lower_regions;
        
        // NHANLT: debug
        //cout << "Added node " << (*sequence)->seq_name << endl;
        //cout << (*sequence)->seq_name << endl;
        //cout << tree.exportTreeString() << ";" << endl;
        
        /*if ((*sequence)->seq_name == "2219")
        {
            //cout << tree.exportTreeString() << ";" << endl;
            string output_file(tree.params->diff_path);
            exportOutput(output_file + "_init.treefile");
            exit(0);
        }*/
    }
    
    // show the runtime for building an initial tree
    auto end = getRealTime();
    cout << " - Time spent on building an initial tree: " << end - start << endl;
    
    // traverse the intial tree from root to re-calculate all likelihoods regarding the latest/final estimated model parameters
    tree.refreshAllLhs();
}

void CMaple::optimizeTree()
{
    string output_file(tree.params->diff_path);
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
    
    // do further optimization on branch lengths (if needed)
    if (tree.params->optimize_branch_length)
        optimizeBranchLengthsOfTree();
}

void CMaple::optimizeTreeTopology(bool short_range_search)
{
    //  NHANLT - DELETE: reset the log file of the depths of pruning and regrafting branches
    // open the tree file
    string output_file(tree.params->diff_path);
    ofstream out = ofstream(output_file + ".statistics.txt");
    out << "Orig_prune \t Orig_Regraft \t Act_prune \t Act_regraft \t Dis_prune_regraft" << endl;
    // close the output file
    out.close();
    
    // record the start time
    auto start = getRealTime();
    int num_tree_improvement = short_range_search ? 1 : tree.params->num_tree_improvement;
    
    for (int i = 0; i < num_tree_improvement; ++i)
    {
        // first, set all nodes outdated
        tree.setAllNodeOutdated();
        
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
    tree.setAllNodeOutdated();
   
    // traverse the tree from root to optimize branch lengths
    PositionType num_improvement = tree.optimizeBranchLengths();
   
    // run improvements only on the nodes that have been affected by some changes in the last round, and so on
    for (int j = 0; j < 20; ++j)
    {
        // stop trying if the improvement is so small
        if (num_improvement < tree.params->thresh_entire_tree_improvement)
        //if (num_improvement == 0)
          //  break;
        
        // traverse the tree from root to optimize branch lengths
        num_improvement = tree.optimizeBranchLengths();
    }

    // show the runtime for optimize the branch lengths
    auto end = getRealTime();
    cout << " - Time spent on optimizing the branch lengths: " << end - start << endl;
}

void CMaple::doInference()
{
    // 1. Build an initial tree
    buildInitialTree();
    
    // 2. Optimize the tree with SPR
    optimizeTree();
}

void CMaple::postInference()
{
    string output_file(tree.params->diff_path);
    output_file += ".treefile";
    exportOutput(output_file);
}

void CMaple::exportOutput(const string &filename)
{
    // open the tree file
    ofstream out = ofstream(filename);
    
    // write tree string into the tree file
    out << tree.exportTreeString(tree.params->export_binary_tree) << ";" << endl;
    
    // close the output file
    out.close();
}

void runCMaple(Params &params)
{
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
