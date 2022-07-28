/*
 *  cmaple.h
 *  Created on: Jan 25, 2022
 *      Author: Nhan Ly-Trong
 */

#include "cmaple.h"

CMaple::CMaple()
{
    cumulative_rate = NULL;
    
    tree = NULL;
}

CMaple::CMaple(Params *n_params)
{
    cumulative_rate = NULL;
    
    tree = new Tree(n_params);
}

CMaple::~CMaple()
{
    if (cumulative_rate)
    {
        delete[] cumulative_rate;
        cumulative_rate = NULL;
    }
    
    if (tree)
    {
        delete tree;
        tree = NULL;
    }
}

void CMaple::loadInput()
{
    // extract sequences (in vectors of mutations) from an input alignment (in PHYLIP or FASTA format)
    if (tree->params->aln_path)
    {
        // record the starting time
        auto start = getRealTime();
        
        tree->aln->extractDiffFile(tree->params);
        
        // record the end time and show the runtime
        auto end = getRealTime();
        cout << "The input alignment is converted into DIFF format at " << tree->params->diff_path << endl;
        cout << " - Converting time: " << end-start << endl;
    }
    // or read sequences (in vectors of mutations) from a DIFF file
    else
    {
        if (tree->params->only_extract_diff)
            outError("To export a Diff file, please supple an alignment via --aln <ALIGNMENT>");
        
        // only want to reconstruc the aln file from the Diff file
        if (tree->params->output_aln)
            tree->aln->reconstructAln(tree->params->diff_path, tree->params->output_aln);
        // otherwise, read the Diff file
        else
            tree->aln->readDiff(tree->params->diff_path, tree->params->ref_path);
    }
}

void CMaple::preInference()
{
    // validate input
    ASSERT(tree->aln->ref_seq.size() > 0);
    if (tree->aln->size() < 3)
        outError("The number of input sequences must be at least 3! Please check and try again!");
    
    // check whether output file is already exists
    string output_file(tree->params->diff_path);
    output_file += ".treefile";
    if (!tree->params->redo_inference && fileExists(output_file))
        outError("File " + output_file + " already exists. Use `-redo` option if you really want to redo the analysis and overwrite all output files.\n");
    
    // sort sequences by their distances to the reference sequence
    tree->aln->sortSeqsByDistances(tree->params->hamming_weight);
    
    // extract related info (freqs, log_freqs) of the ref sequence
    tree->model->extractRefInfo(tree->aln->ref_seq, tree->aln->num_states);
    
    // init the mutation matrix from a model name
    tree->model->initMutationMat(tree->params->model_name, tree->aln->num_states);
    
    // compute cumulative rates of the ref sequence
    tree->model->computeCumulativeRate(cumulative_rate, cumulative_base, tree->aln);
    
    // compute the default initial branch length
    default_blength = 1.0 / tree->aln->ref_seq.size();
    min_blength = tree->params->min_blength_factor * default_blength;
    max_blength = tree->params->max_blength_factor * default_blength;
    min_blength_mid = tree->params->min_blength_mid_factor * default_blength;
    
    // compute thresholds for approximations
    threshold_prob2 = tree->params->threshold_prob * tree->params->threshold_prob;
    threshold_prob4 = threshold_prob2 * threshold_prob2;
}

void CMaple::buildInitialTree()
{
    // record the start time
    auto start = getRealTime();
    
    // dummy variables
    Alignment* aln = tree->aln;
    Model* model = tree->model;
    StateType num_states = aln->num_states;
    
    // place the root node
    Sequence* root_sequence = tree->aln->at(0);
    Node* root = new Node(root_sequence->seq_name);
    tree->root = root;
    root->partial_lh = root_sequence->getLowerLhVector(aln->ref_seq.size(), num_states, aln->seq_type);
    root->total_lh = root->partial_lh->computeTotalLhAtRoot(num_states, model);
    
    // iteratively place other samples (sequences)
    for (PositionType i = 1; i < aln->size(); i++)
    {
        Sequence* sequence = aln->at(i);
        
        // get the lower likelihood vector of the current sequence
        Regions* lower_regions = sequence->getLowerLhVector(aln->ref_seq.size(), num_states, aln->seq_type);
        
        // update the mutation matrix from empirical number of mutations observed from the recent sequences
        if (i % tree->params->mutation_update_period == 0)
            tree->model->updateMutationMatEmpirical(cumulative_rate, cumulative_base, aln);
        
        // NHANLT: debug
       /* if (sequence->seq_name == "614")
            cout << "debug" <<endl;*/
        
        // seek a position for new sample placement
        Node *selected_node, *best_child;
        double best_lh_diff, best_up_lh_diff, best_down_lh_diff;
        bool is_mid_branch;
        tree->seekPlacement(tree->root, sequence->seq_name, lower_regions, selected_node, best_lh_diff, is_mid_branch, best_up_lh_diff, best_down_lh_diff, best_child, cumulative_rate, default_blength, min_blength_mid);
        
        // if new sample is not less informative than existing nodes (~selected_node != NULL) -> place the new sample in the existing tree
        if (selected_node)
            tree->placeNewSample(selected_node, lower_regions, sequence->seq_name, best_lh_diff, is_mid_branch, best_up_lh_diff, best_down_lh_diff, best_child, cumulative_rate, cumulative_base, default_blength, max_blength, min_blength);
        
        // NHANLT: debug
        /*cout << "Added node " << sequence->seq_name << endl;
        cout << tree->exportTreeString() << ";" << endl;*/
        
        // don't delete lower_lh_seq as it is used as the lower lh regions of the newly adding tip
    }
    
    // show the runtime for building an initial tree
    auto end = getRealTime();
    cout << " - Time spent on building an initial tree: " << end - start << endl;
}

void CMaple::doInference()
{
    // 1. Build an initial tree
    buildInitialTree();
    
    // 2. Optimize the tree with SPR
}

void CMaple::postInference()
{
    // open the tree file
    string output_file(tree->params->diff_path);
    output_file += ".treefile";
    ofstream out = ofstream(output_file);
    
    // write tree string into the tree file
    out << tree->exportTreeString() << ";" << endl;
    
    // close the output file
    out.close();
}

void CMaple::tmpTestingMethod()
{
    // test some new methods
    
    cout << endl;
}

void runCMaple(Params &params)
{
    auto start = getRealTime();
    CMaple cmaple = CMaple(&params);
    
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
    
    // just test new method
    cmaple.tmpTestingMethod();
    
    // show runtime
    auto end = getRealTime();
    cout << "Runtime: " << end - start << "s" << endl;
}
