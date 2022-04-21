/*
 *  cmaple.h
 *  Created on: Jan 25, 2022
 *      Author: Nhan Ly-Trong
 */

#include "cmaple.h"

CMaple::CMaple()
{
    cumulative_rate = NULL;
    pseu_mutation_count = NULL;
    
    tree = NULL;
    default_blength = 0;
}

CMaple::CMaple(Params *n_params)
{
    cumulative_rate = NULL;
    pseu_mutation_count = NULL;
    
    tree = new Tree(n_params);
    default_blength = 0;
}

CMaple::~CMaple()
{
    if (cumulative_rate)
    {
        delete[] cumulative_rate;
        cumulative_rate = NULL;
    }
    
    if (pseu_mutation_count)
    {
        delete[] pseu_mutation_count;
        pseu_mutation_count = NULL;
    }
    
    if (tree)
    {
        delete tree;
        tree = NULL;
    }
}

void CMaple::extractDiffFile()
{
    // record the starting time
    auto start = getRealTime();
   
    // read input sequences
    ASSERT(tree->params->aln_path);
    StrVector sequences;
    StrVector seq_names;
    tree->aln->readSequences(tree->params->aln_path, sequences, seq_names);
    
    // validate the input sequences
    if (sequences.size() == 0)
        outError("Empty input sequences. Please check and try again!");
    // make sure all sequences have the same length
    if (detectInputFile(tree->params->aln_path) == IN_FASTA)
        for (PositionType i = 0; i < sequences.size(); i++)
        {
            if (sequences[i].length() != sequences[0].length())
                outError("Sequence " + seq_names[i] + " has a different length compared to the first sequence.");
        }
    
    // generate reference sequence from the input sequences
    string ref_sequence;
    // read the reference sequence from file (if the user supplies it)
    if (tree->params->ref_path)
        ref_sequence = tree->aln->readRef(tree->params->ref_path, tree->params->only_extract_diff);
    else
        ref_sequence = tree->aln->generateRef(sequences, seq_names, tree->params->only_extract_diff);
    
    // prepare output (Diff) file
    // init diff_path if it's null
    if (!tree->params->diff_path)
    {
        string diff_path_str(tree->params->aln_path);
        diff_path_str += ".diff";
        
        tree->params->diff_path = new char[diff_path_str.length() + 1];
        strcpy(tree->params->diff_path, diff_path_str.c_str());
    }
    // check whether the Diff file already exists
    string diff_file(tree->params->diff_path);
    if (!tree->params->redo_inference && fileExists(diff_file))
        outError("File " + diff_file + " already exists. Use `-redo` option if you want overwrite it.\n");
    
    // open the output file
    ofstream out = ofstream(tree->params->diff_path);
    
    // write reference sequence to the output file
    out << ">" << REF_NAME << endl;
    out << ref_sequence << endl;
    
    // extract and write mutations of each sequence to file
    tree->aln->extractMutations(sequences, seq_names, ref_sequence, out, tree->params->only_extract_diff);
    
    // close the output file
    out.close();
    
    // record the end time and show the runtime
    auto end = getRealTime();
    cout << "The input alignment is converted into DIFF format at " << tree->params->diff_path << endl;
    cout << " - Converting time: " << end-start << endl;
}

void CMaple::loadInput()
{
    // extract sequences (in vectors of mutations) from an input alignment (in PHYLIP or FASTA format)
    if (tree->params->aln_path)
        extractDiffFile();
    // or read sequences (in vectors of mutations) from a DIFF file
    else
    {
        if (tree->params->only_extract_diff)
            outError("To export a Diff file, please supple an alignment via --aln <ALIGNMENT>");
        
        tree->aln->readDiff(tree->params->diff_path, tree->params->ref_path);
    }
}


void CMaple::computeThresholds()
{
    threshold_prob2 = tree->params->threshold_prob * tree->params->threshold_prob;
    threshold_prob4 = threshold_prob2 * threshold_prob2;
}

void CMaple::computeCumulativeRate()
{
    PositionType sequence_length = tree->aln->ref_seq.size();
    ASSERT(sequence_length > 0);
    cumulative_rate = new double[sequence_length];
    
    cumulative_rate[0] = tree->model->mutation_mat[tree->aln->ref_seq[0] * (tree->aln->num_states + 1)];
    
    for (PositionType i = 1; i < sequence_length; i++)
        cumulative_rate[i] = cumulative_rate[i - 1] + tree->model->mutation_mat[tree->aln->ref_seq[i] * (tree->aln->num_states + 1)];
        
}

void CMaple::checkRedoInference()
{
    string output_file(tree->params->diff_path);
    output_file += ".treefile";
    
    if (!tree->params->redo_inference && fileExists(output_file))
        outError("File " + output_file + " already exists. Use `-redo` option if you really want to redo the analysis and overwrite all output files.\n");
}

void CMaple::preInference()
{
    // validate input
    ASSERT(tree->aln->ref_seq.size() > 0);
    if (tree->aln->size() < 3)
        outError("The number of input sequences must be at least 3! Please check and try again!");
    
    // check whether output file is already exists
    checkRedoInference();
    
    // sort sequences by their distances to the reference sequence
    tree->aln->sortSeqsByDistances(tree->params->hamming_weight);
    
    // extract related info (freqs, log_freqs) of the ref sequence
    tree->model->extractRefInfo(tree->aln->ref_seq, tree->aln->num_states);
    
    // init the mutation matrix from a model name
    tree->model->initMutationMat(tree->params->model_name, tree->aln->num_states, pseu_mutation_count);
    
    // compute cummulative rates of the ref sequence
    computeCumulativeRate();
    
    // compute the default initial branch length
    default_blength = 1.0 / tree->aln->ref_seq.size();
    
    // compute thresholds for approximations
    computeThresholds();
}

void CMaple::buildInitialTree()
{
    // Place the root node
    Node* root = new Node(tree->aln->at(0), tree->aln);
    root->overall_lh_seq = tree->getOverallLhSeqAtRoot(root->lower_lh_seq, 0, threshold_prob4);
    tree->root = root;
    
    // Iteratively place other samples (sequences)
    
}

void CMaple::doInference()
{
    // 1. Build an initial tree
    buildInitialTree();
    
}

void CMaple::tmpTestingMethod()
{
    // do something
    
    cout << endl;
}

void runCMaple(Params &params)
{
    auto start = getRealTime();
    CMaple cmaple = CMaple(&params);
    
    // load input data
    cmaple.loadInput();
    
    // terminate if the user only wants to export a Diff file from an alignment
    if (params.only_extract_diff)
        return;
    
    // prepare for the inference
    cmaple.preInference();
    
    // inference trees and model params
    cmaple.doInference();
    
    // just test new method
    cmaple.tmpTestingMethod();
    
    // show runtime
    auto end = getRealTime();
    cout << "Runtime: " << end - start << "s" << endl;
}
