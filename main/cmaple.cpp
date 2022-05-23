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

void CMaple::examinePosition(Node* node, double parent_lh, int failure_count, string seq_name, vector<Region*> regions, Node* &selected_node, double &best_lh_diff , bool &is_mid_branch, double &best_up_lh_diff, double &best_down_lh_diff, Node* &best_child, bool &adjust_blen)
{
    // if the current node is a leaf AND the new sample/sequence is strictly less informative than the current node
    // -> add the new sequence into the list of minor sequences of the current node + stop seeking the placement
    if ((!node->relative) && (tree->aln->compareSequences(node->lower_lh_seq, regions) == 1))
    {
        node->less_info_seqs.push_back(seq_name);
        selected_node = node;
        best_lh_diff = 1.0;
        is_mid_branch = false;
        best_up_lh_diff = -DBL_MAX;
        best_down_lh_diff = -DBL_MAX;
        best_child = NULL;
        adjust_blen = false;
        return;
    }

    bool tmp_adjust_blen = false;
    double lh_diff_mid_branch = 0;
    double lh_diff_at_node = 0;
    
    // try first placing as a descendant of the mid-branch point of the branch above the current node
    /*if (node->parent && node->length > threshold_prob2)
    {
        lh_diff_mid_branch = calculatePlacementCost(node->overall_lh_seq, regions, default_blength);
        
        // try also placing with a longer new terminal branch, which can be useful if the sample has many new mutations.
        if (params->fastLK_blen_adjustment > 1)
        {
            double new_lh_diff = calculatePlacementCost(visit_node.node->genome_entries_total_lh_up, genomelist.genome_entries, blength * params->fastLK_blen_adjustment);
            
            if (new_lh_diff > lh_diff_mid_branch)
            {
                lh_diff_mid_branch = new_lh_diff;
                tmp_adjust_blen = true;
            }
            else
                tmp_adjust_blen = false;
        }
            
        if (lh_diff_mid_branch > best_lh_diff)
        {
            adjust_blen = tmp_adjust_blen;
            best_lh_diff = lh_diff_mid_branch;
            selected_node = visit_node.node;
            visit_node.failure_count = 0;
            is_mid_branch = true;
        }
    }
    else
        lh_diff_mid_branch = -DBL_MAX;

    // try to place as descendant of the current node (this is skipped if the node has top branch length 0 and so is part of a polytomy).
    if (visit_node.node->branch_length > thresh_prob2)
    {
        lh_diff_at_node = calculatePlacementCost(visit_node.node->genome_entries_total_lh, genomelist.genome_entries, blength);
        
        // try also placing with a longer new terminal branch, which can be useful if the sample has many new mutations.
        if (params->fastLK_blen_adjustment > 1)
        {
            double new_lh_diff = calculatePlacementCost(visit_node.node->genome_entries_total_lh, genomelist.genome_entries, blength * params->fastLK_blen_adjustment);
            
            if (new_lh_diff > lh_diff_at_node)
            {
                lh_diff_at_node = new_lh_diff;
                tmp_adjust_blen = true;
            }
            else
                tmp_adjust_blen = false;
        }
            
        if (lh_diff_at_node > best_lh_diff)
        {
            adjust_blen = tmp_adjust_blen;
            best_lh_diff = lh_diff_at_node;
            selected_node = visit_node.node;
            visit_node.failure_count = 0;
            is_mid_branch = false;
            best_up_lh_diff = lh_diff_mid_branch;
        }
        else if (lh_diff_mid_branch >= (best_lh_diff - params->fastLK_threshold_prob))
        {
            best_up_lh_diff = visit_node.parent_lh;
            best_down_lh_diff = lh_diff_at_node;
            selected_node = visit_node.node;
        }
        // placement at current node is considered failed if placement likelihood is not improved by a certain margin compared to best placement so far for the nodes above it.
        else if (lh_diff_at_node < (visit_node.parent_lh - 1.0))
            visit_node.failure_count++;
    }
    else
        lh_diff_at_node = visit_node.parent_lh;
     */

    // keep trying to place at children nodes, unless the number of attempts has reaches the failure limit
    /*if (failure_count < tree->params->failure_limit || lh_diff_at_node > (best_lh_diff - tree->params->fastLK_loglh_thresh))
        for (Neighbor neighbor: visit_node.node->neighbors)
            queue_visit_nodes.push(VisitingNode(neighbor.node, lh_diff_at_node, visit_node.failure_count));*/
}

void CMaple::seekPlacement(Node* start_node, string seq_name, vector<Region*> regions, Node* &selected_node, double &best_lh_diff , bool &is_mid_branch, double &best_up_lh_diff, double &best_down_lh_diff, Node* &best_child, bool &adjust_blen)
{
    // init variables
    selected_node = start_node;
    best_lh_diff = -DBL_MAX;
    is_mid_branch = false;
    best_up_lh_diff = -DBL_MAX;
    best_down_lh_diff = -DBL_MAX;
    best_child = NULL;
    adjust_blen = false;
    
    // recursively examine positions for placing the new sample
    examinePosition(start_node, -DBL_MAX, 0, seq_name, regions, selected_node, best_lh_diff, is_mid_branch, best_up_lh_diff, best_down_lh_diff, best_child, adjust_blen);

    /*// exploration of the tree is finished, and we are left with the node found so far with the best appending likelihood cost. Now we explore placement just below this node for more fine-grained placement within its descendant branches.
    double best_child_lh = -DBL_MAX;
    best_child = NULL;
    
    if (!is_mid_branch)
    {
        // current node might be part of a polytomy (represented by 0 branch lengths) so we want to explore all the children of the current node to find out if the best placement is actually in any of the branches below the current node.
        queue<Node*> queue_nodes;
        for (Neighbor neighbor: selected_node->neighbors)
            queue_nodes.push(neighbor.node);
        
        while (!queue_nodes.empty())
        {
            Node* node = queue_nodes.front();
            queue_nodes.pop();

            if (node->branch_length <= thresh_prob2)
            {
                for (Neighbor neighbor:node->neighbors)
                    queue_nodes.push(neighbor.node);
            }
            else
            {
                // now try to place on the current branch below the best node, at an height above the mid-branch.
                double new_blength = node->branch_length / 2;
                
                double new_best_lh_mid_branch = -DBL_MAX;
                int branch_deep = -1;
                vector<GenomeEntry> genome_entries_total_lh_mid_branch = node->genome_entries_total_lh_up;

                while (true)
                {
                    double new_lh_mid_branch = calculatePlacementCost(genome_entries_total_lh_mid_branch, genomelist.genome_entries, blength);
                    
                    if (params->fastLK_blen_adjustment > 1)
                    {
                        double tmp_new_lh_mid_branch = calculatePlacementCost(genome_entries_total_lh_mid_branch, genomelist.genome_entries, blength * params->fastLK_blen_adjustment);
                        
                        if (tmp_new_lh_mid_branch > new_lh_mid_branch)
                            new_lh_mid_branch = tmp_new_lh_mid_branch;
                    }
                    
                    if (new_lh_mid_branch > new_best_lh_mid_branch)
                        new_best_lh_mid_branch = new_lh_mid_branch;
                    else
                        break;
                    
                    new_blength /= params->fastLK_blength_split_factor;
                    
                    if (new_blength <= blength / (params->fastLK_blength_split_factor + params->fastLK_threshold_prob))
                        break;
                    
                    branch_deep++;
                    
                    genome_entries_total_lh_mid_branch = node->mid_branch_genome_entries[branch_deep];
                }
                
                if (new_best_lh_mid_branch > best_child_lh)
                {
                    best_child_lh = new_best_lh_mid_branch;
                    best_child = node;
                }
            }
        }
    }*/
}

void CMaple::updateMutationMatFromEmpiricalCount()
{
    StateType num_states = tree->aln->num_states;
    
    // backup the current mutation matrix
    double* tmp_mutation_mat = new double[num_states * num_states];
    memcpy(tmp_mutation_mat, tree->model->mutation_mat, num_states * num_states * sizeof(double));
    
    // update the mutation matrix regarding the pseu_mutation_count
    tree->model->updateMutationMat(pseu_mutation_count, num_states);
    
    // update cumulative_rate if the mutation matrix changes more than a threshold
    double change_thresh = 1e-3;
    double update = false;
    for (StateType j = 0; j < tree->aln->num_states; j++)
    {
        int index = j * (num_states + 1);
        if (fabs(tmp_mutation_mat[index] - tree->model->mutation_mat[index]) > change_thresh)
        {
            update = true;
            break;
        }
    }
    
    // update the cumulative_rate
    if (update)
        computeCumulativeRate();
    
    // delete tmp_mutation_mat
    delete[] tmp_mutation_mat;
}

void CMaple::addRoot()
{
    Sequence* root_sequence = tree->aln->at(0);
    Node* root = new Node(0, root_sequence->seq_name);
    root->lower_lh_seq = tree->getLowerLhSeqAtTip(root_sequence->mutations, 0);
    root->overall_lh = tree->getOverallLhSeqAtRoot(root->lower_lh_seq, 0, threshold_prob4);
    tree->root = root;
}

void CMaple::buildInitialTree()
{
    // record the start time
    auto start = getRealTime();
    
    // dummy variables
    Alignment* aln = tree->aln;
    Model* model = tree->model;
    StateType num_states = aln->num_states;
    
    // Place the root node
    addRoot();
    
    // Iteratively place other samples (sequences)
    for (PositionType i = 1; i < aln->size(); i++)
    {
        Sequence* sequence = aln->at(i);
        // get the lower likelihood sequence
        vector<Region*> lower_lh_seq = tree->getLowerLhSeqAtTip(sequence->mutations, 0);
        
        // update the mutation matrix from empirical number of mutations observed from the recent sequences
        if (i % tree->params->mutation_update_period == 0)
            updateMutationMatFromEmpiricalCount();
        
        // seek a position for new sample placement
        Node *selected_node, *best_child;
        double best_lh_diff, best_up_lh_diff, best_down_lh_diff;
        bool is_mid_branch, adjust_blen;
        seekPlacement(tree->root, sequence->seq_name, lower_lh_seq, selected_node, best_lh_diff, is_mid_branch, best_up_lh_diff, best_down_lh_diff, best_child, adjust_blen);
    }
    
    // show the runtime for building an initial tree
    auto end = getRealTime();
    cout << " - Time spent on building an initial tree: " << end - start << endl;
}

void CMaple::doInference()
{
    // 1. Build an initial tree
    buildInitialTree();
    
}

void CMaple::tmpTestingMethod()
{
    // do something
    int id = 0;
    tree = new Tree();
    Node* node = new Node(id++, "T0");
    
    tree->root = node;
    
    // add T6
    Node* new_node = new Node(id++, "T6");
    tree->addNode(new_node, node, id);
    
    // add T5
    Node* node_in_queue = node;
    node = new_node;
    new_node = new Node(id++, "T5");
    tree->addNode(new_node, node, id, false);
    
    // add T4
    node = node_in_queue;
    new_node = new Node(id++, "T4");
    tree->addNode(new_node, node, id);
    
    // add T2
    new_node = new Node(id++, "T2");
    tree->addNode(new_node, node, id);
    
    // add T1
    node_in_queue = new_node;
    new_node = new Node(id++, "T1");
    tree->addNode(new_node, node, id);
    
    // add T3
    node = node_in_queue;
    new_node = new Node(id++, "T3");
    tree->addNode(new_node, node, id);
    
    // add T7 at root
    new_node = new Node(id++, "T7");
    new_node->neighbor = tree->root;
    tree->root->neighbor = new_node;
    tree->root = new_node;
    
    traverseTree();
    
    cout << endl;
}

void CMaple::traverseTree(Node* node)
{
    // init starting node from root
    if (!node)
        node = tree->root;
    
    // do something at the current node when traversing tree
    cout << node->id << endl;
    
    // move to its neighbor
    node = node->neighbor;
    
    // do something with its neighbor
    if (node->isLeave())
        cout << "Leave " << node->seq_name << endl;
    else
        cout << node->id << endl;
        
    Node* relative;
    FOR_RELATIVE(node, relative)
        traverseTree(relative);
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
    //cmaple.doInference();
    
    // just test new method
    cmaple.tmpTestingMethod();
    
    // show runtime
    auto end = getRealTime();
    cout << "Runtime: " << end - start << "s" << endl;
}
