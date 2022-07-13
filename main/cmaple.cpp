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
    
    // init cumulative_base
    cumulative_base.resize(sequence_length);
    
    cumulative_rate[0] = tree->model->mutation_mat[tree->aln->ref_seq[0] * (tree->aln->num_states + 1)];
    cumulative_base[0].resize(tree->aln->num_states, 0);
    cumulative_base[0][tree->aln->ref_seq[0]] = 1;
    
    for (PositionType i = 1; i < sequence_length; i++)
    {
        StateType state = tree->aln->ref_seq[i];
        cumulative_rate[i] = cumulative_rate[i - 1] + tree->model->mutation_mat[state * (tree->aln->num_states + 1)];
        
        cumulative_base[i] =  cumulative_base[i -1];
        cumulative_base[i][state] = cumulative_base[i - 1][state] + 1;
    }
        
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
    min_blength = tree->params->min_blength_factor * default_blength;
    max_blength = tree->params->max_blength_factor * default_blength;
    
    // compute thresholds for approximations
    computeThresholds();
}

double CMaple::calculatePlacementCost(Regions* parent_regions, Regions* child_regions, double blength)
{
    // init dummy variables
    double lh_cost = 0;
    PositionType seq1_index = -1;
    PositionType seq2_index = -1;
    PositionType pos = 0;
    double total_factor = 1;
    StateType num_states = tree->aln->num_states;
    Region *seq1_region, *seq2_region;
    PositionType seq1_end = -1, seq2_end = -1;
    PositionType length;
    double minimum_carry_over = DBL_MIN * 1e50;
    if (blength < 0) blength = 0;
    double total_blength = blength;
    
    while (pos < tree->aln->ref_seq.size())
    {
        // get the next shared segment in the two sequences
        getNextSharedSegment(pos, parent_regions, child_regions, seq1_index, seq2_index, seq1_region, seq2_region, seq1_end, seq2_end, length);
        
        // 1. e1.type = N || e2.type = N
        if (seq2_region->type == TYPE_N || seq1_region->type == TYPE_N)
        {
            pos += length;
            continue;
        }
        else
        {
            // total_blength will be here the total length from the root or from the upper node, down to the down node.
            if (seq1_region->plength_from_root >= 0)
                total_blength = seq1_region->plength_from_root + (blength >= 0 ? blength : 0);
            else if (seq1_region->plength_observation >= 0)
                total_blength = seq1_region->plength_observation + (blength >= 0 ? blength : 0);
            else
                total_blength = blength;
                
            if (seq2_region->plength_observation >= 0)
                total_blength = (total_blength > 0 ? total_blength : 0) + seq2_region->plength_observation;
            
            // 2. e1.type = R
            if (seq1_region->type == TYPE_R)
            {
                // 2.1. e1.type = R and e2.type = R
                if (seq2_region->type == TYPE_R)
                {
                    if (seq1_region->plength_from_root >= 0)
                        total_blength += seq1_region->plength_observation;
                    
                    if (total_blength > 0)
                        lh_cost += total_blength * (cumulative_rate[pos + length - 1] - (pos == 0 ? 0 : cumulative_rate[pos - 1]));
                }
                // 2.2. e1.type = R and e2.type = O
                else if (seq2_region->type == TYPE_O)
                {
                    double tot = 0;
                    StateType seq1_state = tree->aln->ref_seq[pos];
                    if (seq1_region->plength_from_root >= 0)
                    {
                        for (StateType i = 0; i < num_states; i++)
                        {
                            double tot2;
                            StateType mutation_index = i * num_states + seq1_state;
                            if (seq1_state == i)
                                tot2 = tree->model->root_freqs[i] * (1.0 + tree->model->mutation_mat[mutation_index] * seq1_region->plength_observation);
                            else
                                tot2 = tree->model->root_freqs[i] * (tree->model->mutation_mat[mutation_index] * seq1_region->plength_observation);
                                
                            double tot3 = 0;
                            if (total_blength > 0)
                            {
                                for (StateType j = 0; j < num_states; j++)
                                    tot3 += tree->model->mutation_mat[i * num_states + j] * seq2_region->likelihood[j];
                            }
                            tot += tot2 * (seq2_region->likelihood[i] + total_blength * tot3);

                        }
                        
                        tot /= tree->model->root_freqs[seq1_state];
                    }
                    else
                    {
                        if (total_blength > 0)
                        {
                            for (StateType j = 0; j < num_states; j++)
                                tot += tree->model->mutation_mat[seq1_state * num_states + j] * seq2_region->likelihood[j];
                            
                            tot *= total_blength;
                        }
                        
                        tot += seq2_region->likelihood[seq1_state];
                    }
                        
                    total_factor *= tot;
                }
                // 2.3. e1.type = R and e2.type = A/C/G/T
                else
                {
                    StateType seq1_state = tree->aln->ref_seq[pos];
                    StateType seq2_state = seq2_region->type;
                    
                    if (seq1_region->plength_from_root >= 0)
                    {
                        if (total_blength > 0)
                            total_factor *= ((tree->model->root_freqs[seq1_state] * tree->model->mutation_mat[seq1_state * num_states + seq2_state] * total_blength * (1.0 + tree->model->mutation_mat[seq1_state * (num_states + 1)] * seq1_region->plength_observation) + tree->model->root_freqs[seq2_state] * tree->model->mutation_mat[seq2_state * num_states + seq1_state] * seq1_region->plength_observation * (1.0 + tree->model->mutation_mat[seq2_state * (num_states + 1)] * total_blength)) / tree->model->root_freqs[seq1_state]);
                        else
                            total_factor *= ((tree->model->root_freqs[seq2_state] * tree->model->mutation_mat[seq2_state * num_states + seq1_state] * seq1_region->plength_observation) / tree->model->root_freqs[seq1_state]);
                    }
                    else
                    {
                        if (total_blength > 0)
                            total_factor *= tree->model->mutation_mat[seq1_state * num_states + seq2_state] * total_blength;
                        else
                            return -DBL_MAX;
                    }
                }
            }
            // 3. e1.type = O
            else if (seq1_region->type == TYPE_O)
            {
                double blength13 = blength;
                if (seq1_region->plength_observation >= 0)
                {
                    blength13 = seq1_region->plength_observation;
                    if (blength > 0)
                        blength13 += blength;
                }
                    
                // 3.1. e1.type = O and e2.type = O
                if (seq2_region->type == TYPE_O)
                {
                    double tot = 0;
                    
                    for (StateType i = 0; i < num_states; i++)
                    {
                        double tot2 = 0;
                        
                        for (StateType j = 0; j < num_states; j++)
                            if (seq2_region->likelihood[j] > 0.1)
                                tot2 += tree->model->mutation_mat[i * num_states + j];
                        tot2 *= blength13;
                        
                        if (seq2_region->likelihood[i] > 0.1)
                            tot2 += 1;
                        
                        tot += tot2 * seq1_region->likelihood[i];
                    }
                        
                    total_factor *= tot;
                }
                // 3.2. e1.type = O and e2.type = R or A/C/G/T
                else
                {
                    StateType seq2_state = seq2_region->type;
                    if (seq2_state == TYPE_R)
                        seq2_state = tree->aln->ref_seq[pos];
                    
                    double tot2 = 0;
                    if (total_blength > 0)
                    {
                        for (StateType j = 0; j < num_states; j++)
                            tot2 += tree->model->mutation_mat[j * num_states + seq2_state] * seq1_region->likelihood[j];
                    }
                    
                    total_factor *= seq1_region->likelihood[seq2_state] + blength13 * tot2;
                }
            }
            // 4. e1.type = A/C/G/T
            else
            {
                int seq1_state = seq1_region->type;
                
                // 4.1. e1.type =  e2.type
                if (seq1_region->type == seq2_region->type)
                {
                    if (seq1_region->plength_from_root >= 0)
                        total_blength += seq1_region->plength_observation;
                    
                    if (total_blength > 0)
                        lh_cost += tree->model->mutation_mat[seq1_state * (num_states + 1)] * total_blength;
                }
                // e1.type = A/C/G/T and e2.type = O/A/C/G/T
                else
                {
                    // 4.2. e1.type = A/C/G/T and e2.type = O
                    if (seq2_region->type == TYPE_O)
                    {
                        double tot = 0.0;
                        
                        if (seq1_region->plength_from_root >= 0)
                        {
                            for (StateType i = 0; i < num_states; i++)
                            {
                                double tot2;
                                StateType mutation_index = i * num_states + seq1_state;
                                if (seq1_state == i)
                                    tot2 = tree->model->root_freqs[i] * (1.0 + tree->model->mutation_mat[mutation_index] * seq1_region->plength_observation);
                                else
                                    tot2 = tree->model->root_freqs[i] * (tree->model->mutation_mat[mutation_index] * seq1_region->plength_observation);
                                    
                                double tot3 = 0;
                                for (StateType j = 0; j < num_states; j++)
                                    tot3 += tree->model->mutation_mat[i * num_states + j] * seq2_region->likelihood[j];
                                tot += tot2 * (seq2_region->likelihood[i] + total_blength * tot3);
                            }
                            total_factor *= (tot / tree->model->root_freqs[seq1_state]);
                        }
                        else
                        {
                            for (StateType j = 0; j < num_states; j++)
                                tot += tree->model->mutation_mat[seq1_state * num_states + j] * seq2_region->likelihood[j];
                            
                            tot *= total_blength;
                            tot += seq2_region->likelihood[seq1_state];
                            
                            total_factor *= tot;
                        }
                    }
                    // 4.3. e1.type = A/C/G/T and e2.type = R or A/C/G/T
                    else
                    {
                        StateType seq2_state = seq2_region->type;
                        if (seq2_state == TYPE_R)
                            seq2_state = tree->aln->ref_seq[pos];
                        
                        if (seq1_region->plength_from_root >= 0)
                        {
                            if (total_blength > 0)
                                total_factor *= ((tree->model->root_freqs[seq1_state] * tree->model->mutation_mat[seq1_state * num_states + seq2_state] * total_blength * (1.0 + tree->model->mutation_mat[seq1_state * (num_states + 1)] * seq1_region->plength_observation) + tree->model->root_freqs[seq2_state] * tree->model->mutation_mat[seq2_state * num_states + seq1_state] * seq1_region->plength_observation * (1.0 + tree->model->mutation_mat[seq2_state * (num_states + 1)] * total_blength)) / tree->model->root_freqs[seq1_state]);
                            else
                                total_factor *= ((tree->model->root_freqs[seq2_state] * tree->model->mutation_mat[seq2_state * num_states + seq1_state] * seq1_region->plength_observation) / tree->model->root_freqs[seq1_state]);
                        }
                        else
                        {
                            if (total_blength > 0)
                                total_factor *= tree->model->mutation_mat[seq1_state * num_states + seq2_state] * total_blength;
                            else
                                return -DBL_MAX;
                        }
                    
                    }
                        
                }
            }
        }
         
        // approximately update lh_cost and total_factor
        if (total_factor <= minimum_carry_over)
        {
            if (total_factor < DBL_MIN)
                return -DBL_MAX;
            lh_cost += log(total_factor);
            total_factor = 1.0;
        }
        pos += length;
    }
    
    return lh_cost + log(total_factor);
}

void CMaple::seekPlacement(Node* start_node, string seq_name, Regions* regions, Node* &selected_node, double &best_lh_diff , bool &is_mid_branch, double &best_up_lh_diff, double &best_down_lh_diff, Node* &best_child, bool &adjust_blen)
{
    // init variables
    selected_node = start_node;
    best_lh_diff = -DBL_MAX;
    is_mid_branch = false;
    best_up_lh_diff = -DBL_MAX;
    best_down_lh_diff = -DBL_MAX;
    best_child = NULL;
    double lh_diff_mid_branch = 0;
    double lh_diff_at_node = 0;
    queue<Node*> node_queue;
    start_node->double_attributes[LH_DIFF] = -DBL_MAX;
    start_node->double_attributes[FAILURE_COUNT] = 0;
    node_queue.push(start_node);
    
    // recursively examine positions for placing the new sample
    while (!node_queue.empty())
    {
        Node* current_node = node_queue.front();
        node_queue.pop();
    
        // examine the current position
        // if the current node is a leaf AND the new sample/sequence is strictly less informative than the current node
        // -> add the new sequence into the list of minor sequences of the current node + stop seeking the placement
        if ((!current_node->next) && (compareSequences(current_node->partial_lh, regions) == 1))
        {
            current_node->less_info_seqs.push_back(seq_name);
            selected_node = current_node;
            best_lh_diff = 1.0;
            is_mid_branch = false;
            best_up_lh_diff = -DBL_MAX;
            best_down_lh_diff = -DBL_MAX;
            best_child = NULL;
            return;
        }
        
        // try first placing as a descendant of the mid-branch point of the branch above the current node
        if (current_node != tree->root && current_node->length >0)
        {
           
            Regions* total_mid_branch_regions = new Regions();
            mergeLhUpDown(total_mid_branch_regions, getPartialLhAtNode(current_node->neighbor), current_node->length / 2, getPartialLhAtNode(current_node), current_node->length / 2);
            lh_diff_mid_branch = calculatePlacementCost(total_mid_branch_regions, regions, default_blength);
            
            if (lh_diff_mid_branch > best_lh_diff)
            {
                best_lh_diff = lh_diff_mid_branch;
                selected_node = current_node;
                current_node->double_attributes[FAILURE_COUNT] = 0;
                is_mid_branch = true;
            }
            
            // delete total_mid_branch_regions
            delete total_mid_branch_regions;
        }
        else
            lh_diff_mid_branch = -DBL_MAX;

        // try to place as descendant of the current node (this is skipped if the node has top branch length 0 and so is part of a polytomy).
        if (current_node == tree->root || current_node->length > 0)
        {
            // compute the placement cost to place the new sample
            lh_diff_at_node = calculatePlacementCost(current_node->total_lh, regions, default_blength);
                
            if (lh_diff_at_node > best_lh_diff)
            {
                best_lh_diff = lh_diff_at_node;
                selected_node = current_node;
                current_node->double_attributes[FAILURE_COUNT] = 0;
                is_mid_branch = false;
                best_up_lh_diff = lh_diff_mid_branch;
            }
            else if (lh_diff_mid_branch >= (best_lh_diff - tree->params->threshold_prob))
            {
                best_up_lh_diff = current_node->double_attributes[LH_DIFF];
                best_down_lh_diff = lh_diff_at_node;
                selected_node = current_node;
            }
            // placement at current node is considered failed if placement likelihood is not improved by a certain margin compared to best placement so far for the nodes above it.
            else if (lh_diff_at_node < (current_node->double_attributes[LH_DIFF] - tree->params->thresh_log_lh_failure))
                current_node->double_attributes[FAILURE_COUNT] = current_node->double_attributes[FAILURE_COUNT] + 1;
        }
        else
            lh_diff_at_node = current_node->double_attributes[LH_DIFF];
        
        // keep trying to place at children nodes, unless the number of attempts has reaches the failure limit
        if ((tree->params->strict_stop_seeking_placement
             && current_node->double_attributes[FAILURE_COUNT]  < tree->params->failure_limit
             && lh_diff_at_node > (best_lh_diff - tree->params->thresh_log_lh_subtree_explore))
            || (!tree->params->strict_stop_seeking_placement
                && (current_node->double_attributes[FAILURE_COUNT]  < tree->params->failure_limit
                    || lh_diff_at_node > (best_lh_diff - tree->params->thresh_log_lh_subtree_explore))))
        {
            Node* neighbor_node;
            FOR_NEIGHBOR(current_node, neighbor_node)
            {
                neighbor_node->double_attributes[LH_DIFF] = lh_diff_at_node;
                neighbor_node->double_attributes[FAILURE_COUNT] = current_node->double_attributes[FAILURE_COUNT];
                node_queue.push(neighbor_node);
            }
        }
        
        // remove double_attributes of the current node
        current_node->double_attributes.erase(FAILURE_COUNT);
        current_node->double_attributes.erase(LH_DIFF);
    }

    // exploration of the tree is finished, and we are left with the node found so far with the best appending likelihood cost. Now we explore placement just below this node for more fine-grained placement within its descendant branches.
    double best_child_lh = -DBL_MAX;
    best_child = NULL;
    
    if (!is_mid_branch)
    {
        // current node might be part of a polytomy (represented by 0 branch lengths) so we want to explore all the children of the current node to find out if the best placement is actually in any of the branches below the current node.
        queue<Node*> queue_nodes;
        Node* neighbor_node;
        FOR_NEIGHBOR(selected_node, neighbor_node)
            queue_nodes.push(neighbor_node);
        
        while (!queue_nodes.empty())
        {
            Node* node = queue_nodes.front();
            queue_nodes.pop();

            if (node->length == 0)
            {
                FOR_NEIGHBOR(node, neighbor_node)
                    queue_nodes.push(neighbor_node);
            }
            else
            {
                // now try to place on the current branch below the best node, at an height above the mid-branch.
                double new_blength = node->length / 2;
                double new_best_lh_mid_branch = -DBL_MAX;
                Regions* upper_left_right_regions = node->neighbor->partial_lh;
                Regions* lower_regions = node->partial_lh;

                while (true)
                {
                    Regions* current_node_regions = new Regions();
                    mergeLhUpDown(current_node_regions, upper_left_right_regions, new_blength, lower_regions, node->length - new_blength);
                    
                    double new_lh_mid_branch = calculatePlacementCost(current_node_regions, regions, default_blength);
                    
                    /*if (params->fastLK_blen_adjustment > 1)
                    {
                        double tmp_new_lh_mid_branch = calculatePlacementCost(genome_entries_total_lh_mid_branch, genomelist.genome_entries, blength * params->fastLK_blen_adjustment);
                        
                        if (tmp_new_lh_mid_branch > new_lh_mid_branch)
                            new_lh_mid_branch = tmp_new_lh_mid_branch;
                    }*/
                    
                    if (new_lh_mid_branch > new_best_lh_mid_branch)
                        new_best_lh_mid_branch = new_lh_mid_branch;
                    else
                        break;
                    
                    new_blength /= 2;
                    
                    if (new_blength <= min_blength / 2)
                        break;
                    
                    delete current_node_regions;
                }
                
                if (new_best_lh_mid_branch > best_child_lh)
                {
                    best_child_lh = new_best_lh_mid_branch;
                    best_child = node;
                }
            }
        }
    }
}

void CMaple::mergeLhUpDown(Regions* &merged_regions, Regions* upper_regions, double upper_plength, Regions* lower_regions, double lower_plength)
{
    // init variables
    PositionType seq1_index = -1;
    PositionType seq2_index = -1;
    PositionType pos = 0;
    StateType num_states = tree->aln->num_states;
    Region *seq1_region, *seq2_region;
    PositionType seq1_end = -1, seq2_end = -1;
    PositionType length;
    
    // init merged_regions
    if (merged_regions) delete merged_regions;
    merged_regions = new Regions();
                
    while (pos <  tree->aln->ref_seq.size())
    {
        // get the next shared segment in the two sequences
        getNextSharedSegment(pos, upper_regions, lower_regions, seq1_index, seq2_index, seq1_region, seq2_region, seq1_end, seq2_end, length);
        
        // seq1_entry = 'N'
        if (seq1_region->type == TYPE_N)
        {
            // seq1_entry = 'N' and seq2_entry = 'N'
            if (seq2_region->type == TYPE_N)
                merged_regions->push_back(new Region(seq1_region->type, pos));
            // seq1_entry = 'N' and seq2_entry = O/R/ACGT
            else
            {
                // seq1_entry = 'N' and seq2_entry = O
                if (seq2_region->type == TYPE_O)
                {
                    double total_blength = lower_plength;
                    if (seq2_region->plength_observation >= 0)
                        total_blength = seq2_region->plength_observation + (lower_plength > 0 ? lower_plength: 0);
                        
                    double* new_lh = new double[num_states];
                    double sum_lh = 0;
                    
                    for (StateType i = 0; i < num_states; i++)
                    {
                        double tot = 0;
                        
                        if (total_blength > 0)
                        {
                            for (StateType j = 0; j < num_states; j++)
                                tot += tree->model->mutation_mat[i * num_states + j] * seq2_region->likelihood[j];
                            
                            tot *= total_blength;
                        }
                        
                        tot += seq2_region->likelihood[i];
                        new_lh[i] = tot * tree->model->root_freqs[i];
                        sum_lh += new_lh[i];
                    }

                    // normalize the new partial likelihood
                    for (StateType i = 0; i < num_states; i++)
                        new_lh[i] /= sum_lh;
                    
                    // add merged region into merged_regions
                    merged_regions->push_back(new Region(TYPE_O, pos, 0, 0, new_lh));
                }
                // seq1_entry = 'N' and seq2_entry = R/ACGT
                else
                {
                    if (seq2_region->plength_observation >= 0)
                    {
                        double new_plength = seq2_region->plength_observation ;
                        if (lower_plength > 0)
                            new_plength += lower_plength;
                        merged_regions->push_back(new Region(seq2_region->type, pos, new_plength, 0));
                    }
                    else
                    {
                        if (lower_plength > 0)
                            merged_regions->push_back(new Region(seq2_region->type, pos, lower_plength, 0));
                        else
                            merged_regions->push_back(new Region(seq2_region->type, pos));
                    }
                    
                }
            }
        }
        // seq2_entry = 'N'
        else if (seq2_region->type == TYPE_N)
        {
            // seq1_entry = 'O' and seq2_entry = N
            if (seq1_region->type == TYPE_O)
            {
                double total_blength = -1;
                
                if (seq1_region->plength_observation >= 0)
                {
                    total_blength = seq1_region->plength_observation;
                    if (upper_plength > 0)
                        total_blength += upper_plength;
                }
                else if (upper_plength > 0)
                    total_blength = upper_plength;
                
                if (total_blength > 0)
                {
                    double* new_lh = new double[num_states];
                    double sum_lh = 0;
                    
                    for (StateType i = 0; i < num_states; i++)
                    {
                        double tot = 0;
                        
                        for (StateType j = 0; j < num_states; j++)
                            tot += tree->model->mutation_mat[j * num_states + i] * seq1_region->likelihood[j];
                        
                        tot *= total_blength;
                        tot += seq1_region->likelihood[i];
                        new_lh[i] = tot;
                        sum_lh += new_lh[i];
                    }

                    // normalize the new partial likelihood
                    for (StateType i = 0; i < num_states; i++)
                        new_lh[i] /= sum_lh;
                    
                    // add merged region into merged_regions
                    merged_regions->push_back(new Region(TYPE_O, pos, 0, 0, new_lh));
                }
                else
                    // add merged region into merged_regions
                    merged_regions->push_back(new Region(seq1_region->type, pos, 0, 0, seq1_region->likelihood));
            }
            // seq2_entry = 'N' and seq1_entry = R/ACGT
            else
            {
                if (seq1_region->plength_from_root >= 0)
                {
                    double plength_from_root = seq1_region->plength_from_root;
                    if (upper_plength > 0)
                        plength_from_root += upper_plength;
                    merged_regions->push_back(new Region(seq1_region->type, pos, seq1_region->plength_observation, plength_from_root));
                }
                else if (seq1_region->plength_observation >= 0)
                {
                    double plength_observation = seq1_region->plength_observation;
                    if (upper_plength > 0)
                        plength_observation += upper_plength;
                    merged_regions->push_back(new Region(seq1_region->type, pos, plength_observation));
                }
                else
                {
                    if (upper_plength > 0)
                        merged_regions->push_back(new Region(seq1_region->type, pos, upper_plength));
                    else
                        merged_regions->push_back(new Region(seq1_region->type, pos));
                }
            }
        }
        // seq1_entry = seq2_entry = R/ACGT
        else if (seq1_region->type == seq2_region->type && (seq1_region->type < num_states || seq1_region->type == TYPE_R))
            merged_regions->push_back(new Region(seq1_region->type, pos));
        // cases where the new genome list entry will likely be of type "O"
        else
        {
            double total_blength_1 = upper_plength;
            if (seq1_region->plength_observation >= 0)
            {
                total_blength_1 = seq1_region->plength_observation;
                if (upper_plength > 0)
                    total_blength_1 += upper_plength;
                
                if (seq1_region->type != TYPE_O && seq1_region->plength_from_root >= 0)
                    total_blength_1 += seq1_region->plength_from_root;
            }
            
            double total_blength_2 = lower_plength;
            if (seq2_region->plength_observation >= 0)
            {
                total_blength_2 = seq2_region->plength_observation;
                if (lower_plength > 0)
                    total_blength_2 += lower_plength;
            }
            
            // due to 0 distance, the entry will be of same type as entry2
            if ((seq2_region->type < num_states || seq2_region->type == TYPE_R) && total_blength_2 < 0)
            {
                if ((seq1_region->type < num_states || seq1_region->type == TYPE_R) && total_blength_1 < 0)
                    //return None
                    outError("Sorry! something went wrong. DEBUG: ((seq2_region->type < num_states || seq2_region->type == TYPE_R) && total_blength_2 == 0) && ((seq1_region->type < num_states || seq1_region->type == TYPE_R) && total_blength_1 == 0)");
                
                merged_regions->push_back(new Region(seq2_region->type, pos));
            }
            // due to 0 distance, the entry will be of same type as entry1
            else if ((seq1_region->type < num_states || seq1_region->type == TYPE_R) && total_blength_1 == 0)
            {
                merged_regions->push_back(new Region(seq1_region->type, pos));
            }
            // seq1_entry = O
            else if (seq1_region->type == TYPE_O)
            {
                double* new_lh = new double[num_states];
                
                // if total_blength_1 > 0 => compute new partial likelihood
                if (total_blength_1 > 0)
                {
                    for (StateType i = 0; i < num_states; i++)
                    {
                        double tot = 0;
                        
                        for (StateType j = 0; j < num_states; j++)
                            tot += tree->model->mutation_mat[j * num_states + i] * seq1_region->likelihood[j];
                        
                        tot *= total_blength_1;
                        tot += seq1_region->likelihood[i];
                        new_lh[i] = tot;
                    }
                }
                // otherwise, clone the partial likelihood from seq1
                else
                    memcpy(new_lh, seq1_region->likelihood, sizeof(double) * num_states);
                
                double sum_new_lh = 0;

                // seq2 = O
                if (seq2_region->type == TYPE_O)
                {
                    for (StateType i = 0; i < num_states; i++)
                    {
                        double tot = 0;
                        
                        if (total_blength_2 > 0)
                        {
                            for (StateType j = 0; j < num_states; j++)
                                tot += tree->model->mutation_mat[i * num_states + j] * seq2_region->likelihood[j];
        
                            tot *= total_blength_2;
                        }
                        tot += seq2_region->likelihood[i];
                        new_lh[i] *= tot;
                        sum_new_lh += new_lh[i];
                    }
                }
                // seq1 = "O" and seq2 = ACGT
                else
                {
                    StateType seq2_state = seq2_region->type;
                    if (seq2_state == TYPE_R)
                        seq2_state = tree->aln->ref_seq[pos];
                    
                    if (total_blength_2 > 0)
                    {
                        for (StateType i = 0; i < num_states; i++)
                        {
                            StateType mutation_index = i * num_states + seq2_state;
                            if (i == seq2_state)
                                new_lh[i] *= (1.0 + tree->model->mutation_mat[mutation_index] * total_blength_2);
                            else
                                new_lh[i] *= (tree->model->mutation_mat[mutation_index] * total_blength_2);
                                
                            sum_new_lh += new_lh[i];
                        }
                    }
                    else
                    {
                        for (StateType i = 0; i < num_states; i++)
                        {
                            if (i != seq2_state)
                                new_lh[i] = 0;
                            
                            sum_new_lh += new_lh[i];
                        }
                    }
                }
                
                // normalize the new partial lh
                if (sum_new_lh == 0)
                    outError("Sum of new partital lh is zero.");
                for (StateType i = 0; i < num_states; i++)
                    new_lh[i] /= sum_new_lh;
                
                StateType new_state = simplify(new_lh, tree->aln->ref_seq[pos]);

                if (new_state == TYPE_O)
                    merged_regions->push_back(new Region(TYPE_O, pos, 0, 0, new_lh));
                else
                    merged_regions->push_back(new Region(new_state, pos));
            }
            // seq1_entry = R/ACGT
            else
            {
                StateType seq1_state = seq1_region->type;
                if (seq1_state == TYPE_R)
                    seq1_state = tree->aln->ref_seq[pos];
                
                double* new_lh = new double[num_states];
                double sum_new_lh = 0;
                
                if (seq1_region->plength_from_root >= 0)
                {
                    double length_to_root = seq1_region->plength_from_root;
                    if (upper_plength > 0)
                        length_to_root += upper_plength;
                    double* root_vec = new double[num_states];
                    memcpy(root_vec, tree->model->root_freqs, sizeof(double) * num_states);
                    
                    for (StateType i = 0; i < num_states; i++)
                    {
                        StateType mutation_index = i * num_states + seq1_state;
                        
                        if (i == seq1_state)
                            root_vec[i] *= (1.0 + tree->model->mutation_mat[mutation_index] * seq1_region->plength_observation);
                        else
                            root_vec[i] *= tree->model->mutation_mat[mutation_index] * seq1_region->plength_observation;
                    }
                        
                    for (StateType i = 0; i < num_states; i++)
                    {
                        double tot = 0;
                        
                        for (StateType j = 0; j < num_states; j++)
                            tot += tree->model->mutation_mat[j * num_states + i] * root_vec[j];
                        
                        tot *= length_to_root;
                        tot += root_vec[i];
                        new_lh[i] = tot;
                    }
                }
                else
                {
                    if (total_blength_1 > 0)
                    {
                        for (StateType i = 0; i < num_states; i++)
                        {
                            StateType mutation_index = seq1_state * num_states + i;
                            
                            if (i == seq1_state)
                                new_lh[i] = 1.0 + tree->model->mutation_mat[mutation_index] * total_blength_1;
                            else
                                new_lh[i] = tree->model->mutation_mat[mutation_index] * total_blength_1;
                        }
                    }
                    else
                    {
                        for (StateType i = 0; i < num_states; i++)
                        {
                            if (i == seq1_state)
                                new_lh[i] = 1;
                            else
                                new_lh[i] = 0;
                        }
                    }
                }
                  
                // seq2 = "O" and seq1 = ACGT
                if (seq2_region->type == TYPE_O)
                {
                    for (StateType i = 0; i < num_states; i++)
                    {
                        double tot = 0;
                        
                        if (total_blength_2 > 0)
                        {
                            for (StateType j = 0; j < num_states; j++)
                                tot += tree->model->mutation_mat[i * num_states + j] * seq2_region->likelihood[j];
                            
                            tot *= total_blength_2;
                        }
                        
                        tot += seq2_region->likelihood[i];
                        new_lh[i] *= tot;
                        sum_new_lh += new_lh[i];
                    }
                    
                    // normalize the new partial lh
                    for (StateType i = 0; i < num_states; i++)
                        new_lh[i] /= sum_new_lh;
                    
                    StateType new_state = simplify(new_lh, tree->aln->ref_seq[pos]);
                    
                    if (new_state == TYPE_O)
                        merged_regions->push_back(new Region(TYPE_O, pos, 0, 0, new_lh));
                    else
                        merged_regions->push_back(new Region(new_state, pos));
                }
                // seq1 = ACGT and different from seq2 = R/ACGT
                else
                {
                    StateType seq2_state = seq2_region->type;
                    
                    if (seq2_state == TYPE_R)
                        seq2_state = tree->aln->ref_seq[pos];
                    
                    if (total_blength_2 > 0)
                    {
                        for (StateType i = 0; i < num_states; i++)
                        {
                            StateType mutation_index = i * num_states + seq2_state;
                            
                            if (i == seq2_state)
                                new_lh[i] *= 1.0 + tree->model->mutation_mat[mutation_index] * total_blength_2;
                            else
                                new_lh[i] *= tree->model->mutation_mat[mutation_index] * total_blength_2;
                            
                            sum_new_lh += new_lh[i];
                        }
                    }
                    else
                    {
                        for (StateType i = 0; i < num_states; i++)
                        {
                            if (i != seq2_state)
                                new_lh[i] = 0;
                                
                            sum_new_lh += new_lh[i];
                        }
                    }
                    
                    // normalize the new partial lh
                    for (StateType i = 0; i < num_states; i++)
                        new_lh[i] /= sum_new_lh;

                    // add new region into the merged regions
                    merged_regions->push_back(new Region(TYPE_O, pos, 0, 0, new_lh));
                }
            }
        }

        // update pos
        pos += length;
    }

    // try to merge consecutive and similar 'R' regions together
    mergeRegionR(merged_regions);
}

double CMaple::mergeLhTwoLower(Regions* &merged_regions, Regions* regions1, double plength1, Regions* regions2, double plength2, bool return_log_lh)
{
    // init variables
    double log_lh = 0;
    PositionType seq1_index = -1;
    PositionType seq2_index = -1;
    PositionType pos = 0;
    StateType num_states = tree->aln->num_states;
    Region *seq1_region, *seq2_region;
    PositionType seq1_end = -1, seq2_end = -1;
    PositionType length;
    
    // init merged_regions
    if (merged_regions) delete merged_regions;
    merged_regions = new Regions();
                
    while (pos < tree->aln->ref_seq.size())
    {
        // get the next shared segment in the two sequences
        getNextSharedSegment(pos, regions1, regions2, seq1_index, seq2_index, seq1_region, seq2_region, seq1_end, seq2_end, length);
        
        // seq1_entry = 'N'
        if (seq1_region->type == TYPE_N)
        {
            // seq1_entry = 'N' and seq2_entry = 'N'
            if (seq2_region->type == TYPE_N)
                merged_regions->push_back(new Region(seq1_region->type, pos));
            // seq1_entry = 'N' and seq2_entry = O/R/ACGT
            else
            {
                // seq1_entry = 'N' and seq2_entry = O
                if (seq2_region->type == TYPE_O)
                {
                    // add merged region into merged_regions
                    Region* new_region = new Region(seq2_region, num_states);
                    new_region->position = pos;
                    if (seq2_region->plength_observation >= 0)
                    {
                        if (plength2 > 0)
                            new_region->plength_observation += plength2;
                    }
                    else
                    {
                        if (plength2 > 0)
                            new_region->plength_observation = plength2;
                    }
                    merged_regions->push_back(new_region);
                }
                // seq1_entry = 'N' and seq2_entry = R/ACGT
                else
                {
                    if (seq2_region->plength_observation >= 0)
                    {
                        double new_plength = seq2_region->plength_observation;
                        if (plength2 > 0)
                            new_plength += plength2;
                        merged_regions->push_back(new Region(seq2_region->type, pos, new_plength));
                    }
                    else
                    {
                        if (plength2 > 0)
                            merged_regions->push_back(new Region(seq2_region->type, pos, plength2));
                        else
                            merged_regions->push_back(new Region(seq2_region->type, pos));
                    }
                }
            }
        }
        // seq2_entry = 'N'
        else if (seq2_region->type == TYPE_N)
        {
            // seq1_entry = 'O' and seq2_entry = N
            if (seq1_region->type == TYPE_O)
            {
                // add merged region into merged_regions
                Region* new_region = new Region(seq1_region, num_states);
                new_region->position = pos;
                if (seq1_region->plength_observation >= 0)
                {
                    if (plength1 > 0)
                        new_region->plength_observation += plength1;
                }
                else
                {
                    if (plength1 > 0)
                        new_region->plength_observation = plength1;
                }
                merged_regions->push_back(new_region);
            }
            // seq1_entry = 'N' and seq2_entry = R/ACGT
            else
            {
                if (seq1_region->plength_observation >= 0)
                {
                    double new_plength = seq1_region->plength_observation;
                    if (plength1 > 0)
                        new_plength += plength1;
                    merged_regions->push_back(new Region(seq1_region->type, pos, new_plength));
                }
                else
                {
                    if (plength1 > 0)
                        merged_regions->push_back(new Region(seq1_region->type, pos, plength1));
                    else
                        merged_regions->push_back(new Region(seq1_region->type, pos));
                        
                }
            }
        }
        // neither seq1_entry nor seq2_entry = N
        else
        {
            double total_blength_1 = plength1;
            if (seq1_region->plength_observation >= 0)
            {
                total_blength_1 = seq1_region->plength_observation;
                if (plength1 > 0)
                    total_blength_1 += plength1;
            }
            
            double total_blength_2 = plength2;
            if (seq2_region->plength_observation >= 0)
            {
                total_blength_2 = seq2_region->plength_observation;
                if (plength2 > 0)
                    total_blength_2 += plength2;
            }
            
            // seq1_entry and seq2_entry are identical seq1_entry = R/ACGT
            if (seq1_region->type == seq2_region->type && (seq1_region->type == TYPE_R || seq1_region->type < num_states))
            {
                merged_regions->push_back(new Region(seq1_region->type, pos));
                
                if (return_log_lh)
                {
                    if (seq1_region->type == TYPE_R)
                        log_lh += (total_blength_1 + total_blength_2) * (cumulative_rate[pos + length - 1] - (pos == 0 ? 0 : cumulative_rate[pos - 1]));
                    else
                        log_lh += tree->model->mutation_mat[seq1_region->type * (num_states + 1)] * (total_blength_1 + total_blength_2);
                }
            }
            // #0 distance between different nucleotides: merge is not possible
            else if (total_blength_1 == 0 && total_blength_2 == 0 && (seq1_region->type == TYPE_R || seq1_region->type < num_states) && (seq2_region->type == TYPE_R || seq2_region->type < num_states))
            {
                outError("#0 distance between different nucleotides: merge is not possible");
                /*if (return_log_lh)
                    return None, float("-inf")
                else
                    return None*/
            }
            // seq1_entry = O
            else if (seq1_region->type == TYPE_O)
            {
                double* new_lh = new double[num_states];
                double sum_lh = 0;
                
                if (total_blength_1 > 0)
                {
                    for (StateType i = 0; i < num_states; i++)
                    {
                        double tot = 0;
                        
                        for (StateType j = 0; j < num_states; j++)
                            tot += tree->model->mutation_mat[i * num_states + j] * seq1_region->likelihood[j];
                        
                        tot *= total_blength_1;
                        tot += seq1_region->likelihood[i];
                        new_lh[i] = tot;
                    }
                }
                // otherwise, clone the partial likelihood from seq1
                else
                    memcpy(new_lh, seq1_region->likelihood, sizeof(double) * num_states);

                // seq1_entry = O and seq2_entry = O
                if (seq2_region->type == TYPE_O)
                {
                    for (StateType i = 0; i < num_states; i++)
                    {
                        double tot = 0;
                        
                        if (total_blength_2 > 0)
                        {
                            for (StateType j = 0; j < num_states; j++)
                                tot += tree->model->mutation_mat[i * num_states + j] * seq2_region->likelihood[j];
                            
                            tot *= total_blength_2;
                        }
                        tot += seq2_region->likelihood[i];
                        new_lh[i] *= tot;
                        sum_lh += new_lh[i];
                    }
                    
                    if (sum_lh == 0)
                        outError("Sum of partial likelihood is zero");
                        
                    // normalize new partial lh
                    for (StateType i = 0; i < num_states; i++)
                        new_lh[i] /= sum_lh;
                    
                    StateType new_state = simplify(new_lh, tree->aln->ref_seq[pos]);

                    if (new_state == TYPE_O)
                        merged_regions->push_back(new Region(new_state, pos, 0, 0, new_lh));
                    else
                        merged_regions->push_back(new Region(new_state, pos));
                    
                    if (return_log_lh)
                        log_lh += log(sum_lh);
                }
                // seq1_entry = O and seq2_entry = ACGT
                else
                {
                    StateType seq2_state = seq2_region->type;
                    if (seq2_state == TYPE_R)
                        seq2_state = tree->aln->ref_seq[pos];
                    
                    if (total_blength_2 > 0)
                    {
                        for (StateType i = 0; i < num_states; i++)
                        {
                            StateType mutation_index = i * num_states + seq2_state;
                            
                            if (seq2_state == i)
                                new_lh[i] *= (1 + tree->model->mutation_mat[mutation_index] * total_blength_2);
                            else
                                new_lh[i] *= (tree->model->mutation_mat[mutation_index] * total_blength_2);
                            
                            sum_lh += new_lh[i];
                        }
                        
                        // normalize new partial lh
                        for (StateType i = 0; i < num_states; i++)
                            new_lh[i] /= sum_lh;
                        
                        StateType new_state = simplify(new_lh, tree->aln->ref_seq[pos]);

                        if (new_state == TYPE_O)
                            merged_regions->push_back(new Region(new_state, pos, 0, 0, new_lh));
                        else
                            merged_regions->push_back(new Region(new_state, pos));
                        
                        if (return_log_lh)
                            log_lh += log(sum_lh);
                    }
                    else
                    {
                        if (new_lh[seq2_state] == 0)
                        {
                            outError("new_lh[seq2_state] == 0");
                            /*if returnLK:
                                return None, float("-inf")
                            else:
                                return None*/
                        }
                        
                        merged_regions->push_back(new Region(seq2_region->type, pos));
                    
                        if (return_log_lh)
                            log_lh += log(new_lh[seq2_state]);
                    }
                }
            }
            // seq1_entry = R/ACGT
            else
            {
                StateType seq1_state = seq1_region->type;
                if (seq1_state == TYPE_R)
                    seq1_state = tree->aln->ref_seq[pos];
                
                double* new_lh = new double[num_states];
                double sum_lh = 0;
                
                if (total_blength_1 > 0)
                {
                    for (StateType i = 0; i < num_states; i++)
                    {
                        StateType mutation_index = i * num_states + seq1_state;
                        
                        if (seq1_state == i)
                            new_lh[i] = 1 + tree->model->mutation_mat[mutation_index] * total_blength_1;
                        else
                            new_lh[i] = tree->model->mutation_mat[mutation_index] * total_blength_1;
                    }
                }
                else
                {
                    for (StateType i = 0; i < num_states; i++)
                    {
                        if (seq1_state == i)
                            new_lh[i] = 1;
                        else
                            new_lh[i] = 0;
                    }
                }

                // seq1_entry = ACGT and seq2_entry = O
                if (seq2_region->type == TYPE_O)
                {
                    for (StateType i = 0; i < num_states; i++)
                    {
                        double tot = 0;
                        
                        if (total_blength_2 > 0)
                        {
                            for (StateType j = 0; j < num_states; j++)
                                tot += tree->model->mutation_mat[i * num_states + j] * seq2_region->likelihood[j];
                            
                            tot *= total_blength_2;
                        }
                        tot += seq2_region->likelihood[i];
                        new_lh[i] *= tot;
                        sum_lh += new_lh[i];
                    }
                    
                    if (sum_lh == 0)
                        outError("Sum of partial likelihood is zero");
                        
                    // normalize new partial lh
                    for (StateType i = 0; i < num_states; i++)
                        new_lh[i] /= sum_lh;
                    
                    StateType new_state = simplify(new_lh, tree->aln->ref_seq[pos]);

                    if (new_state == TYPE_O)
                        merged_regions->push_back(new Region(new_state, pos, 0, 0, new_lh));
                    else
                        merged_regions->push_back(new Region(new_state, pos));
                    
                    if (return_log_lh)
                        log_lh += log(sum_lh);
                }
                // seq1_entry = ACGT and seq2_entry = R/ACGT
                else
                {
                    StateType seq2_state = seq2_region->type;
                    if (seq2_state == TYPE_R)
                        seq2_state = tree->aln->ref_seq[pos];
                    
                    if (total_blength_2 > 0)
                    {
                        for (StateType i = 0; i < num_states; i++)
                        {
                            StateType mutation_index = i * num_states + seq2_state;
                            
                            if (seq2_state == i)
                                new_lh[i] *= (1 + tree->model->mutation_mat[mutation_index] * total_blength_2);
                            else
                                new_lh[i] *= (tree->model->mutation_mat[mutation_index] * total_blength_2);
                            
                            sum_lh += new_lh[i];
                        }
                        
                        // normalize new partial lh
                        for (StateType i = 0; i < num_states; i++)
                            new_lh[i] /= sum_lh;
                        
                        StateType new_state = simplify(new_lh, tree->aln->ref_seq[pos]);

                        if (new_state == TYPE_O)
                            merged_regions->push_back(new Region(new_state, pos, 0, 0, new_lh));
                        else
                            merged_regions->push_back(new Region(new_state, pos));
                        
                        if (return_log_lh)
                            log_lh += log(sum_lh);
                    }
                    else
                    {
                        merged_regions->push_back(new Region(seq2_region->type, pos));
                    
                        if (return_log_lh)
                            log_lh += log(new_lh[seq2_state]);
                    }
                }
            }
        }
        // update pos
        pos += length;
    }

    // try to merge consecutive and similar 'R' regions together
    mergeRegionR(merged_regions);
    
    return log_lh;
}

StateType CMaple::simplify(double* &partial_lh, StateType ref_state)
{
    ASSERT(partial_lh);
    StateType num_states = tree->aln->num_states;
    double max_prob = 0;
    double max_index = 0;
    StateType high_prob_count = 0;
    
    for (StateType i = 0; i < num_states; i++)
    {
        if (partial_lh[i] > max_prob)
        {
            max_prob = partial_lh[i];
            max_index = i;
        }
        
        if (partial_lh[i] > tree->params->threshold_prob)
            high_prob_count++;
    }
    
    // return new state
    if (high_prob_count == 1)
    {
        if (max_index == ref_state)
            return TYPE_R;
        else
            return max_index;
    }
    else
        return TYPE_O;
}

void CMaple::mergeRegionR(Regions* &regions)
{
    Regions* new_regions = new Regions();
    
    Region* region = regions->getRegion(0);
    
    for (PositionType i = 0; i < regions->size() - 1; i++)
    {
        Region* next_region = regions->getRegion(i + 1);
        
        // merge two consecutive regions if they are both 'R' and "similar" to each other (regarding plength_observation, and plength_from_root)
        if (!(region->type == TYPE_R && next_region->type == TYPE_R && fabs(region->plength_observation - next_region->plength_observation) < tree->params->threshold_prob && fabs(region->plength_from_root - next_region->plength_from_root) < tree->params->threshold_prob))
        {
            // add the current region into new_regions
            new_regions->push_back(new Region(region, tree->aln->num_states, true));
            
            // move to the next entry
            region = regions->getRegion(i + 1);
        }
    }
    
    // add the last region into new_regions
    new_regions->push_back(new Region(region, tree->aln->num_states, true));
    
    // update regions
    regions->copyRegions(new_regions, tree->aln->num_states);
}

void CMaple::getLowerLhSeqAtTip(Regions* &regions, vector<Mutation*> mutations, double blength)
{
    // init regions
    if (regions)
        delete regions;
    regions = new Regions();
    
    PositionType sequence_length = tree->aln->ref_seq.size();
    PositionType pos = 0;
    
    for (Mutation* mutation: mutations)
    {
        // insert Region of type R (if necessary)
        if (mutation->position > pos)
            regions->push_back(new Region(TYPE_R, pos, tree->aln->seq_type, tree->aln->num_states));
        
        // convert the current mutation
        pos = mutation->position + mutation->getLength();
        regions->push_back(new Region(mutation, tree->aln->seq_type, tree->aln->num_states));
    }
    
    // insert the last Region of type R (if necessary)
    if (pos < sequence_length)
        regions->push_back(new Region(TYPE_R, pos, tree->aln->seq_type, tree->aln->num_states));
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
        StateType index = j * (num_states + 1);
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

void CMaple::updatePesudoCount(Regions* regions1, Regions* regions2)
{
    if (tree->model->model_name != "JC")
    {
        // init variables
        PositionType seq1_index = -1;
        PositionType seq2_index = -1;
        PositionType pos = 0;
        StateType num_states = tree->aln->num_states;
        Region *seq1_region, *seq2_region;
        PositionType seq1_end = -1, seq2_end = -1;
        PositionType length;
                    
        while (pos < tree->aln->ref_seq.size())
        {
            // get the next shared segment in the two sequences
            getNextSharedSegment(pos, regions1, regions2, seq1_index, seq2_index, seq1_region, seq2_region, seq1_end, seq2_end, length);
        
            if (seq1_region->type != seq2_region->type && (seq1_region->type < num_states || seq1_region->type == TYPE_R) && (seq2_region->type < num_states || seq2_region->type == TYPE_R))
            {
                if (seq1_region->type == TYPE_R)
                    pseu_mutation_count[tree->aln->ref_seq[pos] * num_states + seq2_region->type] += 1;
                else if (seq2_region->type == TYPE_R)
                    pseu_mutation_count[seq1_region->type * num_states + tree->aln->ref_seq[pos]] += 1;
                else
                    pseu_mutation_count[seq1_region->type * num_states + seq2_region->type] += 1;
            }

            // update pos
            pos += length;
        }
    }
}

void CMaple::addRoot()
{
    Sequence* root_sequence = tree->aln->at(0);
    Node* root = new Node(root_sequence->seq_name);
    getLowerLhSeqAtTip(root->partial_lh, root_sequence->mutations, 0);
    root->total_lh = computeTotalLhForRoot(root->partial_lh, -1);
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
        Regions* lower_lh_seq = new Regions();
        getLowerLhSeqAtTip(lower_lh_seq, sequence->mutations, 0);
        
        // update the mutation matrix from empirical number of mutations observed from the recent sequences
        if (i % tree->params->mutation_update_period == 0)
            updateMutationMatFromEmpiricalCount();
        
        // seek a position for new sample placement
        Node *selected_node, *best_child;
        double best_lh_diff, best_up_lh_diff, best_down_lh_diff;
        bool is_mid_branch, adjust_blen;
        seekPlacement(tree->root, sequence->seq_name, lower_lh_seq, selected_node, best_lh_diff, is_mid_branch, best_up_lh_diff, best_down_lh_diff, best_child, adjust_blen);
        
        // place the new sample in the existing tree
        if (best_lh_diff < 0.5)
            placeNewSample(selected_node, lower_lh_seq, sequence->seq_name, best_lh_diff, is_mid_branch, best_up_lh_diff, best_down_lh_diff, best_child, adjust_blen);
        
        cout << tree->exportTreeString() << endl;
        
        // don't delete lower_lh_seq as it is used as the lower lh regions of the newly adding tip
    }
    
    // show the runtime for building an initial tree
    auto end = getRealTime();
    cout << " - Time spent on building an initial tree: " << end - start << endl;
}

double CMaple::computeLhAtRoot(Regions* regions)
{
    double log_lh = 0;
    double log_factor = 1;
    StateType num_states = tree->aln->num_states;
    
    for (PositionType region_index = 0; region_index < regions->size(); region_index++)
    {
        Region* region = regions->getRegion(region_index);
        
        // type R
        if (region->type == TYPE_R)
        {
            PositionType start_pos = region->position;
            PositionType end_pos = (region_index == regions->size() - 1 ? tree->aln->ref_seq.size() - 1 : regions->getRegion(region_index + 1)->position - 1);
            
            for (StateType i = 0; i < num_states; i++)
                log_lh += tree->model->root_log_freqs[i] * (cumulative_base[end_pos][i] - (start_pos == 0 ? 0 : cumulative_base[start_pos - 1][i]));
        }
        // type ACGT
        else if (region->type < num_states)
            log_lh += tree->model->root_log_freqs[region->type];
        else if (region->type == TYPE_O)
        {
            double tot = 0;
            for (StateType i = 0; i < num_states; i++)
                tot += tree->model->root_freqs[i] * region->likelihood[i];
                log_factor *= tot;
        }
    }

    log_lh += log(log_factor);
    return log_lh;
}

Regions* CMaple::getPartialLhAtNode(Node* node)
{
    // compute partial_lh from next nodes
    if (!node->partial_lh)
    {
        // init partial_lh
        node->partial_lh = new Regions();
        
        // node must be an internal
        if (node->next)
        {
            // if node is a top node -> partial_lh is the lower lh regions
            if (node->double_attributes.find(IS_TOP_NODE) != node->double_attributes.end())
            {
                Node* next_node_1 = node->next;
                Regions* regions1 = getPartialLhAtNode(next_node_1->neighbor);
                Node* next_node_2 = next_node_1->next;
                Regions* regions2 = getPartialLhAtNode(next_node_2->neighbor);
                
                // compute partial_lh
                mergeLhTwoLower(node->partial_lh, regions1, next_node_1->length, regions2, next_node_2->length);
            }
            // otherwise -> partial_lh is the upper left/right regions
            else
            {
                Regions* upper_regions, *lower_regions;
                double upper_blength, lower_blength;
                
                Node* next_node_1 = node->next;
                Node* next_node_2 = next_node_1->next;
                
                if (next_node_1->double_attributes.find(IS_TOP_NODE) != next_node_1->double_attributes.end())
                {
                    upper_regions = getPartialLhAtNode(next_node_1->neighbor);
                    upper_blength = next_node_1->length;
                    lower_regions = getPartialLhAtNode(next_node_2->neighbor);
                    lower_blength = next_node_2->length;
                }
                else
                {
                    upper_regions = getPartialLhAtNode(next_node_2->neighbor);
                    upper_blength = next_node_2->length;
                    lower_regions = getPartialLhAtNode(next_node_1->neighbor);
                    lower_blength = next_node_1->length;
                }
                
                // compute partial_lh
                mergeLhUpDown(node->partial_lh, upper_regions, upper_blength, lower_regions, lower_blength);
            }
        }
        // otherwise, partial_lh must be already computed
        else
            outError("Something went wrong! Lower likelihood regions has not been computed at tip!");
    }
    
    // return
    return node->partial_lh;
}

Regions* CMaple::computeTotalLhAtNode(Node* node, bool update)
{
    Regions* new_regions;
    
    // if node is root
    if (node == tree->root)
        new_regions = computeTotalLhForRoot(getPartialLhAtNode(node), -1);
    // if not is normal nodes
    else
    {
        new_regions = new Regions();
        mergeLhUpDown(new_regions, getPartialLhAtNode(node->neighbor), node->length, getPartialLhAtNode(node), -1);
    }
    
    // update if necessary
    if (update)
    {
        if (node->total_lh) delete node->total_lh;
        node->total_lh = new_regions;
    }
    
    return new_regions;
}

Regions* CMaple::computeTotalLhForRoot(Regions* lower_regions, double blength)
{
    Regions* total_lh = new Regions();
    StateType num_states = tree->aln->num_states;
    
    for (Region* region : (*lower_regions))
    {
        // type N
        if (region->type == TYPE_N)
        {
            Region* new_region = new Region(region, num_states, false);
            total_lh->push_back(new_region);
        }
        else
            // type O
            if (region->type == TYPE_O)
            {
                double total_plength_observation = blength;
                if (region->plength_observation >= 0)
                {
                    total_plength_observation = region->plength_observation;
                    if (blength > 0)
                        total_plength_observation += blength;
                }
                double* new_partial_likelihood = new double[num_states];
                double sum_partial_likelihood = 0;
                
                for (StateType i = 0; i < num_states; i++)
                {
                    double tot = 0.0;
                    
                    if (total_plength_observation > 0)
                    {
                        for (StateType j = 0; j < num_states; j++)
                            tot += tree->model->mutation_mat[i * num_states + j] * region->likelihood[j];
                        tot *= total_plength_observation;
                    }
                    
                    tot += region->likelihood[i];
                    new_partial_likelihood[i] = tot * tree->model->root_freqs[i];
                    sum_partial_likelihood += new_partial_likelihood[i];
                }
                
                // normalize partial likelihood
                for (StateType i = 0; i < num_states; i++)
                    new_partial_likelihood[i] /= sum_partial_likelihood;
                
                // add new region to the total_lh_regions
                Region* new_region = new Region(region, num_states, false);
                new_region->likelihood = new_partial_likelihood;
                total_lh->push_back(new_region);
            }
            // other types: R or A/C/G/T
            else
            {
                // add new region to the total_lh_regions
                Region* new_region = new Region(region, num_states, false);
                
                if (new_region->plength_observation >= 0)
                {
                    if (blength > 0)
                        new_region->plength_observation += blength;

                    new_region->plength_from_root = 0;
                }
                else if (blength > 0)
                    new_region->plength_observation = blength;
                
                total_lh->push_back(new_region);
            }
    }
    
    return total_lh;
}

bool CMaple::areDiffRegions(Regions* regions1, Regions* regions2)
{
    if (!regions2 || regions2->size() == 0)
        return true;
    
    // init variables
    PositionType seq1_index = -1;
    PositionType seq2_index = -1;
    PositionType pos = 0;
    StateType num_states = tree->aln->num_states;
    Region *seq1_region, *seq2_region;
    PositionType seq1_end = -1, seq2_end = -1;
    PositionType length;
                
    while (pos < tree->aln->ref_seq.size())
    {
        // get the next shared segment in the two sequences
        getNextSharedSegment(pos, regions1, regions2, seq1_index, seq2_index, seq1_region, seq2_region, seq1_end, seq2_end, length);
        
        // compare each pair of regions
        if (seq1_region->type != seq2_region->type)
            return true;
            
        if (seq1_region->type < num_states || seq1_region->type == TYPE_R)
        {
            if (fabs(seq1_region->plength_from_root - seq2_region->plength_from_root) > tree->params->threshold_prob
                ||fabs(seq1_region->plength_observation - seq2_region->plength_observation) > tree->params->threshold_prob)
                return true;
        }
        
        if (seq1_region->type == TYPE_O)
        {
            if (fabs(seq1_region->plength_observation - seq2_region->plength_observation) > tree->params->threshold_prob)
                return true;
            
            for (StateType i = 0; i < num_states; i++)
            {
                double diff = fabs(seq1_region->likelihood[i] - seq2_region->likelihood[i]);
                
                if (diff > 0)
                {
                    if ((seq1_region->likelihood[i] == 0) || (seq2_region->likelihood[i] == 0))
                        return true;
                    
                    if (diff > tree->params->thresh_diff_update
                        || (diff > tree->params->threshold_prob
                            && ((diff / seq1_region->likelihood[i] > tree->params->thresh_diff_update)
                                || (diff / seq2_region->likelihood[i] > tree->params->thresh_diff_update))))
                        return true;
                }
            }
        }
        
        // update pos
        pos += length;
    }
    
    return false;
}

void CMaple::updateZeroBlength(queue<Node*> &node_queue, Node* node)
{
    // get the top node in the phylo-node
    Node* top_node = node->getTopNode();
    ASSERT(top_node);
    Regions* upper_left_right_regions = getPartialLhAtNode(top_node->neighbor);
    
    double best_lh = calculatePlacementCost(upper_left_right_regions, getPartialLhAtNode(top_node), default_blength);
    double best_length = default_blength;
    
    while (best_length > min_blength)
    {
        double new_blength = best_length/2;
        double new_lh = calculatePlacementCost(upper_left_right_regions, getPartialLhAtNode(top_node), new_blength);
        
        if (new_lh > best_lh)
        {
            best_lh = new_lh;
            best_length = new_blength;
        }
        else
            break;
    }
    
    if (best_length > 0.7 * default_blength)
    {
        while (best_length < max_blength)
        {
            double new_blength = best_length * 2;
            double new_lh = calculatePlacementCost(upper_left_right_regions, getPartialLhAtNode(top_node), new_blength);
            if (new_lh > best_lh)
            {
                best_lh = new_lh;
                best_length = new_blength;
            }
            else
                break;
        }
    }
    
    // update best_length
    top_node->length = best_length;
    top_node->neighbor->length = best_length;
    
    // add current node and its parent to node_queue to for updating partials further from these nodes
    top_node->dirty = true;
    top_node->neighbor->dirty = true;
    node_queue.push(top_node);
    node_queue.push(top_node->neighbor);
}

void CMaple::updatePartialLh(queue<Node*> &node_queue)
{
    StateType num_states = tree->aln->num_states;
    
    while (!node_queue.empty())
    {
        Node* node = node_queue.front();
        node_queue.pop();
        
        bool update_blength = false;
        node->dirty = true;
        
        Regions* parent_upper_regions = NULL;
        if (tree->root != node->getTopNode())
            parent_upper_regions = getPartialLhAtNode(node->neighbor);
            
        // change in likelihoods is coming from parent node
        if (node->double_attributes.find(IS_TOP_NODE) != node->double_attributes.end())
        {
            // if necessary, update the total probabilities at the mid node.
            if (node->length > 0)
            {
                Regions* regions = new Regions();
                mergeLhUpDown(regions, parent_upper_regions, node->length / 2, getPartialLhAtNode(node), node->length / 2);
                
                if (!regions)
                {
                    if (node->length > 1e-100)
                        outError("inside updatePartialLh(), from parent: should not have happened since node->length > 0");
                    /*updateZeroBlength(nodeList,node,mutMatrix)
                    update_blength=True*/
                }
                else
                {
                    // do nothing as I ignore mid-branch points
                    /*node.probVectTotUp=newTot
                    if node.dist>=2*minBLenForMidNode:
                        createFurtherMidNodes(node,parent_upper_regions)*/
                }
                
                // if necessary, update the total probability vector.
                if (!update_blength)
                {
                    computeTotalLhAtNode(node);
                    if (node->total_lh->size() == 0)
                    {
                        outError("inside updatePartialLh(), from parent 2: should not have happened since node->length > 0");
                        
                        /*updateZeroBlength(nodeList,node,mutMatrix)
                        update_blength=True
                        exit()*/
                    }
                }
                
                // delete regions
                delete regions;
            }
            
            // at valid internal node, update upLeft and upRight, and if necessary add children to node_queue.
            if (node->next && !update_blength)
            {
                Node* next_node_1 = node->next;
                Node* next_node_2 = next_node_1->next;
                
                Regions* upper_left_right_regions_1 = new Regions();
                Regions* upper_left_right_regions_2 = new Regions();
                mergeLhUpDown(upper_left_right_regions_1, parent_upper_regions, node->length, getPartialLhAtNode(next_node_1->neighbor), next_node_1->length);
                
                if (upper_left_right_regions_1->size() == 0)
                {
                    if (node->length == 0 && next_node_1->length == 0)
                        updateZeroBlength(node_queue, node);
                    else
                        outError("Strange: None vector from non-zero distances in updatePartialLh() from parent direction.");
                    
                    update_blength = true;
                }
                
                if (!update_blength)
                {
                    mergeLhUpDown(upper_left_right_regions_2, parent_upper_regions, node->length, getPartialLhAtNode(next_node_2->neighbor), next_node_2->length);
                    
                    if (upper_left_right_regions_2->size() == 0)
                    {
                        if (node->length == 0 && next_node_2->length == 0)
                            updateZeroBlength(node_queue, node);
                        else
                            outError("Strange: None vector from non-zero distances in updatePartialLh() from parent direction, child0.");
                    }
                    
                    update_blength = true;
                }
                
                if (!update_blength)
                {
                    if (areDiffRegions(getPartialLhAtNode(next_node_1), upper_left_right_regions_2))
                    {
                        next_node_1->partial_lh->copyRegions(upper_left_right_regions_2, tree->aln->num_states);
                        node_queue.push(next_node_1->neighbor);
                    }
                    
                    if (areDiffRegions(getPartialLhAtNode(next_node_2), upper_left_right_regions_1))
                    {
                        next_node_2->partial_lh->copyRegions(upper_left_right_regions_1, tree->aln->num_states);
                        node_queue.push(next_node_2->neighbor);
                    }
                }
                
                // delete upper_left_right_regions_1, upper_left_right_regions_2
                delete upper_left_right_regions_1;
                delete upper_left_right_regions_2;
            }
        }
        // otherwise, change in likelihoods is coming from a child.
        else
        {
            Node* top_node;
            Node* other_next_node;
            Node* next_node;
            FOR_NEXT(node, next_node)
            {
                if (next_node->double_attributes.find(IS_TOP_NODE) != next_node->double_attributes.end())
                    top_node = next_node;
                else
                    other_next_node = next_node;
            }
            
            ASSERT(top_node && other_next_node);
            
            double this_node_distance = node->length;
            double other_next_node_distance = other_next_node->length;
            
            Regions* this_node_lower_regions = getPartialLhAtNode(node->neighbor);
            Regions* this_node_upper_left_right_regions = getPartialLhAtNode(node);
            Regions* next_node_upper_left_right_regions = getPartialLhAtNode(other_next_node);
            
            // update lower likelihoods
            Regions* merged_two_lower_regions = new Regions();
            Regions* old_lower_regions;
            mergeLhTwoLower(merged_two_lower_regions, getPartialLhAtNode(other_next_node->neighbor), other_next_node_distance, this_node_lower_regions, this_node_distance);
            
            if (!merged_two_lower_regions || merged_two_lower_regions->size() == 0)
            {
                if (this_node_distance == 0 && other_next_node_distance == 0)
                {
                    updateZeroBlength(node_queue, node->neighbor);
                    update_blength = true;
                }
                else
                {
                    outError("Strange: None vector from non-zero distances in updatePartialLh() from child direction.");
                    /*print(this_node_distance)
                    print(this_node_lower_regions)
                    print(other_next_node_distance)
                    print(otherChildVect)
                    exit()*/
                }
            }
            else
            {
                old_lower_regions = new Regions();
                old_lower_regions->copyRegions(getPartialLhAtNode(top_node), num_states);
                top_node->partial_lh->copyRegions(merged_two_lower_regions, num_states);
            }

            // delete merged_two_lower_regions
            delete merged_two_lower_regions;
            
            // update total likelihood
            if (!update_blength)
            {
                if (top_node->length > 0 || top_node == tree->root)
                {
                    Regions* new_total_lh_regions = computeTotalLhAtNode(top_node, false);
                    
                    if (!new_total_lh_regions && top_node->length == 0)
                    {
                        updateZeroBlength(node_queue, top_node);
                        update_blength = true;
                    }
                    else if (!new_total_lh_regions)
                        outError("Strange: None vector from non-zero distances in updatePartialLh() from child direction while doing overall likelihood.");
                    else
                    {
                        // init top_node->total_lh
                        if (top_node->total_lh) delete top_node->total_lh;
                        top_node->total_lh = new Regions();
                        
                        top_node->total_lh->copyRegions(new_total_lh_regions, num_states);
                    }
                    
                    // delete new_total_lh_regions
                    delete new_total_lh_regions;
                }
            }
            
            /*#update total mid-branches likelihood
            if not update_blength:
                if node.dist and node.up!=None:
                    newTot=mergeLhUpDown(parent_upper_regions,node.dist/2,node.probVect,node.dist/2,mutMatrix)
                    if newTot==None:
                        updateZeroBlength(nodeList,node,mutMatrix)
                        update_blength=True
                        print("inside updatePartials(), from child: should not have happened since node.dist>0")
                    else:
                        node.probVectTotUp=newTot
                        if node.dist>=2*minBLenForMidNode:
                            createFurtherMidNodes(node,parent_upper_regions)
                */
            
            if (!update_blength)
            {
                // update likelihoods at parent node
                if (areDiffRegions(getPartialLhAtNode(top_node), old_lower_regions))
                {
                    if (tree->root != top_node)
                        node_queue.push(top_node->neighbor);
                }

                // update likelihoods at sibling node
                Regions* new_upper_regions = new Regions();
                if (tree->root != top_node)
                    mergeLhUpDown(new_upper_regions, parent_upper_regions, top_node->length, this_node_lower_regions, this_node_distance);
                else
                {
                    delete new_upper_regions;
                    new_upper_regions = computeTotalLhForRoot(getPartialLhAtNode(node->neighbor), this_node_distance);
                }
                
                if (!new_upper_regions || new_upper_regions->size() == 0)
                {
                    if (top_node->length == 0 && this_node_distance == 0)
                    {
                        updateZeroBlength(node_queue, top_node);
                        update_blength = true;
                    }
                    else
                        outError("Strange: None vector from non-zero distances in updatePartialLh() from child direction, new_upper_regions.");
                }
                else
                {
                    if (areDiffRegions(next_node_upper_left_right_regions, new_upper_regions))
                    {
                        other_next_node->partial_lh->copyRegions(new_upper_regions, num_states);
                        node_queue.push(other_next_node->neighbor);
                    }
                }
                    
                    // delete new_upper_regions
                    delete new_upper_regions;
                }
                    
            // delete old_lower_regions
            if (!old_lower_regions) delete old_lower_regions;
        }
    }
}

void CMaple::placeNewSample(Node* selected_node, Regions* sample, string seq_name, double best_lh_diff , bool is_mid_branch, double best_up_lh_diff, double best_down_lh_diff, Node* best_child, bool adjust_blen)
{
    // dummy variables
    double best_child_lh;
    double best_child_split;
    double best_parent_lh;
    double best_parent_split;
    Regions* best_parent_regions = new Regions();
    Regions* best_child_regions = new Regions();
    double best_root_blength;
    StateType num_states = tree->aln->num_states;
    
    // try to place the new sample as a descendant of a mid-branch point
    if (is_mid_branch)
    {
        Regions* upper_left_right_regions = getPartialLhAtNode(selected_node->neighbor);
        double best_split = 0.5;
        double best_split_lh = best_lh_diff;
        double new_split = 0.25;
        mergeLhUpDown(best_child_regions, getPartialLhAtNode(selected_node->neighbor), selected_node->length / 2, getPartialLhAtNode(selected_node), selected_node->length / 2);
        Regions* new_parent_regions = new Regions();

        // try different positions on the existing branch
        while (new_split * selected_node->length > min_blength)
        {
            mergeLhUpDown(new_parent_regions, upper_left_right_regions, selected_node->length * new_split, getPartialLhAtNode(selected_node),  selected_node->length * (1 - new_split));
            double placement_cost = calculatePlacementCost(new_parent_regions, sample, default_blength);
            
            if (placement_cost>best_split_lh)
            {
                best_split_lh = placement_cost;
                best_split = new_split;
                best_child_regions->copyRegions(new_parent_regions, num_states);
            }
            else
                break;
            new_split = best_split / 2;
        }
        
        if (best_split > 0.49)
        {
            new_split = 0.25;
            while (new_split * selected_node->length > min_blength)
            {
                mergeLhUpDown(new_parent_regions, upper_left_right_regions, selected_node->length * (1.0 - new_split), getPartialLhAtNode(selected_node),selected_node->length * new_split);
                
                double placement_cost = calculatePlacementCost(new_parent_regions, sample, default_blength);
                if (placement_cost > best_split_lh)
                {
                    best_split_lh = placement_cost;
                    best_split = new_split;
                    best_child_regions->copyRegions(new_parent_regions, num_states);
                }
                else
                    break;
                
                new_split = best_split / 2;
            }
            if (best_split < 0.49)
                best_split = 1.0 - best_split;
        }
        
        // delete new_parent_regions
        delete new_parent_regions;
        
        // now try different lengths for the new branch
        double new_branch_lh = best_split_lh;
        double best_blength = default_blength;
        while (best_blength > min_blength)
        {
            double new_blength = best_blength / 2;
            double placement_cost = calculatePlacementCost(best_child_regions, sample, new_blength);
            
            if (placement_cost > new_branch_lh)
            {
                new_branch_lh = placement_cost;
                best_blength = new_blength;
            }
            else
                break;
        }
        if (best_blength > 0.7 * default_blength)
        {
            while (best_blength < max_blength)
            {
                double new_blength = best_blength * 2;
                double placement_cost = calculatePlacementCost(best_child_regions, sample, new_blength);
                if (placement_cost > new_branch_lh)
                {
                    new_branch_lh = placement_cost;
                    best_blength = new_blength;
                }
                else
                    break;
            }
        }
        if (best_blength < min_blength)
        {
            double zero_branch_lh = calculatePlacementCost(best_child_regions, sample, -1);
            if (zero_branch_lh > new_branch_lh)
                best_blength = -1;
        }
        
        // create new internal node and append child to it
        Node* new_internal_node = new Node(true);
        Node* next_node_1 = new Node();
        Node* next_node_2 = new Node();
        Node* new_sample_node = new Node(seq_name);
        
        new_internal_node->next = next_node_1;
        next_node_1->next = next_node_2;
        next_node_2->next = new_internal_node;
        
        new_internal_node->neighbor = selected_node->neighbor;
        selected_node->neighbor->neighbor = new_internal_node;
        double top_distance = selected_node->length * best_split;
        new_internal_node->length = top_distance;
        new_internal_node->neighbor->length = top_distance;
        
        selected_node->neighbor = next_node_2;
        next_node_2->neighbor = selected_node;
        double down_distance = selected_node->length * (1 - best_split);
        selected_node->length = down_distance;
        selected_node->neighbor->length = down_distance;
        
        new_sample_node->neighbor = next_node_1;
        next_node_1->neighbor = new_sample_node;
        new_sample_node->length = best_blength;
        new_sample_node->neighbor->length = best_blength;
        
        new_sample_node->partial_lh = sample;
        next_node_1->partial_lh = new Regions();
        next_node_1->partial_lh->copyRegions(best_child_regions, num_states);
        mergeLhUpDown(next_node_2->partial_lh, upper_left_right_regions, new_internal_node->length, sample, best_blength);
        mergeLhTwoLower(new_internal_node->partial_lh, getPartialLhAtNode(selected_node), selected_node->length, sample, best_blength);
        //newInternalNode.probVectTotUp=mergeLhUpDown(upper_left_right_regions,distTop/2,newInternalNode.probVect,distTop/2,mutMatrix)
        computeTotalLhAtNode(new_internal_node);
        //mergeLhUpDown(new_internal_node->total_lh, best_child_regions, 0, sample, best_blength);
        
        if (!new_internal_node->total_lh || new_internal_node->total_lh->size() == 0)
            outError("Problem, None vector when placing sample, below node");
        
        /*if distTop>=2*min_blengthForMidNode:
         createFurtherMidNodes(newInternalNode,upper_left_right_regions)*/
        if (best_blength > 0)
        {
            computeTotalLhAtNode(new_sample_node);
            // mergeLhUpDown(new_sample_node->total_lh, best_child_regions, best_blength, sample, 0);
            /*newInternalNode.children[1].probVectTotUp=mergeLhUpDown(best_child_regions,best_blength/2,sample,best_blength/2,mutMatrix)
            if best_blength>=2*min_blengthForMidNode:
                createFurtherMidNodes(newInternalNode.children[1],best_child_regions)*/
        }
        
        // update pseudo_count
        updatePesudoCount(best_child_regions, sample);

        // iteratively traverse the tree to update partials from the current node
        queue<Node*> node_queue;
        node_queue.push(selected_node);
        node_queue.push(new_internal_node->neighbor);
        updatePartialLh(node_queue);
    }
    // otherwise, best lk so far is for appending directly to existing node
    else
    {
        // place the new sample as a descendant of an existing node
        if (best_child)
        {
            double best_split = 0.5;
            double best_split_lh = best_down_lh_diff;
            Regions* upper_left_right_regions = getPartialLhAtNode(best_child->neighbor);
            Regions* lower_regions = getPartialLhAtNode(best_child);
            best_child_regions = new Regions();
            best_child_regions->copyRegions(upper_left_right_regions, num_states);
            double new_split = 0.25;
            
            while (new_split * best_child->length > min_blength)
            {
                Regions* new_parent_regions = new Regions();
                mergeLhUpDown(new_parent_regions, upper_left_right_regions, best_child->length * new_split, lower_regions, best_child->length * (1 - new_split));
                
                double placement_cost = calculatePlacementCost(new_parent_regions, sample, default_blength);
                
                if (placement_cost > best_split_lh)
                {
                    best_split_lh = placement_cost;
                    best_split = new_split;
                    best_child_regions->copyRegions(new_parent_regions, num_states);
                }
                else
                    break;
                
                new_split = best_split / 2;
                
                // delete new_parent_regions
                delete new_parent_regions;
            }
            
            best_child_lh = best_split_lh;
            best_child_split = best_split;
        }
        else
            best_child_lh = -DBL_MAX;
        
        // if node is root, try to place as sibling of the current root.
        double old_root_lh;
        if (tree->root == selected_node)
        {
            old_root_lh = computeLhAtRoot(getPartialLhAtNode(selected_node));
            double new_root_lh;
            Regions* merged_root_sample_regions = new Regions();
            Regions* lower_regions = getPartialLhAtNode(selected_node);
            
            // merge 2 lower vector into one
            new_root_lh = mergeLhTwoLower(merged_root_sample_regions, lower_regions, default_blength, sample, default_blength, true);
            
            new_root_lh += computeLhAtRoot(merged_root_sample_regions);
            best_parent_lh = new_root_lh - old_root_lh;
            best_root_blength = default_blength;
            best_parent_regions->copyRegions(merged_root_sample_regions, num_states);
            double new_blength = 0.5 * default_blength;
            
            while (new_blength > min_blength)
            {
                // merge 2 lower vector into one
                new_root_lh = mergeLhTwoLower(merged_root_sample_regions, lower_regions, new_blength, sample, default_blength, true);
                
                new_root_lh += computeLhAtRoot(merged_root_sample_regions);
                double diff_root_lh = new_root_lh-old_root_lh;
                if (diff_root_lh > best_parent_lh)
                {
                    best_parent_lh = diff_root_lh;
                    best_root_blength = new_blength;
                    best_parent_regions->copyRegions(merged_root_sample_regions, num_states);
                }
                else
                    break;
                
                new_blength = best_root_blength / 2;
            }
            
            // delete merged_root_sample_regions
            delete merged_root_sample_regions;
        }
        // selected_node is not root
        else
        {
            double best_split = 0.5;
            double best_split_lh = best_up_lh_diff;
            Regions* upper_left_right_regions = getPartialLhAtNode(selected_node->neighbor);
            Regions* lower_regions = getPartialLhAtNode(selected_node);
            mergeLhUpDown(best_parent_regions, upper_left_right_regions, selected_node->length / 2, lower_regions, selected_node->length / 2);
            double new_split = 0.25;
            
            Regions* new_parent_regions = new Regions();
            while (new_split * selected_node->length > min_blength)
            {
                mergeLhUpDown(new_parent_regions, upper_left_right_regions, selected_node->length * (1 - new_split), lower_regions, selected_node->length * new_split);
                
                double placement_cost = calculatePlacementCost(new_parent_regions, sample, default_blength);
                
                if (placement_cost > best_split_lh)
                {
                    best_split_lh = placement_cost;
                    best_split = new_split;
                    best_parent_regions->copyRegions(new_parent_regions, num_states);
                }
                else
                    break;
                
                new_split = best_split / 2;
            }
            // delete new_parent_regions
            delete new_parent_regions;
            
            best_parent_lh = best_split_lh;
            best_parent_split = best_split;
        }
        
        // if the best placement is below the selected_node => add an internal node below the selected_node
        if (best_child_lh >= best_parent_lh && best_child_lh >= best_lh_diff)
        {
            Regions* upper_left_right_regions = getPartialLhAtNode(best_child->neighbor);

            double new_branch_length_lh = best_child_lh;
            double best_length = default_blength;
            
            while (best_length > min_blength)
            {
                double new_blength = best_length / 2;
                double placement_cost = calculatePlacementCost(upper_left_right_regions, sample, new_blength);
                
                if (placement_cost > placement_cost)
                {
                    new_branch_length_lh = placement_cost;
                    best_length = new_blength;
                }
                else
                    break;
            }
            
            if (best_length > 0.7 * default_blength)
            {
                while (best_length < max_blength)
                {
                    double new_blength = best_length * 2;
                    double placement_cost = calculatePlacementCost(best_child_regions, sample, new_blength);
                    if (placement_cost > new_branch_length_lh)
                    {
                        new_branch_length_lh = placement_cost;
                        best_length = new_blength;
                    }
                    else
                        break;
                }
            }
            
            if (best_length < min_blength)
            {
                double tmp_lh = calculatePlacementCost(best_child_regions, sample, -1);
                
                if (tmp_lh > new_branch_length_lh)
                    best_length = -1;
            }
            
            // create new internal node and append child to it
            Node* new_internal_node = new Node(true);
            Node* next_node_1 = new Node();
            Node* next_node_2 = new Node();
            Node* new_sample_node = new Node(seq_name);
            
            new_internal_node->next = next_node_1;
            next_node_1->next = next_node_2;
            next_node_2->next = new_internal_node;
            
            new_internal_node->neighbor = best_child->neighbor;
            best_child->neighbor->neighbor = new_internal_node;
            double top_distance = best_child->length * best_child_split;
            new_internal_node->length = top_distance;
            new_internal_node->neighbor->length = top_distance;
            
            best_child->neighbor = next_node_2;
            next_node_2->neighbor = best_child;
            double down_distance = best_child->length * (1 - best_child_split);
            best_child->length = down_distance;
            best_child->neighbor->length = down_distance;
            
            new_sample_node->neighbor = next_node_1;
            next_node_1->neighbor = new_sample_node;
            new_sample_node->length = best_length;
            new_sample_node->neighbor->length = best_length;
            
            new_sample_node->partial_lh = sample;
            next_node_1->partial_lh = new Regions();
            next_node_1->partial_lh->copyRegions(best_child_regions, num_states);
            mergeLhUpDown(next_node_2->partial_lh, upper_left_right_regions, new_internal_node->length, sample, best_length);
            mergeLhTwoLower(new_internal_node->partial_lh, getPartialLhAtNode(best_child), best_child->length, sample, best_length);
            //new_internal_node.probVectTotUp=mergeLhUpDown(this_node_upper_left_right_regions,top_distance/2,new_internal_node.probVect,top_distance/2,mutMatrix)
            //mergeLhUpDown(new_internal_node->total_lh, best_child_regions, 0, sample, best_length);
            computeTotalLhAtNode(new_internal_node);
            
            /*if (top_distance >= 2 * min_blengthForMidNode)
                createFurtherMidNodes(new_internal_node,this_node_upper_left_right_regions)*/
            if (best_length > 0)
            {
                computeTotalLhAtNode(new_sample_node);
                //mergeLhUpDown(new_sample_node->total_lh, best_child_regions, best_length, sample, 0);
                // new_sample_node.probVectTotUp=mergeLhUpDown(best_child_regions,best_length/2,sample,best_length/2,mutMatrix)
                /*if best_length>=2*min_blengthForMidNode:
                    createFurtherMidNodes(new_sample_node,best_child_regions)*/
            }
            
            // update pseudo_count
            updatePesudoCount(best_child_regions, sample);

            // iteratively traverse the tree to update partials from the current node
            queue<Node*> node_queue;
            node_queue.push(best_child);
            node_queue.push(new_internal_node->neighbor);
            updatePartialLh(node_queue);
        }
        // otherwise, add new parent to the selected_node
        else
        {
            // new parent is actually part of a polytomy since best placement is exactly at the node
            if (best_lh_diff >= best_parent_lh)
            {
                best_root_blength = -1;
                best_parent_split = -1;
                best_parent_lh = best_lh_diff;
                best_parent_regions->copyRegions(selected_node->total_lh, num_states);
                if (selected_node == tree->root)
                {
                    delete best_parent_regions;
                    best_parent_regions = new Regions();
                    
                    mergeLhTwoLower(best_parent_regions, getPartialLhAtNode(selected_node), -1, sample, default_blength);
                }
            }

            // add parent to the root
            if (selected_node == tree->root)
            {
                // now try different lengths for right branch
                double best_length2 = default_blength;
                Regions* new_root_lower_regions = new Regions();
                double new_root_lh = 0;
                
                while (best_length2 > min_blength)
                {
                    double new_blength = best_length2 / 2;
                    
                    new_root_lh = mergeLhTwoLower(new_root_lower_regions, getPartialLhAtNode(selected_node), best_root_blength, sample, new_blength, true);
                    
                    new_root_lh += computeLhAtRoot(new_root_lower_regions);
                    
                    double root_lh_diff = new_root_lh - old_root_lh;
                    if (root_lh_diff > best_parent_lh)
                    {
                        best_parent_lh = root_lh_diff;
                        best_length2 = new_blength;
                        best_parent_regions->copyRegions(new_root_lower_regions, num_states);
                    }
                    else
                        break;
                }
                
                if (best_length2 > 0.7 * default_blength)
                {
                    while (best_length2 < max_blength)
                    {
                        double new_blength = best_length2 * 2;
                        new_root_lh = mergeLhTwoLower(new_root_lower_regions, getPartialLhAtNode(selected_node), best_root_blength, sample, new_blength, true);
                        new_root_lh += computeLhAtRoot(new_root_lower_regions);
                        double root_lh_diff = new_root_lh - old_root_lh;
                        
                        if (root_lh_diff > best_parent_lh)
                        {
                            best_parent_lh = root_lh_diff;
                            best_length2 = new_blength;
                            best_parent_regions->copyRegions(new_root_lower_regions, num_states);
                        }
                        else
                            break;
                    }
                }
                
                // try with length zero
                if (best_length2 < min_blength)
                {
                    new_root_lh = mergeLhTwoLower(new_root_lower_regions, getPartialLhAtNode(selected_node), best_root_blength, sample, -1, true);
                    new_root_lh += computeLhAtRoot(new_root_lower_regions);
                    double root_lh_diff = new_root_lh - old_root_lh;
                    if (root_lh_diff > best_parent_lh)
                    {
                        best_length2 = -1;
                        best_parent_lh = root_lh_diff;
                        best_parent_regions->copyRegions(new_root_lower_regions, num_states);
                    }
                }

                // delete new_root_lower_regions
                delete new_root_lower_regions;
                
                // add new root node into tree
                Node* new_root = new Node(true);
                Node* next_node_1 = new Node();
                Node* next_node_2 = new Node();
                Node* new_sample_node = new Node(seq_name);
                
                new_root->next = next_node_1;
                next_node_1->next = next_node_2;
                next_node_2->next = new_root;
                
                // attach the left child
                selected_node->neighbor = next_node_2;
                next_node_2->neighbor = selected_node;
                selected_node->length = best_root_blength;
                selected_node->neighbor->length = best_root_blength;
                
                if (best_root_blength == 0)
                {
                    delete selected_node->total_lh;
                    selected_node->total_lh = NULL;
                    
                    //selected_node.probVectTotUp=None
                    //selected_node.furtherMidNodes=None
                }
                
                // attach the right child
                new_sample_node->neighbor = next_node_1;
                next_node_1->neighbor = new_sample_node;
                new_sample_node->length = best_length2;
                new_sample_node->neighbor->length = best_length2;
                
                new_root->partial_lh = new Regions();
                new_root->partial_lh->copyRegions(best_parent_regions, num_states);
                new_root->total_lh = computeTotalLhForRoot(getPartialLhAtNode(new_root), -1);

                next_node_1->partial_lh = computeTotalLhForRoot(getPartialLhAtNode(selected_node), best_root_blength);
                next_node_2->partial_lh = computeTotalLhForRoot(sample, best_length2);
                
                new_sample_node->partial_lh = sample;
                
                if (!new_root->total_lh || new_root->total_lh->size() == 0)
                {
                    outError("Problem, None vector when placing sample, new root");
                    /*print(merged_root_sample_regions)
                    print(node.probVect)
                    print(sample)
                    print(best_length2)
                    print(best_root_blength)*/
                }
                
                if (best_length2 > 0)
                {
                    computeTotalLhAtNode(new_sample_node);
                    //mergeLhUpDown(new_sample_node->total_lh, getPartialLhAtNode(new_sample_node->neighbor), best_length2, sample, 0);
                    
                    /*new_sample_node.probVectTotUp=mergeLhUpDown(new_root.probVectUpLeft,best_length2/2,sample,best_length2/2,mutMatrix)
                    if best_length2>=2*min_blengthForMidNode:
                        createFurtherMidNodes(new_root.children[1],new_root.probVectUpLeft)*/
                }
                
                // update tree->root;
                tree->root = new_root;
                
                // iteratively traverse the tree to update partials from the current node
                queue<Node*> node_queue;
                node_queue.push(selected_node);
                updatePartialLh(node_queue);
            }
            //add parent to non-root node
            else
            {
                Regions* upper_left_right_regions = getPartialLhAtNode(selected_node->neighbor);
                
                // now try different lengths for the new branch
                double new_branch_length_lh = best_parent_lh;
                double best_length = default_blength;
                while (best_length > min_blength)
                {
                    double new_blength = best_length / 2;
                    double placement_cost = calculatePlacementCost(best_parent_regions, sample, new_blength);
                    
                    if (placement_cost > new_branch_length_lh)
                    {
                        new_branch_length_lh = placement_cost;
                        best_length = new_blength;
                    }
                    else
                        break;
                }
                
                if (best_length > 0.7 * default_blength)
                {
                    while (best_length < max_blength)
                    {
                        double new_blength = best_length * 2;
                        double placement_cost = calculatePlacementCost(best_parent_regions, sample, new_blength);
                        
                        if (placement_cost > new_branch_length_lh)
                        {
                            new_branch_length_lh = placement_cost;
                            best_length = new_blength;
                        }
                        else
                            break;
                    }
                }
                
                // try with length zero
                if (best_length < min_blength)
                {
                    double placement_cost = calculatePlacementCost(best_parent_regions, sample, -1);
                    
                    if (placement_cost > new_branch_length_lh)
                        best_length = -1;
                }
                
                // now create new internal node and append child to it
                Node* new_internal_node = new Node(true);
                Node* next_node_1 = new Node();
                Node* next_node_2 = new Node();
                Node* new_sample_node = new Node(seq_name);
                
                new_internal_node->next = next_node_1;
                next_node_1->next = next_node_2;
                next_node_2->next = new_internal_node;
                
                double down_distance = selected_node->length * best_parent_split;
                double top_distance = selected_node->length * (1.0 - best_parent_split);
                if (best_parent_split < 0)
                {
                    delete selected_node->total_lh;
                    selected_node->total_lh = NULL;
                    down_distance = 0;
                    //down_distance = false;
                    top_distance = selected_node->length;
                    
                    /*node.probVectTotUp=None
                    node.furtherMidNodes=None*/
                }
                
                // attach to the parent node
                new_internal_node->neighbor = selected_node->neighbor;
                selected_node->neighbor->neighbor = new_internal_node;
                new_internal_node->length = top_distance;
                new_internal_node->neighbor->length = top_distance;
                
                // attach to the right child
                selected_node->neighbor = next_node_2;
                next_node_2->neighbor = selected_node;
                selected_node->length = down_distance;
                selected_node->neighbor->length = down_distance;
                
                // attach to the left child
                new_sample_node->neighbor = next_node_1;
                next_node_1->neighbor = new_sample_node;
                new_sample_node->length = best_length;
                new_sample_node->neighbor->length = best_length;
                
                new_sample_node->partial_lh = sample;
                next_node_1->partial_lh = new Regions();
                next_node_1->partial_lh->copyRegions(best_parent_regions, num_states);
                mergeLhUpDown(next_node_2->partial_lh, upper_left_right_regions, new_internal_node->length, sample, best_length);
                mergeLhTwoLower(new_internal_node->partial_lh, getPartialLhAtNode(selected_node), selected_node->length, sample, best_length);
                //new_internal_node.probVectTotUp=mergeLhUpDown(this_node_upper_left_right_regions,top_distance/2,new_internal_node.probVect,top_distance/2,mutMatrix)
                computeTotalLhAtNode(new_internal_node);
                //mergeLhUpDown(new_internal_node->total_lh, best_parent_regions, -1, sample, best_length);
                
                if (!new_internal_node->total_lh || new_internal_node->total_lh->size() == 0)
                {
                    outError("Problem, None vector when placing sample, new parent");
                    /*print(best_parent_regions)
                    print(upper_left_right_regions)
                    print(sample)
                    print(best_length)
                    print(top_distance)
                    print(down_distance)*/
                }
                
                /*if (top_distance >= 2 * min_blengthForMidNode)
                    createFurtherMidNodes(new_internal_node,this_node_upper_left_right_regions)*/
                if (best_length > 0)
                {
                    computeTotalLhAtNode(new_sample_node);
                    //mergeLhUpDown(new_sample_node->total_lh, best_parent_regions, best_length, sample, 0);
                    // new_sample_node.probVectTotUp=mergeLhUpDown(best_parent_regions,best_length/2,sample,best_length/2,mutMatrix)
                    /*if best_length>=2*min_blengthForMidNode:
                        createFurtherMidNodes(new_sample_node,best_parent_regions)*/
                }
                
                // update pseudo_count
                updatePesudoCount(best_parent_regions, sample);

                // iteratively traverse the tree to update partials from the current node
                queue<Node*> node_queue;
                node_queue.push(selected_node);
                node_queue.push(new_internal_node->neighbor);
                updatePartialLh(node_queue);
            }
        }
    }
    
    // delete best_parent_regions and best_child_regions
    if (best_parent_regions)
        delete best_parent_regions;
    if (best_child_regions)
        delete best_child_regions;
}

void CMaple::doInference()
{
    // 1. Build an initial tree
    buildInitialTree();
    
}

void CMaple::tmpTestingMethod()
{
    // do something
   /* int id = 0;
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
    tree->root = new_node;*/
    
    cout << tree->exportTreeString();
    
    cout << endl;
}

int CMaple::compareSequences(Regions* sequence1, Regions* sequence2)
{
    PositionType seq_length = tree->aln->ref_seq.size();
    StateType num_states = tree->aln->num_states;
    ASSERT(seq_length > 0);
    
    // init dummy variables
    PositionType seq1_index = -1;
    PositionType seq2_index = -1;
    bool seq1_more_info = false;
    bool seq2_more_info = false;
    PositionType pos = 0;
    Region *seq1_region, *seq2_region;
    PositionType seq1_end = -1, seq2_end = -1;
    PositionType length;

    while (pos < seq_length && (!seq1_more_info || !seq2_more_info))
    {
        // get the next shared segment in the two sequences
        getNextSharedSegment(pos, sequence1, sequence2, seq1_index, seq2_index, seq1_region, seq2_region, seq1_end, seq2_end, length);
        
        // The two regions have different types from each other
        if (seq1_region->type != seq2_region->type)
        {
            if (seq1_region->type == TYPE_N)
                seq2_more_info = true;
            else
                if (seq2_region->type == TYPE_N)
                    seq1_more_info = true;
                else if (seq1_region->type == TYPE_O)
                        seq2_more_info = true;
                    else
                        if (seq2_region->type == TYPE_O)
                            seq1_more_info = true;
                        else
                        {
                            seq1_more_info = true;
                            seq2_more_info = true;
                        }
        }
        // Both regions are type O
        else if (seq1_region->type == TYPE_O)
        {
            for (StateType i = 0; i < num_states; i++)
            {
                if (seq1_region->likelihood[i] > seq2_region->likelihood[i] + 1e-2)
                    seq2_more_info = true;
                else if (seq2_region->likelihood[i] > seq1_region->likelihood[i] + 1e-2)
                        seq1_more_info = true;
            }
        }

        // update pos
        pos += length;
    }

    // return result
    if (seq1_more_info)
        if (seq2_more_info)
            return 0;
        else
            return 1;
    else
        if (seq2_more_info)
            return -1;
        else
            return 1;
    
    return 0;
}

void CMaple::move2NextRegion(Regions* sequence, PositionType region_index, Region* &region, PositionType &current_pos, PositionType &end_pos)
{
    ASSERT(region_index < sequence->size());
    
    // get the current region
    region = sequence->getRegion(region_index);
    
    // get the current position and end position
    current_pos = region->position;
    PositionType length = (region_index < sequence->size() - 1 ? sequence->getRegion(region_index + 1)->position : tree->aln->ref_seq.size()) - current_pos;
    end_pos = current_pos + length - 1;
}

void CMaple::getNextSharedSegment(PositionType current_pos, Regions* sequence1, Regions* sequence2, PositionType &seq1_index, PositionType &seq2_index, Region* &seq1_region, Region* &seq2_region, PositionType &seq1_end_pos, PositionType &seq2_end_pos, PositionType &length)
{
    PositionType seq1_pos, seq2_pos, end_pos;
    
    // move to the next region in sequence 1
    if (current_pos > seq1_end_pos)
    {
        seq1_index++;
        move2NextRegion(sequence1, seq1_index, seq1_region, seq1_pos, seq1_end_pos);
    }
    
    // move to the next region in sequence 2
    if (current_pos > seq2_end_pos)
    {
        seq2_index++;
        move2NextRegion(sequence2, seq2_index, seq2_region, seq2_pos, seq2_end_pos);
    }
    
    // compute the end_pos for the shared segment
    end_pos = seq1_end_pos < seq2_end_pos ? seq1_end_pos : seq2_end_pos;
    length = end_pos + 1 - current_pos;
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
