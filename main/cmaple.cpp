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
    for (PositionType i = 0; i < sequence_length; i++)
        cumulative_base[i].resize(tree->aln->num_states, 0);
    
    cumulative_rate[0] = tree->model->mutation_mat[tree->aln->ref_seq[0] * (tree->aln->num_states + 1)];
    cumulative_base[0][tree->aln->ref_seq[0]] = 1;
    
    for (PositionType i = 1; i < sequence_length; i++)
    {
        StateType state = tree->aln->ref_seq[i];
        cumulative_rate[i] = cumulative_rate[i - 1] + tree->model->mutation_mat[state * (tree->aln->num_states + 1)];
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
    
    // compute thresholds for approximations
    computeThresholds();
}

double CMaple::calculatePlacementCost(Regions* parent_regions, Regions* child_regions, double blength)
{
    // init dummy variables
    double lh_cost = 0;
    PositionType seq1_index = -1;
    PositionType seq2_index = -1;
    PositionType pos = 1;
    double total_factor = 1;
    StateType num_states = tree->aln->num_states;
    Region *seq1_region, *seq2_region;
    PositionType seq1_pos, seq2_pos;
    PositionType seq1_end = -1, seq2_end = -1;
    PositionType length;
    double minimum_carry_over = DBL_MIN * 1e50;
    double total_blength = blength;
    
    while (pos < tree->aln->ref_seq.size())
    {
        // get the next shared segment in the two sequences
        getNextSharedSegment(pos, parent_regions, child_regions, seq1_index, seq2_index, seq1_region, seq2_region, seq1_end, seq2_end, length);
        
        // 1. e1.type = N || e2.type = N
        if (seq2_region->type == TYPE_N || seq1_region->type == TYPE_N)
            continue;
        else
        {
            if (seq1_region->plength_from_root > 0)
                total_blength += seq1_region->plength_from_root;
            else
                total_blength += seq1_region->plength_observation;
            
            if (seq2_region->plength_from_root > 0)
                total_blength += seq2_region->plength_from_root;
            else
                total_blength += seq2_region->plength_observation;
            
            // 2. e1.type = R
            if (seq1_region->type == TYPE_R)
            {
                // 2.1. e1.type = R and e2.type = R
                if (seq2_region->type == TYPE_R)
                {
                    total_blength = seq1_region->plength_observation;
                    if (total_blength > 0)
                        lh_cost += total_blength * (cumulative_rate[pos + length - 1] - cumulative_rate[pos]);
                }
                // 2.2. e1.type = R and e2.type = O
                else if (seq2_region->type == TYPE_O)
                {
                    double tot = 0;
                    StateType seq1_state = tree->aln->ref_seq[pos];
                    if (seq1_region->plength_from_root > 0)
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
                    
                    if (seq1_region->plength_from_root > 0)
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
                // 3.1. e1.type = O and e2.type = O
                if (seq2_region->type == TYPE_O)
                {
                    double tot = 0;
                    
                    for (StateType i = 0; i < num_states; i++)
                    {
                        double tot2 = 0;
                        if (total_blength > 0)
                        {
                            for (StateType j = 0; j < num_states; j++)
                                tot2 += tree->model->mutation_mat[i * num_states + j] * seq2_region->likelihood[j];
                        }
                        tot += seq1_region->likelihood[i] * (seq2_region->likelihood[i] + total_blength * tot2);
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
                    
                    total_factor *= seq1_region->likelihood[seq2_state] + total_blength * tot2;
                }
            }
            // 4. e1.type = A/C/G/T
            else
            {
                int seq1_state = seq1_region->type;
                
                // 4.1. e1.type =  e2.type
                if (seq1_region->type == seq2_region->type)
                {
                    if (seq1_region->plength_from_root > 0)
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
                        
                        if (seq1_region->plength_from_root > 0)
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
                        
                        if (seq1_region->plength_from_root > 0)
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
            lh_diff_mid_branch = -DBL_MAX;*/

        // try to place as descendant of the current node (this is skipped if the node has top branch length 0 and so is part of a polytomy).
        if (current_node->length > 0)
        {
            // compute the placement cost to place the new sample
            lh_diff_at_node = calculatePlacementCost(current_node->total_lh, regions, default_blength);
                
            if (lh_diff_at_node > best_lh_diff)
            {
                best_lh_diff = lh_diff_at_node;
                selected_node = current_node;
                current_node->double_attributes[FAILURE_COUNT] = 0;
                is_mid_branch = false;
                //best_up_lh_diff = lh_diff_mid_branch;
            }
            /*else if (lh_diff_mid_branch >= (best_lh_diff - params->fastLK_threshold_prob))
            {
                best_up_lh_diff = visit_node.parent_lh;
                best_down_lh_diff = lh_diff_at_node;
                selected_node = visit_node.node;
            }*/
            // placement at current node is considered failed if placement likelihood is not improved by a certain margin compared to best placement so far for the nodes above it.
            else if (lh_diff_at_node < (current_node->double_attributes[LH_DIFF] - tree->params->thresh_log_lh_failure))
                current_node->double_attributes[FAILURE_COUNT] = current_node->double_attributes[FAILURE_COUNT] + 1;
        }
        else
            lh_diff_at_node = current_node->double_attributes[LH_DIFF];
        
        // keep trying to place at children nodes, unless the number of attempts has reaches the failure limit
        if (tree->params->strict_stop_seeking_placement)
        {
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

                while (true)
                {
                    Regions* current_node_regions = new Regions();
                    mergeLhUpDown(current_node_regions, node->neighbor->partial_lh, new_blength, node->partial_lh, node->length - new_blength);
                    
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
    PositionType pos = 1;
    StateType num_states = tree->aln->num_states;
    Region *seq1_region, *seq2_region;
    PositionType seq1_pos, seq2_pos;
    PositionType seq1_end = -1, seq2_end = -1;
    PositionType length;
                
    while (pos <=  tree->aln->ref_seq.size())
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
                    double total_blength = lower_plength + seq2_region->plength_observation;
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
                    merged_regions->push_back(new Region(seq2_region->type, pos, seq2_region->plength_observation + lower_plength));
            }
        }
        // seq2_entry = 'N'
        else if (seq2_region->type == TYPE_N)
        {
            // seq1_entry = 'O' and seq2_entry = N
            if (seq1_region->type == TYPE_O)
            {
                double total_blength = upper_plength + seq1_region->plength_observation;
                
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
            // seq1_entry = 'N' and seq2_entry = R/ACGT
            else
            {
                if (seq1_region->plength_from_root > 0)
                    merged_regions->push_back(new Region(seq1_region->type, pos, seq1_region->plength_observation, seq1_region->plength_from_root + upper_plength));
                else
                    merged_regions->push_back(new Region(seq1_region->type, pos, seq1_region->plength_observation + upper_plength));
            }
        }
        // seq1_entry = seq2_entry = R/ACGT
        else if (seq1_region->type == seq2_region->type && (seq1_region->type < num_states || seq1_region->type == TYPE_R))
            merged_regions->push_back(new Region(seq1_region->type, pos));
        // cases where the new genome list entry will likely be of type "O"
        else
        {
            double total_blength_1 = upper_plength + seq1_region->plength_observation;
            if (seq1_region->type != TYPE_O)
                total_blength_1 += seq1_region->plength_from_root;
            
            double total_blength_2 = lower_plength + seq2_region->plength_observation;
            
            // due to 0 distance, the entry will be of same type as entry2
            if ((seq2_region->type < num_states || seq2_region->type == TYPE_R) && total_blength_2 == 0)
            {
                if ((seq1_region->type < num_states || seq1_region->type == TYPE_R) && total_blength_1 == 0)
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
                
                if (seq1_region->plength_from_root > 0)
                {
                    double length_to_root = seq1_region->plength_from_root + upper_plength;
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

void CMaple::mergeLhTwoLower(Regions* &merged_regions, Regions* regions1, double plength1, Regions* regions2, double plength2, double log_lh, bool return_log_lh)
{
    // init variables
    PositionType seq1_index = -1;
    PositionType seq2_index = -1;
    PositionType pos = 1;
    StateType num_states = tree->aln->num_states;
    Region *seq1_region, *seq2_region;
    PositionType seq1_pos, seq2_pos;
    PositionType seq1_end = -1, seq2_end = -1;
    PositionType length;
                
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
                    new_region->plength_observation += plength2;
                    merged_regions->push_back(new_region);
                }
                // seq1_entry = 'N' and seq2_entry = R/ACGT
                else
                    merged_regions->push_back(new Region(seq2_region->type, pos, seq2_region->plength_observation + plength2));
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
                new_region->plength_observation += plength1;
                merged_regions->push_back(new_region);
            }
            // seq1_entry = 'N' and seq2_entry = R/ACGT
            else
                merged_regions->push_back(new Region(seq1_region->type, pos, seq1_region->plength_observation + plength1));
        }
        // neither seq1_entry nor seq2_entry = N
        else
        {
            double total_blength_1 = plength1 + seq1_region->plength_observation;
            double total_blength_2 = plength2 + seq2_region->plength_observation;
            
            // seq1_entry and seq2_entry are identical seq1_entry = R/ACGT
            if (seq1_region->type == seq2_region->type && (seq1_region->type == TYPE_R || seq1_region->type < num_states))
            {
                merged_regions->push_back(new Region(seq1_region->type, pos));
                
                if (return_log_lh)
                {
                    if (seq1_region->type == TYPE_R)
                        log_lh += (total_blength_1 + total_blength_2) * (cumulative_rate[pos + length - 1] - cumulative_rate[pos]);
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

void CMaple::addRoot()
{
    Sequence* root_sequence = tree->aln->at(0);
    Node* root = new Node(0, root_sequence->seq_name);
    getLowerLhSeqAtTip(root->partial_lh, root_sequence->mutations, 0);
    root->computeTotalLhForRoot(tree->model, tree->aln->num_states);
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
        {
            Node* new_root = placeNewSample(selected_node, lower_lh_seq, sequence->seq_name, best_lh_diff, is_mid_branch, best_up_lh_diff, best_down_lh_diff, best_child, adjust_blen);
            if (new_root)
                tree->root = new_root;
        }
        
        // delete lower_lh_seq
        delete lower_lh_seq;
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
                log_lh += tree->model->root_log_freqs[i] * (cumulative_base[end_pos][i] - cumulative_base[start_pos][i]);
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

Node* CMaple::placeNewSample(Node* selected_node, Regions* sample, string seq_name, double best_lh_diff , bool is_mid_branch, double best_up_lh_diff, double best_down_lh_diff, Node* best_child, bool adjust_blen)
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
    
    // temporarily ignore placing the new sample as a descendant of a mid-branch point
    if (is_mid_branch)
        return NULL;
    
    // place the new sample as a descendant of an existing node
    if (best_child)
    {
        double best_split = 0.5;
        double best_split_lh = best_child_lh;
        best_child_regions = new Regions();
        mergeLhUpDown(best_child_regions, best_child->neighbor->partial_lh, best_child->length / 2, best_child->partial_lh, best_child->length / 2);
        double new_split = 0.25;
        
        while (new_split * best_child->length > min_blength)
        {
            Regions* new_parent_regions = new Regions();
            mergeLhUpDown(new_parent_regions, best_child->neighbor->partial_lh, best_child->length * new_split, best_child->partial_lh, best_child->length * (1 - new_split));
            
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
    if (tree->root == selected_node)
    {
        double old_root_lh = computeLhAtRoot(selected_node->partial_lh);
        double new_root_lh;
        Regions* merged_root_sample_regions = new Regions();
        
        // merge 2 lower vector into one
        mergeLhTwoLower(merged_root_sample_regions, selected_node->partial_lh, default_blength, sample, default_blength, new_root_lh, true);
        
        new_root_lh += computeLhAtRoot(merged_root_sample_regions);
        best_parent_lh = new_root_lh - old_root_lh;
        best_root_blength = default_blength;
        best_parent_regions->copyRegions(merged_root_sample_regions, num_states);
        double new_blength = 0.5 * default_blength;
        
        while (new_blength > min_blength)
        {
            // merge 2 lower vector into one
            mergeLhTwoLower(merged_root_sample_regions, selected_node->partial_lh, new_blength, sample, default_blength, new_root_lh, true);
            
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
        mergeLhUpDown(best_parent_regions, selected_node->neighbor->partial_lh, best_child->length / 2, selected_node->partial_lh, best_child->length / 2);
        double new_split = 0.25;
        
        while (new_split * selected_node->length > min_blength)
        {
            Regions* new_parent_regions = new Regions();
            mergeLhUpDown(new_parent_regions, best_child->neighbor->partial_lh, selected_node->length * (1 - new_split), selected_node->partial_lh, best_child->length * new_split);
            
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
            
            // delete new_parent_regions
            delete new_parent_regions;
        }
        
        best_parent_lh = best_split_lh;
        best_parent_split = best_split;
    }
    
    /*// if the best placement is below the selected_node => add an internal node below the selected_node
    if bestChildLK>=best_parent_lh and bestChildLK>=newChildLK
    {
        if bestDownNode==bestDownNode.up.children[0]:
            child=0
            vectUp=bestDownNode.up.probVectUpRight
        else:
            child=1
            vectUp=bestDownNode.up.probVectUpLeft

        LK1=bestChildLK
        bestLen=oneMutBLen*factor
        while bestLen>minBLen:
            new_blengthen=bestLen/2
            placement_cost=appendProb(best_child_regions,newPartials,new_blengthen,mutMatrix)
            if placement_cost>LK1:
                LK1=placement_cost
                bestLen=new_blengthen
            else:
                break
        if bestLen>0.7*oneMutBLen*factor:
            while bestLen<maxBLen:
                new_blengthen=bestLen*2
                placement_cost=appendProb(best_child_regions,newPartials,new_blengthen,mutMatrix)
                if placement_cost>LK1:
                    LK1=placement_cost
                    bestLen=new_blengthen
                else:
                    break
        if bestLen<minBLen:
            LK0=appendProb(best_child_regions,newPartials,False,mutMatrix)
            if LK0>LK1:
                bestLen=False
        #now create new internal node and append child to it
        newInternalNode=Tree()
        bestDownNode.up.children[child]=newInternalNode
        newInternalNode.up=bestDownNode.up
        distBottom=bestDownNode.dist*(1.0-best_child_split)
        distTop=bestDownNode.dist*best_child_split
        bestDownNode.up=newInternalNode
        bestDownNode.dist=distBottom
        newInternalNode.add_child(bestDownNode)
        newNode=Tree(name=sample,dist=bestLen)
        newNode.minorSequences=[]
        newNode.up=newInternalNode
        newInternalNode.add_child(newNode)
        newInternalNode.dist=distTop
        newInternalNode.children[1].probVect=newPartials
        newInternalNode.probVectUpRight=mergeVectorsUpDown(vectUp,newInternalNode.dist,newPartials,bestLen,mutMatrix)
        newInternalNode.probVectUpLeft=best_child_regions
        newInternalNode.probVect=mergeVectors(bestDownNode.probVect,bestDownNode.dist,newPartials,bestLen,mutMatrix)
        newInternalNode.probVectTotUp=mergeVectorsUpDown(vectUp,distTop/2,newInternalNode.probVect,distTop/2,mutMatrix)
        newInternalNode.probVectTot=mergeVectorsUpDown(best_child_regions,False,newPartials,bestLen,mutMatrix)
        if newInternalNode.probVectTot==None:
                print("Problem, None vector when placing sample, below node")
                print(best_child_regions)
                print(vectUp)
                print(newPartials)
                print(bestLen)
                print(distTop)
                print(distBottom)
        if distTop>=2*minBLenForMidNode:
            createFurtherMidNodes(newInternalNode,vectUp)
        if bestLen:
            newNode.probVectTot=mergeVectorsUpDown(best_child_regions,bestLen,newPartials,False,mutMatrix)
            newNode.probVectTotUp=mergeVectorsUpDown(best_child_regions,bestLen/2,newPartials,bestLen/2,mutMatrix)
            if bestLen>=2*minBLenForMidNode:
                createFurtherMidNodes(newNode,best_child_regions)
        updatePesudoCounts(best_child_regions,newPartials,pseudoMutCounts)
        if verbose:
            print("new internal node added to tree")
            print(newInternalNode.probVect)
            print(newInternalNode.probVectUpRight)
            print(newInternalNode.probVectUpLeft)
            print(newInternalNode.probVectTot)
        #updatePartialsFromTop(bestDownNode,newInternalNode.probVectUpRight,mutMatrix)
        #updatePartialsFromBottom(newInternalNode.up,newInternalNode.probVect,child,newInternalNode,mutMatrix)
        nodeList=[(bestDownNode,2),(newInternalNode.up,child)]
        updatePartials(nodeList,mutMatrix)
    }
    // otherwise, add new parent to the selected_node
    else
    {
        #new parent is actually part of a polytomy since best placement is exactly at the node
        if newChildLK>=best_parent_lh:
            best_root_blength=False
            best_parent_split=False
            best_parent_lh=newChildLK
            best_parent_regions->copyRegions(node.probVectTot, num_states);
            if node.up==None:
                best_parent_regions,new_root_lh = mergeVectors(node.probVect,False,newPartials,oneMutBLen*factor,mutMatrix,returnLK=True)

        #add parent to the root
        if node.up==None:
            #now try different lengths for right branch
            bestLen2=oneMutBLen*factor
            while bestLen2>minBLen:
                new_blengthen=bestLen2/2
                newProbVectRoot,newProbRoot = mergeVectors(node.probVect,best_root_blength,newPartials,new_blengthen,mutMatrix,returnLK=True)
                newProbRoot+= computeLhAtRoot(newProbVectRoot)
                LKdiffRoot=newProbRoot-old_root_lh
                if LKdiffRoot>best_parent_lh:
                    best_parent_lh=LKdiffRoot
                    bestLen2=new_blengthen
                    best_parent_regions->copyRegions(newProbVectRoot, num_states);
                else:
                    break
            if bestLen2>0.7*oneMutBLen*factor:
                while bestLen2<maxBLen:
                    new_blengthen=bestLen2*2
                    newProbVectRoot,newProbRoot = mergeVectors(node.probVect,best_root_blength,newPartials,new_blengthen,mutMatrix,returnLK=True)
                    newProbRoot+= computeLhAtRoot(newProbVectRoot)
                    LKdiffRoot=newProbRoot-old_root_lh
                    if LKdiffRoot>best_parent_lh:
                        best_parent_lh=LKdiffRoot
                        bestLen2=new_blengthen
                        best_parent_regions->copyRegions(newProbVectRoot, num_states);
                    else:
                        break
            if bestLen2<minBLen:
                newProbVectRoot,newProbRoot = mergeVectors(node.probVect,best_root_blength,newPartials,False,mutMatrix,returnLK=True)
                newProbRoot+= computeLhAtRoot(newProbVectRoot)
                LK0=newProbRoot-old_root_lh
                if LK0>best_parent_lh:
                    bestLen2=False
                    best_parent_lh=LK0
                    best_parent_regions->copyRegions(newProbVectRoot, num_states);

            newRoot=Tree()
            newRoot.probVect=best_parent_regions
            newRoot.probVectTot=rootVector(best_parent_regions,False,mutMatrix)
            newRoot.probVectUpRight=rootVector(newPartials,bestLen2,mutMatrix)
            newRoot.probVectUpLeft=rootVector(node.probVect,best_root_blength,mutMatrix)
            if newRoot.probVectTot==None:
                print("Problem, None vector when placing sample, new root")
                print(merged_root_sample_regions)
                print(node.probVect)
                print(newPartials)
                print(bestLen2)
                print(best_root_blength)
            node.up=newRoot
            node.dist=best_root_blength
            if not best_root_blength:
                node.probVectTot=None
                node.probVectTotUp=None
                node.furtherMidNodes=None
            newRoot.add_child(node)
            newNode=Tree(name=sample,dist=bestLen2)
            newNode.minorSequences=[]
            newNode.up=newRoot
            newRoot.add_child(newNode)
            newNode.probVect=newPartials
            if bestLen2:
                newNode.probVectTot=mergeVectorsUpDown(newRoot.probVectUpLeft,bestLen2,newPartials,False,mutMatrix)
                newNode.probVectTotUp=mergeVectorsUpDown(newRoot.probVectUpLeft,bestLen2/2,newPartials,bestLen2/2,mutMatrix)
                if bestLen2>=2*minBLenForMidNode:
                    createFurtherMidNodes(newRoot.children[1],newRoot.probVectUpLeft)
            if verbose:
                print("new root added to tree")
                print(newRoot.probVect)
                print(newRoot.children[0].probVect)
                print(newNode.probVect)
            #updatePartialsFromTop(node,newRoot.probVectUpRight,mutMatrix)
            nodeList=[(node,2)]
            updatePartials(nodeList,mutMatrix)
            return newRoot

        #add parent to non-root node
        else:
            if node==node.up.children[0]:
                child=0
                vectUp=node.up.probVectUpRight
            else:
                child=1
                vectUp=node.up.probVectUpLeft
            
            #now try different lengths for the new branch
            LK1=best_parent_lh
            bestLen=oneMutBLen*factor
            while bestLen>minBLen:
                new_blengthen=bestLen/2
                placement_cost=appendProb(best_parent_regions,newPartials,new_blengthen,mutMatrix)
                if placement_cost>LK1:
                    LK1=placement_cost
                    bestLen=new_blengthen
                else:
                    break
            if bestLen>0.7*oneMutBLen*factor:
                while bestLen<maxBLen:
                    new_blengthen=bestLen*2
                    placement_cost=appendProb(best_parent_regions,newPartials,new_blengthen,mutMatrix)
                    if placement_cost>LK1:
                        LK1=placement_cost
                        bestLen=new_blengthen
                    else:
                        break
            if bestLen<minBLen:
                LK0=appendProb(best_parent_regions,newPartials,False,mutMatrix)
                if LK0>LK1:
                    bestLen=False
            #now create new internal node and append child to it
            newInternalNode=Tree()
            node.up.children[child]=newInternalNode
            newInternalNode.up=node.up
            if best_parent_split:
                distBottom=node.dist*best_parent_split
                distTop=node.dist*(1.0-best_parent_split)
            else:
                distBottom=False
                distTop=node.dist
                node.probVectTot=None
                node.probVectTotUp=None
                node.furtherMidNodes=None
            node.dist=distBottom
            node.up=newInternalNode
            newInternalNode.add_child(node)
            newNode=Tree(name=sample,dist=bestLen)
            newNode.minorSequences=[]
            newNode.up=newInternalNode
            newInternalNode.add_child(newNode)
            newInternalNode.dist=distTop
            newInternalNode.children[1].probVect=newPartials
            newInternalNode.probVectUpRight=mergeVectorsUpDown(vectUp,newInternalNode.dist,newPartials,bestLen,mutMatrix)
            newInternalNode.probVectUpLeft=best_parent_regions
            newInternalNode.probVect=mergeVectors(node.probVect,node.dist,newPartials,bestLen,mutMatrix)
            newInternalNode.probVectTotUp=mergeVectorsUpDown(vectUp,distTop/2,newInternalNode.probVect,distTop/2,mutMatrix)
            newInternalNode.probVectTot=mergeVectorsUpDown(best_parent_regions,False,newPartials,bestLen,mutMatrix)
            if newInternalNode.probVectTot==None:
                print("Problem, None vector when placing sample, new parent")
                print(best_parent_regions)
                print(vectUp)
                print(newPartials)
                print(bestLen)
                print(distTop)
                print(distBottom)
            if distTop>=2*minBLenForMidNode:
                createFurtherMidNodes(newInternalNode,vectUp)
            if bestLen:
                newNode.probVectTot=mergeVectorsUpDown(best_parent_regions,bestLen,newPartials,False,mutMatrix)
                newNode.probVectTotUp=mergeVectorsUpDown(best_parent_regions,bestLen/2,newPartials,bestLen/2,mutMatrix)
                if bestLen>=2*minBLenForMidNode:
                    createFurtherMidNodes(newNode,best_parent_regions)
            updatePesudoCounts(best_parent_regions,newPartials,pseudoMutCounts)
            if verbose:
                print("new internal node added to tree")
                print(newInternalNode.probVect)
                print(newInternalNode.probVectUpRight)
                print(newInternalNode.probVectUpLeft)
                print(newInternalNode.probVectTot)
            #updatePartialsFromTop(node,newInternalNode.probVectUpRight,mutMatrix)
            #updatePartialsFromBottom(newInternalNode.up,newInternalNode.probVect,child,newInternalNode,mutMatrix)
            nodeList=[(node,2),(newInternalNode.up,child)]
            updatePartials(nodeList,mutMatrix)
    }*/
    
    // delete best_parent_regions and best_child_regions
    if (best_parent_regions)
        delete best_parent_regions;
    if (best_child_regions)
        delete best_child_regions;
        
    return NULL;
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
        
    Node* next;
    FOR_NEXT(node, next)
        traverseTree(next);
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
