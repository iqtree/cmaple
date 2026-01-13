#include "model_dna_rate_variation.h"
#include "../tree/tree.h"
#include "../tree/phylonode.h"

using namespace cmaple;

ModelDNARateVariation::ModelDNARateVariation( const cmaple::ModelBase::SubModel sub_model, PositionType _genome_size, 
                                                bool _use_site_rates, cmaple::RealNumType _wt_pseudocount, std::string _rates_filename)
    : ModelDNA(sub_model) {
    
    genome_size = _genome_size;
    use_site_rates = _use_site_rates;
    mat_size = row_index[num_states_];
    waiting_time_pseudocount = _wt_pseudocount;
    rates_filename = _rates_filename;

    mutation_matrices = new RealNumType[mat_size * genome_size]();
    transposed_mutation_matrices = new RealNumType[mat_size * genome_size]();
    diagonal_mutation_matrices = new RealNumType[num_states_ * genome_size]();
    freqi_freqj_Qijs = new RealNumType[mat_size * genome_size]();
    freqj_transposedijs = new RealNumType[mat_size * genome_size]();

    if(use_site_rates) {
        rates = new cmaple::RealNumType[genome_size]();
    }
    if(rates_filename.length() > 0) {
        this->readRatesFile();
    }
}

ModelDNARateVariation::~ModelDNARateVariation() { 
    delete[] mutation_matrices;
    delete[] transposed_mutation_matrices;
    delete[] diagonal_mutation_matrices;
    delete[] freqi_freqj_Qijs;
    delete[] freqj_transposedijs;
    if(use_site_rates) {
        delete[] rates;
    }
}

void ModelDNARateVariation::printMatrix(const RealNumType* matrix, std::ostream* out_stream) {
    for(int j = 0; j < num_states_; j++) {
        std::string line = "|";
        for(int k = 0; k < num_states_; k++) {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(5) << matrix[row_index[j] + k];
            line += "\t" + oss.str();
            //line += "\t" + convertDoubleToString(matrix[row_index[j] + k], 5);
        }
        line += "\t|";
        (*out_stream) << line << std::endl;
    }
}

void ModelDNARateVariation::printCountsAndWaitingTimes(const RealNumType* counts, const RealNumType* waiting_times, std::ostream* out_stream) {
    for(int j = 0; j < num_states_; j++) {
        std::string line = "|";
        for(int k = 0; k < num_states_; k++) {
            line += "\t" + convertDoubleToString(counts[row_index[j] + k]);
        }
        line += "\t|\t" + convertDoubleToString(waiting_times[j]) + "\t|";
        (*out_stream) << line << std::endl;
    }
}

bool ModelDNARateVariation::updateMutationMatEmpirical() {
    if(rates_filename.length() > 0) {
        return true;
    }
    if(rates_estimated) {
        std::cout << "[ModelDNARateVariation] Warning: Overwriting estimated rate matrices with single empirical mutation matrix." << std::endl;
    }
    bool val = ModelDNA::updateMutationMatEmpirical();
    for(int i = 0; i < genome_size; i++) {
        for(int j = 0; j < mat_size; j++) {
            mutation_matrices[i * mat_size + j] = mutation_mat[j];
            transposed_mutation_matrices[i * mat_size + j] = transposed_mut_mat[j];
            freqi_freqj_Qijs[i * mat_size + j] = freqi_freqj_qij[j];
            freqj_transposedijs[i * mat_size + j] = freq_j_transposed_ij[j];
        }
        for(int j = 0; j < num_states_; j++) {
            diagonal_mutation_matrices[i * num_states_ + j] = diagonal_mut_mat[j];
        }
    }
    return val;
}

void ModelDNARateVariation::estimateRates(cmaple::Tree* tree) {
    rates_estimated = true;
    if(use_site_rates) {
        estimateRatePerSite(tree);

    } else {
        if(rates_filename.size() == 0) {
            RealNumType old_LK = -std::numeric_limits<double>::infinity();
            RealNumType new_LK = tree->computeLh();
            int num_steps = 0;
            while(new_LK - old_LK > 1 && num_steps < 20) {
                estimateRatesPerSitePerEntry(tree);
                old_LK = new_LK;
                new_LK = tree->computeLh();
            }  
        }
    }

    // Write out rate matrices to file
    if(cmaple::verbose_mode > VB_MIN) 
    {
        const std::string prefix = tree->params->output_prefix.length() ? 
            tree->params->output_prefix : tree->params->aln_path;
        //std::cout << "Writing rate matrices to file " << prefix << ".rateMatrices.txt" << std::endl;
        std::ofstream out_file(prefix + ".rateMatrices.txt");
        out_file << "Rate matrix for all sites: " << std::endl;
        printMatrix(getOriginalRateMatrix(), &out_file);
        for(int i = 0; i < genome_size; i++) {
            out_file << "Position: " << i << std::endl;
            if(use_site_rates) {
                out_file << "Rate: " << rates[i] << std::endl;
            }
            out_file << "Rate Matrix: " << std::endl;
            printMatrix(getMutationMatrix(i), &out_file);
            out_file << std::endl;

        }
        out_file.close();
    } 
}

void ModelDNARateVariation::estimateRatePerSite(cmaple::Tree* tree){
    std::cout << "Estimating mutation rate per site..." << std::endl;
    RealNumType* waiting_times = new RealNumType[num_states_ * genome_size];
    RealNumType* num_substitutions = new RealNumType[genome_size];
    for(int i = 0; i < genome_size; i++) {
        for(int j = 0; j < num_states_; j++) {
            waiting_times[i * num_states_ + j] = 0;
        }
        num_substitutions[i] = 0;
    }

    std::stack<Index> node_stack;
    const PhyloNode& root = tree->nodes[tree->root_vector_index];
    node_stack.push(root.getNeighborIndex(RIGHT));
    node_stack.push(root.getNeighborIndex(LEFT));
    while(!node_stack.empty()) {

        Index index = node_stack.top();
        node_stack.pop();
        PhyloNode& node = tree->nodes[index.getVectorIndex()];
        RealNumType blength = node.getUpperLength();
        //std::cout << "blength: " << blength  << std::endl;

        if (node.isInternal()) {
            node_stack.push(node.getNeighborIndex(RIGHT));
            node_stack.push(node.getNeighborIndex(LEFT));
        }

        if(blength <= 0.) {
            continue;
        }

        Index parent_index = node.getNeighborIndex(TOP);
        PhyloNode& parent_node = tree->nodes[parent_index.getVectorIndex()];
        // const std::unique_ptr<SeqRegions>& parent_regions = parent_node.getPartialLh(parent_index.getMiniIndex());
        
        // use the same local ref as the current/child node
        // 0. extract the mutations at the selected node
        std::unique_ptr<SeqRegions>& selected_node_mutations =
            tree->node_mutations[index.getVectorIndex()];
        // 1. create a new upper_lr_regions that integrate the mutations, if any
        std::unique_ptr<SeqRegions> mut_integrated_upper_lr_regions =
            (selected_node_mutations && selected_node_mutations->size())
            ? parent_node.getPartialLh(parent_index.getMiniIndex())
              ->integrateMutations<4>(selected_node_mutations, tree->aln)
            : nullptr;
        // 2. create the pointer that points to the appropriate upper_lr_regions
        const std::unique_ptr<SeqRegions>* upper_lr_regions_ptr =
            (selected_node_mutations && selected_node_mutations->size())
            ? &(mut_integrated_upper_lr_regions)
            : &(parent_node.getPartialLh(parent_index.getMiniIndex()));
        // 3. create a reference from that pointer
        auto& parent_regions = *upper_lr_regions_ptr;
        
        const std::unique_ptr<SeqRegions>& child_regions = node.getPartialLh(TOP);

        PositionType pos = 0;
        const SeqRegions& seq1_regions = *parent_regions;
        const SeqRegions& seq2_regions = *child_regions;
        size_t iseq1 = 0;
        size_t iseq2 = 0;

        while(pos < genome_size) {
            PositionType end_pos;
            SeqRegions::getNextSharedSegment(pos, seq1_regions, seq2_regions, iseq1, iseq2, end_pos);
            const auto* seq1_region = &seq1_regions[iseq1];
            const auto* seq2_region = &seq2_regions[iseq2];

            if(seq1_region->type == TYPE_R && seq2_region->type == TYPE_R) {
                // both states are type REF
                for(int i = pos; i <= end_pos; i++) {
                    int state = tree->aln->ref_seq[static_cast<std::vector<cmaple::StateType>::size_type>(i)];
                    waiting_times[i * num_states_ + state] += blength;
                }
            }  else if(seq1_region->type == seq2_region->type && seq1_region->type < TYPE_R) {
                // both states are equal but not of type REF
                 for(int i = pos; i <= end_pos; i++) {
                    waiting_times[i * num_states_ + seq1_region->type] += blength;
                }               
            } else if(seq1_region->type <= TYPE_R && seq2_region->type <= TYPE_R) {
                // both states are not equal
                  for(int i = pos; i <= end_pos; i++) {
                    num_substitutions[i] += 1;
                }                 
            }
            pos = end_pos + 1;
        }
    }

    RealNumType rate_count = 0;
    for(int i = 0; i < genome_size; i++) {
        if(num_substitutions[i] == 0) {
            rates[i] = 0.001;
        } else {
            RealNumType expected_rate_no_substitution = 0;
            for(int j = 0; j < num_states_; j++) {
                RealNumType summand = waiting_times[i * num_states_ + j] * abs(diagonal_mut_mat[j]);
                expected_rate_no_substitution += summand;
            }
            if(expected_rate_no_substitution <= 0.01) {
                rates[i] = 1.;
            } else {
                rates[i] = num_substitutions[i] / expected_rate_no_substitution; 
            }
        }
        rate_count += rates[i];
    }

    // RealNumType average_rate = rate_count / genome_size;
    RealNumType inverse_average_rate = genome_size / rate_count;
    for(int i = 0; i < genome_size; i++) {
        // rates[i] /= average_rate;
        rates[i] *= inverse_average_rate;
        rates[i] = std::min(100.0, std::max(0.0001, rates[i]));
        for(int stateA = 0; stateA < num_states_; stateA++) {
            RealNumType row_sum = 0;
            for(int stateB = 0; stateB < num_states_; stateB++) {
                if(stateA != stateB) {
                    RealNumType val = mutation_matrices[i * mat_size + (stateB + row_index[stateA])] * rates[i];
                    mutation_matrices[i * mat_size + (stateB + row_index[stateA])] = val;
                    transposed_mutation_matrices[i * mat_size + (stateA + row_index[stateB])] = val;
                    freqi_freqj_Qijs[i * mat_size + (stateB + row_index[stateA])] = root_freqs[stateA] * inverse_root_freqs[stateB] * val;
                    row_sum += val;
                }
            }
            mutation_matrices[i * mat_size + (stateA + row_index[stateA])] = -row_sum;
            transposed_mutation_matrices[i * mat_size + (stateA + row_index[stateA])] = -row_sum;
            diagonal_mutation_matrices[i * num_states_ + stateA] = -row_sum;
            freqi_freqj_Qijs[i * mat_size + (stateA + row_index[stateA])] = -row_sum;

            // pre-compute matrix to speedup
            const RealNumType* transposed_mut_mat_row = getTransposedMutationMatrixRow(stateA, i);
            RealNumType* freqj_transposedijs_row = freqj_transposedijs + (i * mat_size) + row_index[stateA];
            setVecByProduct<4>(freqj_transposedijs_row, root_freqs, transposed_mut_mat_row);
        
        }
    }

    delete[] waiting_times;
    delete[] num_substitutions;
}

void ModelDNARateVariation::estimateRatesPerSitePerEntry(cmaple::Tree* tree) {

    RealNumType* C = new RealNumType[genome_size * mat_size];
    RealNumType* W = new RealNumType[genome_size * num_states_];
    for(int i = 0; i < genome_size; i++) {
        for(int j = 0; j < num_states_; j++) {
            W[i * num_states_ + j] = 0;
            for(int k = 0; k < num_states_; k++) {
                C[i * num_states_ + row_index[j] + k] = 0;
            }
        }
    }
    std::stack<Index> node_stack;
    const PhyloNode& root = tree->nodes[tree->root_vector_index];
    node_stack.push(root.getNeighborIndex(RIGHT));
    node_stack.push(root.getNeighborIndex(LEFT));
    while(!node_stack.empty()) {

        Index index = node_stack.top();
        node_stack.pop();
        PhyloNode& node = tree->nodes[index.getVectorIndex()];
        RealNumType blength = node.getUpperLength();

        if (node.isInternal()) {
            node_stack.push(node.getNeighborIndex(RIGHT));
            node_stack.push(node.getNeighborIndex(LEFT));
        }

        if(blength <= 0.) {
            continue;
        }

        Index parent_index = node.getNeighborIndex(TOP);
        PhyloNode& parent_node = tree->nodes[parent_index.getVectorIndex()];
        // const std::unique_ptr<SeqRegions>& parent_regions = parent_node.getPartialLh(parent_index.getMiniIndex());
        
        // use the same local ref as the current/child node
        // 0. extract the mutations at the selected node
        std::unique_ptr<SeqRegions>& selected_node_mutations =
            tree->node_mutations[index.getVectorIndex()];
        // 1. create a new upper_lr_regions that integrate the mutations, if any
        std::unique_ptr<SeqRegions> mut_integrated_upper_lr_regions =
            (selected_node_mutations && selected_node_mutations->size())
            ? parent_node.getPartialLh(parent_index.getMiniIndex())
              ->integrateMutations<4>(selected_node_mutations, tree->aln)
            : nullptr;
        // 2. create the pointer that points to the appropriate upper_lr_regions
        const std::unique_ptr<SeqRegions>* upper_lr_regions_ptr =
            (selected_node_mutations && selected_node_mutations->size())
            ? &(mut_integrated_upper_lr_regions)
            : &(parent_node.getPartialLh(parent_index.getMiniIndex()));
        // 3. create a reference from that pointer
        auto& parent_regions = *upper_lr_regions_ptr;
        
        const std::unique_ptr<SeqRegions>& child_regions = node.getPartialLh(TOP);

        PositionType pos = 0;
        const SeqRegions& seqP_regions = *parent_regions;
        const SeqRegions& seqC_regions = *child_regions;
        size_t iseq1 = 0;
        size_t iseq2 = 0;

        while(pos < genome_size) {
            PositionType end_pos;
            SeqRegions::getNextSharedSegment(pos, seqP_regions, seqC_regions, iseq1, iseq2, end_pos);
            const auto* seqP_region = &seqP_regions[iseq1];
            const auto* seqC_region = &seqC_regions[iseq2];

            // if the child of this branch does not observe its state directly then 
            // skip this branch.
            if(seqC_region->plength_observation2node > 0) {
                pos = end_pos + 1;
                continue;
            }

            // distance to last observation or root if last observation was across the root.
            RealNumType branch_length_to_observation = blength;
            if(seqP_region->plength_observation2node > 0 && seqP_region->plength_observation2root <= 0) {
                branch_length_to_observation = blength + seqP_region->plength_observation2node;
            }
            else if(seqP_region->plength_observation2root > 0) {
                branch_length_to_observation = blength + seqP_region->plength_observation2root;
            }

            if(seqP_region->type == TYPE_R && seqC_region->type == TYPE_R) {
                // both states are type REF
                for(int i = pos; i <= end_pos; i++) {
                    StateType state = tree->aln->ref_seq[static_cast<std::vector<cmaple::StateType>::size_type>(i)];
                    W[i * num_states_ + state] += branch_length_to_observation;
                }
            }  else if(seqP_region->type == seqC_region->type && seqP_region->type < TYPE_R) {
                // both states are equal but not of type REF
                 for(int i = pos; i <= end_pos; i++) {
                    W[i * num_states_ + seqP_region->type] += branch_length_to_observation;
                }               
            } else if(seqP_region->type <= TYPE_R && seqC_region->type <= TYPE_R) {
                //states are not equal
                StateType stateA = seqP_region->type;
                StateType stateB = seqC_region->type;
                if(seqP_region->type == TYPE_R) {
                    // stateA = tree->aln->ref_seq[static_cast<std::vector<cmaple::StateType>::size_type>(end_pos)];
                    stateA = seqC_region->prev_state;
                }
                if(seqC_region->type == TYPE_R) {
                    // stateB = tree->aln->ref_seq[static_cast<std::vector<cmaple::StateType>::size_type>(end_pos)];
                    stateB = seqP_region->prev_state;
                }
                // Case 1: Last observation was this side of the root node
                if(seqP_region->plength_observation2root <= 0) {
                    for(int i = pos; i <= end_pos; i++) {
                        W[i * num_states_ + stateA] += branch_length_to_observation/2;
                        W[i * num_states_ + stateB] += branch_length_to_observation/2;
                        C[i * mat_size + stateB + row_index[stateA]] += 1;
                    }
                } else {
                    // Case 2: Last observation was the other side of the root.
                    // In this case there are two further cases - the mutation happened either side of the root.
                    // We calculate the relative likelihood of each case and use this to weight waiting times etc.
                    RealNumType dist_to_root = seqP_region->plength_observation2root + blength;
                    RealNumType dist_to_observed = seqP_region->plength_observation2node;
                    updateCountsAndWaitingTimesAcrossRoot(pos, end_pos, stateA, stateB, dist_to_root, dist_to_observed, W, C);
                }              
            } else if(seqP_region->type <= TYPE_R && seqC_region->type == TYPE_O) {
                StateType stateA = seqP_region->type;
                if(seqP_region->type == TYPE_R) {
                    // stateA = tree->aln->ref_seq[static_cast<std::vector<cmaple::StateType>::size_type>(end_pos)];
                    stateA = seqC_region->prev_state;
                }

                // Calculate a weight vector giving the relative probabilities of observing
                // each state at the O node.
                std::vector<RealNumType> weight_vector(num_states_);
                RealNumType sum = 0.0;
                for(StateType stateB = 0; stateB < num_states_; stateB++) {
                    RealNumType likelihoodB = std::round(seqC_region->getLH(stateB) * 1000) / 1000.0;
                    if(stateB != stateA) {
                        RealNumType prob = likelihoodB * branch_length_to_observation * getMutationMatrixEntry(stateA, stateB, end_pos);
                        weight_vector[stateB] += prob;
                        sum += prob;
                    } else {
                        RealNumType prob = likelihoodB * (1 - branch_length_to_observation * getMutationMatrixEntry(stateB, stateB, end_pos));
                        weight_vector[stateB] += prob;
                        sum += prob; 
                    }
                }
                // Normalise weight vector 
                normalize_arr(weight_vector.data(), num_states_, sum);

                // Case 1: Last observation was this side of the root node
                if(seqP_region->plength_observation2root < 0) {
                    for(StateType stateB = 0; stateB < num_states_; stateB++) {
                        RealNumType prob = weight_vector[stateB];
                        if(stateB != stateA) {
                            C[end_pos * mat_size + stateB + row_index[stateA]] += prob;

                            W[end_pos * num_states_ + stateA] += prob * branch_length_to_observation/2;
                            W[end_pos * num_states_ + stateB] += prob * branch_length_to_observation/2;
                        } else {
                            W[end_pos * num_states_ + stateA] += prob * branch_length_to_observation;
                        }
                    }
                } else {
                    // Case 2: Last observation was the other side of the root.
                    RealNumType dist_to_root = seqP_region->plength_observation2root + blength;
                    RealNumType dist_to_observed = seqP_region->plength_observation2node;
                     for(StateType stateB = 0; stateB < num_states_; stateB++) {
                        RealNumType prob = weight_vector[stateB];
                        updateCountsAndWaitingTimesAcrossRoot(pos, end_pos, stateA, stateB, dist_to_root, dist_to_observed, W, C, prob);
                     }
                }
            } else if(seqP_region->type == TYPE_O && seqC_region->type <= TYPE_R) {
                StateType stateB = seqC_region->type;
                if(seqC_region->type == TYPE_R) {
                    // stateB = tree->aln->ref_seq[static_cast<std::vector<cmaple::StateType>::size_type>(end_pos)];
                    stateB = seqP_region->prev_state;
                }

                // Calculate a weight vector giving the relative probabilities of observing
                // each state at the O node.
                std::vector<RealNumType> weight_vector(num_states_);
                RealNumType sum = 0.0;
                for(StateType stateA = 0; stateA < num_states_; stateA++) {
                    RealNumType likelihoodA = std::round(seqP_region->getLH(stateA) * 1000)/1000.0;
                    if(stateA != stateB) {
                        RealNumType prob = likelihoodA * branch_length_to_observation * getMutationMatrixEntry(stateA, stateB, end_pos);
                        weight_vector[stateA] += prob;
                        sum += prob;
                    } else {
                        RealNumType prob = likelihoodA * (1 - branch_length_to_observation * getMutationMatrixEntry(stateA, stateA, end_pos));
                        weight_vector[stateA] += prob;
                        sum += prob;
                    }
                }
                // Normalise weight vector 
                normalize_arr(weight_vector.data(), num_states_, sum);

                // Case 1: Last observation was this side of the root node 
                if(seqP_region->plength_observation2root < 0) {
                    for(StateType stateA = 0; stateA < num_states_; stateA++) {
                        RealNumType prob = weight_vector[stateA];
                        if(stateB != stateA) {
                            C[end_pos * mat_size + stateB + row_index[stateA]] += prob;

                            W[end_pos * num_states_ + stateA] += prob * branch_length_to_observation/2;
                            W[end_pos * num_states_ + stateB] += prob * branch_length_to_observation/2;
                        } else {
                            W[end_pos * num_states_ + stateA] += prob * branch_length_to_observation;
                        }
                    }
                } else {
                    // Case 2: Last observation was the other side of the root.
                    RealNumType dist_to_root = seqP_region->plength_observation2root + blength;
                    RealNumType dist_to_observed = seqP_region->plength_observation2node;
                    for(StateType stateA = 0; stateA < num_states_; stateA++) {
                        RealNumType prob = weight_vector[stateA];
                        updateCountsAndWaitingTimesAcrossRoot(pos, end_pos, stateA, stateB, dist_to_root, dist_to_observed, W, C, prob);
                    }                    
                }
            } else if(seqP_region->type == TYPE_O && seqC_region->type == TYPE_O) {
                // Calculate a weight vector giving the relative probabilities of observing
                // each state at each of the O nodes.
                std::vector<RealNumType> weight_vector(num_states_ * num_states_);
                RealNumType sum = 0.0;
                for(StateType stateA = 0; stateA < num_states_; stateA++) {
                    RealNumType likelihoodA = std::round(seqP_region->getLH(stateA) * 1000)/1000.0;
                    for(StateType stateB = 0; stateB < num_states_; stateB++) {
                        RealNumType likelihoodB = std::round(seqC_region->getLH(stateB)*1000)/1000.0;
                        if(stateA != stateB) {
                            RealNumType prob = likelihoodA * likelihoodB * branch_length_to_observation * getMutationMatrixEntry(stateA, stateB, end_pos);
                            weight_vector[row_index[stateA] + stateB] += prob;
                            sum += prob;
                        } else {
                            RealNumType prob = likelihoodA * likelihoodB * (1 - branch_length_to_observation * getMutationMatrixEntry(stateA, stateA, end_pos));
                            weight_vector[row_index[stateA] + stateB] += prob;
                            sum += prob;
                        }
                    }
                }
                // Normalise weight vector 
                normalize_arr(weight_vector.data(), num_states_*num_states_, sum);

                // Case 1: Last observation was this side of the root node
                if(seqP_region->plength_observation2root < 0) {
                    for(StateType stateA = 0; stateA < num_states_; stateA++) {
                        for(StateType stateB = 0; stateB < num_states_; stateB++) {
                            RealNumType prob = weight_vector[row_index[stateA] + stateB];
                            if(stateB != stateA) {
                                C[end_pos * mat_size + stateB + row_index[stateA]] += prob;

                                W[end_pos * num_states_ + stateA] +=  prob * branch_length_to_observation/2;
                                W[end_pos * num_states_ + stateB] +=  prob * branch_length_to_observation/2;
                            } else {
                                W[end_pos * num_states_ + stateA] +=  prob * branch_length_to_observation;
                            }
                        }
                    }
                } else {
                     // Case 2: Last observation was the other side of the root.
                    RealNumType dist_to_root = seqP_region->plength_observation2root + blength;
                    RealNumType dist_to_observed = seqP_region->plength_observation2node;
                    for(StateType stateA = 0; stateA < num_states_; stateA++) {
                        for(StateType stateB = 0; stateB < num_states_; stateB++) {
                            RealNumType prob = weight_vector[row_index[stateA] + stateB];
                            updateCountsAndWaitingTimesAcrossRoot(pos, end_pos, stateA, stateB, dist_to_root, dist_to_observed, W, C, prob);
                        }
                    }                
                }
            } 
            pos = end_pos + 1;
        }
    }

    if(cmaple::verbose_mode > VB_MIN) 
    {
        const std::string prefix = tree->params->output_prefix.length() ? tree->params->output_prefix : tree->params->aln_path;
        //std::cout << "Writing rate matrices to file " << prefix << ".rateMatrices.txt" << std::endl;
        std::ofstream out_file(prefix + ".countMatrices.txt");
        for(int i = 0; i < genome_size; i++) {
            out_file << "Position: " << i << std::endl;
            out_file << "Count Matrix: " << std::endl;
            printCountsAndWaitingTimes(C + (i * mat_size), W + (i * num_states_), &out_file);
            out_file << std::endl;
        }
        out_file.close();
    }

    // Get genome-wide average mutation counts and waiting times
    RealNumType* global_counts = new RealNumType[mat_size];
    RealNumType* global_waiting_times = new RealNumType[num_states_];
    for(int i = 0; i < num_states_; i++) {
        global_waiting_times[i] = 0;
        for(int j = 0; j < num_states_; j++) {
            global_counts[row_index[i] + j] = 0;
        }
    }

    for(int i = 0; i < genome_size; i++) {
        for(int j = 0; j < num_states_; j++) {
            global_waiting_times[j] += W[i * num_states_ + j];
            for(int k = 0; k < num_states_; k++) {
                global_counts[row_index[j] + k] += C[i * mat_size + row_index[j] + k];
            }
        }        
    }

    for(int j = 0; j < num_states_; j++) {
        global_waiting_times[j] /= genome_size;
        for(int k = 0; k < num_states_; k++) {
            global_counts[row_index[j] + k] /= genome_size;
        }
    }

    // Add pseudocounts
    for(int i = 0; i < genome_size; i++) {
        for(int j = 0; j < num_states_; j++) {
            // Add pseudocounts to waiting_times
            W[i * num_states_ + j] += waiting_time_pseudocount;
            for(int k = 0; k < num_states_; k++) {
                // Add pseudocount of average rate across genome * waitingTime pseudocount for counts
                C[i * mat_size + row_index[j] + k] += global_counts[row_index[j] + k] * waiting_time_pseudocount / global_waiting_times[j];
            }
        }
    }

    /*if(cmaple::verbose_mode > VB_MIN) 
    {
        RealNumType* reference_freqs = new RealNumType[num_states_];
        for(int j = 0; j < num_states_; j++) {
            reference_freqs[j] = 0;
        }
        for(int i = 0; i < genome_size; i++) {
            reference_freqs[tree->aln->ref_seq[i]]++;
        }
        for(int j = 0; j < num_states_; j++) {
            reference_freqs[j] /= genome_size;
        }

        std::cout << "Genome-wide average waiting times:\t\t";
        for(int j = 0; j < num_states_; j++) {
            std::cout << global_waiting_times[j] << "\t" ;
        }
        std::cout << std::endl;
        std::cout << "Reference nucleotide frequencies:\t\t";
        for(int j = 0; j < num_states_; j++) {
            std::cout << reference_freqs[j] << "\t" ;
        }
        std::cout << std::endl;
        delete[] reference_freqs;
    }*/

    delete[] global_counts;
    delete[] global_waiting_times;

    RealNumType total_rate = 0;
    // Update mutation matrices with new rate estimation
    for(int i = 0; i < genome_size; i++) {
        RealNumType* Ci = C + (i * mat_size);
        RealNumType* Wi = W + (i * num_states_);
        StateType ref_state = tree->aln->ref_seq[static_cast<std::vector<cmaple::StateType>::size_type>(i)];

        for(int stateA = 0; stateA < num_states_; stateA++) {
            for(int stateB = 0; stateB < num_states_; stateB++) {
                if(stateA != stateB) { 
                    RealNumType new_rate = Ci[stateB + row_index[stateA]] / Wi[stateA];                
                    mutation_matrices[i * mat_size + (stateB + row_index[stateA])] = new_rate;

                    // Approximate total rate by considering rates from reference nucleotide
                    if(ref_state == stateA) {
                        total_rate += new_rate;
                    }
                }
            }
        }
    } 

    // Normalise entries of mutation matrices
    total_rate /= genome_size;
    //RealNumType average_rate = total_rate / genome_size;
    for(int i = 0; i < genome_size; i++) {
        for(int stateA = 0; stateA < num_states_; stateA++) {
            RealNumType row_sum = 0;
            for(int stateB = 0; stateB < num_states_; stateB++) {
                if(stateA != stateB) {
                    RealNumType val = mutation_matrices[i * mat_size + (stateB + row_index[stateA])];
                    //val /= average_rate;
                    val /= total_rate;
                    val = std::min(250.0, std::max(0.001, val));

                    mutation_matrices[i * mat_size + (stateB + row_index[stateA])] = val;
                    transposed_mutation_matrices[i * mat_size + (stateA + row_index[stateB])] = val;
                    freqi_freqj_Qijs[i * mat_size + (stateB + row_index[stateA])] = root_freqs[stateA] * inverse_root_freqs[stateB] * val;

                    row_sum += val;
                } 
            }
            mutation_matrices[i * mat_size + (stateA + row_index[stateA])] = -row_sum;
            transposed_mutation_matrices[i * mat_size + (stateA + row_index[stateA])] = -row_sum;
            diagonal_mutation_matrices[i * num_states_ + stateA] = -row_sum;
            freqi_freqj_Qijs[i * mat_size + (stateA + row_index[stateA])] = -row_sum;

            // pre-compute matrix to speedup
            const RealNumType* transposed_mut_mat_row = getTransposedMutationMatrixRow(stateA, i);
            RealNumType* freqj_transposedijs_row = freqj_transposedijs + (i * mat_size) + row_index[stateA];
            setVecByProduct<4>(freqj_transposedijs_row, root_freqs, transposed_mut_mat_row);
        
        }
    } 

    // Clean-up
    delete[] C;
    delete[] W;
}

void ModelDNARateVariation::updateCountsAndWaitingTimesAcrossRoot( PositionType start, PositionType end, 
                                                StateType parent_state, StateType child_state,
                                                RealNumType dist_to_root, RealNumType dist_to_observed,
                                                RealNumType* waiting_times, RealNumType* counts,
                                                RealNumType weight)
{
    if(parent_state != child_state) {
        for(int i = start; i <= end; i++) {
            RealNumType p_root_is_state_parent = root_freqs[parent_state] * getMutationMatrixEntry(parent_state, child_state, i) * dist_to_root;
            RealNumType p_root_is_state_child = root_freqs[child_state] * getMutationMatrixEntry(child_state, parent_state, i) * dist_to_observed;
            RealNumType relative_root_is_state_parent = p_root_is_state_parent / (p_root_is_state_parent + p_root_is_state_child);
            waiting_times[i * num_states_ +  parent_state] += weight * relative_root_is_state_parent * dist_to_root/2;
            waiting_times[i * num_states_ + child_state] += weight * relative_root_is_state_parent * dist_to_root/2;
            counts[i * mat_size + child_state + row_index[parent_state]] += weight * relative_root_is_state_parent;

            RealNumType relative_root_is_state_child = 1 - relative_root_is_state_parent;
            waiting_times[i * num_states_ + child_state] += weight * relative_root_is_state_child * dist_to_root;
        }
    } else {
        for(int i = start; i <= end; i++) {
            waiting_times[i * num_states_ + child_state] += weight * dist_to_root;
        }      
    }
}

void ModelDNARateVariation::setAllMatricesToDefault() {
    for(int i = 0; i < genome_size; i++) {
        for(int stateA = 0; stateA < num_states_; stateA++) {
            RealNumType row_sum = 0;
            for(int stateB = 0; stateB < num_states_; stateB++) {
                mutation_matrices[i * mat_size + (stateB + row_index[stateA])] = mutation_mat[stateB + row_index[stateA]];
                transposed_mutation_matrices[i * mat_size + (stateB + row_index[stateA])] = transposed_mut_mat[stateB + row_index[stateA]];
                freqi_freqj_Qijs[i * mat_size + (stateB + row_index[stateA])] = freqi_freqj_qij[stateB + row_index[stateA]];
            }

            diagonal_mutation_matrices[i * num_states_ + stateA] = diagonal_mut_mat[stateA];

            // pre-compute matrix to speedup
            const RealNumType* transposed_mut_mat_row = getTransposedMutationMatrixRow(stateA, i);
            RealNumType* freqj_transposedijs_row = freqj_transposedijs + (i * mat_size) + row_index[stateA];
            setVecByProduct<4>(freqj_transposedijs_row, root_freqs, transposed_mut_mat_row);
        
        }
    }  
}

void ModelDNARateVariation::setMatrixAtPosition(RealNumType* matrix, PositionType i) {
    for(int stateA = 0; stateA < num_states_; stateA++) {
        diagonal_mutation_matrices[i * num_states_ + stateA] = matrix[stateA + row_index[stateA]];
        for(int stateB = 0; stateB < num_states_; stateB++) {
            mutation_matrices[i * mat_size + (stateB + row_index[stateA])] = matrix[stateB + row_index[stateA]];
            transposed_mutation_matrices[i * mat_size + (stateB + row_index[stateA])] = matrix[stateA + row_index[stateB]];
        }
    }
}

void ModelDNARateVariation::readRatesFile() {
    std::ifstream infile(rates_filename);
    std::string line;
    if (infile.is_open()) {
        PositionType genome_position = 0;
        while (std::getline(infile, line)) {
            std::stringstream ss(line);
            std::string field;
            std::vector<std::string> fields;
            while (std::getline(ss, field, ' ')) {
                // Add the word to the vector
                fields.push_back(field);
            }

            if(fields.size() != 12) {
                std::cerr << "Error: Rate matrix file not in expected format." << std::endl;
                std::cerr << "Expected exactly 12 entries." << std::endl;
                continue;
            }
            RealNumType* rate_matrix = mutation_matrices + (genome_position * mat_size);
            
            // A row
            rate_matrix[1] = std::stof(fields[0]);
            rate_matrix[2] = std::stof(fields[1]);
            rate_matrix[3] = std::stof(fields[2]);

            // C row
            rate_matrix[4] = std::stof(fields[3]);
            rate_matrix[6] = std::stof(fields[4]);
            rate_matrix[7] = std::stof(fields[5]);

            // G row
            rate_matrix[8] = std::stof(fields[6]);
            rate_matrix[9] = std::stof(fields[7]);
            rate_matrix[11] = std::stof(fields[8]);

             // T row
            rate_matrix[12] = std::stof(fields[9]);
            rate_matrix[13] = std::stof(fields[10]);
            rate_matrix[14] = std::stof(fields[11]);

            genome_position++;
        }
        infile.close();

        if(genome_position != genome_size) {
            std::cerr << "Error: Rate matrix file contained unexpected number of lines." << std::endl;
            return;
        }

        for(int i = 0; i < genome_size; i++) {
            for(int stateA = 0; stateA < num_states_; stateA++) {
                RealNumType row_sum = 0;
                for(int stateB = 0; stateB < num_states_; stateB++) {
                    if(stateA != stateB) {
                        RealNumType val = mutation_matrices[i * mat_size + (stateB + row_index[stateA])];
                        mutation_matrices[i * mat_size + (stateB + row_index[stateA])] = val;
                        transposed_mutation_matrices[i * mat_size + (stateA + row_index[stateB])] = val;
                        freqi_freqj_Qijs[i * mat_size + (stateB + row_index[stateA])] = root_freqs[stateA] * inverse_root_freqs[stateB] * val;

                        row_sum += val;
                    } 
                }
                mutation_matrices[i * mat_size + (stateA + row_index[stateA])] = -row_sum;
                transposed_mutation_matrices[i * mat_size + (stateA + row_index[stateA])] = -row_sum;
                diagonal_mutation_matrices[i * num_states_ + stateA] = -row_sum;
                freqi_freqj_Qijs[i * mat_size + (stateA + row_index[stateA])] = -row_sum;

                // pre-compute matrix to speedup
                const RealNumType* transposed_mut_mat_row = getTransposedMutationMatrixRow(stateA, i);
                RealNumType* freqj_transposedijs_row = freqj_transposedijs + (i * mat_size) + row_index[stateA];
                setVecByProduct<4>(freqj_transposedijs_row, root_freqs, transposed_mut_mat_row);      
            }
        } 
    }
    else {
        std::cerr << "Unable to open rate matrix file " << rates_filename << std::endl;
    }
}
