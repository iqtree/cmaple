#include "model_dna_rate_variation.h"
#include "../tree/tree.h"
#include "../tree/phylonode.h"

using namespace cmaple;

ModelDNARateVariation::ModelDNARateVariation(const cmaple::ModelBase::SubModel sub_model, PositionType _genomeSize, bool _useSiteRates)
    : ModelDNA(sub_model) {
    
    genomeSize = _genomeSize;
    useSiteRates = _useSiteRates;
    matSize = row_index[num_states_];

    mutationMatrices = new RealNumType[matSize * genomeSize]();
    transposedMutationMatrices = new RealNumType[matSize * genomeSize]();
    diagonalMutationMatrices = new RealNumType[num_states_ * genomeSize]();
    freqiFreqjQijs = new RealNumType[matSize * genomeSize]();
    freqjTransposedijs = new RealNumType[matSize * genomeSize]();

    if(useSiteRates) {
        rates = new cmaple::RealNumType[genomeSize]();
    }
}

ModelDNARateVariation::~ModelDNARateVariation() { 
    delete[] mutationMatrices;
    delete[] transposedMutationMatrices;
    delete[] diagonalMutationMatrices;
    delete[] freqiFreqjQijs;
    delete[] freqjTransposedijs;
    if(useSiteRates) {
        delete[] rates;
    }
}

void ModelDNARateVariation::printMatrix(const RealNumType* matrix, std::ostream* outStream) {
    for(int j = 0; j < num_states_; j++) {
        std::string line = "|";
        for(int k = 0; k < num_states_; k++) {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(5) << matrix[row_index[j] + k];
            line += "\t" + oss.str();
            //line += "\t" + convertDoubleToString(matrix[row_index[j] + k], 5);
        }
        line += "\t|";
        (*outStream) << line << std::endl;
    }
}

void ModelDNARateVariation::printCountsAndWaitingTimes(const RealNumType* counts, const RealNumType* waitingTImes, std::ostream* outStream) {
    for(int j = 0; j < num_states_; j++) {
        std::string line = "|";
        for(int k = 0; k < num_states_; k++) {
            line += "\t" + convertDoubleToString(counts[row_index[j] + k]);
        }
        line += "\t|\t" + convertDoubleToString(waitingTImes[j]) + "\t|";
        (*outStream) << line << std::endl;
    }
}

bool ModelDNARateVariation::updateMutationMatEmpirical() {
    bool val = ModelDNA::updateMutationMatEmpirical();
    for(int i = 0; i < genomeSize; i++) {
        for(int j = 0; j < matSize; j++) {
            mutationMatrices[i * matSize + j] = mutation_mat[j];
            transposedMutationMatrices[i * matSize + j] = transposed_mut_mat[j];
            freqiFreqjQijs[i * matSize + j] = freqi_freqj_qij[j];
            freqjTransposedijs[i * matSize + j] = freq_j_transposed_ij[j];
        }
        //std::cout << std::endl;
        for(int j = 0; j < num_states_; j++) {
            diagonalMutationMatrices[i * num_states_ + j] = diagonal_mut_mat[j];
        }
    }
    return val;
}

void ModelDNARateVariation::estimateRates(cmaple::Tree* tree) {
    if(useSiteRates) {
        estimateRatePerSite(tree);

        //todo: write rates to file

    } else {
        RealNumType oldLK = -std::numeric_limits<double>::infinity();
        RealNumType newLK = tree->computeLh();
        int numSteps = 0;
        while(newLK - oldLK > 1 && numSteps < 20) {
            estimateRatesPerSitePerEntry(tree);
            oldLK = newLK;
            newLK = tree->computeLh();
        }  
    }

    // Write out rate matrices to file
    if(cmaple::verbose_mode > VB_MIN) 
    {
        const std::string prefix = tree->params->output_prefix.length() ? 
            tree->params->output_prefix : tree->params->aln_path;
        //std::cout << "Writing rate matrices to file " << prefix << ".rateMatrices.txt" << std::endl;
        std::ofstream outFile(prefix + ".rateMatrices.txt");
        outFile << "Rate matrix for all sites: " << std::endl;
        printMatrix(getOriginalRateMatrix(), &outFile);
        for(int i = 0; i < genomeSize; i++) {
            outFile << "Position: " << i << std::endl;
            if(useSiteRates) {
                outFile << "Rate: " << rates[i] << std::endl;
            }
            outFile << "Rate Matrix: " << std::endl;
            printMatrix(getMutationMatrix(i), &outFile);
            outFile << std::endl;
        }
        outFile.close();
    } 
}

void ModelDNARateVariation::estimateRatePerSite(cmaple::Tree* tree){
    std::cout << "Estimating mutation rate per site..." << std::endl;
    RealNumType** totals = new RealNumType*[num_states_];
    for(int j = 0; j < num_states_; j++) {
        totals[j] = new RealNumType[genomeSize];
    }
    RealNumType* numSubstitutions = new RealNumType[genomeSize];
    for(int i = 0; i < genomeSize; i++) {
        for(int j = 0; j < num_states_; j++) {
            totals[j][i] = 0;
        }
        numSubstitutions[i] = 0;
    }

    std::stack<Index> nodeStack;
    const PhyloNode& root = tree->nodes[tree->root_vector_index];
    nodeStack.push(root.getNeighborIndex(RIGHT));
    nodeStack.push(root.getNeighborIndex(LEFT));
    while(!nodeStack.empty()) {

        Index index = nodeStack.top();
        nodeStack.pop();
        PhyloNode& node = tree->nodes[index.getVectorIndex()];
        RealNumType blength = node.getUpperLength();
        //std::cout << "blength: " << blength  << std::endl;

        if (node.isInternal()) {
            nodeStack.push(node.getNeighborIndex(RIGHT));
            nodeStack.push(node.getNeighborIndex(LEFT));
        }

        if(blength <= 0.) {
            continue;
        }

        Index parent_index = node.getNeighborIndex(TOP);
        PhyloNode& parent_node = tree->nodes[parent_index.getVectorIndex()];
        const std::unique_ptr<SeqRegions>& parentRegions = parent_node.getPartialLh(parent_index.getMiniIndex());
        const std::unique_ptr<SeqRegions>& childRegions = node.getPartialLh(TOP);

        PositionType pos = 0;
        const SeqRegions& seq1_regions = *parentRegions;
        const SeqRegions& seq2_regions = *childRegions;
        size_t iseq1 = 0;
        size_t iseq2 = 0;

        while(pos < genomeSize) {
            PositionType end_pos;
            SeqRegions::getNextSharedSegment(pos, seq1_regions, seq2_regions, iseq1, iseq2, end_pos);
            const auto* seq1_region = &seq1_regions[iseq1];
            const auto* seq2_region = &seq2_regions[iseq2];

            if(seq1_region->type == TYPE_R && seq2_region->type == TYPE_R) {
                // both states are type REF
                for(int i = pos; i <= end_pos; i++) {
                    int state = tree->aln->ref_seq[static_cast<std::vector<cmaple::StateType>::size_type>(i)];
                    totals[state][i] += blength;
                }
            }  else if(seq1_region->type == seq2_region->type && seq1_region->type < TYPE_R) {
                // both states are equal but not of type REF
                 for(int i = pos; i <= end_pos; i++) {
                    totals[seq1_region->type][i] += blength;
                }               
            } else if(seq1_region->type <= TYPE_R && seq2_region->type <= TYPE_R) {
                // both states are not equal
                  for(int i = pos; i <= end_pos; i++) {
                    numSubstitutions[i] += 1;
                }                 
            }
            pos = end_pos + 1;
        }
    }

    RealNumType rateCount = 0;
    for(int i = 0; i < genomeSize; i++) {
        if(numSubstitutions[i] == 0) {
            rates[i] = 0.001;
        } else {
            RealNumType expectedRateNoSubstitution = 0;
            for(int j = 0; j < num_states_; j++) {
                RealNumType summand = totals[j][i] * abs(diagonal_mut_mat[j]);
                expectedRateNoSubstitution += summand;
            }
            if(expectedRateNoSubstitution <= 0.01) {
                rates[i] = 1.;
            } else {
                rates[i] = numSubstitutions[i] / expectedRateNoSubstitution;
                rates[i] = MIN(100.0, MAX(0.0001, rates[i]));  
            }
        }
        rateCount += rates[i];
    }

    RealNumType averageRate = rateCount / genomeSize;
    for(int i = 0; i < genomeSize; i++) {
        rates[i] /= averageRate;
        //std::cout << rates[i] << " ";
        for(int stateA = 0; stateA < num_states_; stateA++) {
            RealNumType rowSum = 0;
            for(int stateB = 0; stateB < num_states_; stateB++) {
                if(stateA != stateB) {
                    RealNumType val = mutationMatrices[i * matSize + (stateB + row_index[stateA])] * rates[i];
                    mutationMatrices[i * matSize + (stateB + row_index[stateA])] = val;
                    transposedMutationMatrices[i * matSize + (stateA + row_index[stateB])] = val;
                    freqiFreqjQijs[i * matSize + (stateB + row_index[stateA])] = root_freqs[stateA] * inverse_root_freqs[stateB] * val;
                    rowSum += val;
                }
            }
            mutationMatrices[i * matSize + (stateA + row_index[stateA])] = -rowSum;
            transposedMutationMatrices[i * matSize + (stateA + row_index[stateA])] = -rowSum;
            diagonalMutationMatrices[i * num_states_ + stateA] = -rowSum;
            freqiFreqjQijs[i * matSize + (stateA + row_index[stateA])] = -rowSum;

            // pre-compute matrix to speedup
            const RealNumType* transposed_mut_mat_row = getTransposedMutationMatrixRow(stateA, i);
            RealNumType* freqjTransposedijsRow = freqjTransposedijs + (i * matSize) + row_index[stateA];
            setVecByProduct<4>(freqjTransposedijsRow, root_freqs, transposed_mut_mat_row);
        
        }
    }

    for(int j = 0; j < num_states_; j++) {
        delete[] totals[j];
    }
    delete[] totals;
    delete[] numSubstitutions;
}

void ModelDNARateVariation::estimateRatesPerSitePerEntry(cmaple::Tree* tree) {

    // Possibly better to keep this memory allocated and reuse
    // since this will be called several times?
    RealNumType** C = new RealNumType*[genomeSize];
    RealNumType** W = new RealNumType*[genomeSize];
    for(int i = 0; i < genomeSize; i++) {
        C[i] = new RealNumType[matSize];
        W[i] = new RealNumType[num_states_];
        for(int j = 0; j < num_states_; j++) {
            W[i][j] = 0.0001;
            for(int k = 0; k < num_states_; k++) {
                C[i][row_index[j] + k] = 0;
            }
        }
    }

    std::stack<Index> nodeStack;
    const PhyloNode& root = tree->nodes[tree->root_vector_index];
    nodeStack.push(root.getNeighborIndex(RIGHT));
    nodeStack.push(root.getNeighborIndex(LEFT));
    while(!nodeStack.empty()) {

        Index index = nodeStack.top();
        nodeStack.pop();
        PhyloNode& node = tree->nodes[index.getVectorIndex()];
        RealNumType blength = node.getUpperLength();

        if (node.isInternal()) {
            nodeStack.push(node.getNeighborIndex(RIGHT));
            nodeStack.push(node.getNeighborIndex(LEFT));
        }

        if(blength <= 0.) {
            continue;
        }

        Index parent_index = node.getNeighborIndex(TOP);
        PhyloNode& parent_node = tree->nodes[parent_index.getVectorIndex()];
        const std::unique_ptr<SeqRegions>& parentRegions = parent_node.getPartialLh(parent_index.getMiniIndex());
        const std::unique_ptr<SeqRegions>& childRegions = node.getPartialLh(TOP);

        PositionType pos = 0;
        const SeqRegions& seqP_regions = *parentRegions;
        const SeqRegions& seqC_regions = *childRegions;
        size_t iseq1 = 0;
        size_t iseq2 = 0;

        while(pos < genomeSize) {
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
            RealNumType branchLengthToObservation = blength;
            if(seqP_region->plength_observation2node > 0 && seqP_region->plength_observation2root <= 0) {
                branchLengthToObservation = blength + seqP_region->plength_observation2node;
            }
            else if(seqP_region->plength_observation2root > 0) {
                branchLengthToObservation = blength + seqP_region->plength_observation2root;
            }

            if(seqP_region->type == TYPE_R && seqC_region->type == TYPE_R) {
                // both states are type REF
                for(int i = pos; i <= end_pos; i++) {
                    StateType state = tree->aln->ref_seq[static_cast<std::vector<cmaple::StateType>::size_type>(i)];
                    W[i][state] += branchLengthToObservation;
                }
            }  else if(seqP_region->type == seqC_region->type && seqP_region->type < TYPE_R) {
                // both states are equal but not of type REF
                 for(int i = pos; i <= end_pos; i++) {
                    W[i][seqP_region->type] += branchLengthToObservation;
                }               
            } else if(seqP_region->type <= TYPE_R && seqC_region->type <= TYPE_R) {
                //states are not equal
                StateType stateA = seqP_region->type;
                StateType stateB = seqC_region->type;
                if(seqP_region->type == TYPE_R) {
                    stateA = tree->aln->ref_seq[static_cast<std::vector<cmaple::StateType>::size_type>(end_pos)];
                }
                if(seqC_region->type == TYPE_R) {
                    stateB = tree->aln->ref_seq[static_cast<std::vector<cmaple::StateType>::size_type>(end_pos)];
                }
                // Case 1: Last observation was this side of the root node
                if(seqP_region->plength_observation2root <= 0) {
                    for(int i = pos; i <= end_pos; i++) {
                        W[i][stateA] += branchLengthToObservation/2;
                        W[i][stateB] += branchLengthToObservation/2;
                        C[i][stateB + row_index[stateA]] += 1;
                    }
                } else {
                    // Case 2: Last observation was the other side of the root.
                    // In this case there are two further cases - the mutation happened either side of the root.
                    // We calculate the relative likelihood of each case and use this to weight waiting times etc.
                    RealNumType distToRoot = seqP_region->plength_observation2root + blength;
                    RealNumType distToObserved = seqP_region->plength_observation2node;
                    updateCountsAndWaitingTimesAcrossRoot(pos, end_pos, stateA, stateB, distToRoot, distToObserved, W, C);
                }              
            } else if(seqP_region->type <= TYPE_R && seqC_region->type == TYPE_O) {
                StateType stateA = seqP_region->type;
                if(seqP_region->type == TYPE_R) {
                    stateA = tree->aln->ref_seq[static_cast<std::vector<cmaple::StateType>::size_type>(end_pos)];
                }
                // Case 1: Last observation was this side of the root node
                if(seqP_region->plength_observation2root <= 0) {
                    for(StateType stateB = 0; stateB < num_states_; stateB++) {
                        RealNumType prob = seqC_region->getLH(stateB);
                        if(stateB != stateA) {
                            C[end_pos][stateB + row_index[stateA]] += prob;
                            W[end_pos][stateA] += prob * branchLengthToObservation/2;
                            W[end_pos][stateB] += prob * branchLengthToObservation/2;
                        } else {
                            W[end_pos][stateA] += prob * branchLengthToObservation;
                        }
                    }
                } else {
                    // Case 2: Last observation was the other side of the root.
                    RealNumType distToRoot = seqP_region->plength_observation2root + blength;
                    RealNumType distToObserved = seqP_region->plength_observation2node;
                     for(StateType stateB = 0; stateB < num_states_; stateB++) {
                        RealNumType prob = seqC_region->getLH(stateB);
                        updateCountsAndWaitingTimesAcrossRoot(pos, end_pos, stateA, stateB, distToRoot, distToObserved, W, C, prob);
                     }
                }
            } else if(seqP_region->type == TYPE_O && seqC_region->type <= TYPE_R) {
                StateType stateB = seqC_region->type;
                if(seqC_region->type == TYPE_R) {
                    stateB = tree->aln->ref_seq[static_cast<std::vector<cmaple::StateType>::size_type>(end_pos)];
                } 
                // Case 1: Last observation was this side of the root node
                if(seqP_region->plength_observation2root <= 0) {
                    for(StateType stateA = 0; stateA < num_states_; stateA++) {
                        RealNumType prob = seqP_region->getLH(stateA);
                        if(stateB != stateA) {
                            C[end_pos][stateB + row_index[stateA]] += prob;
                            W[end_pos][stateA] += prob * branchLengthToObservation/2;
                            W[end_pos][stateB] += prob * branchLengthToObservation/2;
                        } else {
                            W[end_pos][stateA] += prob * branchLengthToObservation;
                        }
                    }
                } else {
                    // Case 2: Last observation was the other side of the root.
                    RealNumType distToRoot = seqP_region->plength_observation2root + blength;
                    RealNumType distToObserved = seqP_region->plength_observation2node;
                    for(StateType stateA = 0; stateA < num_states_; stateA++) {
                        RealNumType prob = seqP_region->getLH(stateA);
                        updateCountsAndWaitingTimesAcrossRoot(pos, end_pos, stateA, stateB, distToRoot, distToObserved, W, C, prob);
                    }                    
                }
            } else if(seqP_region->type == TYPE_O && seqC_region->type == TYPE_O) {
                // Case 1: Last observation was this side of the root node
                if(seqP_region->plength_observation2root <= 0) {
                    for(StateType stateA = 0; stateA < num_states_; stateA++) {
                        RealNumType probA = seqP_region->getLH(stateA);
                        for(StateType stateB = 0; stateB < num_states_; stateB++) {
                            RealNumType probB = seqC_region->getLH(stateB);
                            if(stateB != stateA) {
                                C[end_pos][stateB + row_index[stateA]] += probA * probB;
                                W[end_pos][stateA] +=  probA * probB * branchLengthToObservation/2;
                                W[end_pos][stateB] +=  probA * probB * branchLengthToObservation/2;
                            } else {
                                W[end_pos][stateA] +=  probA * probB * branchLengthToObservation;
                            }
                        }
                    }
                } else {
                     // Case 2: Last observation was the other side of the root.
                    RealNumType distToRoot = seqP_region->plength_observation2root + blength;
                    RealNumType distToObserved = seqP_region->plength_observation2node;
                    for(StateType stateA = 0; stateA < num_states_; stateA++) {
                        RealNumType probA = seqP_region->getLH(stateA);
                        for(StateType stateB = 0; stateB < num_states_; stateB++) {
                            RealNumType probB = seqC_region->getLH(stateB);
                            updateCountsAndWaitingTimesAcrossRoot(pos, end_pos, stateA, stateB, distToRoot, distToObserved, W, C, probA * probB);
                        }
                    }                
                }
            } 
            pos = end_pos + 1;
        }
    }

    RealNumType totalRate = 0;

    // Update mutation matrices with new rate estimation
    for(int i = 0; i < genomeSize; i++) {
        RealNumType* Ci = C[i];
        RealNumType* Wi = W[i];
        for(int stateA = 0; stateA < num_states_; stateA++) {
            for(int stateB = 0; stateB < num_states_; stateB++) {
                if(stateA != stateB) { 
                    RealNumType newRate;     
                    if(Ci[stateB + row_index[stateA]] == 0) {
                        newRate = 0.001;
                    }  else if (Wi[stateA] <= 0.01) {
                        newRate = 1.0;
                    } else {
                        newRate = Ci[stateB + row_index[stateA]] / Wi[stateA];
                        newRate = MIN(100.0, MAX(0.0001, newRate));                 
                    }
                    mutationMatrices[i * matSize + (stateB + row_index[stateA])] = newRate;
                    totalRate += newRate;
                }
            }
        }
    }

    if(cmaple::verbose_mode > VB_MIN) 
    {
        const std::string prefix = tree->params->output_prefix.length() ? tree->params->output_prefix : tree->params->aln_path;
        //std::cout << "Writing rate matrices to file " << prefix << ".rateMatrices.txt" << std::endl;
        std::ofstream outFile(prefix + ".countMatrices.txt");
        for(int i = 0; i < genomeSize; i++) {
            outFile << "Position: " << i << std::endl;
            outFile << "Count Matrix: " << std::endl;
            printCountsAndWaitingTimes(C[i], W[i], &outFile);
            outFile << std::endl;
        }
        outFile.close();
    }  

    // Normalise entries of mutation matrices
    RealNumType averageRate = totalRate / genomeSize;
    for(int i = 0; i < genomeSize; i++) {
        for(int stateA = 0; stateA < num_states_; stateA++) {
            RealNumType rowSum = 0;
            for(int stateB = 0; stateB < num_states_; stateB++) {
                if(stateA != stateB) {
                    RealNumType val = mutationMatrices[i * matSize + (stateB + row_index[stateA])];
                    val /= averageRate;

                    mutationMatrices[i * matSize + (stateB + row_index[stateA])] = val;
                    transposedMutationMatrices[i * matSize + (stateA + row_index[stateB])] = val;
                    freqiFreqjQijs[i * matSize + (stateB + row_index[stateA])] = root_freqs[stateA] * inverse_root_freqs[stateB] * val;

                    rowSum += val;
                } 
            }
            mutationMatrices[i * matSize + (stateA + row_index[stateA])] = -rowSum;
            transposedMutationMatrices[i * matSize + (stateA + row_index[stateA])] = -rowSum;
            diagonalMutationMatrices[i * num_states_ + stateA] = -rowSum;
            freqiFreqjQijs[i * matSize + (stateA + row_index[stateA])] = -rowSum;

            // pre-compute matrix to speedup
            const RealNumType* transposed_mut_mat_row = getTransposedMutationMatrixRow(stateA, i);
            RealNumType* freqjTransposedijsRow = freqjTransposedijs + (i * matSize) + row_index[stateA];
            setVecByProduct<4>(freqjTransposedijsRow, root_freqs, transposed_mut_mat_row);
        
        }
    } 

    // Clean-up
    for(int i = 0; i < genomeSize; i++) {
        delete[] C[i];
        delete[] W[i];
    }
    delete[] C;
    delete[] W;

}

void ModelDNARateVariation::updateCountsAndWaitingTimesAcrossRoot( PositionType start, PositionType end, 
                                                StateType parentState, StateType childState,
                                                RealNumType distToRoot, RealNumType distToObserved,
                                                RealNumType** waitingTimes, RealNumType** counts,
                                                RealNumType weight)
{
    if(parentState != childState) {
        for(int i = start; i <= end; i++) {
            RealNumType pRootIsStateParent = root_freqs[parentState] * getMutationMatrixEntry(parentState, childState, i) * distToRoot;
            RealNumType pRootIsStateChild = root_freqs[childState] * getMutationMatrixEntry(childState, parentState, i) * distToObserved;
            RealNumType relativeRootIsStateParent = pRootIsStateParent / (pRootIsStateParent + pRootIsStateChild);
            waitingTimes[i][parentState] += weight * relativeRootIsStateParent * distToRoot/2;
            waitingTimes[i][childState] += weight * relativeRootIsStateParent * distToRoot/2;
            counts[i][childState + row_index[parentState]] += weight * relativeRootIsStateParent;

            RealNumType relativeRootIsStateChild = 1 - relativeRootIsStateParent;
            waitingTimes[i][childState] += weight * relativeRootIsStateChild * distToRoot;
        }
    } else {
        for(int i = start; i <= end; i++) {
            waitingTimes[i][childState] += weight * distToRoot;
        }      
    }
}