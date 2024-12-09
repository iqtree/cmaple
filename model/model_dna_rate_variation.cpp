#include "model_dna_rate_variation.h"
#include "../tree/tree.h"
#include "../tree/phylonode.h"

using namespace cmaple;

ModelDNARateVariation::ModelDNARateVariation(const cmaple::ModelBase::SubModel sub_model, PositionType _genomeSize)
    : ModelDNA(sub_model) {
    
    genomeSize = _genomeSize;
    matSize = row_index[num_states_];

    rates = new cmaple::RealNumType[genomeSize]();
    mutationMatrices = new RealNumType[matSize * genomeSize]();
    transposedMutationMatrices = new RealNumType[matSize * genomeSize]();
    diagonalMutationMatrices = new RealNumType[num_states_ * genomeSize]();
}

ModelDNARateVariation::~ModelDNARateVariation() { 
    delete[] rates; 
    delete[] mutationMatrices;
    delete[] transposedMutationMatrices;
    delete[] diagonalMutationMatrices;
}

void ModelDNARateVariation::printMatrix(const RealNumType* matrix, std::ostream* outStream) {
    for(int j = 0; j < num_states_; j++) {
        std::string line = "|";
        for(int k = 0; k < num_states_; k++) {
            line += "\t" + convertDoubleToString(matrix[row_index[j] + k]);
        }
        line += "\t|";
        (*outStream) << line << std::endl;
    }
}

bool ModelDNARateVariation::updateMutationMatEmpirical() {
    bool val = ModelDNA::updateMutationMatEmpirical();
    for(int i = 0; i < genomeSize; i++) {
        for(int j = 0; j < matSize; j++) {
            mutationMatrices[i * matSize + j] = mutation_mat[j];
            transposedMutationMatrices[i * matSize + j] = transposed_mut_mat[j];
        }
        //std::cout << std::endl;
        for(int j = 0; j < num_states_; j++) {
            diagonalMutationMatrices[i * num_states_ + j] = diagonal_mut_mat[j];
        }
    }
    return val;
}

void ModelDNARateVariation::estimateRatesPerSite(cmaple::Tree* tree){
    std::cout << "Estimating per site rates" << std::endl;
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

        Index parent_index = node.getNeighborIndex(TOP);
        PhyloNode& parent_node = tree->nodes[parent_index.getVectorIndex()];
        const std::unique_ptr<SeqRegions>& parentRegions = parent_node.getPartialLh(parent_index.getMiniIndex());
        const std::unique_ptr<SeqRegions>& childRegions = node.getPartialLh(TOP);

        if (node.isInternal()) {
            nodeStack.push(node.getNeighborIndex(RIGHT));
            nodeStack.push(node.getNeighborIndex(LEFT));
        }

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
            //std::cout << "blength: " << blength  << " ";

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
                RealNumType summand = totals[j][i] * diagonal_mut_mat[j];
                expectedRateNoSubstitution += summand;
            }
            if(expectedRateNoSubstitution == 0) {
                rates[i] = 1.;
            } else {
                rates[i] = numSubstitutions[i] / expectedRateNoSubstitution;
            }
        }
        rateCount += rates[i];
    }

    RealNumType averageRate = rateCount / genomeSize;
    for(int i = 0; i < genomeSize; i++) {
        rates[i] /= averageRate;
        //std::cout << rates[i] << " ";
        for(int j = 0; j < matSize; j++) {
            mutationMatrices[i * matSize + j] *= rates[i];
            transposedMutationMatrices[i * matSize + j] *= rates[i];
            //std::cout << mutationMatrices[i * matSize + j] << " ";
        }
        //std::cout << std::endl;
        for(int j = 0; j < num_states_; j++) {
            diagonalMutationMatrices[i * num_states_ + j] *= rates[i];
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
            W[i][j] = 0.001;
            for(int k = 0; k < num_states_; k++) {
                C[i][row_index[j] + k] = 0.001;
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

        if(blength == 0.) {
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

    RealNumType* totalRates = new RealNumType[matSize];
    for(int j  = 0; j < matSize; j++) {
        totalRates[j] = 0;
    }

    // Update mutation matrices with new rate estimation
    // TODO: Update other matrices.
    for(int i = 0; i < genomeSize; i++) {
        for(int stateA = 0; stateA < num_states_; stateA++) {
            for(int stateB = 0; stateB < num_states_; stateB++) {
                if(stateA != stateB) {             
                    RealNumType newRate = C[i][stateB + row_index[stateA]] / W[i][stateA];          
                    mutationMatrices[i * matSize + (stateB + row_index[stateA])] = newRate;
                    totalRates[stateB + row_index[stateA]] += newRate;
                }
            }
        }
    }

    /*
    std::cout << "Counts matrix: " << std::endl;
    this->printMatrix(C[1000], &std::cout);
    std::cout << "Waiting times: " << std::endl;
    for(int i = 0; i < num_states_; i ++) {
        std::cout << convertDoubleToString(W[1000][i]) << "\t";
    }
    std::cout << std::endl << std::endl;
    */

    // Normalise entries of mutation matrices
   for(int j  = 0; j < matSize; j++) {
        totalRates[j] /= genomeSize;
    }
    for(int i = 0; i < genomeSize; i++) {
        for(int stateA = 0; stateA < num_states_; stateA++) {
            RealNumType rowSum = 0;
            for(int stateB = 0; stateB < num_states_; stateB++) {
                if(stateA != stateB) {
                    RealNumType val = mutationMatrices[i * matSize + (stateB + row_index[stateA])];
                    val /= totalRates[stateB + row_index[stateA]];
                    val = MIN(100.0, MAX(0.001, val));  
                    mutationMatrices[i * matSize + (stateB + row_index[stateA])] = val;
                    transposedMutationMatrices[i * matSize + (stateA + row_index[stateB])] = val;
                    rowSum += val;
                }
            }
            mutationMatrices[i * matSize + (stateA + row_index[stateA])] = -rowSum;
            transposedMutationMatrices[i * matSize + (stateA + row_index[stateA])] = -rowSum;
            diagonalMutationMatrices[i * num_states_ + stateA] = -rowSum;
        }
    } 

    // Clean-up
    delete[] totalRates;
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