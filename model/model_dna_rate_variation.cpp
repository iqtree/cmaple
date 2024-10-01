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

void ModelDNARateVariation::printMatrix(RealNumType* matrix) {
    for(int j = 0; j < num_states_; j++) {
        std::string line = "|";
        for(int k = 0; k < num_states_; k++) {
            line += "\t" + std::to_string(matrix[row_index[j] + k]);
        }
        line += "\t|";
        std::cout << line << std::endl;
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

void ModelDNARateVariation::estimateRates(cmaple::Tree* tree){
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

        const Index parent_index = node.getNeighborIndex(TOP);
        std::unique_ptr<SeqRegions>& parentRegions = tree->getPartialLhAtNode(parent_index); 
        std::unique_ptr<SeqRegions>& childRegions = node.getPartialLh(TOP);

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

void ModelDNARateVariation::estimateRatesWithEM(cmaple::Tree* tree) {

    // Possibly better to keep this memory allocated and reuse
    // since this will be called several times.
    
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

        const Index parent_index = node.getNeighborIndex(TOP);
        std::unique_ptr<SeqRegions>& parentRegions = tree->getPartialLhAtNode(parent_index); 
        std::unique_ptr<SeqRegions>& childRegions = node.getPartialLh(TOP);

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

            if(seq1_region->type == TYPE_R && seq2_region->type == TYPE_R) {
                // both states are type REF
                for(int i = pos; i <= end_pos; i++) {
                    auto state = tree->aln->ref_seq[static_cast<std::vector<cmaple::StateType>::size_type>(end_pos)];
                    W[i][state] += blength;
                }
            }  else if(seq1_region->type == seq2_region->type && seq1_region->type < TYPE_R) {
                // both states are equal but not of type REF
                 for(int i = pos; i <= end_pos; i++) {
                    W[i][seq1_region->type] += blength;
                }               
            } else if(seq1_region->type <= TYPE_R && seq2_region->type <= TYPE_R) {
                // both states are not equal
                auto stateA = seq1_region->type;
                auto stateB = seq2_region->type;
                  for(int i = pos; i <= end_pos; i++) {
                    W[i][stateA] += blength/2;
                    W[i][stateB] += blength/2;
                    C[i][stateB + row_index[stateA]] += 1;
                }                 
            }
            pos = end_pos + 1;
        }
    }

    RealNumType* totalRates = new RealNumType[matSize];
    for(int j  = 0; j < matSize; j++) {
        totalRates[j] = 0;
    }

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

    //Normalise:
    for(int j  = 0; j < matSize; j++) {
        totalRates[j] /= genomeSize;
    }
    for(int i = 0; i < genomeSize; i++) {
        for(int stateA = 0; stateA < num_states_; stateA++) {
            RealNumType rowSum = 0;
            for(int stateB = 0; stateB < num_states_; stateB++) {
                if(stateA != stateB) {
                    mutationMatrices[i * matSize + (stateB + row_index[stateA])] /= totalRates[stateB + row_index[stateA]];
                    mutationMatrices[i * matSize + (stateB + row_index[stateA])] = MIN(100.0, MAX(0.001, mutationMatrices[i * matSize + (stateB + row_index[stateA])]));  
                    rowSum += mutationMatrices[i * matSize + (stateB + row_index[stateA])];
                }
            }
            mutationMatrices[i * matSize + (stateA + row_index[stateA])] = -rowSum;
        }
    }

    delete[] totalRates;
    for(int i = 0; i < genomeSize; i++) {
        delete[] C[i];
        delete[] W[i];
    }
    delete[] C;
    delete[] W;
}