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
                auto state = tree->aln->ref_seq[static_cast<std::vector<cmaple::StateType>::size_type>(end_pos)];
                for(int i = pos; i <= end_pos; i++) {
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
            rates[i] = 0.01;
        } else {
            RealNumType expectedRateNoSubstitution = 0;
            for(int j = 0; j < num_states_; j++) {
                RealNumType summand = totals[j][i] * diagonal_mut_mat[j];
                expectedRateNoSubstitution += summand;
            }
            if(expectedRateNoSubstitution == 0) {
                rates[i] = 0.01;
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
}