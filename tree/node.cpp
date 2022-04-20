#include "node.h"

Node::Node()
{
    id = -1;
    seq_name = "";
    length = 0;
    circle = NULL;
    next = NULL;
}

Node::Node(Sequence* sequence, Alignment* aln)
{
    id = -1;
    seq_name = sequence->seq_name;
    length = 0;
    circle = NULL;
    next = NULL;
    
    // convert vector of mutations into vector of regions
    lower_lh_seq = sequence->getRegionsFromMutations(aln->ref_seq.size(), aln->seq_type, aln->num_states);
}

Node::~Node()
{
    // delete lower_lh_seq
    for (Region* region:lower_lh_seq)
        delete region;
    lower_lh_seq.resize(0);
    
    // delete overall_lh_seq
    for (Region* region:overall_lh_seq)
        delete region;
    overall_lh_seq.resize(0);
        
    // delete other neighbors
    // NHANLT: TODO
}
