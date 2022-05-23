#include "node.h"

Node::Node()
{
    id = -1;
    seq_name = "";
    length = 0;
    relative = NULL;
    neighbor = NULL;
}

Node::Node(int n_id, string n_seq_name)
{
    id = n_id;
    seq_name = n_seq_name;
    length = 0;
    relative = NULL;
    neighbor = NULL;
}

Node::~Node()
{
    // delete lower_lh_seq
    for (Region* region:lower_lh_seq)
        delete region;
    lower_lh_seq.resize(0);
    
    // delete overall_lh_seq
    for (Region* region:overall_lh)
        delete region;
    overall_lh.resize(0);
        
    // delete other neighbors
    // NHANLT: TODO
}

bool Node::isLeave()
{
    return relative == NULL;
}
