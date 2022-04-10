#include "node.h"

Node::Node()
{
    id = -1;
    seq_name = "";
    length = 0;
    parent = NULL;
}

Node::Node(string n_seq_name, Sequence sequence)
{
    id = -1;
    seq_name = n_seq_name;
    length = 0;
    parent = NULL;
    
    // convert sequence (vector of mutations) into genome list (vector of regions)
    // NHANLT: TODO
}
