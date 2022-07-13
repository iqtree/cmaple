#include "node.h"

Node::Node(bool is_top_node)
{
    id = -1;
    seq_name = "";
    length = 0;
    if (is_top_node)
        double_attributes[IS_TOP_NODE] = 1;
    next = NULL;
    neighbor = NULL;
    partial_lh = NULL;
    total_lh = NULL;
    dirty = false;
}

Node::Node(int n_id, string n_seq_name)
{
    id = n_id;
    seq_name = n_seq_name;
    length = 0;
    next = NULL;
    neighbor = NULL;
    partial_lh = NULL;
    total_lh = NULL;
    dirty = false;
}

Node::Node(string n_seq_name)
{
    id = -1;
    seq_name = n_seq_name;
    length = 0;
    double_attributes[IS_TOP_NODE] = 1;
    next = NULL;
    neighbor = NULL;
    partial_lh = NULL;
    total_lh = NULL;
    dirty = false;
}

Node::~Node()
{
    // delete other neighbors
    // NHANLT: TODO
}

bool Node::isLeave()
{
    return next == NULL;
}

Node* Node::getTopNode()
{
    Node* next_node;
    Node* node = this;
    
    if (node->double_attributes.find(IS_TOP_NODE) != node->double_attributes.end())
        return node;
    
    FOR_NEXT(node, next_node)
    {
        if (next_node->double_attributes.find(IS_TOP_NODE) != next_node->double_attributes.end())
            return next_node;
    }
    
    return NULL;
}

string Node::exportString()
{
    if (isLeave())
    {
        // without minor sequences -> simply return node's name and its branch length
        if (less_info_seqs.size() == 0)
            return seq_name + ":" + convertDoubleToString(length);
        // with minor sequences -> return minor sequences' names with zero branch lengths
        else
        {
            string output = "(" + seq_name + ":0";
            for (string minor_seq_name : less_info_seqs)
                output += "," + minor_seq_name + ":0";
            output += "):" + convertDoubleToString(length);
            
            return output;
        }
    }
    
    return "";
}
