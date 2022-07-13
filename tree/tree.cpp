#include "tree.h"

Tree::Tree()
{
    params = NULL;
    aln = new Alignment();
    model = new Model();
    root = NULL;
}

Tree::Tree(Params *n_params, Node* n_root)
{
    params = n_params;
    aln = new Alignment();
    model = new Model();
    root = n_root;
}

Tree::~Tree()
{
    if (aln)
    {
        delete aln;
        aln = NULL;
    }
    
    if (model)
    {
        delete model;
        model = NULL;
    }
    
    // browse tree to delete all nodes
    // NHANLT: TODO
}

void Tree::addNode(Node* new_node, Node* current_node, int &new_id, bool insert_to_right)
{
    Node* new_node_1 = new Node(new_id++);
    Node* new_node_2 = new Node(new_id++);
    Node* new_node_3 = new Node(new_id++);
    new_node_1->next = new_node_2;
    new_node_2->next = new_node_3;
    new_node_3->next = new_node_1;
    
    new_node_1->neighbor = current_node->neighbor;
    if (current_node == root)
        root = new_node_1;
    else
        new_node_1->neighbor->neighbor = new_node_1;
    
    Node *left_tip, *right_tip;
    if (insert_to_right)
    {
        left_tip = current_node;
        right_tip = new_node;
    }
    else
    {
        left_tip = new_node;
        right_tip = current_node;
    }
    
    new_node_2->neighbor = right_tip;
    right_tip->neighbor = new_node_2;
    
    new_node_3->neighbor = left_tip;
    left_tip->neighbor = new_node_3;
}

string Tree::exportTreeString(Node* node)
{
    // init starting node from root
    if (!node)
        node = root;
    
    // move to its neighbor
    if (node->neighbor)
        node = node->neighbor;
    
    // do something with its neighbor
    if (node->isLeave())
        return node->exportString();
        
    string output = "(";
    bool add_comma = false;
    Node* next;
    FOR_NEXT(node, next)
    {
        if (!add_comma)
            add_comma = true;
        else
            output += ",";
        output += exportTreeString(next);
    }
    output += "):" + convertDoubleToString(node->length);
    
    return output;
}
