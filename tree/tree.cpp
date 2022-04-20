#include "tree.h"

Tree::Tree()
{
    root = NULL;
}

Tree::Tree(Node* node)
{
    root = node;
}

Tree::~Tree()
{
    // browse tree to delete all nodes
    // NHANLT: TODO
}
