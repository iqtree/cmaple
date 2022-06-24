#include "alignment/regions.h"
#include "alignment/sequence.h"
#include "alignment/alignment.h"
#include "model/model.h"

#ifndef NODE_H
#define NODE_H

class Node {
public:
    // node's id
    int id;
    
    // sequence name
    string seq_name;
    
    // likelihood of the data arrived from other next nodes
    Regions* partial_lh = NULL;
    
    // total likelihood at the current node
    Regions* total_lh = NULL;
    
    // length of branch connecting to parent
    double length;
    
    // next node in the circle of neighbors. For tips, circle = NULL
    Node* next;
    
    // next node in the phylo tree
    Node* neighbor;
    
    // vector of less informative sequences
    vector<string> less_info_seqs;

    // flexible string attributes
    map<string,string> str_attributes;
    
    // flexible double attributes
    map<string,double> double_attributes;
    
    /**
        constructor
     */
    Node();
    
    /**
        constructor
        @param n_seq_name the sequence name
     */
    Node(int n_id = -1, string n_seq_name = "");
    
    /**
        deconstructor
     */
    ~Node();
    
    /**
        TRUE if this node is a leaf
     */
    bool isLeave();
    
    /**
        compute total lh for root node
     */
    void computeTotalLhForRoot(Model* model, StateType num_states, double blength = 0);
};

#endif

#define FOR_NEXT(node, next_node) \
for(next_node = node->next; next_node && next_node != node; next_node = next_node->next)

#define FOR_NEIGHBOR(node, neighbor_node) \
if (node->next) \
for(neighbor_node = node->next->neighbor; neighbor_node && neighbor_node != node->neighbor; neighbor_node = neighbor_node->neighbor->next->neighbor)
