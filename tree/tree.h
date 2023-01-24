#include "node.h"
#include "alignment/alignment.h"
#include "model/model.h"
#include <optional>

#ifndef TREE_H
#define TREE_H

/** The tree structure */
class Tree {
private:
    /**
        Pointer  to updatePartialLh method
     */
    typedef void (Tree::*UpdatePartialLhPointerType)(std::stack<Node*> &, RealNumType*, RealNumType, RealNumType, RealNumType);
    UpdatePartialLhPointerType updatePartialLhPointer;
    
    /**
        Template of updatePartialLh
     */
    template <const StateType num_states>
    void updatePartialLhTemplate(std::stack<Node*> &node_stack, RealNumType* cumulative_rate, RealNumType default_blength, RealNumType min_blength, RealNumType max_blength);
    
    /**
        Pointer  to updatePartialLh method
     */
    typedef RealNumType (Tree::*CalculatePlacementCostType)(RealNumType*, const SeqRegions* const, const SeqRegions* const, RealNumType);
    CalculatePlacementCostType calculateSamplePlacementCostPointer;
    
    /**
        Template of calculateSamplePlacementCost
     */
    template <const StateType num_states>
    RealNumType calculateSamplePlacementCostTemplate(RealNumType* cumulative_rate, const SeqRegions* const parent_regions, const SeqRegions* const child_regions, RealNumType blength);
    
    /**
        Pointer  to updatePartialLh method
     */
    CalculatePlacementCostType calculateSubTreePlacementCostPointer;
    
    /**
        Template of calculateSubTreePlacementCost
     */
    template <const StateType num_states>
    RealNumType calculateSubTreePlacementCostTemplate(RealNumType* cumulative_rate, const SeqRegions* const parent_regions, const SeqRegions* const child_regions, RealNumType blength);
    
    /**
        Traverse the intial tree from root to re-calculate all lower likelihoods regarding the latest/final estimated model parameters
     */
    void refreshAllLowerLhs(RealNumType *cumulative_rate, RealNumType default_blength, RealNumType max_blength, RealNumType min_blength);
    
    /**
        Traverse the intial tree from root to re-calculate all non-lower likelihoods regarding the latest/final estimated model parameters
     */
    void refreshAllNonLowerLhs(RealNumType *cumulative_rate, RealNumType default_blength, RealNumType max_blength, RealNumType min_blength);
    
    /**
        Try to improve a subtree rooted at node with SPR moves
        @return total improvement
     */
    RealNumType improveSubTree(Node* node, bool short_range_search, RealNumType *cumulative_rate, std::vector< std::vector<PositionType> > &cumulative_base, RealNumType default_blength, RealNumType max_blength, RealNumType min_blength, RealNumType min_blength_mid);
    
    /**
       Calculate derivative starting from coefficients.
       @return derivative
    */
    RealNumType calculateDerivative(std::vector<RealNumType> &coefficient_vec, RealNumType delta_t);

    /**
       Examine placing a sample at a mid-branch point
    */
    void examineSamplePlacementMidBranch(Node* &selected_node, RealNumType &best_lh_diff, bool& is_mid_branch, RealNumType& lh_diff_mid_branch, RealNumType* cumulative_rate, TraversingNode& current_extended_node, const SeqRegions* const sample_regions, const RealNumType default_blength);
    
    /**
       Examine placing a sample as a descendant of an existing node
    */
    void examineSamplePlacementAtNode(Node* &selected_node, RealNumType &best_lh_diff, bool& is_mid_branch, RealNumType& lh_diff_at_node, RealNumType& lh_diff_mid_branch, RealNumType &best_up_lh_diff, RealNumType &best_down_lh_diff, Node* &best_child, RealNumType* cumulative_rate, TraversingNode& current_extended_node, const SeqRegions* const sample_regions, const RealNumType default_blength);\
    
    /**
       Traverse downwards polytomy for more fine-grained placement
    */
    void finetuneSamplePlacementAtNode(const Node* const selected_node, RealNumType &best_down_lh_diff, Node* &best_child, RealNumType* cumulative_rate, const SeqRegions* const sample_regions, const RealNumType default_blength, const RealNumType min_blength_mid);
    
    /**
       Add start nodes for seeking a placement for a subtree
    */
    void addStartingNodes(const Node* const node, Node* const other_child_node, RealNumType* cumulative_rate, const RealNumType threshold_prob, SeqRegions* &parent_upper_lr_regions, const RealNumType best_lh_diff, std::stack<UpdatingNode*> &node_stack);
    
    /**
       Examine placing a subtree at a mid-branch point
    */
    bool examineSubtreePlacementMidBranch(Node* &best_node, RealNumType &best_lh_diff, bool& is_mid_branch, RealNumType& lh_diff_at_node, RealNumType& lh_diff_mid_branch, RealNumType &best_up_lh_diff, RealNumType &best_down_lh_diff, RealNumType* cumulative_rate, UpdatingNode* const updating_node, const SeqRegions* const subtree_regions, const RealNumType threshold_prob, const RealNumType removed_blength, Node* const top_node, SeqRegions* &bottom_regions);
    
    /**
       Examine placing a subtree as a descendant of an existing node
    */
    bool examineSubTreePlacementAtNode(Node* &best_node, RealNumType &best_lh_diff, bool& is_mid_branch, RealNumType& lh_diff_at_node, RealNumType& lh_diff_mid_branch, RealNumType &best_up_lh_diff, RealNumType &best_down_lh_diff, RealNumType* cumulative_rate, UpdatingNode* const updating_node, const SeqRegions* const subtree_regions, const RealNumType threshold_prob, const RealNumType removed_blength, Node* const top_node);
    
    /**
       Add a child node for downwards traversal when seeking a new subtree placement
    */
    void addChildSeekSubtreePlacement(Node* const child_1, Node* const child_2, const RealNumType& lh_diff_at_node, RealNumType* cumulative_rate, UpdatingNode* const updating_node, std::stack<UpdatingNode*>& node_stack, const RealNumType threshold_prob);
    
    /**
       Add neighbor nodes (parent/sibling) for traversal when seeking a new subtree placement
    */
    bool addNeighborsSeekSubtreePlacement(Node* const top_node, Node* const other_child, SeqRegions* &parent_upper_lr_regions, SeqRegions* &bottom_regions, const RealNumType& lh_diff_at_node, RealNumType* cumulative_rate, UpdatingNode* const updating_node, std::stack<UpdatingNode*>& node_stack, const RealNumType threshold_prob);
    
    /**
        Place a new sample at a mid-branch point
     */
    void placeNewSampleMidBranch(Node* const selected_node, SeqRegions* const sample, const std::string &seq_name, const RealNumType best_lh_diff, RealNumType* cumulative_rate, const RealNumType default_blength, const RealNumType max_blength, const RealNumType min_blength);
    
    /**
        Place a new sample as a descendant of a node
     */
    void placeNewSampleAtNode(Node* const selected_node, SeqRegions* const sample, const std::string &seq_name, const RealNumType best_lh_diff, const RealNumType best_up_lh_diff, const RealNumType best_down_lh_diff, Node* const best_child, RealNumType* cumulative_rate, std::vector< std::vector<PositionType> > &cumulative_base, const RealNumType default_blength, const RealNumType max_blength, const RealNumType min_blength);
    
    /**
        Check whether we can obtain a higher likelihood with a shorter length for an existing branch
     */
    bool tryShorterBranch(const RealNumType current_blength, SeqRegions* &best_child_regions, const SeqRegions* const sample, const SeqRegions* const upper_left_right_regions, const SeqRegions* const lower_regions, RealNumType* cumulative_rate, RealNumType &best_split_lh, RealNumType &best_branch_length_split, const RealNumType default_blength, const RealNumType min_blength, const bool tryFirstBranch);
    
    /**
        Check whether we can obtain a higher likelihood with a shorter length for the new branch
     */
    bool tryShorterNewBranch(const SeqRegions* const best_child_regions, const SeqRegions* const sample, RealNumType &best_blength, RealNumType &new_branch_lh, RealNumType* cumulative_rate, const RealNumType min_blength);
    
    /**
        Check whether we can obtain a higher likelihood with a longer length for the new branch
     */
    void tryLongerNewBranch(const SeqRegions* const best_child_regions, const SeqRegions* const sample, RealNumType &best_blength, RealNumType &new_branch_lh, RealNumType* cumulative_rate, const RealNumType max_blength);
    
    /**
        Estimate the length for a new branch
     */
    void estimateLengthNewBranch(const RealNumType best_split_lh, const SeqRegions* const best_child_regions, const SeqRegions* const sample, RealNumType &best_blength, RealNumType* cumulative_rate, const RealNumType min_blength, const RealNumType max_blength);
    
public:
    /**
        Program parameters
     */
    std::optional<Params> params;
    
    /**
        Alignment
     */
    Alignment aln;
    
    /**
        Evolutionary model
     */
    Model model;
    
    /**
        Root node of the tree
     */
    Node* root = nullptr;
    
    /**
        Constructor
     */
    Tree() = default;
    
    /**
        Constructor
    */
    Tree(Params&& params, Node* n_root = nullptr);
    
    /**
        Destructor
     */
    ~Tree();
    
    /**
        Setup function pointers
     */
    void setupFunctionPointers();
    
    /**
        Export tree std::string in Newick format
     */
    std::string exportTreeString(bool binary = false, Node* node = NULL);
    
    /**
        Increase the length of a 0-length branch (connecting this node to its parent) to resolve the inconsistency when updating regions in updatePartialLh()
     */
    void updateZeroBlength(Node* node, std::stack<Node*> &node_stack, RealNumType threshold_prob, RealNumType* cumulative_rate, RealNumType default_blength, RealNumType min_blength, RealNumType max_blength);
    
    /**
        Iteratively update partial_lh starting from the nodes in node_stack
        @param node_stack stack of nodes;
     */
    void updatePartialLh(std::stack<Node*> &node_stack, RealNumType* cumulative_rate, RealNumType default_blength, RealNumType min_blength, RealNumType max_blength);
    
    /**
        Seek a position for a sample placement starting at the start_node
     */
    void seekSamplePlacement(Node* start_node, const std::string &seq_name, SeqRegions* sample_regions, Node* &selected_node, RealNumType &best_lh_diff , bool &is_mid_branch, RealNumType &best_up_lh_diff, RealNumType &best_down_lh_diff, Node* &best_child, RealNumType* cumulative_rate, RealNumType default_blength, RealNumType min_blength_mid);
    
    /**
        Seek a position for placing a subtree/sample starting at the start_node
     */
    void seekSubTreePlacement(Node* &best_node, RealNumType &best_lh_diff, bool &is_mid_branch, RealNumType &best_up_lh_diff, RealNumType &best_down_lh_diff, Node* &best_child, bool short_range_search, Node* start_node, RealNumType &removed_blength, RealNumType* cumulative_rate, RealNumType default_blength, RealNumType min_blength_mid, bool search_subtree_placement = true, SeqRegions* sample_regions = NULL);

    /**
        Place a new sample on the tree
     */
    void placeNewSample(Node* selected_node, SeqRegions* sample, const std::string &seq_name, RealNumType best_lh_diff, bool is_mid_branch, RealNumType best_up_lh_diff, RealNumType best_down_lh_diff, Node* best_child, RealNumType* cumulative_rate, std::vector< std::vector<PositionType> > &cumulative_base, RealNumType default_blength, RealNumType max_blength, RealNumType min_blength);
    
    /**
        Place a subtree
     */
    void placeSubtree(Node* selected_node, Node* subtree, SeqRegions* subtree_regions, bool is_mid_branch, RealNumType branch_length, RealNumType new_lh, RealNumType* cumulative_rate, std::vector< std::vector<PositionType> > &cumulative_base, RealNumType default_blength, RealNumType max_blength, RealNumType min_blength, RealNumType min_blength_mid);
    
    /**
        Apply SPR move
        pruning a subtree then regrafting it to a new position
     */
    void applySPR(Node* subtree, Node* best_node, bool is_mid_branch, RealNumType branch_length, RealNumType best_lh_diff, RealNumType* cumulative_rate, std::vector< std::vector<PositionType> > &cumulative_base, RealNumType default_blength, RealNumType max_blength, RealNumType min_blength, RealNumType min_blength_mid);
    
    /**
        Traverse the intial tree from root to re-calculate all likelihoods regarding the latest/final estimated model parameters
     */
    void refreshAllLhs(RealNumType *cumulative_rate, RealNumType default_blength, RealNumType max_blength, RealNumType min_blength);
    
    /**
        Traverse the tree to set the outdated flag (which is used to prevent traversing the same part of the tree multiple times) of all nodes to TRUE
     */
    void setAllNodeOutdated();
    
    /**
        Try to improve the entire tree with SPR moves
        @return total improvement
     */
    RealNumType improveEntireTree(bool short_range_search, RealNumType *cumulative_rate, std::vector< std::vector<PositionType> > &cumulative_base, RealNumType default_blength, RealNumType max_blength, RealNumType min_blength, RealNumType min_blength_mid);
    
    /**
        Try to optimize branch lengths of the tree
        @return num of improvements
     */
    PositionType optimizeBranchLengths(RealNumType *cumulative_rate, RealNumType default_blength, RealNumType max_blength, RealNumType min_blength, RealNumType min_blength_sensitivity);
    
    /**
        Estimate the length of a branch using the derivative of the likelihood cost function wrt the branch length
     */
    RealNumType estimateBranchLength(const SeqRegions* const parent_regions, const SeqRegions* const child_regions, const RealNumType* const cumulative_rate, RealNumType min_blength_sensitivity);
    
    /**
        Calculate the placement cost of a sample
        @param child_regions: vector of regions of the new sample
     */
    RealNumType calculateSamplePlacementCost(RealNumType* cumulative_rate, const SeqRegions* const parent_regions, const SeqRegions* const child_regions, RealNumType blength);
    
    /**
        Calculate the placement cost of a subtree
        @param child_regions: vector of regions of the new sample
     */
    RealNumType calculateSubTreePlacementCost(RealNumType* cumulative_rate, const SeqRegions* const parent_regions, const SeqRegions* const child_regions, RealNumType blength);

};

#endif
