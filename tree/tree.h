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
        Setup function pointers
     */
    void setupFunctionPointers();
    
    /**
        Setup function pointers
     */
    void setupBlengthThresh();
    
    /**
        Pointer  to updatePartialLh method
     */
    typedef void (Tree::*UpdatePartialLhPointerType)(std::stack<Node*>&);
    UpdatePartialLhPointerType updatePartialLhPointer;
    
    /**
        Template of updatePartialLh
     */
    template <const StateType num_states>
    void updatePartialLhTemplate(std::stack<Node*> &node_stack);
    
    /**
        Pointer  to updatePartialLh method
     */
    typedef RealNumType (Tree::*CalculatePlacementCostType)(const SeqRegions* const, const SeqRegions* const, RealNumType);
    CalculatePlacementCostType calculateSamplePlacementCostPointer;
    
    /**
        Template of calculateSamplePlacementCost
     */
    template <const StateType num_states>
    RealNumType calculateSamplePlacementCostTemplate(const SeqRegions* const parent_regions, const SeqRegions* const child_regions, RealNumType blength);
    
    /**
        Pointer  to updatePartialLh method
     */
    CalculatePlacementCostType calculateSubTreePlacementCostPointer;
    
    /**
        Template of calculateSubTreePlacementCost
     */
    template <const StateType num_states>
    RealNumType calculateSubTreePlacementCostTemplate(const SeqRegions* const parent_regions, const SeqRegions* const child_regions, RealNumType blength);
    
    /**
        Traverse the intial tree from root to re-calculate all lower likelihoods regarding the latest/final estimated model parameters
     */
    void refreshAllLowerLhs();
    
    /**
        Traverse the intial tree from root to re-calculate all non-lower likelihoods regarding the latest/final estimated model parameters
     */
    void refreshAllNonLowerLhs();
    
    /**
        Try to improve a subtree rooted at node with SPR moves
        @return total improvement
     */
    RealNumType improveSubTree(Node* node, bool short_range_search);
    
    /**
       Calculate derivative starting from coefficients.
       @return derivative
    */
    RealNumType calculateDerivative(const std::vector<RealNumType> &coefficient_vec, const RealNumType delta_t);

    /**
       Examine placing a sample at a mid-branch point
    */
    void examineSamplePlacementMidBranch(Node* &selected_node, RealNumType &best_lh_diff, bool& is_mid_branch, RealNumType& lh_diff_mid_branch, TraversingNode& current_extended_node, const SeqRegions* const sample_regions);
    
    /**
       Examine placing a sample as a descendant of an existing node
    */
    void examineSamplePlacementAtNode(Node* &selected_node, RealNumType &best_lh_diff, bool& is_mid_branch, RealNumType& lh_diff_at_node, RealNumType& lh_diff_mid_branch, RealNumType &best_up_lh_diff, RealNumType &best_down_lh_diff, Node* &best_child, TraversingNode& current_extended_node, const SeqRegions* const sample_regions);
    
    /**
       Traverse downwards polytomy for more fine-grained placement
    */
    void finetuneSamplePlacementAtNode(const Node* const selected_node, RealNumType &best_down_lh_diff, Node* &best_child, const SeqRegions* const sample_regions);
    
    /**
       Add start nodes for seeking a placement for a subtree
    */
    void addStartingNodes(const Node* const node, Node* const other_child_node, const RealNumType threshold_prob, SeqRegions* &parent_upper_lr_regions, const RealNumType best_lh_diff, std::stack<UpdatingNode*> &node_stack);
    
    /**
       Examine placing a subtree at a mid-branch point
    */
    bool examineSubtreePlacementMidBranch(Node* &best_node, RealNumType &best_lh_diff, bool& is_mid_branch, RealNumType& lh_diff_at_node, RealNumType& lh_diff_mid_branch, RealNumType &best_up_lh_diff, RealNumType &best_down_lh_diff, UpdatingNode* const updating_node, const SeqRegions* const subtree_regions, const RealNumType threshold_prob, const RealNumType removed_blength, Node* const top_node, SeqRegions* &bottom_regions);
    
    /**
       Examine placing a subtree as a descendant of an existing node
    */
    bool examineSubTreePlacementAtNode(Node* &best_node, RealNumType &best_lh_diff, bool& is_mid_branch, RealNumType& lh_diff_at_node, RealNumType& lh_diff_mid_branch, RealNumType &best_up_lh_diff, RealNumType &best_down_lh_diff, UpdatingNode* const updating_node, const SeqRegions* const subtree_regions, const RealNumType threshold_prob, const RealNumType removed_blength, Node* const top_node);
    
    /**
       Add a child node for downwards traversal when seeking a new subtree placement
    */
    void addChildSeekSubtreePlacement(Node* const child_1, Node* const child_2, const RealNumType& lh_diff_at_node, UpdatingNode* const updating_node, std::stack<UpdatingNode*>& node_stack, const RealNumType threshold_prob);
    
    /**
       Add neighbor nodes (parent/sibling) for traversal when seeking a new subtree placement
    */
    bool addNeighborsSeekSubtreePlacement(Node* const top_node, Node* const other_child, SeqRegions* &parent_upper_lr_regions, SeqRegions* &bottom_regions, const RealNumType& lh_diff_at_node, UpdatingNode* const updating_node, std::stack<UpdatingNode*>& node_stack, const RealNumType threshold_prob);
    
    /**
        Check whether we can obtain a higher likelihood with a shorter length for an existing branch
     */
    template <RealNumType(Tree::*calculatePlacementCost)(const SeqRegions* const, const SeqRegions* const, RealNumType)>
    bool tryShorterBranch(const RealNumType current_blength, SeqRegions* &best_child_regions, const SeqRegions* const sample, const SeqRegions* const upper_left_right_regions, const SeqRegions* const lower_regions, RealNumType &best_split_lh, RealNumType &best_branch_length_split, const RealNumType new_branch_length, const bool try_first_branch);
    
    /**
        Check whether we can obtain a higher likelihood with a shorter length at root
     */
    void tryShorterBranchAtRoot(const SeqRegions* const sample, const SeqRegions* const lower_regions, SeqRegions* &best_parent_regions, RealNumType &best_root_blength, RealNumType &best_parent_lh, const RealNumType fixed_blength);
    
    /**
        Check whether we can obtain a higher likelihood with a shorter length for the new branch at root
     */
    bool tryShorterNewBranchAtRoot(const SeqRegions* const sample, const SeqRegions* const lower_regions, SeqRegions* &best_parent_regions, RealNumType &best_root_blength, RealNumType &best_parent_lh, const RealNumType fixed_blength);
    
    /**
        Check whether we can obtain a higher likelihood with a longer length for the new branch at root
     */
    bool tryLongerNewBranchAtRoot(const SeqRegions* const sample, const SeqRegions* const lower_regions, SeqRegions* &best_parent_regions, RealNumType &best_length, RealNumType &best_parent_lh, const RealNumType fixed_blength);
    
    /**
        Estimate the length for a new branch at root
     */
    void estimateLengthNewBranchAtRoot(const SeqRegions* const sample, const SeqRegions* const lower_regions, SeqRegions* &best_parent_regions, RealNumType &best_length, RealNumType &best_parent_lh, const RealNumType fixed_blength, const RealNumType short_blength_thresh, const bool optional_check);
    
    /**
        Check whether we can obtain a higher likelihood with a shorter length for the new branch
     */
    template <RealNumType(Tree::*calculatePlacementCost)(const SeqRegions* const, const SeqRegions* const, RealNumType)>
    bool tryShorterNewBranch(const SeqRegions* const best_child_regions, const SeqRegions* const sample, RealNumType &best_blength, RealNumType &new_branch_lh, const RealNumType short_blength_thresh);
    
    /**
        Check whether we can obtain a higher likelihood with a longer length for the new branch
     */
    template <RealNumType(Tree::*calculatePlacementCost)(const SeqRegions* const, const SeqRegions* const, RealNumType)>
    void tryLongerNewBranch(const SeqRegions* const best_child_regions, const SeqRegions* const sample, RealNumType &best_blength, RealNumType &new_branch_lh, const RealNumType long_blength_thresh);
    
    /**
        Estimate the length for a new branch
     */
    template <RealNumType(Tree::*calculatePlacementCost)(const SeqRegions* const, const SeqRegions* const, RealNumType)>
    void estimateLengthNewBranch(const RealNumType best_split_lh, const SeqRegions* const best_child_regions, const SeqRegions* const sample, RealNumType &best_blength, const RealNumType long_blength_thresh, const RealNumType short_blength_thresh, const bool optional_check);
    
    /**
        Connect a new sample to a branch
     */
    void connectNewSample2Branch(SeqRegions* const sample, const std::string &seq_name, Node* const sibling_node, const RealNumType top_distance, const RealNumType down_distance, const RealNumType best_blength, SeqRegions* &best_child_regions, const SeqRegions* const upper_left_right_regions);
    
    /**
        Connect a new sample to root
     */
    void connectNewSample2Root(SeqRegions* const sample, const std::string &seq_name, Node* const sibling_node, const RealNumType best_root_blength, const RealNumType best_length2, SeqRegions* &best_parent_regions);
    
    /**
        Place a subtree as a descendant of a node
     */
    void placeSubTreeAtNode(Node* const selected_node, Node* const subtree, const SeqRegions* const subtree_regions, const RealNumType new_branch_length, const RealNumType new_lh);
    
    /**
        Place a subtree at a mid-branch point
     */
    void placeSubTreeMidBranch(Node* const selected_node, Node* const subtree, const SeqRegions* const subtree_regions, const RealNumType new_branch_length, const RealNumType new_lh);
    
    /**
        Connect a subtree to a branch
     */
    template<void (Tree::*updateRegionsSubTree)(Node* const, Node* const, Node* const, Node* const, SeqRegions* &, const SeqRegions* const, const  SeqRegions* const, const SeqRegions* const, RealNumType&)>
    void connectSubTree2Branch(const SeqRegions* const subtree_regions, const SeqRegions* const lower_regions, Node* const subtree, Node* const sibling_node, const RealNumType top_distance, const RealNumType down_distance, RealNumType &best_blength, SeqRegions* &best_child_regions, const SeqRegions* const upper_left_right_regions);
    
    /**
        Connect a subtree to root
     */
    void connectSubTree2Root(Node* const subtree, const SeqRegions* const subtree_regions, const SeqRegions* const lower_regions, Node* const sibling_node, const RealNumType best_root_blength, const RealNumType best_length2, SeqRegions* &best_parent_regions);
    
    /**
        Update next_node_1->partial_lh and new_internal_node->partial_lh after placing a subtree in common cases (e.g., at a mid-branch point, under a node)
     */
    void updateRegionsPlaceSubTree(Node* const subtree, Node* const next_node_1, Node* const sibling_node, Node* const new_internal_node, SeqRegions* &best_child_regions, const SeqRegions* const subtree_regions, const  SeqRegions* const upper_left_right_regions, const SeqRegions* const lower_regions, RealNumType &best_length);
    
    /**
        Update next_node_1->partial_lh and new_internal_node->partial_lh after placing a subtree in other cases (i.e., above a node)
     */
    void updateRegionsPlaceSubTreeAbove(Node* const subtree, Node* const next_node_1, Node* const sibling_node, Node* const new_internal_node, SeqRegions* &best_child_regions, const SeqRegions* const subtree_regions, const  SeqRegions* const upper_left_right_regions, const SeqRegions* const lower_regions, RealNumType &best_length);
    
    /**
        Handle polytomy when placing a subtree
     */
    void handlePolytomyPlaceSubTree(Node* const selected_node, const SeqRegions* const subtree_regions, const RealNumType new_branch_length, RealNumType &best_down_lh_diff, Node* &best_child, RealNumType &best_child_blength_split, SeqRegions* &best_child_regions);
    
    /**
        Update likelihood at mid-branch point
     */
    void updateMidBranchLh(Node* const node, const SeqRegions* const parent_upper_regions, std::stack<Node*> &node_stack, bool &update_blength);
    
    /**
        Compute Upper Left/Right regions at a node, updating the top branch length if neccessary
     */
    SeqRegions* computeUpperLeftRightRegions(Node* const next_node, Node* const node, const SeqRegions* const parent_upper_regions, std::stack<Node*> &node_stack, bool &update_blength);
    
    /**
        Update the PartialLh (seqregions) at a node if the new one is different from the current one
     */
    void updateNewPartialIfDifferent(Node* const next_node, SeqRegions* &upper_left_right_regions, std::stack<Node*> &node_stack, const PositionType seq_length);
    
    /**
        Handle cases when the new seqregions is null/empty: (1) update the branch length; or (2) return an error message
     */
    void inline handleNullNewRegions(Node* const node_update_zero_blength, const bool do_update_zeroblength, std::stack<Node*> &node_stack, bool &update_blength, const std::string err_msg)
    {
        if (do_update_zeroblength)
        {
            updateZeroBlength(node_update_zero_blength, node_stack, params->threshold_prob);
            update_blength = true;
        }
        else
            outError(err_msg);
    }
    
    /**
        Replace a seqregions by a new one
     */
    inline void replacePartialLH(SeqRegions* &old_regions, SeqRegions* &new_regions)
    {
        if (old_regions) delete old_regions;
        old_regions = new_regions;
        new_regions = NULL;
    }
    
    /**
        Update partial_lh comming from the parent node
     */
    void updatePartialLhFromParent(PhyloNode& node, const Index index, std::stack<Node*> &node_stack, const SeqRegions* const parent_upper_regions, const PositionType seq_length);
    
    /**
        Update partial_lh comming from the children
     */
    void updatePartialLhFromChildren(PhyloNode& node, const Index index, std::stack<Node*> &node_stack, const SeqRegions* const parent_upper_regions, const bool is_non_root, const PositionType seq_length);
    
    /**
        Compute the mid-branch region for a node/branch
     */
    inline void computeMidBranchRegions(Node* const node, SeqRegions* &regions_2_update, const SeqRegions &parent_upper_lr_lh)
    {
        /*SeqRegions* lower_lh = node->getPartialLhAtNode(aln, model, params->threshold_prob);
        RealNumType half_branch_length = node->length * 0.5;
        parent_upper_lr_lh.mergeUpperLower(regions_2_update, half_branch_length, *lower_lh, half_branch_length, aln, model, params->threshold_prob);*/
    }
    
    /**
        Refresh all non-lowerlhs traversing from a parent node
     */
    void refreshNonLowerLhsFromParent(Node* &node, Node* &last_node);
    
    /**
        Refresh upper left/right regions
     */
    void refreshUpperLR(Node* const node, Node* const next_node, SeqRegions* &replaced_regions, const SeqRegions &parent_upper_lr_lh);
    
    /**
        Calculate coefficients when merging R with O to estimate a branch length
     */
    void estimateBlength_R_O(const SeqRegion* const seq1_region, const SeqRegion* const seq2_region, const RealNumType total_blength, const PositionType end_pos, RealNumType &coefficient, std::vector<RealNumType> &coefficient_vec);
    
    /**
        Calculate coefficients when merging R with ACGT to estimate a branch length
     */
    void estimateBlength_R_ACGT(const SeqRegion* const seq1_region, const StateType seq2_state, const RealNumType total_blength, const PositionType end_pos, std::vector<RealNumType> &coefficient_vec);
    
    /**
        Calculate coefficients when merging O with X(i.e., O, R, ACGT) to estimate a branch length
     */
    void estimateBlength_O_X(const SeqRegion* const seq1_region, const SeqRegion* const seq2_region, const RealNumType total_blength, const PositionType end_pos, RealNumType &coefficient, std::vector<RealNumType> &coefficient_vec);
    
    /**
        Calculate coefficients when merging ACGT with O to estimate a branch length
     */
    void estimateBlength_ACGT_O(const SeqRegion* const seq1_region, const SeqRegion* const seq2_region, const RealNumType total_blength, RealNumType &coefficient, std::vector<RealNumType> &coefficient_vec);
    
    /**
        Calculate coefficients when merging ACGT with R/ACGT to estimate a branch length
     */
    void estimateBlength_ACGT_RACGT(const SeqRegion* const seq1_region, const SeqRegion* const seq2_region, const RealNumType total_blength, const PositionType end_pos, std::vector<RealNumType> &coefficient_vec);
    
    /**
        Estimate a branch length from coefficients
     */
    RealNumType estimateBlengthFromCoeffs(RealNumType &coefficient, const std::vector<RealNumType> coefficient_vec);
    
    /**
        Handle branch length changed when improve a subtree
     */
    void handleBlengthChanged(Node* const node, const RealNumType best_blength);
    
    /**
        Optimize a branch length before seeking an SPR move for a subtree
     */
    void optimizeBlengthBeforeSeekingSPR(Node* const node, RealNumType &best_blength, RealNumType &best_lh, bool &blength_changed, const SeqRegions* const parent_upper_lr_lh, const SeqRegions* const lower_lh);
    
    /**
        Check and apply SPR move
     */
    void checkAndApplySPR(const RealNumType best_lh_diff, const RealNumType best_blength, const RealNumType best_lh, Node* const node, Node* const best_node, Node* const parent_node, const bool is_mid_node, RealNumType& total_improvement, bool& topology_updated);
    
public:
    /*
     Branch length thresholds
     */
    RealNumType default_blength, min_blength, max_blength, min_blength_mid, min_blength_sensitivity, half_max_blength, half_min_blength_mid, double_min_blength;
    
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
    
    /*
        Vector of phylonodes
     */
    std::vector<PhyloNode> nodes;
    
    /*
        (vector) Index of root in the vector of phylonodes
     */
    uint32_t root_vector_index;
    
    /**
        Constructor
     */
    Tree() = default;
    
    /**
        Constructor
    */
    Tree(Params&& params);
    
    /**
        Destructor
     */
    ~Tree();
    
    /**
        Setup tree parameters/thresholds
     */
    void setup();
    
    /**
        Export tree std::string in Newick format
     */
    std::string exportTreeString(bool binary = false, Node* node = NULL);
    
    /**
        Increase the length of a 0-length branch (connecting this node to its parent) to resolve the inconsistency when updating regions in updatePartialLh()
     */
    void updateZeroBlength(Node* node, std::stack<Node*> &node_stack, RealNumType threshold_prob);
    
    /**
        Iteratively update partial_lh starting from the nodes in node_stack
        @param node_stack stack of nodes;
     */
    void updatePartialLh(std::stack<Node*> &node_stack);
    
    /**
        Seek a position for a sample placement starting at the start_node
     */
    void seekSamplePlacement(Node* start_node, const std::string &seq_name, SeqRegions* sample_regions, Node* &selected_node, RealNumType &best_lh_diff , bool &is_mid_branch, RealNumType &best_up_lh_diff, RealNumType &best_down_lh_diff, Node* &best_child);
    
    /**
        Seek a position for placing a subtree/sample starting at the start_node
     */
    void seekSubTreePlacement(Node* &best_node, RealNumType &best_lh_diff, bool &is_mid_branch, RealNumType &best_up_lh_diff, RealNumType &best_down_lh_diff, Node* &best_child, bool short_range_search, Node* start_node, RealNumType &removed_blength, bool search_subtree_placement = true, SeqRegions* sample_regions = NULL);
    
    /**
        Place a new sample at a mid-branch point
     */
    void placeNewSampleMidBranch(Node* const selected_node, SeqRegions* const sample, const std::string &seq_name, const RealNumType best_lh_diff);
    
    /**
        Place a new sample as a descendant of a node
     */
    void placeNewSampleAtNode(Node* const selected_node, SeqRegions* const sample, const std::string &seq_name, const RealNumType best_lh_diff, const RealNumType best_up_lh_diff, const RealNumType best_down_lh_diff, Node* const best_child);
    
    /**
        Apply SPR move
        pruning a subtree then regrafting it to a new position
     */
    void applySPR(Node* const subtree, Node* const best_node, const bool is_mid_branch, const RealNumType branch_length, const RealNumType best_lh_diff);
    
    /**
        Traverse the intial tree from root to re-calculate all likelihoods regarding the latest/final estimated model parameters
     */
    void refreshAllLhs();
    
    /**
        Traverse the tree to set the outdated flag (which is used to prevent traversing the same part of the tree multiple times) of all nodes to TRUE
     */
    void setAllNodeOutdated();
    
    /**
        Try to improve the entire tree with SPR moves
        @return total improvement
     */
    RealNumType improveEntireTree(bool short_range_search);
    
    /**
        Try to optimize branch lengths of the tree
        @return num of improvements
     */
    PositionType optimizeBranchLengths();
    
    /**
        Estimate the length of a branch using the derivative of the likelihood cost function wrt the branch length
     */
    RealNumType estimateBranchLength(const SeqRegions* const parent_regions, const SeqRegions* const child_regions);
    
    /**
        Calculate the placement cost of a sample
        @param child_regions: vector of regions of the new sample
     */
    RealNumType calculateSamplePlacementCost(const SeqRegions* const parent_regions, const SeqRegions* const child_regions, RealNumType blength);
    
    /**
        Calculate the placement cost of a subtree
        @param child_regions: vector of regions of the new sample
     */
    RealNumType calculateSubTreePlacementCost(const SeqRegions* const parent_regions, const SeqRegions* const child_regions, RealNumType blength);
};

#endif
