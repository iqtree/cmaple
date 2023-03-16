#include "updatingnode.h"
#include "alignment/alignment.h"
#include "model/model.h"
#include <optional>
#ifdef _OPENMP
    #include <omp.h>
#endif

#pragma once

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
    typedef void (Tree::*UpdatePartialLhPointerType)(std::stack<Index>&);
    UpdatePartialLhPointerType updatePartialLhPointer;
    
    /**
        Template of updatePartialLh
     */
    template <const StateType num_states>
    void updatePartialLhTemplate(std::stack<Index> &node_stack);
    
    /**
        Pointer  to updatePartialLh method
     */
    typedef RealNumType (Tree::*CalculatePlacementCostType)(const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, const RealNumType);
    CalculatePlacementCostType calculateSamplePlacementCostPointer;
    
    /**
        Template of calculateSamplePlacementCost
     */
    template <const StateType num_states>
    RealNumType calculateSamplePlacementCostTemplate(const std::unique_ptr<SeqRegions>& parent_regions, const std::unique_ptr<SeqRegions>& child_regions, const RealNumType blength);
    
    /**
        Pointer  to updatePartialLh method
     */
    CalculatePlacementCostType calculateSubTreePlacementCostPointer;
    
    /**
        Template of calculateSubTreePlacementCost
     */
    template <const StateType num_states>
    RealNumType calculateSubTreePlacementCostTemplate(const std::unique_ptr<SeqRegions>& parent_regions, const std::unique_ptr<SeqRegions>& child_regions, const RealNumType blength);
    
    /**
        Traverse the intial tree from root to re-calculate all non-lower likelihoods regarding the latest/final estimated model parameters
     */
    void refreshAllNonLowerLhs();
    
    /**
        Try to improve a subtree rooted at node with SPR moves
        @return total improvement
     */
    RealNumType improveSubTree(const Index index, PhyloNode& node, bool short_range_search);
    
    /**
       Calculate derivative starting from coefficients.
       @return derivative
    */
    RealNumType calculateDerivative(const std::vector<RealNumType> &coefficient_vec, const RealNumType delta_t);

    /**
       Examine placing a sample at a mid-branch point
    */
    void examineSamplePlacementMidBranch(Index& selected_node_index, const std::unique_ptr<SeqRegions>& mid_branch_lh, RealNumType &best_lh_diff, bool& is_mid_branch, RealNumType& lh_diff_mid_branch, TraversingNode& current_extended_node, const std::unique_ptr<SeqRegions>& sample_regions);
    
    /**
       Examine placing a sample as a descendant of an existing node
    */
    void examineSamplePlacementAtNode(Index& selected_node_index, const std::unique_ptr<SeqRegions>& total_lh, RealNumType &best_lh_diff, bool& is_mid_branch, RealNumType& lh_diff_at_node, RealNumType& lh_diff_mid_branch, RealNumType &best_up_lh_diff, RealNumType &best_down_lh_diff, Index& best_child_index, TraversingNode& current_extended_node, const std::unique_ptr<SeqRegions>& sample_regions);
    
    /**
       Traverse downwards polytomy for more fine-grained placement
    */
    void finetuneSamplePlacementAtNode(const PhyloNode& selected_node, RealNumType &best_down_lh_diff, Index& best_child_index, const std::unique_ptr<SeqRegions>& sample_regions);
    
    /**
       Add start nodes for seeking a placement for a subtree
    */
    void addStartingNodes(const Index& node_index, PhyloNode& node, const Index& other_child_node_index, const RealNumType best_lh_diff, std::stack<std::unique_ptr<UpdatingNode>>& node_stack);
    
    /**
       Examine placing a subtree at a mid-branch point
    */
    bool examineSubtreePlacementMidBranch(Index& best_node_index, PhyloNode& current_node, RealNumType& best_lh_diff, bool& is_mid_branch, RealNumType& lh_diff_at_node, RealNumType& lh_diff_mid_branch, RealNumType& best_up_lh_diff, RealNumType& best_down_lh_diff, std::unique_ptr<UpdatingNode>& updating_node, const std::unique_ptr<SeqRegions>& subtree_regions, const RealNumType threshold_prob, const RealNumType removed_blength, const Index top_node_index, std::unique_ptr<SeqRegions>& bottom_regions);
    
    /**
       Examine placing a subtree as a descendant of an existing node
    */
    bool examineSubTreePlacementAtNode(Index& best_node_index, PhyloNode& current_node, RealNumType &best_lh_diff, bool& is_mid_branch, RealNumType& lh_diff_at_node, RealNumType& lh_diff_mid_branch, RealNumType &best_up_lh_diff, RealNumType &best_down_lh_diff, std::unique_ptr<UpdatingNode>& updating_node, const std::unique_ptr<SeqRegions>& subtree_regions, const RealNumType threshold_prob, const RealNumType removed_blength, const Index top_node_index);
    
    /**
       Add a child node for downwards traversal when seeking a new subtree placement
    */
    void addChildSeekSubtreePlacement(const Index child_1_index, const Index child_2_index, PhyloNode& child_1, PhyloNode& child_2, const RealNumType& lh_diff_at_node, const std::unique_ptr<UpdatingNode>& updating_node, std::stack<std::unique_ptr<UpdatingNode>>& node_stack, const RealNumType threshold_prob);
    
    /**
       Add neighbor nodes (parent/sibling) for traversal when seeking a new subtree placement
    */
    bool addNeighborsSeekSubtreePlacement(PhyloNode& current_node, const Index other_child_index, std::unique_ptr<SeqRegions>&& bottom_regions, const RealNumType& lh_diff_at_node, const std::unique_ptr<UpdatingNode>& updating_node, std::stack<std::unique_ptr<UpdatingNode>>& node_stack, const RealNumType threshold_prob);
    
    /**
        Check whether we can obtain a higher likelihood with a shorter length for an existing branch
     */
    template <RealNumType(Tree::*calculatePlacementCost)(const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, const RealNumType)>
    bool tryShorterBranch(const RealNumType current_blength, std::unique_ptr<SeqRegions>& best_child_regions, const std::unique_ptr<SeqRegions>& sample, const std::unique_ptr<SeqRegions>& upper_left_right_regions, const std::unique_ptr<SeqRegions>& lower_regions, RealNumType &best_split_lh, RealNumType &best_branch_length_split, const RealNumType new_branch_length, const bool try_first_branch);
    
    /**
        Check whether we can obtain a higher likelihood with a shorter length at root
     */
    void tryShorterBranchAtRoot(const std::unique_ptr<SeqRegions>& sample, const std::unique_ptr<SeqRegions>& lower_regions, std::unique_ptr<SeqRegions>& best_parent_regions, RealNumType &best_root_blength, RealNumType &best_parent_lh, const RealNumType fixed_blength);
    
    /**
        Check whether we can obtain a higher likelihood with a shorter length for the new branch at root
     */
    bool tryShorterNewBranchAtRoot(const std::unique_ptr<SeqRegions>& sample, const std::unique_ptr<SeqRegions>& lower_regions, std::unique_ptr<SeqRegions>&best_parent_regions, RealNumType &best_length, RealNumType &best_parent_lh, const RealNumType fixed_blength);
    
    /**
        Check whether we can obtain a higher likelihood with a longer length for the new branch at root
     */
    bool tryLongerNewBranchAtRoot(const std::unique_ptr<SeqRegions>& sample, const std::unique_ptr<SeqRegions>& lower_regions, std::unique_ptr<SeqRegions>& best_parent_regions, RealNumType &best_length, RealNumType &best_parent_lh, const RealNumType fixed_blength);
    
    /**
        Estimate the length for a new branch at root
     */
    void estimateLengthNewBranchAtRoot(const std::unique_ptr<SeqRegions>& sample, const std::unique_ptr<SeqRegions>& lower_regions, std::unique_ptr<SeqRegions>& best_parent_regions, RealNumType &best_length, RealNumType &best_parent_lh, const RealNumType fixed_blength, const RealNumType short_blength_thresh, const bool optional_check);
    
    /**
        Check whether we can obtain a higher likelihood with a shorter length for the new branch
     */
    template <RealNumType(Tree::*calculatePlacementCost)(const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, const RealNumType)>
    bool tryShorterNewBranch(const std::unique_ptr<SeqRegions>& best_child_regions, const std::unique_ptr<SeqRegions>& sample, RealNumType &best_blength, RealNumType &new_branch_lh, const RealNumType short_blength_thresh);
    
    /**
        Check whether we can obtain a higher likelihood with a longer length for the new branch
     */
    template <RealNumType(Tree::*calculatePlacementCost)(const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, const RealNumType)>
    void tryLongerNewBranch(const std::unique_ptr<SeqRegions>& best_child_regions, const std::unique_ptr<SeqRegions>& sample, RealNumType &best_blength, RealNumType &new_branch_lh, const RealNumType long_blength_thresh);
    
    /**
        Estimate the length for a new branch
     */
    template <RealNumType(Tree::*calculatePlacementCost)(const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, const RealNumType)>
    void estimateLengthNewBranch(const RealNumType best_split_lh, const std::unique_ptr<SeqRegions>& best_child_regions, const std::unique_ptr<SeqRegions>& sample, RealNumType &best_blength, const RealNumType long_blength_thresh, const RealNumType short_blength_thresh, const bool optional_check);
    
    /**
        Connect a new sample to a branch
     */
    void connectNewSample2Branch(std::unique_ptr<SeqRegions>& sample, const NumSeqsType seq_name_index, const Index sibling_node_index, PhyloNode& sibling_node, const RealNumType top_distance, const RealNumType down_distance, const RealNumType best_blength, std::unique_ptr<SeqRegions>& best_child_regions, const std::unique_ptr<SeqRegions>& upper_left_right_regions);
    
    /**
        Connect a new sample to root
     */
    void connectNewSample2Root(std::unique_ptr<SeqRegions>& sample, const NumSeqsType seq_name_index, const Index sibling_node_index, PhyloNode& sibling_node, const RealNumType best_root_blength, const RealNumType best_length2, std::unique_ptr<SeqRegions>& best_parent_regions);
    
    /**
        Place a subtree as a descendant of a node
     */
    void placeSubTreeAtNode(const Index selected_node_index, const Index subtree_index, PhyloNode& subtree, const std::unique_ptr<SeqRegions>& subtree_regions, const RealNumType new_branch_length, const RealNumType new_lh);
    
    /**
        Place a subtree at a mid-branch point
     */
    void placeSubTreeMidBranch(const Index selected_node_index, const Index subtree_index, PhyloNode& subtree, const std::unique_ptr<SeqRegions>& subtree_regions, const RealNumType new_branch_length, const RealNumType new_lh);
    
    /**
        Connect a subtree to a branch
     */
    template<void (Tree::*updateRegionsSubTree)(PhyloNode&, PhyloNode&, PhyloNode&, std::unique_ptr<SeqRegions>&&, const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, RealNumType&)>
    void connectSubTree2Branch(const std::unique_ptr<SeqRegions>& subtree_regions, const std::unique_ptr<SeqRegions>& lower_regions, const Index subtree_index, PhyloNode& subtree, const Index sibling_node_index, PhyloNode& sibling_node, const RealNumType top_distance, const RealNumType down_distance, RealNumType &best_blength, std::unique_ptr<SeqRegions>&& best_child_regions, const std::unique_ptr<SeqRegions>& upper_left_right_regions);
    
    /**
        Connect a subtree to root
     */
    void connectSubTree2Root(const Index subtree_index, PhyloNode& subtree, const std::unique_ptr<SeqRegions>& subtree_regions, const std::unique_ptr<SeqRegions>& lower_regions, const Index sibling_node_index, PhyloNode& sibling_node, const RealNumType best_root_blength, const RealNumType best_length2, std::unique_ptr<SeqRegions>&& best_parent_regions);
    
    /**
        Update next_node_1->partial_lh and new_internal_node->partial_lh after placing a subtree in common cases (e.g., at a mid-branch point, under a node)
     */
    void updateRegionsPlaceSubTree(PhyloNode& subtree, PhyloNode& sibling_node, PhyloNode& internal, std::unique_ptr<SeqRegions>&& best_child_regions, const std::unique_ptr<SeqRegions>& subtree_regions, const std::unique_ptr<SeqRegions>& upper_left_right_regions, const std::unique_ptr<SeqRegions>& lower_regions, RealNumType& best_blength);
    
    /**
        Update next_node_1->partial_lh and new_internal_node->partial_lh after placing a subtree in other cases (i.e., above a node)
     */
    void updateRegionsPlaceSubTreeAbove(PhyloNode& subtree, PhyloNode& sibling_node, PhyloNode& internal, std::unique_ptr<SeqRegions>&& best_child_regions, const std::unique_ptr<SeqRegions>& subtree_regions, const std::unique_ptr<SeqRegions>& upper_left_right_regions, const std::unique_ptr<SeqRegions>& lower_regions, RealNumType& best_blength);
    
    /**
        Handle polytomy when placing a subtree
     */
    void handlePolytomyPlaceSubTree(const Index selected_node_index, PhyloNode& selected_node, const std::unique_ptr<SeqRegions>& subtree_regions, const RealNumType new_branch_length, RealNumType& best_down_lh_diff, Index& best_child_index, RealNumType& best_child_blength_split, std::unique_ptr<SeqRegions>& best_child_regions);
    
    /**
        Update likelihood at mid-branch point
     */
    void updateMidBranchLh(const Index node_index, PhyloNode& node, const std::unique_ptr<SeqRegions>& parent_upper_regions, std::stack<Index> &node_stack, bool &update_blength);
    
    /**
        Compute Upper Left/Right regions at a node, updating the top branch length if neccessary
     */
    std::unique_ptr<SeqRegions> computeUpperLeftRightRegions(const Index node_index, PhyloNode& node, const MiniIndex next_node_mini, const std::unique_ptr<SeqRegions>& parent_upper_regions, std::stack<Index> &node_stack, bool &update_blength);
    
    /**
        Update the PartialLh (seqregions) at a node if the new one is different from the current one
     */
    void updateNewPartialIfDifferent(PhyloNode& node, const MiniIndex next_node_mini, std::unique_ptr<SeqRegions>& upper_left_right_regions, std::stack<Index> &node_stack, const PositionType seq_length);
    
    /**
        Handle cases when the new seqregions is null/empty: (1) update the branch length; or (2) return an error message
     */
    void inline handleNullNewRegions(const Index index, PhyloNode& node, const bool do_update_zeroblength, std::stack<Index> &node_stack, bool &update_blength, const std::string err_msg)
    {
        if (do_update_zeroblength)
        {
            updateZeroBlength(index, node, node_stack);
            update_blength = true;
        }
        else
            outError(err_msg);
    }
    
    /**
        Update partial_lh comming from the parent node
     */
    void updatePartialLhFromParent(const Index index, PhyloNode& node, std::stack<Index> &node_stack, const std::unique_ptr<SeqRegions>& parent_upper_regions, const PositionType seq_length);
    
    /**
        Update partial_lh comming from the children
     */
    void updatePartialLhFromChildren(const Index index, PhyloNode& node, std::stack<Index> &node_stack, const std::unique_ptr<SeqRegions>& parent_upper_regions, const bool is_non_root, const PositionType seq_length);
    
    /**
        Compute the mid-branch region for a node/branch
     */
    inline void computeMidBranchRegions(PhyloNode& node, std::unique_ptr<SeqRegions>& regions_2_update, const SeqRegions &parent_upper_lr_lh)
    {
        std::unique_ptr<SeqRegions>& lower_lh = node.getPartialLh(TOP);
        RealNumType half_branch_length = node.getUpperLength() * 0.5;
        parent_upper_lr_lh.mergeUpperLower(regions_2_update, half_branch_length, *lower_lh, half_branch_length, aln, model, params->threshold_prob);
    }
    
    /**
        Refresh all non-lowerlhs traversing from a parent node
     */
    void refreshNonLowerLhsFromParent(Index& node_index, Index& last_node_index);
    
    /**
        Refresh upper left/right regions
     */
    void refreshUpperLR(const Index node_index, PhyloNode& node, const Index neighbor_index, std::unique_ptr<SeqRegions>& replaced_regions, const SeqRegions& parent_upper_lr_lh);
    
    /**
        Calculate coefficients when merging R with O to estimate a branch length
     */
    void estimateBlength_R_O(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType total_blength, const PositionType end_pos, RealNumType &coefficient, std::vector<RealNumType> &coefficient_vec);
    
    /**
        Calculate coefficients when merging R with ACGT to estimate a branch length
     */
    void estimateBlength_R_ACGT(const SeqRegion& seq1_region, const StateType seq2_state, const RealNumType total_blength, const PositionType end_pos, std::vector<RealNumType> &coefficient_vec);
    
    /**
        Calculate coefficients when merging O with X(i.e., O, R, ACGT) to estimate a branch length
     */
    void estimateBlength_O_X(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType total_blength, const PositionType end_pos, RealNumType &coefficient, std::vector<RealNumType> &coefficient_vec);
    
    /**
        Calculate coefficients when merging ACGT with O to estimate a branch length
     */
    void estimateBlength_ACGT_O(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType total_blength, RealNumType &coefficient, std::vector<RealNumType> &coefficient_vec);
    
    /**
        Calculate coefficients when merging ACGT with R/ACGT to estimate a branch length
     */
    void estimateBlength_ACGT_RACGT(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType total_blength, const PositionType end_pos, std::vector<RealNumType> &coefficient_vec);
    
    /**
        Estimate a branch length from coefficients
     */
    RealNumType estimateBlengthFromCoeffs(RealNumType &coefficient, const std::vector<RealNumType> coefficient_vec);
    
    /**
        Handle branch length changed when improve a subtree
     */
    void handleBlengthChanged(PhyloNode& node, const Index node_index, const RealNumType best_blength);
    
    /**
        Optimize a branch length before seeking an SPR move for a subtree
     */
    void optimizeBlengthBeforeSeekingSPR(PhyloNode& node, RealNumType &best_blength, RealNumType &best_lh, bool &blength_changed, const std::unique_ptr<SeqRegions>& parent_upper_lr_lh, const std::unique_ptr<SeqRegions>& lower_lh);
    
    /**
        Check and apply SPR move
     */
    void checkAndApplySPR(const RealNumType best_lh_diff, const RealNumType best_blength, const RealNumType best_lh, const Index node_index, PhyloNode& node, const Index best_node_index, const Index parent_node_index, const bool is_mid_node, RealNumType& total_improvement, bool& topology_updated);
    
    /**
        Create a new internal phylonode
     */
    inline void createAnInternalNode()
    {
        nodes.emplace_back(InternalNode());
    }
    
    /**
        Create a new leaf phylonode
     */
    inline void createALeafNode(const NumSeqsType new_seq_name_index)
    {
        nodes.emplace_back(LeafNode(new_seq_name_index)); //(PhyloNode(std::move(LeafNode(new_seq_name_index))));
    }
    
    /**
        Get partial_lh at a node by its index
     */
    std::unique_ptr<SeqRegions>& getPartialLhAtNode(const Index index);
    
    /**
        Calculate the likelihood of an NNI neighbor
     */
    bool calculateNNILh(std::stack<Index>& node_stack_aLRT, RealNumType& lh_diff, PhyloNode& current_node, PhyloNode& child_1, PhyloNode& child_2, PhyloNode& sibling, PhyloNode& parent, const Index parent_index, RealNumType& lh_at_root);
    
    /**
        Calculate the likelihood of an NNI neighbor on the branch connecting to root
     */
    bool calculateNNILhRoot(std::stack<Index>& node_stack_aLRT, RealNumType& lh_diff, std::unique_ptr<SeqRegions>& parent_new_lower_lh, const RealNumType& child_2_new_blength, PhyloNode& current_node, PhyloNode& child_1, PhyloNode& child_2, PhyloNode& sibling, PhyloNode& parent, const Index parent_index, RealNumType& lh_at_root);
    
    /**
        Calculate the likelihood of an NNI neighbor on the branch connecting to a non-root node
     */
    bool calculateNNILhNonRoot(std::stack<Index>& node_stack_aLRT, RealNumType& lh_diff, std::unique_ptr<SeqRegions>& parent_new_lower_lh, const RealNumType& child_2_new_blength, PhyloNode& current_node, PhyloNode& child_1, PhyloNode& child_2, PhyloNode& sibling, PhyloNode& parent, const Index parent_index, RealNumType& lh_at_root);
    
    /**
        Replace the current ML Tree by an NNI neighbor on a branch connecting to root
     */
    void replaceMLTreebyNNIRoot(std::stack<Index>& node_stack_aLRT, RealNumType& lh_diff, PhyloNode& current_node, PhyloNode& child_1, PhyloNode& child_2, PhyloNode& sibling, PhyloNode& parent, RealNumType& lh_at_root, const RealNumType child_1_best_blength, const RealNumType child_2_best_blength, const RealNumType sibling_best_blength, const RealNumType parent_best_blength);
    
    /**
        Replace the current ML Tree by an NNI neighbor on a branch connecting to a non-root node
     */
    void replaceMLTreebyNNINonRoot(std::stack<Index>& node_stack_aLRT, RealNumType& lh_diff, PhyloNode& current_node, PhyloNode& child_1, PhyloNode& child_2, PhyloNode& sibling, PhyloNode& parent, RealNumType& lh_at_root, const RealNumType child_1_best_blength, const RealNumType child_2_best_blength, const RealNumType sibling_best_blength, const RealNumType parent_best_blength, const RealNumType new_parent_best_blength);
    
    /**
        Traverse downward to update the upper_left/right_region until the changes is insignificant
     */
    void updateUpperLR(std::stack<Index>& node_stack);
    
    /**
        Calculate aLRT for each internal branches
     */
    void calculate_aRLT();
    
    /**
        Perform a DFS to calculate the Site-lh-contribution
     */
    RealNumType calculateSiteLhs(std::vector<RealNumType>& site_lh_contributions, std::vector<RealNumType>& site_lh_root);
    
    /**
        Calculate aLRT-SH for each internal branches
     */
    void calculate_aRLT_SH(std::vector<RealNumType>& site_lh_contributions, std::vector<RealNumType>& site_lh_root, const RealNumType& LT1);
    
    /**
        Count aLRT-SH for an internal branch
     */
    PositionType count_aRLT_SH_branch(std::vector<RealNumType>& site_lh_contributions, std::vector<RealNumType>& site_lh_root, PhyloNode& node, const RealNumType& LT1);
    
    /**
        Calculate the site-lh differences  between an NNI neighbor on the branch connecting to root and the ML tree
     */
    void calSiteLhDiffRoot(std::vector<RealNumType>& site_lh_diff, std::vector<RealNumType>& site_lh_root_diff, const std::vector<RealNumType>& site_lh_root, std::unique_ptr<SeqRegions>& parent_new_lower_lh, const RealNumType& child_2_new_blength, PhyloNode& current_node, PhyloNode& child_1, PhyloNode& child_2, PhyloNode& sibling, PhyloNode& parent, const Index parent_index);
    
    /**
        Calculate the site-lh differences  between an NNI neighbor on the branch connecting to a non-root node and the ML tree
     */
    void calSiteLhDiffNonRoot(std::vector<RealNumType>& site_lh_diff, std::vector<RealNumType>& site_lh_root_diff, const std::vector<RealNumType>& site_lh_root, std::unique_ptr<SeqRegions>& parent_new_lower_lh, const RealNumType& child_2_new_blength, PhyloNode& current_node, PhyloNode& child_1, PhyloNode& child_2, PhyloNode& sibling, PhyloNode& parent, const Index parent_index);
    
    /**
        Calculate the site-lh differences  between an NNI neighbor and the ML tree
     */
    void calSiteLhDiff(std::vector<RealNumType>& site_lh_diff, std::vector<RealNumType>& site_lh_root_diff, const std::vector<RealNumType>& site_lh_root, PhyloNode& current_node, PhyloNode& child_1, PhyloNode& child_2, PhyloNode& sibling, PhyloNode& parent, const Index parent_index);
    
    /**
        Read the next character from the treefile
     */
    const char readNextChar(std::istream& in, PositionType& in_line, PositionType& in_column, const char& current_ch = 0) const;
    
    /**
        Read string from tree file to create new nodes
     */
    NumSeqsType parseFile(std::istream &infile, char& ch, RealNumType& branch_len, PositionType& in_line, PositionType& in_column, const std::map<std::string, NumSeqsType>& map_seqname_index);
    
    /**
        Collapse leaves with zero-branch-lengths into a vector of less-info-seqs of a leaf
     */
    void collapseAllZeroLeave();
    
    /**
        Collapse one zero-branch-length leaf into its sibling's vector of less-info-seqs
     */
    void collapseOneZeroLeaf(PhyloNode& node, Index& node_index, PhyloNode& neighbor_1, const Index neighbor_1_index, PhyloNode& neighbor_2);
    
    /**
        Update the pesudocount of the model based on the sequence of a leaf
     */
    void updatePesudoCountModel(PhyloNode& node, const Index node_index, const Index parent_index);
    
    /**
        Expand the tree by placing one less-info-seq
     */
    void expandTreeByOneLessInfoSeq(PhyloNode& node, const Index node_index, const Index parent_index);
    
    /**
        Carefully update blength of a node when replacing the ML tree by an NNI neighbor
        Expand the new tree by adding one less-info -seq of the current node (if neccessary) to make sure we compute aLRT for all non-zero internal branches
     */
    void updateBlengthReplaceMLTree(std::stack<Index>& node_stack_aLRT, RealNumType& lh_diff, PhyloNode& node, const Index node_index, const RealNumType best_blength);
    
    /**
        Expand the new tree by adding one less-info -seq of the current node after replacing the ML tree by an NNI neighbor to make sure we compute aLRT for all non-zero internal branches
     */
    void addLessInfoSeqReplacingMLTree(std::stack<Index>& node_stack_aLRT, RealNumType& lh_diff, PhyloNode& node, const Index node_index, const Index parent_index);
    
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
        Vector of likelihood contributions of internal nodes
     */
    std::vector<NodeLh> node_lhs;
    
    /*
        (vector) Index of root in the vector of phylonodes
     */
    NumSeqsType root_vector_index;
    
    /**
        Constructor
     */
    Tree() = default;
    
    /**
        Constructor
    */
    Tree(Params&& params);
    
    /**
        Setup tree parameters/thresholds
     */
    void setup();
    
    /**
        Export tree std::string in Newick format
     */
    const std::string exportTreeString(const bool binary, const NumSeqsType node_vec_index, const bool show_branch_supports);
    
    /**
        Increase the length of a 0-length branch (connecting this node to its parent) to resolve the inconsistency when updating regions in updatePartialLh()
     */
    void updateZeroBlength(const Index index, PhyloNode& node, std::stack<Index> &node_stack);
    
    /**
        Iteratively update partial_lh starting from the nodes in node_stack
        @param node_stack stack of nodes;
     */
    void updatePartialLh(std::stack<Index> &node_stack);
    
    /**
        Seek a position for a sample placement starting at the start_node
     */
    void seekSamplePlacement(const Index start_node_index, const NumSeqsType seq_name_index, const std::unique_ptr<SeqRegions>& sample_regions, Index& selected_node_index, RealNumType &best_lh_diff, bool &is_mid_branch, RealNumType &best_up_lh_diff, RealNumType &best_down_lh_diff, Index& best_child_index);
    
    /**
        Seek a position for placing a subtree/sample starting at the start_node
     */
    void seekSubTreePlacement(Index& best_node_index, RealNumType &best_lh_diff, bool &is_mid_branch, RealNumType &best_up_lh_diff, RealNumType &best_down_lh_diff, Index& best_child_index, const bool short_range_search, const Index child_node_index, RealNumType &removed_blength); //, bool search_subtree_placement = true, SeqRegions* sample_regions = NULL);
    
    /**
        Place a new sample at a mid-branch point
     */
    void placeNewSampleMidBranch(const Index& selected_node_index, std::unique_ptr<SeqRegions>& sample, const NumSeqsType seq_name_index, const RealNumType best_lh_diff);
    
    /**
        Place a new sample as a descendant of a node
     */
    void placeNewSampleAtNode(const Index selected_node_index, std::unique_ptr<SeqRegions>& sample, const NumSeqsType seq_name_index, const RealNumType best_lh_diff, const RealNumType best_up_lh_diff, const RealNumType best_down_lh_diff, const Index best_child_index);
    
    /**
        Apply SPR move
        pruning a subtree then regrafting it to a new position
     */
    void applySPR(const Index subtree_index, PhyloNode& subtree, const Index best_node_index, const bool is_mid_branch, const RealNumType branch_length, const RealNumType best_lh_diff);
    
    /**
        Traverse the intial tree from root to re-calculate all likelihoods regarding the latest/final estimated model parameters
     */
    void refreshAllLhs(bool avoid_using_upper_lr_lhs = false);
    
    /**
        Traverse the tree to reset the outdated/spr_applied flags (which is used to prevent traversing the same part of the tree multiple times) of all nodes
        outdated = true;
        spr_applied = false;
     */
    void resetSPRFlags(const bool reset_outdated);
    
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
    RealNumType estimateBranchLength(const std::unique_ptr<SeqRegions>& parent_regions, const std::unique_ptr<SeqRegions>& child_regions);
    
    /**
        Estimate the length of a branch and check whether the new branch is different from the current one
     */
    RealNumType estimateBranchLengthWithCheck(const std::unique_ptr<SeqRegions>& upper_lr_regions, const std::unique_ptr<SeqRegions>& lower_regions, const RealNumType current_blength);
    
    /**
        Calculate the placement cost of a sample
        @param child_regions: vector of regions of the new sample
     */
    RealNumType calculateSamplePlacementCost(const std::unique_ptr<SeqRegions>& parent_regions, const std::unique_ptr<SeqRegions>& child_regions, const RealNumType blength);
    
    /**
        Calculate the placement cost of a subtree
        @param child_regions: vector of regions of the new sample
     */
    RealNumType calculateSubTreePlacementCost(const std::unique_ptr<SeqRegions>& parent_regions, const std::unique_ptr<SeqRegions>& child_regions, const RealNumType blength);
    
    /**
        Update lower lh of a node
     */
    void updateLowerLh(RealNumType& total_lh, std::unique_ptr<SeqRegions>& new_lower_lh, PhyloNode& node, const std::unique_ptr<SeqRegions>& lower_lh_1, const std::unique_ptr<SeqRegions>& lower_lh_2, const Index neighbor_1_index, PhyloNode& neighbor_1, const Index neighbor_2_index, PhyloNode& neighbor_2, const PositionType& seq_length);
    
    /**
        Update lower lh of a node but avoid using UpperLeft/Right lhs to update zero-blength
        This function is called after reading a tree from an input file, thus, UpperLeft/Right lhs have not yet been computed
     */
    void updateLowerLhAvoidUsingUpperLRLh(RealNumType& total_lh, std::unique_ptr<SeqRegions>& new_lower_lh, PhyloNode& node, const std::unique_ptr<SeqRegions>& lower_lh_1, const std::unique_ptr<SeqRegions>& lower_lh_2, const Index neighbor_1_index, PhyloNode& neighbor_1, const Index neighbor_2_index, PhyloNode& neighbor_2, const PositionType& seq_length);
    
    /**
        compute the likelihood contribution of (the upper branch of) a node
     */
    void computeLhContribution(RealNumType& total_lh, std::unique_ptr<SeqRegions>& new_lower_lh, PhyloNode& node, const std::unique_ptr<SeqRegions>& lower_lh_1, const std::unique_ptr<SeqRegions>& lower_lh_2, const Index neighbor_1_index, PhyloNode& neighbor_1, const Index neighbor_2_index, PhyloNode& neighbor_2, const PositionType& seq_length);
    
    /**
        Calculate the likelihood of the tree
     */
    RealNumType calculateTreeLh();
    
    /**
        Employ Depth First Search to do a task at internal nodes
     */
    template <void(Tree::*task)(RealNumType&, std::unique_ptr<SeqRegions>&, PhyloNode&, const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, const Index, PhyloNode&, const Index, PhyloNode&, const PositionType&)>
    RealNumType performDFS();
    
    /**
        Calculate branch supports
     */
    void calculateBranchSupports();
    
    /**
        Read an input tree from a file
     */
    void readTree(const char* input_treefile);
    
    /**
        Update model parameters from an alignment and a tree
     */
    void updateModelParams();
    
    /**
        Show model parameters
     */
    inline void showModelParams()
    {
        std::cout << model.exportString(aln) << std::endl;
    }
    
    /**
        Employ Depth First Search to do a task at leaves
     */
    template <void(Tree::*task)(PhyloNode&, const Index, const Index)>
    void performDFSAtLeave();
};
