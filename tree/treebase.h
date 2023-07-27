#include "updatingnode.h"
#include "../alignment/alignmentbase.h"
#include "../model/model_dna.h"
#include "../model/model_aa.h"
#include <optional>
#ifdef _OPENMP
    #include <omp.h>
#endif

#pragma once

namespace cmaple
{
    /** The tree structure */
    class TreeBase {
    private:
        /**
            Pointer  to LoadTree method
         */
        typedef void (TreeBase::*LoadTreePtrType)(std::istream&, const bool);
        LoadTreePtrType loadTreePtr;
        
        /**
            Pointer  to changeModel method
         */
        typedef void (TreeBase::*ChangeModelPtrType)(ModelBase*);
        ChangeModelPtrType changeModelPtr;
        
        /**
            Pointer  to changeAln method
         */
        typedef void (TreeBase::*ChangeAlnPtrType)(AlignmentBase*);
        ChangeAlnPtrType changeAlnPtr;
        
        /**
            Pointer  to doInference method
         */
        typedef void (TreeBase::*DoInferencePtrType)(const TreeSearchType, const bool);
        DoInferencePtrType doInferencePtr;
        
        /**
            Pointer  to calculateLh method
         */
        typedef RealNumType (TreeBase::*CalculateLhPtrType)();
        CalculateLhPtrType calculateLhPtr;
        
        /**
            Pointer  to calculateBranchSupport method
         */
        typedef void (TreeBase::*CalculateBranchSupportPtrType)(const int, const int, const double, const bool);
        CalculateBranchSupportPtrType calculateBranchSupportPtr;
        
        /*! Template of loadTree()
         @param[in] tree_stream A stream of the input tree
         @param[in] fixed_blengths TRUE to keep the input branch lengths unchanged (optional)
         */
        template <const cmaple::StateType  num_states>
        void loadTreeTemplate(std::istream& tree_stream, const bool fixed_blengths);
        
        /*! Template of changeAln()
         */
        template <const cmaple::StateType  num_states>
        void changeAlnTemplate(AlignmentBase* aln);
        
        /*! Template of changeModel()
         */
        template <const cmaple::StateType  num_states>
        void changeModelTemplate(ModelBase* model);
        
        /*! Template of doInference()
         */
        template <const cmaple::StateType  num_states>
        void doInferenceTemplate(const TreeSearchType tree_search_type, const bool shallow_tree_search);
        
        /*! Template of calculateLh()
         */
        template <const cmaple::StateType  num_states>
        RealNumType calculateLhTemplate();
        
        /*! Template of calculateBranchSupport()
         */
        template <const cmaple::StateType  num_states>
        void calculateBranchSupportTemplate(const int num_threads, const int num_replicates, const double epsilon, const bool allow_replacing_ML_tree);
        
        /*! Setup function pointers
         */
        void setupFuncPtrs();
        
        /**
         Setup function pointers
         */
        void setupBlengthThresh();
        
        /*! Build an Initial Tree
         */
        template <const cmaple::StateType  num_states>
        void buildInitialTree(const bool from_input_tree);
        
        /*! Optimize the current tree
         */
        template <const cmaple::StateType  num_states>
        void optimizeTree(const bool from_input_tree, const TreeSearchType tree_search_type, const bool shallow_tree_search);
        
        /*! Optimize the tree topology
         */
        template <const cmaple::StateType  num_states>
        void optimizeTreeTopology(bool short_range_search = false);
        
        /*! Optimize the branch lengths of the current tree
         */
        template <const cmaple::StateType  num_states>
        void optimizeBranchLengthsOfTree();
        
        /**
         Traverse the intial tree from root to re-calculate all non-lower likelihoods regarding the latest/final estimated model parameters
         */
        template <const cmaple::StateType  num_states>
        void refreshAllNonLowerLhs();
        
        /**
         Try to improve a subtree rooted at node with SPR moves
         @return total improvement
         */
        template <const cmaple::StateType  num_states>
        cmaple::RealNumType  improveSubTree(const cmaple::Index  index, PhyloNode& node, bool short_range_search);
        
        /**
         Calculate derivative starting from coefficients.
         @return derivative
         */
        cmaple::RealNumType  calculateDerivative(const std::vector<cmaple::RealNumType > &coefficient_vec, const cmaple::RealNumType  delta_t);
        
        /**
         Examine placing a sample at a mid-branch point
         */
        template <const cmaple::StateType  num_states>
        void examineSamplePlacementMidBranch(cmaple::Index & selected_node_index, const std::unique_ptr<SeqRegions>& mid_branch_lh, cmaple::RealNumType  &best_lh_diff, bool& is_mid_branch, cmaple::RealNumType & lh_diff_mid_branch, TraversingNode& current_extended_node, const std::unique_ptr<SeqRegions>& sample_regions);
        
        /**
         Examine placing a sample as a descendant of an existing node
         */
        template <const cmaple::StateType  num_states>
        void examineSamplePlacementAtNode(cmaple::Index & selected_node_index, const std::unique_ptr<SeqRegions>& total_lh, cmaple::RealNumType  &best_lh_diff, bool& is_mid_branch, cmaple::RealNumType & lh_diff_at_node, cmaple::RealNumType & lh_diff_mid_branch, cmaple::RealNumType  &best_up_lh_diff, cmaple::RealNumType  &best_down_lh_diff, cmaple::Index & best_child_index, TraversingNode& current_extended_node, const std::unique_ptr<SeqRegions>& sample_regions);
        
        /**
         Traverse downwards polytomy for more fine-grained placement
         */
        template <const cmaple::StateType  num_states>
        void finetuneSamplePlacementAtNode(const PhyloNode& selected_node, cmaple::RealNumType  &best_down_lh_diff, cmaple::Index & best_child_index, const std::unique_ptr<SeqRegions>& sample_regions);
        
        /**
         Add start nodes for seeking a placement for a subtree
         */
        template <const cmaple::StateType  num_states>
        void addStartingNodes(const cmaple::Index & node_index, PhyloNode& node, const cmaple::Index & other_child_node_index, const cmaple::RealNumType  best_lh_diff, std::stack<std::unique_ptr<UpdatingNode>>& node_stack);
        
        /**
         Examine placing a subtree at a mid-branch point
         */
        template <const cmaple::StateType  num_states>
        bool examineSubtreePlacementMidBranch(cmaple::Index & best_node_index, PhyloNode& current_node, cmaple::RealNumType & best_lh_diff, bool& is_mid_branch, cmaple::RealNumType & lh_diff_at_node, cmaple::RealNumType & lh_diff_mid_branch, cmaple::RealNumType & best_up_lh_diff, cmaple::RealNumType & best_down_lh_diff, std::unique_ptr<UpdatingNode>& updating_node, const std::unique_ptr<SeqRegions>& subtree_regions, const cmaple::RealNumType  threshold_prob, const cmaple::RealNumType  removed_blength, const cmaple::Index  top_node_index, std::unique_ptr<SeqRegions>& bottom_regions);
        
        /**
         Examine placing a subtree as a descendant of an existing node
         */
        template <const cmaple::StateType  num_states>
        bool examineSubTreePlacementAtNode(cmaple::Index & best_node_index, PhyloNode& current_node, cmaple::RealNumType  &best_lh_diff, bool& is_mid_branch, cmaple::RealNumType & lh_diff_at_node, cmaple::RealNumType & lh_diff_mid_branch, cmaple::RealNumType  &best_up_lh_diff, cmaple::RealNumType  &best_down_lh_diff, std::unique_ptr<UpdatingNode>& updating_node, const std::unique_ptr<SeqRegions>& subtree_regions, const cmaple::RealNumType  threshold_prob, const cmaple::RealNumType  removed_blength, const cmaple::Index  top_node_index);
        
        /**
         Add a child node for downwards traversal when seeking a new subtree placement
         */
        template <const cmaple::StateType  num_states>
        void addChildSeekSubtreePlacement(const cmaple::Index  child_1_index, const cmaple::Index  child_2_index, PhyloNode& child_1, PhyloNode& child_2, const cmaple::RealNumType & lh_diff_at_node, const std::unique_ptr<UpdatingNode>& updating_node, std::stack<std::unique_ptr<UpdatingNode>>& node_stack, const cmaple::RealNumType  threshold_prob);
        
        /**
         Add neighbor nodes (parent/sibling) for traversal when seeking a new subtree placement
         */
        template <const cmaple::StateType  num_states>
        bool addNeighborsSeekSubtreePlacement(PhyloNode& current_node, const cmaple::Index  other_child_index, std::unique_ptr<SeqRegions>&& bottom_regions, const cmaple::RealNumType & lh_diff_at_node, const std::unique_ptr<UpdatingNode>& updating_node, std::stack<std::unique_ptr<UpdatingNode>>& node_stack, const cmaple::RealNumType  threshold_prob);
        
        /**
         Check whether we can obtain a higher likelihood with a shorter length for an existing branch
         */
        template <const cmaple::StateType  num_states, cmaple::RealNumType (TreeBase::*calculatePlacementCost)(const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, const cmaple::RealNumType )>
        bool tryShorterBranch(const cmaple::RealNumType  current_blength, std::unique_ptr<SeqRegions>& best_child_regions, const std::unique_ptr<SeqRegions>& sample, const std::unique_ptr<SeqRegions>& upper_left_right_regions, const std::unique_ptr<SeqRegions>& lower_regions, cmaple::RealNumType  &best_split_lh, cmaple::RealNumType  &best_branch_length_split, const cmaple::RealNumType  new_branch_length, const bool try_first_branch);
        
        /**
         Check whether we can obtain a higher likelihood with a shorter length at root
         */
        template <const cmaple::StateType  num_states>
        void tryShorterBranchAtRoot(const std::unique_ptr<SeqRegions>& sample, const std::unique_ptr<SeqRegions>& lower_regions, std::unique_ptr<SeqRegions>& best_parent_regions, cmaple::RealNumType  &best_root_blength, cmaple::RealNumType  &best_parent_lh, const cmaple::RealNumType  fixed_blength);
        
        /**
         Check whether we can obtain a higher likelihood with a shorter length for the new branch at root
         */
        template <const cmaple::StateType  num_states>
        bool tryShorterNewBranchAtRoot(const std::unique_ptr<SeqRegions>& sample, const std::unique_ptr<SeqRegions>& lower_regions, std::unique_ptr<SeqRegions>&best_parent_regions, cmaple::RealNumType  &best_length, cmaple::RealNumType  &best_parent_lh, const cmaple::RealNumType  fixed_blength);
        
        /**
         Check whether we can obtain a higher likelihood with a longer length for the new branch at root
         */
        template <const cmaple::StateType  num_states>
        bool tryLongerNewBranchAtRoot(const std::unique_ptr<SeqRegions>& sample, const std::unique_ptr<SeqRegions>& lower_regions, std::unique_ptr<SeqRegions>& best_parent_regions, cmaple::RealNumType  &best_length, cmaple::RealNumType  &best_parent_lh, const cmaple::RealNumType  fixed_blength);
        
        /**
         Estimate the length for a new branch at root
         */
        template <const cmaple::StateType  num_states>
        void estimateLengthNewBranchAtRoot(const std::unique_ptr<SeqRegions>& sample, const std::unique_ptr<SeqRegions>& lower_regions, std::unique_ptr<SeqRegions>& best_parent_regions, cmaple::RealNumType  &best_length, cmaple::RealNumType  &best_parent_lh, const cmaple::RealNumType  fixed_blength, const cmaple::RealNumType  short_blength_thresh, const bool optional_check);
        
        /**
         Check whether we can obtain a higher likelihood with a shorter length for the new branch
         */
        template <cmaple::RealNumType (TreeBase::*calculatePlacementCost)(const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, const cmaple::RealNumType )>
        bool tryShorterNewBranch(const std::unique_ptr<SeqRegions>& best_child_regions, const std::unique_ptr<SeqRegions>& sample, cmaple::RealNumType  &best_blength, cmaple::RealNumType  &new_branch_lh, const cmaple::RealNumType  short_blength_thresh);
        
        /**
         Check whether we can obtain a higher likelihood with a longer length for the new branch
         */
        template <cmaple::RealNumType (TreeBase::*calculatePlacementCost)(const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, const cmaple::RealNumType )>
        void tryLongerNewBranch(const std::unique_ptr<SeqRegions>& best_child_regions, const std::unique_ptr<SeqRegions>& sample, cmaple::RealNumType  &best_blength, cmaple::RealNumType  &new_branch_lh, const cmaple::RealNumType  long_blength_thresh);
        
        /**
         Estimate the length for a new branch
         */
        template <cmaple::RealNumType (TreeBase::*calculatePlacementCost)(const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, const cmaple::RealNumType )>
        void estimateLengthNewBranch(const cmaple::RealNumType  best_split_lh, const std::unique_ptr<SeqRegions>& best_child_regions, const std::unique_ptr<SeqRegions>& sample, cmaple::RealNumType  &best_blength, const cmaple::RealNumType  long_blength_thresh, const cmaple::RealNumType  short_blength_thresh, const bool optional_check);
        
        /**
         Connect a new sample to a branch
         */
        template <const cmaple::StateType  num_states>
        void connectNewSample2Branch(std::unique_ptr<SeqRegions>& sample, const cmaple::NumSeqsType  seq_name_index, const cmaple::Index  sibling_node_index, PhyloNode& sibling_node, const cmaple::RealNumType  top_distance, const cmaple::RealNumType  down_distance, const cmaple::RealNumType  best_blength, std::unique_ptr<SeqRegions>& best_child_regions, const std::unique_ptr<SeqRegions>& upper_left_right_regions);
        
        /**
         Connect a new sample to root
         */
        template <const cmaple::StateType  num_states>
        void connectNewSample2Root(std::unique_ptr<SeqRegions>& sample, const cmaple::NumSeqsType  seq_name_index, const cmaple::Index  sibling_node_index, PhyloNode& sibling_node, const cmaple::RealNumType  best_root_blength, const cmaple::RealNumType  best_length2, std::unique_ptr<SeqRegions>& best_parent_regions);
        
        /**
         Place a subtree as a descendant of a node
         */
        template <const cmaple::StateType  num_states>
        void placeSubTreeAtNode(const cmaple::Index  selected_node_index, const cmaple::Index  subtree_index, PhyloNode& subtree, const std::unique_ptr<SeqRegions>& subtree_regions, const cmaple::RealNumType  new_branch_length, const cmaple::RealNumType  new_lh);
        
        /**
         Place a subtree at a mid-branch point
         */
        template <const cmaple::StateType  num_states>
        void placeSubTreeMidBranch(const cmaple::Index  selected_node_index, const cmaple::Index  subtree_index, PhyloNode& subtree, const std::unique_ptr<SeqRegions>& subtree_regions, const cmaple::RealNumType  new_branch_length, const cmaple::RealNumType  new_lh);
        
        /**
         Connect a subtree to a branch
         */
        template<const cmaple::StateType  num_states, void (TreeBase::*updateRegionsSubTree)(PhyloNode&, PhyloNode&, PhyloNode&, std::unique_ptr<SeqRegions>&&, const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, cmaple::RealNumType &)>
        void connectSubTree2Branch(const std::unique_ptr<SeqRegions>& subtree_regions, const std::unique_ptr<SeqRegions>& lower_regions, const cmaple::Index  subtree_index, PhyloNode& subtree, const cmaple::Index  sibling_node_index, PhyloNode& sibling_node, const cmaple::RealNumType  top_distance, const cmaple::RealNumType  down_distance, cmaple::RealNumType  &best_blength, std::unique_ptr<SeqRegions>&& best_child_regions, const std::unique_ptr<SeqRegions>& upper_left_right_regions);
        
        /**
         Connect a subtree to root
         */
        template <const cmaple::StateType  num_states>
        void connectSubTree2Root(const cmaple::Index  subtree_index, PhyloNode& subtree, const std::unique_ptr<SeqRegions>& subtree_regions, const std::unique_ptr<SeqRegions>& lower_regions, const cmaple::Index  sibling_node_index, PhyloNode& sibling_node, const cmaple::RealNumType  best_root_blength, const cmaple::RealNumType  best_length2, std::unique_ptr<SeqRegions>&& best_parent_regions);
        
        /**
         Update next_node_1->partial_lh and new_internal_node->partial_lh after placing a subtree in common cases (e.g., at a mid-branch point, under a node)
         */
        template <const cmaple::StateType  num_states>
        void updateRegionsPlaceSubTree(PhyloNode& subtree, PhyloNode& sibling_node, PhyloNode& internal, std::unique_ptr<SeqRegions>&& best_child_regions, const std::unique_ptr<SeqRegions>& subtree_regions, const std::unique_ptr<SeqRegions>& upper_left_right_regions, const std::unique_ptr<SeqRegions>& lower_regions, cmaple::RealNumType & best_blength);
        
        /**
         Update next_node_1->partial_lh and new_internal_node->partial_lh after placing a subtree in other cases (i.e., above a node)
         */
        template <const cmaple::StateType  num_states>
        void updateRegionsPlaceSubTreeAbove(PhyloNode& subtree, PhyloNode& sibling_node, PhyloNode& internal, std::unique_ptr<SeqRegions>&& best_child_regions, const std::unique_ptr<SeqRegions>& subtree_regions, const std::unique_ptr<SeqRegions>& upper_left_right_regions, const std::unique_ptr<SeqRegions>& lower_regions, cmaple::RealNumType & best_blength);
        
        /**
         Handle polytomy when placing a subtree
         */
        template <const cmaple::StateType  num_states>
        void handlePolytomyPlaceSubTree(const cmaple::Index  selected_node_index, PhyloNode& selected_node, const std::unique_ptr<SeqRegions>& subtree_regions, const cmaple::RealNumType  new_branch_length, cmaple::RealNumType & best_down_lh_diff, cmaple::Index & best_child_index, cmaple::RealNumType & best_child_blength_split, std::unique_ptr<SeqRegions>& best_child_regions);
        
        /**
         Update likelihood at mid-branch point
         */
        template <const cmaple::StateType  num_states>
        void updateMidBranchLh(const cmaple::Index  node_index, PhyloNode& node, const std::unique_ptr<SeqRegions>& parent_upper_regions, std::stack<cmaple::Index > &node_stack, bool &update_blength);
        
        /**
         Compute Upper Left/Right regions at a node, updating the top branch length if neccessary
         */
        template <const cmaple::StateType  num_states>
        std::unique_ptr<SeqRegions> computeUpperLeftRightRegions(const cmaple::Index  node_index, PhyloNode& node, const cmaple::MiniIndex  next_node_mini, const std::unique_ptr<SeqRegions>& parent_upper_regions, std::stack<cmaple::Index > &node_stack, bool &update_blength);
        
        /**
         Update the PartialLh (seqregions) at a node if the new one is different from the current one
         */
        bool updateNewPartialIfDifferent(PhyloNode& node, const cmaple::MiniIndex  next_node_mini, std::unique_ptr<SeqRegions>& upper_left_right_regions, std::stack<cmaple::Index > &node_stack, const cmaple::PositionType  seq_length);
        
        /**
         Handle cases when the new seqregions is null/empty: (1) update the branch length; or (2) return an error message
         */
        template <const cmaple::StateType  num_states>
        void inline handleNullNewRegions(const cmaple::Index  index, PhyloNode& node, const bool do_update_zeroblength, std::stack<cmaple::Index > &node_stack, bool &update_blength, const std::string err_msg)
        {
            if (do_update_zeroblength)
            {
                updateZeroBlength<num_states>(index, node, node_stack);
                update_blength = true;
            }
            else
                cmaple::outError(err_msg);
        }
        
        /**
         Update partial_lh comming from the parent node
         */
        template <const cmaple::StateType  num_states>
        void updatePartialLhFromParent(const cmaple::Index  index, PhyloNode& node, std::stack<cmaple::Index > &node_stack, const std::unique_ptr<SeqRegions>& parent_upper_regions, const cmaple::PositionType  seq_length);
        
        /**
         Update partial_lh comming from the children
         */
        template <const cmaple::StateType  num_states>
        void updatePartialLhFromChildren(const cmaple::Index  index, PhyloNode& node, std::stack<cmaple::Index > &node_stack, const std::unique_ptr<SeqRegions>& parent_upper_regions, const bool is_non_root, const cmaple::PositionType  seq_length);
        
        /**
         Compute the mid-branch region for a node/branch
         */
        template <const cmaple::StateType  num_states>
        inline void computeMidBranchRegions(PhyloNode& node, std::unique_ptr<SeqRegions>& regions_2_update, const SeqRegions &parent_upper_lr_lh)
        {
            std::unique_ptr<SeqRegions>& lower_lh = node.getPartialLh(cmaple::TOP);
            cmaple::RealNumType  half_branch_length = node.getUpperLength() * 0.5;
            parent_upper_lr_lh.mergeUpperLower<num_states>(regions_2_update, half_branch_length, *lower_lh, half_branch_length, aln, model, params->threshold_prob);
        }
        
        /**
         Refresh all non-lowerlhs traversing from a parent node
         */
        template <const cmaple::StateType  num_states>
        void refreshNonLowerLhsFromParent(cmaple::Index & node_index, cmaple::Index & last_node_index);
        
        /**
         Refresh upper left/right regions
         */
        template <const cmaple::StateType  num_states>
        void refreshUpperLR(const cmaple::Index  node_index, PhyloNode& node, const cmaple::Index  neighbor_index, std::unique_ptr<SeqRegions>& replaced_regions, const SeqRegions& parent_upper_lr_lh);
        
        /**
         Calculate coefficients when merging R with O to estimate a branch length
         */
        template <const cmaple::StateType  num_states>
        void estimateBlength_R_O(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const cmaple::RealNumType  total_blength, const cmaple::PositionType  end_pos, cmaple::RealNumType  &coefficient, std::vector<cmaple::RealNumType > &coefficient_vec);
        
        /**
         Calculate coefficients when merging R with ACGT to estimate a branch length
         */
        void estimateBlength_R_ACGT(const SeqRegion& seq1_region, const cmaple::StateType  seq2_state, const cmaple::RealNumType  total_blength, const cmaple::PositionType  end_pos, std::vector<cmaple::RealNumType > &coefficient_vec);
        
        /**
         Calculate coefficients when merging O with X(i.e., O, R, ACGT) to estimate a branch length
         */
        template <const cmaple::StateType  num_states>
        void estimateBlength_O_X(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const cmaple::RealNumType  total_blength, const cmaple::PositionType  end_pos, cmaple::RealNumType  &coefficient, std::vector<cmaple::RealNumType > &coefficient_vec);
        
        /**
         Calculate coefficients when merging ACGT with O to estimate a branch length
         */
        template <const cmaple::StateType  num_states>
        void estimateBlength_ACGT_O(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const cmaple::RealNumType  total_blength, cmaple::RealNumType  &coefficient, std::vector<cmaple::RealNumType > &coefficient_vec);
        
        /**
         Calculate coefficients when merging ACGT with R/ACGT to estimate a branch length
         */
        void estimateBlength_ACGT_RACGT(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const cmaple::RealNumType  total_blength, const cmaple::PositionType  end_pos, std::vector<cmaple::RealNumType > &coefficient_vec);
        
        /**
         Estimate a branch length from coefficients
         */
        cmaple::RealNumType  estimateBlengthFromCoeffs(cmaple::RealNumType  &coefficient, const std::vector<cmaple::RealNumType > coefficient_vec);
        
        /**
         Handle branch length changed when improve a subtree
         */
        template <const cmaple::StateType  num_states>
        void handleBlengthChanged(PhyloNode& node, const cmaple::Index  node_index, const cmaple::RealNumType  best_blength);
        
        /**
         Optimize a branch length before seeking an SPR move for a subtree
         */
        template <const cmaple::StateType  num_states>
        void optimizeBlengthBeforeSeekingSPR(PhyloNode& node, cmaple::RealNumType  &best_blength, cmaple::RealNumType  &best_lh, bool &blength_changed, const std::unique_ptr<SeqRegions>& parent_upper_lr_lh, const std::unique_ptr<SeqRegions>& lower_lh);
        
        /**
         Check and apply SPR move
         */
        template <const cmaple::StateType  num_states>
        void checkAndApplySPR(const cmaple::RealNumType  best_lh_diff, const cmaple::RealNumType  best_blength, const cmaple::RealNumType  best_lh, const cmaple::Index  node_index, PhyloNode& node, const cmaple::Index  best_node_index, const cmaple::Index  parent_node_index, const bool is_mid_node, cmaple::RealNumType & total_improvement, bool& topology_updated);
        
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
        inline void createALeafNode(const cmaple::NumSeqsType  new_seq_name_index)
        {
            nodes.emplace_back(LeafNode(new_seq_name_index)); //(PhyloNode(std::move(LeafNode(new_seq_name_index))));
        }
        
        /**
         Get partial_lh at a node by its index
         */
        std::unique_ptr<SeqRegions>& getPartialLhAtNode(const cmaple::Index  index);
        
        /**
         Calculate the likelihood of an NNI neighbor
         */
        template <const cmaple::StateType  num_states>
        bool calculateNNILh(std::stack<cmaple::Index >& node_stack_aLRT, cmaple::RealNumType & lh_diff, PhyloNode& current_node, PhyloNode& child_1, PhyloNode& child_2, PhyloNode& sibling, PhyloNode& parent, const cmaple::Index  parent_index, cmaple::RealNumType & lh_at_root, const bool allow_replacing_ML_tree);
        
        /**
         Calculate the likelihood of an NNI neighbor on the branch connecting to root
         */
        template <const cmaple::StateType  num_states>
        bool calculateNNILhRoot(std::stack<cmaple::Index >& node_stack_aLRT, cmaple::RealNumType & lh_diff, std::unique_ptr<SeqRegions>& parent_new_lower_lh, const cmaple::RealNumType & child_2_new_blength, PhyloNode& current_node, PhyloNode& child_1, PhyloNode& child_2, PhyloNode& sibling, PhyloNode& parent, const cmaple::Index  parent_index, cmaple::RealNumType & lh_at_root, const bool allow_replacing_ML_tree);
        
        /**
         Calculate the likelihood of an NNI neighbor on the branch connecting to a non-root node
         */
        template <const cmaple::StateType  num_states>
        bool calculateNNILhNonRoot(std::stack<cmaple::Index >& node_stack_aLRT, cmaple::RealNumType & lh_diff, std::unique_ptr<SeqRegions>& parent_new_lower_lh, const cmaple::RealNumType & child_2_new_blength, PhyloNode& current_node, PhyloNode& child_1, PhyloNode& child_2, PhyloNode& sibling, PhyloNode& parent, const cmaple::Index  parent_index, cmaple::RealNumType & lh_at_root, const bool allow_replacing_ML_tree);
        
        /**
         Replace the current ML Tree by an NNI neighbor on a branch connecting to root
         */
        template <const cmaple::StateType  num_states>
        void replaceMLTreebyNNIRoot(std::stack<cmaple::Index >& node_stack_aLRT, cmaple::RealNumType & lh_diff, PhyloNode& current_node, PhyloNode& child_1, PhyloNode& child_2, PhyloNode& sibling, PhyloNode& parent, cmaple::RealNumType & lh_at_root, const cmaple::RealNumType  child_1_best_blength, const cmaple::RealNumType  child_2_best_blength, const cmaple::RealNumType  sibling_best_blength, const cmaple::RealNumType  parent_best_blength);
        
        /**
         Replace the current ML Tree by an NNI neighbor on a branch connecting to a non-root node
         */
        template <const cmaple::StateType  num_states>
        void replaceMLTreebyNNINonRoot(std::stack<cmaple::Index >& node_stack_aLRT, cmaple::RealNumType & lh_diff, PhyloNode& current_node, PhyloNode& child_1, PhyloNode& child_2, PhyloNode& sibling, PhyloNode& parent, cmaple::RealNumType & lh_at_root, const cmaple::RealNumType  child_1_best_blength, const cmaple::RealNumType  child_2_best_blength, const cmaple::RealNumType  sibling_best_blength, const cmaple::RealNumType  parent_best_blength, const cmaple::RealNumType  new_parent_best_blength);
        
        /**
         Traverse downward to update the upper_left/right_region until the changes is insignificant
         */
        template <const cmaple::StateType  num_states>
        void updateUpperLR(std::stack<cmaple::Index >& node_stack, std::stack<cmaple::Index >& node_stack_aLRT);
        
        /**
         Add grand-children of a node to node_stack_aLRT to recompute their aLRT after replacing the ML tree with an NNI neighbor
         */
        void recompute_aLRT_GrandChildren(PhyloNode& parent, std::stack<cmaple::Index >& node_stack_aLRT);
        
        /**
         Calculate aLRT for each internal branches
         @param[in] allow_replacing_ML_tree TRUE to allow replacing the ML tree by a higher likelihood tree found when computing branch supports
         */
        template <const cmaple::StateType  num_states>
        void calculate_aRLT(const bool allow_replacing_ML_tree);
        
        /**
         Perform a DFS to calculate the Site-lh-contribution
         */
        template <const cmaple::StateType  num_states>
        cmaple::RealNumType  calculateSiteLhs(std::vector<cmaple::RealNumType >& site_lh_contributions, std::vector<cmaple::RealNumType >& site_lh_root);
        
        /**
         Calculate aLRT-SH for each internal branches
         */
        template <const cmaple::StateType  num_states>
        void calculate_aRLT_SH(std::vector<cmaple::RealNumType >& site_lh_contributions, std::vector<cmaple::RealNumType >& site_lh_root, const cmaple::RealNumType & LT1);
        
        /**
         Count aLRT-SH for an internal branch
         */
        template <const cmaple::StateType  num_states>
        cmaple::PositionType  count_aRLT_SH_branch(std::vector<cmaple::RealNumType >& site_lh_contributions, std::vector<cmaple::RealNumType >& site_lh_root, PhyloNode& node, const cmaple::RealNumType & LT1);
        
        /**
         Calculate the site-lh differences  between an NNI neighbor on the branch connecting to root and the ML tree
         */
        template <const cmaple::StateType  num_states>
        void calSiteLhDiffRoot(std::vector<cmaple::RealNumType >& site_lh_diff, std::vector<cmaple::RealNumType >& site_lh_root_diff, const std::vector<cmaple::RealNumType >& site_lh_root, std::unique_ptr<SeqRegions>& parent_new_lower_lh, const cmaple::RealNumType & child_2_new_blength, PhyloNode& current_node, PhyloNode& child_1, PhyloNode& child_2, PhyloNode& sibling, PhyloNode& parent, const cmaple::Index  parent_index);
        
        /**
         Calculate the site-lh differences  between an NNI neighbor on the branch connecting to a non-root node and the ML tree
         */
        template <const cmaple::StateType  num_states>
        void calSiteLhDiffNonRoot(std::vector<cmaple::RealNumType >& site_lh_diff, std::vector<cmaple::RealNumType >& site_lh_root_diff, const std::vector<cmaple::RealNumType >& site_lh_root, std::unique_ptr<SeqRegions>& parent_new_lower_lh, const cmaple::RealNumType & child_2_new_blength, PhyloNode& current_node, PhyloNode& child_1, PhyloNode& child_2, PhyloNode& sibling, PhyloNode& parent, const cmaple::Index  parent_index);
        
        /**
         Calculate the site-lh differences  between an NNI neighbor and the ML tree
         */
        template <const cmaple::StateType  num_states>
        void calSiteLhDiff(std::vector<cmaple::RealNumType >& site_lh_diff, std::vector<cmaple::RealNumType >& site_lh_root_diff, const std::vector<cmaple::RealNumType >& site_lh_root, PhyloNode& current_node, PhyloNode& child_1, PhyloNode& child_2, PhyloNode& sibling, PhyloNode& parent, const cmaple::Index  parent_index);
        
        /**
         Read the next character from the treefile
         */
        const char readNextChar(std::istream& in, cmaple::PositionType & in_line, cmaple::PositionType & in_column, const char& current_ch = 0) const;
        
        /**
         Read string from tree file to create new nodes
         */
        cmaple::NumSeqsType  parseFile(std::istream &infile, char& ch, cmaple::RealNumType & branch_len, cmaple::PositionType & in_line, cmaple::PositionType & in_column, const std::map<std::string, cmaple::NumSeqsType >& map_seqname_index, bool& missing_blengths);
        
        /**
         Collapse leaves with zero-branch-lengths into a vector of less-info-seqs of a leaf
         */
        void collapseAllZeroLeave();
        
        /**
         Collapse one zero-branch-length leaf into its sibling's vector of less-info-seqs
         */
        void collapseOneZeroLeaf(PhyloNode& node, cmaple::Index & node_index, PhyloNode& neighbor_1, const cmaple::Index  neighbor_1_index, PhyloNode& neighbor_2);
        
        /**
         Update the pesudocount of the model based on the sequence of a leaf
         */
        template <const cmaple::StateType  num_states>
        void updatePesudoCountModel(PhyloNode& node, const cmaple::Index  node_index, const cmaple::Index  parent_index);
        
        /**
         Expand the tree by placing one less-info-seq
         */
        template <const cmaple::StateType  num_states>
        void expandTreeByOneLessInfoSeq(PhyloNode& node, const cmaple::Index  node_index, const cmaple::Index  parent_index);
        
        /**
         Carefully update blength of a node when replacing the ML tree by an NNI neighbor
         Expand the new tree by adding one less-info -seq of the current node (if neccessary) to make sure we compute aLRT for all non-zero internal branches
         */
        template <const cmaple::StateType  num_states>
        void updateBlengthReplaceMLTree(std::stack<cmaple::Index >& node_stack_aLRT, cmaple::RealNumType & lh_diff, PhyloNode& node, const cmaple::Index  node_index, const cmaple::RealNumType  best_blength);
        
        /**
         Expand the new tree by adding one less-info -seq of the current node after replacing the ML tree by an NNI neighbor to make sure we compute aLRT for all non-zero internal branches
         */
        template <const cmaple::StateType  num_states>
        void addLessInfoSeqReplacingMLTree(std::stack<cmaple::Index >& node_stack_aLRT, cmaple::RealNumType & lh_diff, PhyloNode& node, const cmaple::Index  node_index, const cmaple::Index  parent_index);
        
        /**
         Export Node in Newick format
         */
        const std::string exportNodeString(const bool binary, const cmaple::NumSeqsType  node_vec_index, const bool show_branch_supports);
        
        /**
         Read an input tree from a stream
         @Return TRUE if the tree contains any branch without a length
         */
        bool readTree(std::istream& tree_stream);
        
        /**
         Check if the current tree is complete (i.e., containing all sequences from the alignment)
         @Return TRUE if the tree is complete
         */
        bool isComplete();
        
        /**
         Update model according to alignment data
         */
        void updateModelByAln();
        
        /**
         Update model params and partial likelihoods after loading a tree
         */
        template <const cmaple::StateType  num_states>
        void updateModelLhAfterLoading();
        
        /**
         Initialize a mapping between sequence names and their index in the alignment
         */
        std::map<std::string, NumSeqsType> initMapSeqNameIndex();
        
        /**
         * Re-mark the sequences in the alignment, which already existed in the current tree
         * @param[in] old_aln the Old alignment
         */
        void remarkExistingSeqs(AlignmentBase* old_aln);
        
        /**
         * Find and mark a sequence existed in the tree
         * @param[in] seq_name Name of the sequence
         * @param[in] map_name_index a mapping between sequence names and its index in the alignment
         */
        void markAnExistingSeq(const std::string& seq_name, const std::map<std::string, NumSeqsType>& map_name_index);
        
        // NHANLT: Debug aLRT
        // void log_current(std::stack<cmaple::Index>& node_stack_aLRT);
        
    public:
        /*
         TRUE to keep the branch lengths fixed
         */
        bool fixed_blengths = false;
        
        /*
         Branch length thresholds
         */
        cmaple::RealNumType  default_blength, min_blength, max_blength, min_blength_mid, min_blength_sensitivity, half_max_blength, half_min_blength_mid, double_min_blength;
        
        /**
         Program parameters
         */
        //std::optional<cmaple::Params> params;
        std::unique_ptr<cmaple::Params> params = nullptr;
        
        /**
         Alignment
         */
        AlignmentBase* aln;
        
        /**
         Evolutionary model
         */
        ModelBase* model = nullptr;
        
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
        cmaple::NumSeqsType  root_vector_index;
        
        /**
         Constructor
         */
        TreeBase():params(cmaple::make_unique<cmaple::Params>()), aln(nullptr), model(nullptr), fixed_blengths(false) {
            // bug fixed: don't use the first element to store node_lh because node_lh_index is usigned int -> we use 0 for UNINITIALIZED node_lh
            if (node_lhs.size() == 0)
                node_lhs.emplace_back(0);
        };
        
        // TODO: remove or disable
        /**
         Constructor
         */
        // Tree(cmaple::Params && n_params):params(std::move(n_params)) {
        /*TreeBase(cmaple::Params && n_params):params(cmaple::make_unique<cmaple::Params>(std::move(n_params))),aln(new AlignmentBase()), fixed_blengths(false){
            aln->setSeqType(params->seq_type);
        };*/
        
        /**
         Attach alignment and model
         */
        void attachAlnModel(AlignmentBase* aln, ModelBase* model);
        
        /**
         Change the alignment
         */
        void changeAln(AlignmentBase* aln);
        
        /**
         Change the substitution model
         */
        void changeModel(ModelBase* model);
        
        /*! Load an input tree
         @param[in] tree_stream A stream of the input tree
         @param[in] fixed_blengths TRUE to keep the input branch lengths unchanged (optional)
         */
        void loadTree(std::istream& tree_stream, const bool fixed_blengths = false);
        
        /*! Do the inference
         * @param[in] tree_search_type one of the following tree search:
         * - FAST_TREE_SEARCH: no tree search (placement only).
         * - NORMAL_TREE_SEARCH: only consider pruning branches at newly-added nodes when seeking SPR moves.
         * - MORE_ACCURATE_TREE_SEARCH: consider all nodes when seeking SPR moves.
         * @param[in] shallow_tree_search TRUE ton enable a shallow tree search before a deeper tree search
         */
        void doInference(const TreeSearchType tree_search_type, const bool shallow_tree_search);
        
        /*! Calculate the likelihood of the tree
         */
        RealNumType calculateLh();
        
        /**
         Export tree std::string in Newick format
         */
        const std::string exportTreeString(const bool binary, const bool show_branch_supports);
        
        /**
         Export tree std::string in Newick format
         */
        const std::string exportTreeString(const TreeType n_tree_type, const bool show_branch_supports);
        
        /**
         Increase the length of a 0-length branch (connecting this node to its parent) to resolve the inconsistency when updating regions in updatePartialLh()
         */
        template <const cmaple::StateType  num_states>
        void updateZeroBlength(const cmaple::Index  index, PhyloNode& node, std::stack<cmaple::Index > &node_stack);
        
        /**
         Iteratively update partial_lh starting from the nodes in node_stack
         @param node_stack stack of nodes;
         */
        template <const cmaple::StateType  num_states>
        void updatePartialLh(std::stack<cmaple::Index > &node_stack);
        
        /**
         Seek a position for a sample placement starting at the start_node
         */
        template <const cmaple::StateType  num_states>
        void seekSamplePlacement(const cmaple::Index  start_node_index, const cmaple::NumSeqsType  seq_name_index, const std::unique_ptr<SeqRegions>& sample_regions, cmaple::Index & selected_node_index, cmaple::RealNumType  &best_lh_diff, bool &is_mid_branch, cmaple::RealNumType  &best_up_lh_diff, cmaple::RealNumType  &best_down_lh_diff, cmaple::Index & best_child_index);
        
        /**
         Seek a position for placing a subtree/sample starting at the start_node
         */
        template <const cmaple::StateType  num_states>
        void seekSubTreePlacement(cmaple::Index & best_node_index, cmaple::RealNumType  &best_lh_diff, bool &is_mid_branch, cmaple::RealNumType  &best_up_lh_diff, cmaple::RealNumType  &best_down_lh_diff, cmaple::Index & best_child_index, const bool short_range_search, const cmaple::Index  child_node_index, cmaple::RealNumType  &removed_blength); //, bool search_subtree_placement = true, SeqRegions* sample_regions = NULL);
        
        /**
         Place a new sample at a mid-branch point
         */
        template <const cmaple::StateType  num_states>
        void placeNewSampleMidBranch(const cmaple::Index & selected_node_index, std::unique_ptr<SeqRegions>& sample, const cmaple::NumSeqsType  seq_name_index, const cmaple::RealNumType  best_lh_diff);
        
        /**
         Place a new sample as a descendant of a node
         */
        template <const cmaple::StateType  num_states>
        void placeNewSampleAtNode(const cmaple::Index  selected_node_index, std::unique_ptr<SeqRegions>& sample, const cmaple::NumSeqsType  seq_name_index, const cmaple::RealNumType  best_lh_diff, const cmaple::RealNumType  best_up_lh_diff, const cmaple::RealNumType  best_down_lh_diff, const cmaple::Index  best_child_index);
        
        /**
         Apply SPR move
         pruning a subtree then regrafting it to a new position
         */
        template <const cmaple::StateType  num_states>
        void applySPR(const cmaple::Index  subtree_index, PhyloNode& subtree, const cmaple::Index  best_node_index, const bool is_mid_branch, const cmaple::RealNumType  branch_length, const cmaple::RealNumType  best_lh_diff);
        
        /**
         Traverse the intial tree from root to re-calculate all likelihoods regarding the latest/final estimated model parameters
         */
        template <const cmaple::StateType  num_states>
        void refreshAllLhs(bool avoid_using_upper_lr_lhs = false);
        
        /**
         Reset the SPR flags
         @param n_SPR_applied the new value of SPR_applied
         @param update_outdated TRUE to update outdated
         @param n_outdated the new value of outdated
         */
        void resetSPRFlags(const bool n_SPR_applied, const bool update_outdated, const bool n_outdated);
        
        /**
         Try to improve the entire tree with SPR moves
         @return total improvement
         */
        template <const cmaple::StateType  num_states>
        cmaple::RealNumType  improveEntireTree(bool short_range_search);
        
        /**
         Try to optimize branch lengths of the tree
         @return num of improvements
         */
        template <const cmaple::StateType  num_states>
        cmaple::PositionType  optimizeBranchLengths();
        
        /**
         Estimate the length of a branch using the derivative of the likelihood cost function wrt the branch length
         */
        template <const cmaple::StateType  num_states>
        cmaple::RealNumType  estimateBranchLength(const std::unique_ptr<SeqRegions>& parent_regions, const std::unique_ptr<SeqRegions>& child_regions);
        
        /**
         Estimate the length of a branch and check whether the new branch is different from the current one
         */
        template <const cmaple::StateType  num_states>
        cmaple::RealNumType  estimateBranchLengthWithCheck(const std::unique_ptr<SeqRegions>& upper_lr_regions, const std::unique_ptr<SeqRegions>& lower_regions, const cmaple::RealNumType  current_blength);
        
        /**
         Calculate the placement cost of a sample
         @param child_regions: vector of regions of the new sample
         */
        template <const cmaple::StateType  num_states>
        cmaple::RealNumType  calculateSamplePlacementCost(const std::unique_ptr<SeqRegions>& parent_regions, const std::unique_ptr<SeqRegions>& child_regions, const cmaple::RealNumType  blength);
        
        /**
         Calculate the placement cost of a subtree
         @param child_regions: vector of regions of the new sample
         */
        template <const cmaple::StateType  num_states>
        cmaple::RealNumType  calculateSubTreePlacementCost(const std::unique_ptr<SeqRegions>& parent_regions, const std::unique_ptr<SeqRegions>& child_regions, const cmaple::RealNumType  blength);
        
        /**
         Update lower lh of a node
         */
        template <const cmaple::StateType  num_states>
        void updateLowerLh(cmaple::RealNumType & total_lh, std::unique_ptr<SeqRegions>& new_lower_lh, PhyloNode& node, const std::unique_ptr<SeqRegions>& lower_lh_1, const std::unique_ptr<SeqRegions>& lower_lh_2, const cmaple::Index  neighbor_1_index, PhyloNode& neighbor_1, const cmaple::Index  neighbor_2_index, PhyloNode& neighbor_2, const cmaple::PositionType & seq_length);
        
        /**
         Update lower lh of a node but avoid using UpperLeft/Right lhs to update zero-blength
         This function is called after reading a tree from an input file, thus, UpperLeft/Right lhs have not yet been computed
         */
        template <const cmaple::StateType  num_states>
        void updateLowerLhAvoidUsingUpperLRLh(cmaple::RealNumType & total_lh, std::unique_ptr<SeqRegions>& new_lower_lh, PhyloNode& node, const std::unique_ptr<SeqRegions>& lower_lh_1, const std::unique_ptr<SeqRegions>& lower_lh_2, const cmaple::Index  neighbor_1_index, PhyloNode& neighbor_1, const cmaple::Index  neighbor_2_index, PhyloNode& neighbor_2, const cmaple::PositionType & seq_length);
        
        /**
         compute the likelihood contribution of (the upper branch of) a node
         */
        template <const cmaple::StateType  num_states>
        void computeLhContribution(cmaple::RealNumType & total_lh, std::unique_ptr<SeqRegions>& new_lower_lh, PhyloNode& node, const std::unique_ptr<SeqRegions>& lower_lh_1, const std::unique_ptr<SeqRegions>& lower_lh_2, const cmaple::Index  neighbor_1_index, PhyloNode& neighbor_1, const cmaple::Index  neighbor_2_index, PhyloNode& neighbor_2, const cmaple::PositionType & seq_length);
        
        /**
         Employ Depth First Search to do a task at internal nodes
         */
        template <void(TreeBase::*task)(cmaple::RealNumType &, std::unique_ptr<SeqRegions>&, PhyloNode&, const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, const cmaple::Index , PhyloNode&, const cmaple::Index , PhyloNode&, const cmaple::PositionType &)>
        cmaple::RealNumType  performDFS();
        
        /**
         Calculate branch supports
         * @param[in] num_threads number of threads (optional)
         * @param[in] num_replicates a positive number of replicates (optional)
         * @param[in] epsilon a positive epsilon, which is used to avoid rounding effects, when the best and second best NNI trees have nearly identical site log-likelihood values (see Guindon et al., 2010) (optional)
         * @param[in] allow_replacing_ML_tree TRUE to allow replacing the ML tree by a higher likelihood tree found when computing branch supports (optional)
         */
        void calculateBranchSupport(const int num_threads = 1, const int num_replicates = 1000, const double epsilon = 0.1, const bool allow_replacing_ML_tree = true);
        
        /**
         Update model parameters from an alignment and a tree
         */
        template <const cmaple::StateType  num_states>
        void updateModelParams ();
        
        /**
         Export model parameters in string
         */
        inline cmaple::ModelParams exportModelParams ()
        {
            return model->exportModelParams();
        }
        
        /**
         Employ Depth First Search to do a task at leaves
         */
        template <void(TreeBase::*task)(PhyloNode&, const cmaple::Index , const cmaple::Index )>
        void performDFSAtLeave();
    };
}
