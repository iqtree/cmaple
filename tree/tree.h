#include "updatingnode.h"
#include "../alignment/alignment.h"
#include "../model/model.h"
#ifdef _OPENMP
    #include <omp.h>
#endif

#pragma once

namespace cmaple
{
    /** The structure of a phylogenetic tree */
    class Tree {
    public:
        
        /*!
         * Types of tree search
         */
        enum TreeSearchType {
            FAST_TREE_SEARCH,/*!< No tree search (placement only) */
            NORMAL_TREE_SEARCH, /*!< Only consider pruning branches at newly-added nodes when seeking SPR moves */
            MORE_ACCURATE_TREE_SEARCH, /*!< Consider all nodes when seeking SPR moves */
            UNKNOWN_TREE_SEARCH, /*!< Unknown (not specified) */
        };
        
        /*!
         * Types of trees
         */
        enum TreeType {
            BIN_TREE, /*!< Binary tree */
            MUL_TREE, /*!< Multifurcating tree */
            UNKNOWN_TREE, /*!< Unknown tree type */
        };
        
        // ----------------- BEGIN OF PUBLIC APIs ------------------------------------ //
        /*! \brief Constructor from a stream of a (bifurcating or multifurcating) tree (with/without branch lengths in NEWICK format), which may or may not contain all taxa in the alignment. Model parameters (if not fixed) will be estimated according to the input tree and the alignment.
         * @param[in] aln An alignment
         * @param[in] model A substitution model
         * @param[in] tree_stream A stream of an input tree
         * @param[in] fixed_blengths TRUE to keep the input branch lengths unchanged (optional)
         * @param[in] params an instance of program parameters (optional)
         * @throw std::invalid\_argument If any of the following situation occurs.
         * - the sequence type is unsupported (neither DNA (for nucleotide data) nor AA (for protein data))
         * - the alignment is empty
         * - the model is unknown/unsupported
         * - the tree is in an incorrect format
         *
         * @throw std::logic\_error if any of the following situations occur.
         * - taxa in the tree (is specified) is not found in the alignment
         * - unexpected values/behaviors found during the operations
         *
         * @throw std::bad\_alloc if failing to allocate memory to store the tree
         */
        Tree(Alignment* aln, Model* model, std::istream& tree_stream, const bool fixed_blengths = false, std::unique_ptr<cmaple::Params>&& params = nullptr);
        
        /*! \brief Constructor from an optional (bifurcating or multifurcating) tree (with/without branch lengths in NEWICK format), which may or may not contain all taxa in the alignment. If users specify an input tree, model parameters (if not fixed) will be estimated according to that tree and the alignment.
         * @param[in] aln An alignment
         * @param[in] model A substitution model
         * @param[in] tree_filename Name of a tree file (optinal)
         * @param[in] fixed_blengths TRUE to keep the input branch lengths unchanged (optional)
         * @param[in] params an instance of program parameters (optional)
         * @throw std::invalid\_argument If any of the following situation occurs.
         * - the sequence type is unsupported (neither DNA (for nucleotide data) nor AA (for protein data))
         * - the alignment is empty
         * - the model is unknown/unsupported
         * - the tree (is specified) but in an incorrect format
         *
         * @throw ios::failure if the tree file (is specified)  is not found
         * @throw std::logic\_error if any of the following situations occur.
         * - taxa in the tree (is specified) is not found in the alignment
         * - unexpected values/behaviors found during the operations
         *
         * @throw std::bad\_alloc if failing to allocate memory to store the tree
         */
        Tree(Alignment* aln, Model* model, const std::string& tree_filename = "", const bool fixed_blengths = false, std::unique_ptr<cmaple::Params>&& params = nullptr);
        
        /*! Destructor
         */
        ~Tree();
        
        /*! \brief Load a tree from a stream of a (bifurcating or multifurcating) tree (with/without branch lengths) in NEWICK format, which may or may not contain all taxa in the alignment. Model parameters (if not fixed) will be estimated according to the input tree and the alignment.
         * @param[in] tree_stream A stream of an input tree
         * @param[in] fixed_blengths TRUE to keep the input branch lengths unchanged (optional)
         * @throw std::invalid\_argument if tree is empty or in an incorrect format
         * @throw std::logic\_error if any of the following situations occur.
         * - the attached substitution model is unknown/unsupported
         * - any taxa in the tree is not found in the alignment
         * - unexpected values/behaviors found during the operations
         *
         *@throw std::bad\_alloc if failing to allocate memory to store the tree
         */
        void load(std::istream& tree_stream, const bool fixed_blengths = false);
        
        /*! \brief Load a tree from a (bifurcating or multifurcating) tree (with/without branch lengths) in NEWICK format, which may or may not contain all taxa in the alignment. Model parameters (if not fixed) will be estimated according to the input tree and the alignment.
         * @param[in] tree_filename Name of a tree file
         * @param[in] fixed_blengths TRUE to keep the input branch lengths unchanged (optional)
         * @throw std::invalid\_argument if tree is empty or in an incorrect format
         * @throw ios::failure if the tree file is not found
         * @throw std::logic\_error if any of the following situations occur.
         * - the attached substitution model is unknown/unsupported
         * - any taxa in the tree is not found in the alignment
         * - unexpected values/behaviors found during the operations
         *
         * @throw std::bad\_alloc if failing to allocate memory to store the tree
         */
        void load(const std::string& tree_filename, const bool fixed_blengths = false);
        
        /*! \brief Change the alignment
         * @param[in] aln An alignment
         * @throw std::invalid\_argument If the alignment is empty
         * @throw std::logic\_error if any of the following situations occur.
         * - taxa in the current tree is not found in the new alignment
         * - the sequence type of the new alignment is different from the old one
         * - unexpected values/behaviors found during the operations
         */
        void changeAln(Alignment* aln);
        
        /*! \brief Change the substitution model
         * @param[in] model A substitution model
         * @throw std::invalid\_argument if the model is unknown/unsupported
         * @throw std::logic\_error if any of the following situations occur.
         * - the sequence type of the new model is different from the old one
         * - unexpected values/behaviors found during the operations
         */
        void changeModel(Model* model);
        
        /*! Do placement (using stepwise addition) to build an initial tree. Model parameters (if not fixed) will be estimated during the placement process.
         * - If users didn't supply an input tree or supplied an incomplete tree (which doesn't contain all the taxa in the alignment) when initializing the tree (by Tree() constructor), this function will add new taxa (which are not existed in the input tree) from the alignment to the tree.
         * - If users already supplied a complete tree, this function does nothing.
         *
         * @return a string contains all messages redirected from std::cout (for information and debugging purpuses only). To output the tree in NEWICK format, one could call exportNewick() later
         * @throw std::logic\_error if any of the following situations occur.
         * - the attached substitution model is unknown/unsupported
         * - unexpected values/behaviors found during the operations
         */
        std::string doPlacement();
        
        /*! Apply SPR moves to optimize the tree
         * @param[in] tree_search_type A type of tree search
         * @param[in] shallow_tree_search TRUE to enable a shallow tree search before a deeper tree search
         * @return a string contains all messages redirected from std::cout (for information and debugging purpuses only). To output the tree in NEWICK format, one could call exportNewick() later
         * @throw std::logic\_error if any of the following situations occur.
         * - the tree is empty
         * - the attached substitution model is unknown/unsupported
         * - unexpected values/behaviors found during the operations
         */
        std::string applySPR(const TreeSearchType tree_search_type, const bool shallow_tree_search);
        
        /*! Optimize the branch lengths of the current tree
         * @return a string contains all messages redirected from std::cout (for information and debugging purpuses only). To output the tree in NEWICK format, one could call exportNewick() later
         * @throw std::logic\_error if any of the following situations occur.
         * - the tree is empty
         * - the attached substitution model is unknown/unsupported
         * - unexpected values/behaviors found during the operations
         */
        std::string optimizeBranch();
        
        /*! \brief Infer a phylogenetic tree by executing doPlacement(), applySPR(), and optimizeBranch()
         * - If users didn't supply an input tree or supplied an incomplete tree (which doesn't contain all the taxa in the alignment) when initializing the tree (by Tree() constructor), this function:
         *  + performs placement (i.e., adding missing taxa from the alignment to the tree)
         *  + applies a NORMAL tree search (which does SPR moves only on newly-added nodes)
         *  + optimizes all branch lengths
         * - If users already supplied a complete tree, this function:
         *  + by default, does neither placment nor tree search, but it optimizes all branch lengths.
         *  + If users want to keep the branch lengths fixed, they should set fixed_blengths = true when initializing the tree (by Tree() constructor);
         *  + If users want to use the input tree as a starting tree (then performs SPR moves and optimizes branch lengths), they should set tree_search_type = MORE_ACCURATE
         *
         * @param[in] tree_search_type A type of tree search
         * @param[in] shallow_tree_search TRUE to enable a shallow tree search before a deeper tree search
         * @return a string contains all messages redirected from std::cout (for information and debugging purpuses only). To output the tree in NEWICK format, one could call exportNewick() later
         * @throw std::invalid\_argument if tree\_search\_type is unknown
         * @throw std::logic\_error if any of the following situations occur.
         * - the attached substitution model is unknown/unsupported
         * - unexpected values/behaviors found during the operations
         */
        std::string autoProceedMAPLE(const TreeSearchType tree_search_type = NORMAL_TREE_SEARCH, const bool shallow_tree_search = false);
        
        /*! \brief Compute the log likelihood of the current tree, which may or may not contain all taxa in the alignment
         * @return The log likelihood of the current tree
         * @throw std::logic\_error if any of the following situations occur.
         * - the tree is empty
         * - the attached substitution model is unknown/unsupported
         * - unexpected values/behaviors found during the operations
         */
        RealNumType computeLh();
        
        /*! \brief Compute branch supports ([aLRT-SH](https://academic.oup.com/sysbio/article/59/3/307/1702850)) of the current tree, which may or may not contain all taxa in the alignment
         * @param[in] num_threads The number of threads (optional)
         * @param[in] num_replicates A positive number of replicates (optional)
         * @param[in] epsilon A positive epsilon (optional), which is used to avoid rounding effects, when the best and second best NNI trees have nearly identical site log-likelihood values (see [Guindon et al., 2010](https://academic.oup.com/sysbio/article/59/3/307/1702850))
         * @param[in] allow_replacing_ML_tree TRUE to allow replacing the ML tree by a higher likelihood tree found when computing branch supports (optional)
         * @return A string contains all messages redirected from std::cout (for information and debugging purpuses only). To output the branch supports values, one could call exportNewick(BIN, true) later
         * @throw std::invalid\_argument if any of the following situations occur.
         * - num_threads < 0 or num_threads > the number of CPU cores
         * - num_replicates <= 0
         * - epsilon < 0
         *
         * @throw std::logic\_error if any of the following situations occur.
         * - the tree is empty
         * - the attached substitution model is unknown/unsupported
         * - unexpected values/behaviors found during the operations
         */
        std::string computeBranchSupport(const int num_threads = 1, const int num_replicates = 1000, const double epsilon = 0.1, const bool allow_replacing_ML_tree = true);
        
        /*! \brief Export the phylogenetic tree  to a string in NEWICK format.
         * @param[in] tree_type The type of the output tree (optional): BIN_TREE (bifurcating tree), MUL_TREE (multifurcating tree)
         * @param[in] show_branch_supports TRUE to output the branch supports (aLRT-SH values)
         * @return A tree string in NEWICK format
         * @throw std::invalid\_argument if any of the following situations occur.
         * - n\_tree\_type is unknown
         * - show\_branch\_supports = true but branch support values have yet been computed
         */
        std::string exportNewick(const TreeType tree_type = BIN_TREE, const bool show_branch_supports = false);
        
        // ----------------- END OF PUBLIC APIs ------------------------------------ //
        
        /**
         @private
         TRUE to keep the branch lengths fixed
         */
        bool fixed_blengths = false;
        
        /**
         @private
         Branch length thresholds
         */
        cmaple::RealNumType  default_blength, min_blength, max_blength, min_blength_mid, min_blength_sensitivity, half_max_blength, half_min_blength_mid, double_min_blength;
        
        /**
         @private
         Program parameters
         */
        //std::optional<cmaple::Params> params;
        std::unique_ptr<cmaple::Params> params = nullptr;
        
        /**
         @private
         Alignment
         */
        Alignment* aln = nullptr;
        
        /**
         @private
         Evolutionary model
         */
        ModelBase* model = nullptr;
        
        /**
         @private
         cumulative rates
         */
        cmaple::RealNumType *cumulative_rate = nullptr;
        
        /**
         @private
         cumulative bases
         */
        std::vector< std::vector<cmaple::PositionType> > cumulative_base;
        
        /**
         @private
         Vector of phylonodes
         */
        std::vector<PhyloNode> nodes;
        
        /**
         @private
         Vector of likelihood contributions of internal nodes
         */
        std::vector<NodeLh> node_lhs;
        
        /**
         @private
         (vector) Index of root in the vector of phylonodes
         */
        cmaple::NumSeqsType  root_vector_index;
        
        /**
         @private
         A backup of sequence names attached to the current tree, in cases that users re-read the alignment
         NHANLT: @TODO - avoid this redundant vector (sequence names are store in alignment and here)
         */
        std::vector<std::string> seq_names;
        
        /**
         @private
         a vector denote whether a sequence in the alignment is added to the tree or not
         */
        std::vector<bool> sequence_added;
        
        /*!
         * @private
         * Apply some minor changes (collapsing zero-branch leaves into less-info sequences, re-estimating model parameters) to make the processes of outputting then re-inputting a tree result in a consistent tree
         * @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        void makeTreeInOutConsistent();
    
        /**
         * @private
         * Parse type of tree search from a string
         * @param[in] tree_search_type Tree search type in string
         * @return a TreeSearchType
         */
        static TreeSearchType parseTreeSearchType(const std::string& tree_search_type);

        /**
         * @private
         * Get tree search type in string
         * @param[in] tree_search_type a type of tree search
         * @return a  type of tree search from a string
         */
        static std::string getTreeSearchStr(const TreeSearchType tree_search_type);
        
        /**
         * @private
         * Parse tree type from a string
         * @param tree_type_str a tree type in string
         * @return a TreeType
         */
        static TreeType parseTreeType(const std::string& tree_type_str);
        
    private:
        
        /**
            Pointer  to LoadTree method
         */
        typedef void (Tree::*LoadTreePtrType)(std::istream&, const bool);
        LoadTreePtrType loadTreePtr;
        
        /**
            Pointer  to changeModel method
         */
        typedef void (Tree::*ChangeModelPtrType)(Model*);
        ChangeModelPtrType changeModelPtr;
        
        /**
            Pointer  to changeAln method
         */
        typedef void (Tree::*ChangeAlnPtrType)(Alignment*);
        ChangeAlnPtrType changeAlnPtr;
        
        /**
            Pointer  to doInference method
         */
        typedef std::string (Tree::*DoInferencePtrType)(const TreeSearchType, const bool);
        DoInferencePtrType doInferencePtr;
        
        /**
            Pointer  to doPlacement method
         */
        typedef std::string (Tree::*DoPlacementPtrType)();
        DoPlacementPtrType doPlacementPtr;
        
        /**
            Pointer  to applySPR method
         */
        typedef std::string (Tree::*ApplySPRPtrType)(const TreeSearchType, const bool);
        ApplySPRPtrType applySPRPtr;
        
        /**
            Pointer  to optimizeBranch method
         */
        typedef std::string (Tree::*OptimizeBranchPtrType)();
        OptimizeBranchPtrType optimizeBranchPtr;
        
        /**
            Pointer  to computeLh method
         */
        typedef RealNumType (Tree::*computeLhPtrType)();
        computeLhPtrType computeLhPtr;
        
        /**
            Pointer  to computeBranchSupport method
         */
        typedef std::string (Tree::*computeBranchSupportPtrType)(const int, const int, const double, const bool);
        computeBranchSupportPtrType computeBranchSupportPtr;
        
        typedef void (Tree::*MakeTreeInOutConsistentPtrType)();
        MakeTreeInOutConsistentPtrType makeTreeInOutConsistentPtr;
        
        /*! Template of loadTree()
         @param[in] tree_stream A stream of the input tree
         @param[in] fixed_blengths TRUE to keep the input branch lengths unchanged (optional)
         */
        template <const cmaple::StateType  num_states>
        void loadTreeTemplate(std::istream& tree_stream, const bool fixed_blengths);
        
        /*! Template of changeAln()
         */
        template <const cmaple::StateType  num_states>
        void changeAlnTemplate(Alignment* aln);
        
        /*! Template of changeModel()
         */
        template <const cmaple::StateType  num_states>
        void changeModelTemplate(Model* model);
        
        /*! Template of doInference()
         */
        template <const cmaple::StateType  num_states>
        std::string doInferenceTemplate(const TreeSearchType tree_search_type, const bool shallow_tree_search);
        
        /*! Template of doPlacement()
         */
        template <const cmaple::StateType  num_states>
        std::string doPlacementTemplate();
        
        /*! Template of applySPR()
         */
        template <const cmaple::StateType  num_states>
        std::string applySPRTemplate(const TreeSearchType tree_search_type, const bool shallow_tree_search);
        
        /*! Template of optimizeBranch()
         */
        template <const cmaple::StateType  num_states>
        std::string optimizeBranchTemplate();
        
        /*! Template of computeLh()
         */
        template <const cmaple::StateType  num_states>
        RealNumType computeLhTemplate();
        
        /*! Template of computeBranchSupport()
         */
        template <const cmaple::StateType  num_states>
        std::string computeBranchSupportTemplate(const int num_threads, const int num_replicates, const double epsilon, const bool allow_replacing_ML_tree);
        
        /*! Template of makeTreeInOutConsistent()
         */
        template <const cmaple::StateType  num_states>
        void makeTreeInOutConsistentTemplate();
        
        /*! Setup function pointers
         @throw std::invalid\_argument If the sequence type is unsupported (neither DNA (for nucleotide data) nor AA (for protein data))
         */
        void setupFuncPtrs();
        
        /**
         Setup function pointers
         */
        void setupBlengthThresh();
        
        /*! \brief Initialize tree base instance
         * @param[in] aln An alignment
         * @param[in] model A substitution model
         * @param[in] params an instance of program parameters
         * @throw std::invalid\_argument If any of the following situation occurs.
         * - the sequence type is unsupported (neither DNA (for nucleotide data) nor AA (for protein data))
         * - the alignment is empty
         * - model is unknown/unsupported
         *
         * @throw std::logic\_error if the reference genome is empty
         */
        void initTree(Alignment* aln, Model* model, std::unique_ptr<cmaple::Params>&& params);
        
        /**
         @private
         Compute cumulative rate of the ref genome
         @throw std::logic\_error if the reference genome is empty
         */
        void computeCumulativeRate();
        
        /*! Optimize the tree topology
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void optimizeTreeTopology(bool short_range_search = false);
        
        /**
         Traverse the intial tree from root to re-calculate all non-lower likelihoods regarding the latest/final estimated model parameters
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void refreshAllNonLowerLhs();
        
        /**
         Try to improve a subtree rooted at node with SPR moves
         @return total improvement
         @throw std::logic\_error if unexpected values/behaviors found during the operations
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
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void finetuneSamplePlacementAtNode(const PhyloNode& selected_node, cmaple::RealNumType  &best_down_lh_diff, cmaple::Index & best_child_index, const std::unique_ptr<SeqRegions>& sample_regions);
        
        /**
         Add start nodes for seeking a placement for a subtree
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void addStartingNodes(const cmaple::Index & node_index, PhyloNode& node, const cmaple::Index & other_child_node_index, const cmaple::RealNumType  best_lh_diff, std::stack<std::unique_ptr<UpdatingNode>>& node_stack);
        
        /**
         Examine placing a subtree at a mid-branch point
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        bool examineSubtreePlacementMidBranch(cmaple::Index & best_node_index, PhyloNode& current_node, cmaple::RealNumType & best_lh_diff, bool& is_mid_branch, cmaple::RealNumType & lh_diff_at_node, cmaple::RealNumType & lh_diff_mid_branch, cmaple::RealNumType & best_up_lh_diff, cmaple::RealNumType & best_down_lh_diff, std::unique_ptr<UpdatingNode>& updating_node, const std::unique_ptr<SeqRegions>& subtree_regions, const cmaple::RealNumType  threshold_prob, const cmaple::RealNumType  removed_blength, const cmaple::Index  top_node_index, std::unique_ptr<SeqRegions>& bottom_regions);
        
        /**
         Examine placing a subtree as a descendant of an existing node
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        bool examineSubTreePlacementAtNode(cmaple::Index & best_node_index, PhyloNode& current_node, cmaple::RealNumType  &best_lh_diff, bool& is_mid_branch, cmaple::RealNumType & lh_diff_at_node, cmaple::RealNumType & lh_diff_mid_branch, cmaple::RealNumType  &best_up_lh_diff, cmaple::RealNumType  &best_down_lh_diff, std::unique_ptr<UpdatingNode>& updating_node, const std::unique_ptr<SeqRegions>& subtree_regions, const cmaple::RealNumType  threshold_prob, const cmaple::RealNumType  removed_blength, const cmaple::Index  top_node_index);
        
        /**
         Add a child node for downwards traversal when seeking a new subtree placement
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void addChildSeekSubtreePlacement(const cmaple::Index  child_1_index, const cmaple::Index  child_2_index, PhyloNode& child_1, PhyloNode& child_2, const cmaple::RealNumType & lh_diff_at_node, const std::unique_ptr<UpdatingNode>& updating_node, std::stack<std::unique_ptr<UpdatingNode>>& node_stack, const cmaple::RealNumType  threshold_prob);
        
        /**
         Add neighbor nodes (parent/sibling) for traversal when seeking a new subtree placement
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        bool addNeighborsSeekSubtreePlacement(PhyloNode& current_node, const cmaple::Index  other_child_index, std::unique_ptr<SeqRegions>&& bottom_regions, const cmaple::RealNumType & lh_diff_at_node, const std::unique_ptr<UpdatingNode>& updating_node, std::stack<std::unique_ptr<UpdatingNode>>& node_stack, const cmaple::RealNumType  threshold_prob);
        
        /**
         Check whether we can obtain a higher likelihood with a shorter length for an existing branch
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states, cmaple::RealNumType (Tree::*calculatePlacementCost)(const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, const cmaple::RealNumType )>
        bool tryShorterBranch(const cmaple::RealNumType  current_blength, std::unique_ptr<SeqRegions>& best_child_regions, const std::unique_ptr<SeqRegions>& sample, const std::unique_ptr<SeqRegions>& upper_left_right_regions, const std::unique_ptr<SeqRegions>& lower_regions, cmaple::RealNumType  &best_split_lh, cmaple::RealNumType  &best_branch_length_split, const cmaple::RealNumType  new_branch_length, const bool try_first_branch);
        
        /**
         Check whether we can obtain a higher likelihood with a shorter length at root
         @throw std::logic\_error if unexpected values/behaviors found during the operations
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
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        bool tryLongerNewBranchAtRoot(const std::unique_ptr<SeqRegions>& sample, const std::unique_ptr<SeqRegions>& lower_regions, std::unique_ptr<SeqRegions>& best_parent_regions, cmaple::RealNumType  &best_length, cmaple::RealNumType  &best_parent_lh, const cmaple::RealNumType  fixed_blength);
        
        /**
         Estimate the length for a new branch at root
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void estimateLengthNewBranchAtRoot(const std::unique_ptr<SeqRegions>& sample, const std::unique_ptr<SeqRegions>& lower_regions, std::unique_ptr<SeqRegions>& best_parent_regions, cmaple::RealNumType  &best_length, cmaple::RealNumType  &best_parent_lh, const cmaple::RealNumType  fixed_blength, const cmaple::RealNumType  short_blength_thresh, const bool optional_check);
        
        /**
         Check whether we can obtain a higher likelihood with a shorter length for the new branch
         */
        template <cmaple::RealNumType (Tree::*calculatePlacementCost)(const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, const cmaple::RealNumType )>
        bool tryShorterNewBranch(const std::unique_ptr<SeqRegions>& best_child_regions, const std::unique_ptr<SeqRegions>& sample, cmaple::RealNumType  &best_blength, cmaple::RealNumType  &new_branch_lh, const cmaple::RealNumType  short_blength_thresh);
        
        /**
         Check whether we can obtain a higher likelihood with a longer length for the new branch
         */
        template <cmaple::RealNumType (Tree::*calculatePlacementCost)(const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, const cmaple::RealNumType )>
        void tryLongerNewBranch(const std::unique_ptr<SeqRegions>& best_child_regions, const std::unique_ptr<SeqRegions>& sample, cmaple::RealNumType  &best_blength, cmaple::RealNumType  &new_branch_lh, const cmaple::RealNumType  long_blength_thresh);
        
        /**
         Estimate the length for a new branch
         */
        template <cmaple::RealNumType (Tree::*calculatePlacementCost)(const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, const cmaple::RealNumType )>
        void estimateLengthNewBranch(const cmaple::RealNumType  best_split_lh, const std::unique_ptr<SeqRegions>& best_child_regions, const std::unique_ptr<SeqRegions>& sample, cmaple::RealNumType  &best_blength, const cmaple::RealNumType  long_blength_thresh, const cmaple::RealNumType  short_blength_thresh, const bool optional_check);
        
        /**
         Connect a new sample to a branch
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void connectNewSample2Branch(std::unique_ptr<SeqRegions>& sample, const cmaple::NumSeqsType  seq_name_index, const cmaple::Index  sibling_node_index, PhyloNode& sibling_node, const cmaple::RealNumType  top_distance, const cmaple::RealNumType  down_distance, const cmaple::RealNumType  best_blength, std::unique_ptr<SeqRegions>& best_child_regions, const std::unique_ptr<SeqRegions>& upper_left_right_regions);
        
        /**
         Connect a new sample to root
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void connectNewSample2Root(std::unique_ptr<SeqRegions>& sample, const cmaple::NumSeqsType  seq_name_index, const cmaple::Index  sibling_node_index, PhyloNode& sibling_node, const cmaple::RealNumType  best_root_blength, const cmaple::RealNumType  best_length2, std::unique_ptr<SeqRegions>& best_parent_regions);
        
        /**
         Place a subtree as a descendant of a node
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void placeSubTreeAtNode(const cmaple::Index  selected_node_index, const cmaple::Index  subtree_index, PhyloNode& subtree, const std::unique_ptr<SeqRegions>& subtree_regions, const cmaple::RealNumType  new_branch_length, const cmaple::RealNumType  new_lh);
        
        /**
         Place a subtree at a mid-branch point
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void placeSubTreeMidBranch(const cmaple::Index  selected_node_index, const cmaple::Index  subtree_index, PhyloNode& subtree, const std::unique_ptr<SeqRegions>& subtree_regions, const cmaple::RealNumType  new_branch_length, const cmaple::RealNumType  new_lh);
        
        /**
         Connect a subtree to a branch
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template<const cmaple::StateType  num_states, void (Tree::*updateRegionsSubTree)(PhyloNode&, PhyloNode&, PhyloNode&, std::unique_ptr<SeqRegions>&&, const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, cmaple::RealNumType &)>
        void connectSubTree2Branch(const std::unique_ptr<SeqRegions>& subtree_regions, const std::unique_ptr<SeqRegions>& lower_regions, const cmaple::Index  subtree_index, PhyloNode& subtree, const cmaple::Index  sibling_node_index, PhyloNode& sibling_node, const cmaple::RealNumType  top_distance, const cmaple::RealNumType  down_distance, cmaple::RealNumType  &best_blength, std::unique_ptr<SeqRegions>&& best_child_regions, const std::unique_ptr<SeqRegions>& upper_left_right_regions);
        
        /**
         Connect a subtree to root
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void connectSubTree2Root(const cmaple::Index  subtree_index, PhyloNode& subtree, const std::unique_ptr<SeqRegions>& subtree_regions, const std::unique_ptr<SeqRegions>& lower_regions, const cmaple::Index  sibling_node_index, PhyloNode& sibling_node, const cmaple::RealNumType  best_root_blength, const cmaple::RealNumType  best_length2, std::unique_ptr<SeqRegions>&& best_parent_regions);
        
        /**
         Update next_node_1->partial_lh and new_internal_node->partial_lh after placing a subtree in common cases (e.g., at a mid-branch point, under a node)
         
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void updateRegionsPlaceSubTree(PhyloNode& subtree, PhyloNode& sibling_node, PhyloNode& internal, std::unique_ptr<SeqRegions>&& best_child_regions, const std::unique_ptr<SeqRegions>& subtree_regions, const std::unique_ptr<SeqRegions>& upper_left_right_regions, const std::unique_ptr<SeqRegions>& lower_regions, cmaple::RealNumType & best_blength);
        
        /**
         Update next_node_1->partial_lh and new_internal_node->partial_lh after placing a subtree in other cases (i.e., above a node)
         
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void updateRegionsPlaceSubTreeAbove(PhyloNode& subtree, PhyloNode& sibling_node, PhyloNode& internal, std::unique_ptr<SeqRegions>&& best_child_regions, const std::unique_ptr<SeqRegions>& subtree_regions, const std::unique_ptr<SeqRegions>& upper_left_right_regions, const std::unique_ptr<SeqRegions>& lower_regions, cmaple::RealNumType & best_blength);
        
        /**
         Handle polytomy when placing a subtree
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void handlePolytomyPlaceSubTree(const cmaple::Index  selected_node_index, PhyloNode& selected_node, const std::unique_ptr<SeqRegions>& subtree_regions, const cmaple::RealNumType  new_branch_length, cmaple::RealNumType & best_down_lh_diff, cmaple::Index & best_child_index, cmaple::RealNumType & best_child_blength_split, std::unique_ptr<SeqRegions>& best_child_regions);
        
        /**
         Update likelihood at mid-branch point
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void updateMidBranchLh(const cmaple::Index  node_index, PhyloNode& node, const std::unique_ptr<SeqRegions>& parent_upper_regions, std::stack<cmaple::Index > &node_stack, bool &update_blength);
        
        /**
         Compute Upper Left/Right regions at a node, updating the top branch length if neccessary
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        std::unique_ptr<SeqRegions> computeUpperLeftRightRegions(const cmaple::Index  node_index, PhyloNode& node, const cmaple::MiniIndex  next_node_mini, const std::unique_ptr<SeqRegions>& parent_upper_regions, std::stack<cmaple::Index > &node_stack, bool &update_blength);
        
        /**
         Update the PartialLh (seqregions) at a node if the new one is different from the current one
         */
        bool updateNewPartialIfDifferent(PhyloNode& node, const cmaple::MiniIndex  next_node_mini, std::unique_ptr<SeqRegions>& upper_left_right_regions, std::stack<cmaple::Index > &node_stack, const cmaple::PositionType  seq_length);
        
        /**
         Handle cases when the new seqregions is null/empty: (1) update the branch length; or (2) return an error message
         @throw std::logic\_error if unexpected values/behaviors found during the operations
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
                throw std::logic_error(err_msg);
        }
        
        /**
         Update partial_lh comming from the parent node
         
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void updatePartialLhFromParent(const cmaple::Index  index, PhyloNode& node, std::stack<cmaple::Index > &node_stack, const std::unique_ptr<SeqRegions>& parent_upper_regions, const cmaple::PositionType  seq_length);
        
        /**
         Update partial_lh comming from the children
         
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void updatePartialLhFromChildren(const cmaple::Index  index, PhyloNode& node, std::stack<cmaple::Index > &node_stack, const std::unique_ptr<SeqRegions>& parent_upper_regions, const bool is_non_root, const cmaple::PositionType  seq_length);
        
        /**
         Compute the mid-branch region for a node/branch
         @throw std::logic\_error if unexpected values/behaviors found during the operations
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
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void refreshNonLowerLhsFromParent(cmaple::Index & node_index, cmaple::Index & last_node_index);
        
        /**
         Refresh upper left/right regions
         @throw std::logic\_error if unexpected values/behaviors found during the operations
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
         @throw std::logic\_error if unexpected values/behaviors found during the operations
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
         @throw std::logic\_error if unexpected values/behaviors found during the operations
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
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        bool calculateNNILh(std::stack<cmaple::Index >& node_stack_aLRT, cmaple::RealNumType & lh_diff, PhyloNode& current_node, PhyloNode& child_1, PhyloNode& child_2, PhyloNode& sibling, PhyloNode& parent, const cmaple::Index  parent_index, cmaple::RealNumType & lh_at_root, const bool allow_replacing_ML_tree);
        
        /**
         Calculate the likelihood of an NNI neighbor on the branch connecting to root
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        bool calculateNNILhRoot(std::stack<cmaple::Index >& node_stack_aLRT, cmaple::RealNumType & lh_diff, std::unique_ptr<SeqRegions>& parent_new_lower_lh, const cmaple::RealNumType & child_2_new_blength, PhyloNode& current_node, PhyloNode& child_1, PhyloNode& child_2, PhyloNode& sibling, PhyloNode& parent, const cmaple::Index  parent_index, cmaple::RealNumType & lh_at_root, const bool allow_replacing_ML_tree);
        
        /**
         Calculate the likelihood of an NNI neighbor on the branch connecting to a non-root node
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        bool calculateNNILhNonRoot(std::stack<cmaple::Index >& node_stack_aLRT, cmaple::RealNumType & lh_diff, std::unique_ptr<SeqRegions>& parent_new_lower_lh, const cmaple::RealNumType & child_2_new_blength, PhyloNode& current_node, PhyloNode& child_1, PhyloNode& child_2, PhyloNode& sibling, PhyloNode& parent, const cmaple::Index  parent_index, cmaple::RealNumType & lh_at_root, const bool allow_replacing_ML_tree);
        
        /**
         Replace the current ML Tree by an NNI neighbor on a branch connecting to root
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void replaceMLTreebyNNIRoot(std::stack<cmaple::Index >& node_stack_aLRT, cmaple::RealNumType & lh_diff, PhyloNode& current_node, PhyloNode& child_1, PhyloNode& child_2, PhyloNode& sibling, PhyloNode& parent, cmaple::RealNumType & lh_at_root, const cmaple::RealNumType  child_1_best_blength, const cmaple::RealNumType  child_2_best_blength, const cmaple::RealNumType  sibling_best_blength, const cmaple::RealNumType  parent_best_blength);
        
        /**
         Replace the current ML Tree by an NNI neighbor on a branch connecting to a non-root node
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void replaceMLTreebyNNINonRoot(std::stack<cmaple::Index >& node_stack_aLRT, cmaple::RealNumType & lh_diff, PhyloNode& current_node, PhyloNode& child_1, PhyloNode& child_2, PhyloNode& sibling, PhyloNode& parent, cmaple::RealNumType & lh_at_root, const cmaple::RealNumType  child_1_best_blength, const cmaple::RealNumType  child_2_best_blength, const cmaple::RealNumType  sibling_best_blength, const cmaple::RealNumType  parent_best_blength, const cmaple::RealNumType  new_parent_best_blength);
        
        /**
         Traverse downward to update the upper_left/right_region until the changes is insignificant
         
         @throw std::logic\_error if unexpected values/behaviors found during the operations
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
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void calculate_aRLT(const bool allow_replacing_ML_tree);
        
        /**
         Perform a DFS to calculate the Site-lh-contribution
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        cmaple::RealNumType  calculateSiteLhs(std::vector<cmaple::RealNumType >& site_lh_contributions, std::vector<cmaple::RealNumType >& site_lh_root);
        
        /**
         Calculate aLRT-SH for each internal branches
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void calculate_aRLT_SH(std::vector<cmaple::RealNumType >& site_lh_contributions, std::vector<cmaple::RealNumType >& site_lh_root, const cmaple::RealNumType & LT1);
        
        /**
         Count aLRT-SH for an internal branch
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        cmaple::PositionType  count_aRLT_SH_branch(std::vector<cmaple::RealNumType >& site_lh_contributions, std::vector<cmaple::RealNumType >& site_lh_root, PhyloNode& node, const cmaple::RealNumType & LT1);
        
        /**
         Calculate the site-lh differences  between an NNI neighbor on the branch connecting to root and the ML tree
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void calSiteLhDiffRoot(std::vector<cmaple::RealNumType >& site_lh_diff, std::vector<cmaple::RealNumType >& site_lh_root_diff, const std::vector<cmaple::RealNumType >& site_lh_root, std::unique_ptr<SeqRegions>& parent_new_lower_lh, const cmaple::RealNumType & child_2_new_blength, PhyloNode& current_node, PhyloNode& child_1, PhyloNode& child_2, PhyloNode& sibling, PhyloNode& parent, const cmaple::Index  parent_index);
        
        /**
         Calculate the site-lh differences  between an NNI neighbor on the branch connecting to a non-root node and the ML tree
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void calSiteLhDiffNonRoot(std::vector<cmaple::RealNumType >& site_lh_diff, std::vector<cmaple::RealNumType >& site_lh_root_diff, const std::vector<cmaple::RealNumType >& site_lh_root, std::unique_ptr<SeqRegions>& parent_new_lower_lh, const cmaple::RealNumType & child_2_new_blength, PhyloNode& current_node, PhyloNode& child_1, PhyloNode& child_2, PhyloNode& sibling, PhyloNode& parent, const cmaple::Index  parent_index);
        
        /**
         Calculate the site-lh differences  between an NNI neighbor and the ML tree
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void calSiteLhDiff(std::vector<cmaple::RealNumType >& site_lh_diff, std::vector<cmaple::RealNumType >& site_lh_root_diff, const std::vector<cmaple::RealNumType >& site_lh_root, PhyloNode& current_node, PhyloNode& child_1, PhyloNode& child_2, PhyloNode& sibling, PhyloNode& parent, const cmaple::Index  parent_index);
        
        /**
         Read the next character from the treefile
         */
        const char readNextChar(std::istream& in, cmaple::PositionType & in_line, cmaple::PositionType & in_column, const char& current_ch = 0) const;
        
        /**
         Read string from tree file to create new nodes
         @throw std::logic\_error if any of the following situations occur.
         - taxa in the tree is not found in the  alignment
         - unexpected values/behaviors found during the operations
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
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void expandTreeByOneLessInfoSeq(PhyloNode& node, const cmaple::Index  node_index, const cmaple::Index  parent_index);
        
        /**
         Carefully update blength of a node when replacing the ML tree by an NNI neighbor
         Expand the new tree by adding one less-info -seq of the current node (if neccessary) to make sure we compute aLRT for all non-zero internal branches
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void updateBlengthReplaceMLTree(std::stack<cmaple::Index >& node_stack_aLRT, cmaple::RealNumType & lh_diff, PhyloNode& node, const cmaple::Index  node_index, const cmaple::RealNumType  best_blength);
        
        /**
         Expand the new tree by adding one less-info -seq of the current node after replacing the ML tree by an NNI neighbor to make sure we compute aLRT for all non-zero internal branches
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void addLessInfoSeqReplacingMLTree(std::stack<cmaple::Index >& node_stack_aLRT, cmaple::RealNumType & lh_diff, PhyloNode& node, const cmaple::Index  node_index, const cmaple::Index  parent_index);
        
        /**
         Export Node in Newick format
         @throw std::invalid\_argument if show\_branch\_supports = true but branch support values have yet been computed
         */
        std::string exportNodeString(const bool binary, const cmaple::NumSeqsType  node_vec_index, const bool show_branch_supports);
        
        /**
         Read an input tree from a stream
         @Return TRUE if the tree contains any branch without a length
         @throw std::invalid\_argument if the tree in an incorrect format
         @throw std::logic\_error if any of the following situations occur.
         - any taxa in the tree is not found in the alignment
         - unexpected values/behaviors found during the operations
         
         @throw std::bad\_alloc if failing to allocate memory to store the tree
         */
        bool readTree(std::istream& tree_stream);
        
        /**
         Check if the current tree is complete (i.e., containing all sequences from the alignment)
         @Return TRUE if the tree is complete
         */
        bool isComplete();
        
        /**
         Update model according to alignment data
         @throw std::logic\_error if the reference genome is empty
         */
        void updateModelByAln();
        
        /**
         Update model params and partial likelihoods after loading a tree
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void updateModelLhAfterLoading();
        
        /**
         Initialize a mapping between sequence names and their index in the alignment
         */
        std::map<std::string, NumSeqsType> initMapSeqNameIndex();
        
        /**
         * Re-mark the sequences in the alignment, which already existed in the current tree
         * @throw std::logic\_error if any taxa in the current tree is not found in the new alignment
         */
        void remarkExistingSeqs();
        
        /**
         * Find and mark a sequence existed in the tree
         * @param[in] seq_name Name of the sequence
         * @param[in] map_name_index a mapping between sequence names and its index in the alignment
         * @return the new index of the corresponding sequence in the new alignment
         * @throw std::logic\_error if the taxon named seq\_name is not found the alignment
         */
        NumSeqsType markAnExistingSeq(const std::string& seq_name, const std::map<std::string, NumSeqsType>& map_name_index);
        
        /**
         * Mark all sequences (in the alignment) as not yet added to the current tree
         */
        void resetSeqAdded();
        
        /**
         @private
         Attach alignment and model
         @throw std::invalid\_argument If the sequence type is unsupported (neither DNA (for nucleotide data) nor AA (for protein data))
         @throw std::logic\_error if the reference genome is empty
         */
        void attachAlnModel(Alignment* aln, ModelBase* model);
        
        /**
         @private
         Export tree std::string in Newick format
         @throw std::invalid\_argument if show\_branch\_supports = true but branch support values have yet been computed
         */
        std::string exportNewick(const bool binary, const bool show_branch_supports);
        
        /**
         @private
         Increase the length of a 0-length branch (connecting this node to its parent) to resolve the inconsistency when updating regions in updatePartialLh()
         */
        template <const cmaple::StateType  num_states>
        void updateZeroBlength(const cmaple::Index  index, PhyloNode& node, std::stack<cmaple::Index > &node_stack);
        
        /**
         @private
         Iteratively update partial_lh starting from the nodes in node_stack
         
         @param node_stack stack of nodes;
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void updatePartialLh(std::stack<cmaple::Index > &node_stack);
        
        /**
         @private
         Seek a position for a sample placement starting at the start_node
         
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void seekSamplePlacement(const cmaple::Index  start_node_index, const cmaple::NumSeqsType  seq_name_index, const std::unique_ptr<SeqRegions>& sample_regions, cmaple::Index & selected_node_index, cmaple::RealNumType  &best_lh_diff, bool &is_mid_branch, cmaple::RealNumType  &best_up_lh_diff, cmaple::RealNumType  &best_down_lh_diff, cmaple::Index & best_child_index);
        
        /**
         @private
         Seek a position for placing a subtree/sample starting at the start_node
         
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void seekSubTreePlacement(cmaple::Index & best_node_index, cmaple::RealNumType  &best_lh_diff, bool &is_mid_branch, cmaple::RealNumType  &best_up_lh_diff, cmaple::RealNumType  &best_down_lh_diff, cmaple::Index & best_child_index, const bool short_range_search, const cmaple::Index  child_node_index, cmaple::RealNumType  &removed_blength); //, bool search_subtree_placement = true, SeqRegions* sample_regions = NULL);
        
        /**
         @private
         Place a new sample at a mid-branch point
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void placeNewSampleMidBranch(const cmaple::Index & selected_node_index, std::unique_ptr<SeqRegions>& sample, const cmaple::NumSeqsType  seq_name_index, const cmaple::RealNumType  best_lh_diff);
        
        /**
         @private
         Place a new sample as a descendant of a node
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void placeNewSampleAtNode(const cmaple::Index  selected_node_index, std::unique_ptr<SeqRegions>& sample, const cmaple::NumSeqsType  seq_name_index, const cmaple::RealNumType  best_lh_diff, const cmaple::RealNumType  best_up_lh_diff, const cmaple::RealNumType  best_down_lh_diff, const cmaple::Index  best_child_index);
        
        /**
         @private
         Apply a single SPR move
         pruning a subtree then regrafting it to a new position
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void applyOneSPR(const cmaple::Index  subtree_index, PhyloNode& subtree, const cmaple::Index  best_node_index, const bool is_mid_branch, const cmaple::RealNumType  branch_length, const cmaple::RealNumType  best_lh_diff);
        
        /**
         @private
         Traverse the intial tree from root to re-calculate all likelihoods regarding the latest/final estimated model parameters
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void refreshAllLhs(bool avoid_using_upper_lr_lhs = false);
        
        /**
         @private
         Reset the SPR flags
         @param n_SPR_applied the new value of SPR_applied
         @param update_outdated TRUE to update outdated
         @param n_outdated the new value of outdated
         */
        void resetSPRFlags(const bool n_SPR_applied, const bool update_outdated, const bool n_outdated);
        
        /**
         @private
         Try to improve the entire tree with SPR moves
         @return total improvement
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        cmaple::RealNumType  improveEntireTree(bool short_range_search);
        
        /**
         @private
         Try to optimize branch lengths of the tree by one round of tree traversal 
         @return num of improvements
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        cmaple::PositionType  optimizeBranchIter();
        
        /**
         @private
         Estimate the length of a branch using the derivative of the likelihood cost function wrt the branch length
         */
        template <const cmaple::StateType  num_states>
        cmaple::RealNumType  estimateBranchLength(const std::unique_ptr<SeqRegions>& parent_regions, const std::unique_ptr<SeqRegions>& child_regions);
        
        /**
         @private
         Estimate the length of a branch and check whether the new branch is different from the current one
         */
        template <const cmaple::StateType  num_states>
        cmaple::RealNumType  estimateBranchLengthWithCheck(const std::unique_ptr<SeqRegions>& upper_lr_regions, const std::unique_ptr<SeqRegions>& lower_regions, const cmaple::RealNumType  current_blength);
        
        /**
         @private
         Calculate the placement cost of a sample
         @param child_regions: vector of regions of the new sample
         */
        template <const cmaple::StateType  num_states>
        cmaple::RealNumType  calculateSamplePlacementCost(const std::unique_ptr<SeqRegions>& parent_regions, const std::unique_ptr<SeqRegions>& child_regions, const cmaple::RealNumType  blength);
        
        /**
         @private
         Calculate the placement cost of a subtree
         @param child_regions: vector of regions of the new sample
         */
        template <const cmaple::StateType  num_states>
        cmaple::RealNumType  calculateSubTreePlacementCost(const std::unique_ptr<SeqRegions>& parent_regions, const std::unique_ptr<SeqRegions>& child_regions, const cmaple::RealNumType  blength);
        
        /**
         @private
         Update lower lh of a node
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void updateLowerLh(cmaple::RealNumType & total_lh, std::unique_ptr<SeqRegions>& new_lower_lh, PhyloNode& node, const std::unique_ptr<SeqRegions>& lower_lh_1, const std::unique_ptr<SeqRegions>& lower_lh_2, const cmaple::Index  neighbor_1_index, PhyloNode& neighbor_1, const cmaple::Index  neighbor_2_index, PhyloNode& neighbor_2, const cmaple::PositionType & seq_length);
        
        /**
         @private
         Update lower lh of a node but avoid using UpperLeft/Right lhs to update zero-blength
         This function is called after reading a tree from an input file, thus, UpperLeft/Right lhs have not yet been computed
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void updateLowerLhAvoidUsingUpperLRLh(cmaple::RealNumType & total_lh, std::unique_ptr<SeqRegions>& new_lower_lh, PhyloNode& node, const std::unique_ptr<SeqRegions>& lower_lh_1, const std::unique_ptr<SeqRegions>& lower_lh_2, const cmaple::Index  neighbor_1_index, PhyloNode& neighbor_1, const cmaple::Index  neighbor_2_index, PhyloNode& neighbor_2, const cmaple::PositionType & seq_length);
        
        /**
         @private
         compute the likelihood contribution of (the upper branch of) a node
         @throw std::logic\_error if unexpected values/behaviors found during the operations
         */
        template <const cmaple::StateType  num_states>
        void computeLhContribution(cmaple::RealNumType & total_lh, std::unique_ptr<SeqRegions>& new_lower_lh, PhyloNode& node, const std::unique_ptr<SeqRegions>& lower_lh_1, const std::unique_ptr<SeqRegions>& lower_lh_2, const cmaple::Index  neighbor_1_index, PhyloNode& neighbor_1, const cmaple::Index  neighbor_2_index, PhyloNode& neighbor_2, const cmaple::PositionType & seq_length);
        
        /**
         @private
         Employ Depth First Search to do a task at internal nodes
         */
        template <void(Tree::*task)(cmaple::RealNumType &, std::unique_ptr<SeqRegions>&, PhyloNode&, const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, const cmaple::Index , PhyloNode&, const cmaple::Index , PhyloNode&, const cmaple::PositionType &)>
        cmaple::RealNumType  performDFS();
        
        /**
         @private
         Update model parameters from an alignment and a tree
         @throw std::logic\_error if the reference genome is empty
         */
        template <const cmaple::StateType  num_states>
        void updateModelParams ();
        
        /**
         @private
         Employ Depth First Search to do a task at leaves
         */
        template <void(Tree::*task)(PhyloNode&, const cmaple::Index , const cmaple::Index )>
        void performDFSAtLeave();
        
        // NHANLT: Debug aLRT
        // void log_current(std::stack<cmaple::Index>& node_stack_aLRT);
    };

    /** \brief Customized << operator to output the tree string in a (bifurcating) NEWICK format
     */
    std::ostream& operator<<(std::ostream& out_stream, cmaple::Tree& tree);

    /** \brief Customized >> operator to read the tree from a stream
     */
    std::istream& operator>>(std::istream& in_stream, cmaple::Tree& tree);
}
