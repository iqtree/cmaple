#include "../alignment/alignment.h"
#include "../model/model.h"
#include "updatingnode.h"
#include "rootcandidate.h"
#include "altbranch.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#pragma once

namespace cmaple {
/** The structure of a phylogenetic tree */
class Tree {
 public:
  /*!
   * Types of tree search
   */
  enum TreeSearchType {
    FAST_TREE_SEARCH,   /*!< No tree search (placement only) */
    NORMAL_TREE_SEARCH, /*!< Consider pruning branches only at newly-added nodes
                           when seeking SPR moves */
    EXHAUSTIVE_TREE_SEARCH, /*!< Consider all nodes when seeking SPR moves */
    UNKNOWN_TREE_SEARCH,       /*!< Unknown (not specified) */
  };

  /*!
   * Types of trees
   */
  enum TreeType {
    BIN_TREE,     /*!< Binary tree */
    MUL_TREE,     /*!< Multifurcating tree */
    UNKNOWN_TREE, /*!< Unknown tree type */
  };

  // ----------------- BEGIN OF PUBLIC APIs ------------------------------------
  // //
  /*! \brief Constructor from a stream of a (bifurcating or multifurcating) tree
   * (with/without branch lengths in NEWICK format), which may or may not
   * contain all taxa in the alignment. Model parameters (if not fixed) will be
   * estimated according to the input tree and the alignment.
   * @param[in] aln An alignment
   * @param[in] model A substitution model
   * @param[in] tree_stream A stream of an input tree
   * @param[in] fixed_blengths TRUE to keep the input branch lengths unchanged
   * (optional)
   * @param[in] params an instance of Params which stores all program parameters
   * (optional). Users can build an instance of Params and specify several
   * parameters by using ParamsBuilder
   * @throw std::invalid\_argument If any of the following situation occurs.
   * - the sequence type is unsupported (neither DNA (for nucleotide data) nor
   * AA (for protein data))
   * - the alignment is empty
   * - the model is unknown/unsupported
   * - the tree is in an incorrect format
   *
   * @throw std::logic\_error if any of the following situations occur.
   * - any taxa in the tree (if specified) is not found in the alignment
   * - unexpected values/behaviors found during the operations
   *
   * @throw std::bad\_alloc if failing to allocate memory to store the tree
   */
  Tree(Alignment* aln,
       Model* model,
       std::istream& tree_stream,
       const bool fixed_blengths = false,
       std::unique_ptr<cmaple::Params>&& params = nullptr);

  /*! \brief Constructor from an optional (bifurcating or multifurcating) tree
   * (with/without branch lengths in NEWICK format), which may or may not
   * contain all taxa in the alignment. If users specify an input tree, model
   * parameters (if not fixed) will be estimated according to that tree and the
   * alignment.
   * @param[in] aln An alignment
   * @param[in] model A substitution model
   * @param[in] tree_filename Name of a tree file (optinal)
   * @param[in] fixed_blengths TRUE to keep the input branch lengths unchanged
   * (optional)
   * @param[in] params an instance of Params which stores all program parameters
   * (optional). Users can build an instance of Params and specify several
   * parameters by using ParamsBuilder
   * @throw std::invalid\_argument If any of the following situation occurs.
   * - the sequence type is unsupported (neither DNA (for nucleotide data) nor
   * AA (for protein data))
   * - the alignment is empty
   * - the model is unknown/unsupported
   * - the tree (if specified) is in an incorrect format
   *
   * @throw ios::failure if the tree file (if specified) is not found
   * @throw std::logic\_error if any of the following situations occur.
   * - any taxa in the tree (if specified) is not found in the alignment
   * - unexpected values/behaviors found during the operations
   *
   * @throw std::bad\_alloc if failing to allocate memory to store the tree
   */
  Tree(Alignment* aln,
       Model* model,
       const std::string& tree_filename = "",
       const bool fixed_blengths = false,
       std::unique_ptr<cmaple::Params>&& params = nullptr);

  /*! Destructor
   */
  ~Tree();

  /*! \brief Load a tree from a stream of a (bifurcating or multifurcating) tree
   *(with/without branch lengths) in NEWICK format, which may or may not contain
   *all taxa in the alignment. Model parameters (if not fixed) will be estimated
   *according to the input tree and the alignment.
   * @param[in] tree_stream A stream of an input tree
   * @param[in] fixed_blengths TRUE to keep the input branch lengths unchanged
   *(optional)
   * @throw std::invalid\_argument if the tree is empty or in an incorrect format
   * @throw std::logic\_error if any of the following situations occur.
   * - the attached substitution model is unknown/unsupported
   * - any taxa in the tree is not found in the alignment
   * - unexpected values/behaviors found during the operations
   *
   *@throw std::bad\_alloc if failing to allocate memory to store the tree
   */
  void load(std::istream& tree_stream, const bool fixed_blengths = false);

  /*! \brief Load a tree from a (bifurcating or multifurcating) tree
   * (with/without branch lengths) in NEWICK format, which may or may not
   * contain all taxa in the alignment. Model parameters (if not fixed) will be
   * estimated according to the input tree and the alignment.
   * @param[in] tree_filename Name of a tree file
   * @param[in] fixed_blengths TRUE to keep the input branch lengths unchanged
   * (optional)
   * @throw std::invalid\_argument if the tree is empty or in an incorrect format
   * @throw ios::failure if the tree file is not found
   * @throw std::logic\_error if any of the following situations occur.
   * - the attached substitution model is unknown/unsupported
   * - any taxa in the tree is not found in the alignment
   * - unexpected values/behaviors found during the operations
   *
   * @throw std::bad\_alloc if failing to allocate memory to store the tree
   */
  void load(const std::string& tree_filename,
            const bool fixed_blengths = false);

  /*! \brief Change the alignment
   * @param[in] aln An alignment
   * @throw std::invalid\_argument If the alignment is empty
   * @throw std::logic\_error if any of the following situations occur.
   * - any taxa in the tree is not found in the new alignment
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

  /*! \brief Do placement (using stepwise addition) to build an initial tree. Model
   * parameters (if not fixed) will be estimated during the placement process.
   * - If users didn't supply an input tree or supplied an incomplete tree
   * (which doesn't contain all the taxa in the alignment) when initializing the
   * tree (by Tree() constructor), this function will add new taxa (which are
   * not existed in the input tree) from the alignment to the tree.
   * - If users already supplied a complete tree, this function does nothing.
   *
   * @param[in] num_threads The number of threads
   * @param[out] out_stream The output message stream (optional)
   * @throw std::logic\_error if any of the following situations occur.
   * - the attached substitution model is unknown/unsupported
   * - unexpected values/behaviors found during the operations
   */
  void doPlacement(const int num_threads, std::ostream& out_stream = std::cout);

  /*! \brief Apply SPR moves to optimize the tree.
   * @param[in] num_threads The number of threads
   * @param[in] tree_search_type A type of tree search
   * @param[in] shallow_tree_search TRUE to enable a shallow tree search before
   * a deeper tree search
   * @param[out] out_stream The output message stream (optional)
   * @throw std::logic\_error if any of the following situations occur.
   * - the tree is empty
   * - the attached substitution model is unknown/unsupported
   * - unexpected values/behaviors found during the operations
   */
  void applySPR(const int num_threads, const TreeSearchType tree_search_type,
                       const bool shallow_tree_search, std::ostream& out_stream = std::cout);

  /*! \brief Optimize the branch lengths of the tree
   * @param[out] out_stream The output message stream (optional)
   * @throw std::logic\_error if any of the following situations occur.
   * - the tree is empty
   * - the attached substitution model is unknown/unsupported
   * - unexpected values/behaviors found during the operations
   */
  void optimizeBranch(std::ostream& out_stream = std::cout);

  /*! \brief Infer a maximum likelihood tree by executing doPlacement(), applySPR(),
   * and optimizeBranch()
   * - If users didn't supply an input tree or supplied an incomplete tree
   * (which doesn't contain all the taxa in the alignment) when initializing the
   * tree (by Tree() constructor), this function:
   *  + performs placement (i.e., adding missing taxa from the alignment to the
   * tree)
   *  + applies a NORMAL tree search (which does SPR moves only on newly-added
   * nodes)
   *  + optimizes all branch lengths
   * - If users already supplied a complete tree, this function, by default,
   * does neither placement nor tree search, but it optimizes all branch lengths.
   * If users want to keep the branch lengths fixed, they should set
   * fixed_blengths = true when initializing the tree (by Tree() constructor);
   * - If users want to use the input (complete/incomplete) tree as a starting
   * tree to perform placement (for an incomplete tree), then consider SPR moves
   * on all nodes, and optimize branch lengths, they can set tree_search_type
   * = EXHAUSTIVE
   *
   * @param[in] num_threads The number of threads (optional)
   * @param[in] tree_search_type A type of tree search (optional)
   * @param[in] shallow_tree_search TRUE to enable a shallow tree search before
   * a deeper tree search (optional)
   * @param[out] out_stream The output message stream (optional)
   * @throw std::invalid\_argument if tree\_search\_type is unknown
   * @throw std::logic\_error if any of the following situations occur.
   * - the attached substitution model is unknown/unsupported
   * - unexpected values/behaviors found during the operations
   */
  void infer(
      const int num_threads = 1,
      const TreeSearchType tree_search_type = NORMAL_TREE_SEARCH,
      const bool shallow_tree_search = false, std::ostream& out_stream = std::cout);

  /*! \brief Compute the log likelihood of the current tree, which may or may
   * not contain all taxa in the alignment
   * @return The log likelihood of the current tree
   * @throw std::logic\_error if any of the following situations occur.
   * - the tree is empty
   * - the attached substitution model is unknown/unsupported
   * - unexpected values/behaviors found during the operations
   */
  RealNumType computeLh();

  /*! \brief Compute branch supports
   * ([aLRT-SH](https://academic.oup.com/sysbio/article/59/3/307/1702850)) of
   * the current tree, which may or may not contain all taxa in the alignment
   * @param[in] num_threads The number of threads (optional)
   * @param[in] num_replicates A positive number of replicates (optional)
   * @param[in] epsilon A positive epsilon (optional), which is used to avoid
   * rounding effects, when the best and second best NNI trees have nearly
   * identical site log-likelihood values (see [Guindon et al.,
   * 2010](https://academic.oup.com/sysbio/article/59/3/307/1702850))
   * @param[in] allow_replacing_ML_tree TRUE to allow replacing the ML tree by a
   * higher likelihood tree found when computing branch supports (optional)
   * @param[out] out_stream The output message stream (optional)
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
  void computeBranchSupport(const int num_threads = 1,
                                   const int num_replicates = 1000,
                                   const double epsilon = 0.1,
                                   const bool allow_replacing_ML_tree = true,
                                   std::ostream& out_stream = std::cout);

  /*! \brief Export the phylogenetic tree  to a string in NEWICK format.
   * @param[in] tree_type The type of the output tree (optional): BIN_TREE
   * (bifurcating tree), MUL_TREE (multifurcating tree)
   * @param[in] show_branch_supports TRUE to output the branch supports (aLRT-SH
   * values)
   * @param[in] print_internal_id TRUE to print the id of internal nodes
   * @return A tree string in NEWICK format
   * @throw std::invalid\_argument if any of the following situations occur.
   * - tree\_type is unknown
   */
  std::string exportNewick(const TreeType tree_type = BIN_TREE,
                           const bool print_internal_id = false,
                           const bool show_branch_supports = true);

  // ----------------- END OF PUBLIC APIs ------------------------------------
  // //

  /*! \cond PRIVATE */
  /**
   TRUE to keep the branch lengths fixed
   */
  bool fixed_blengths = false;

  /**
   Branch length thresholds
   */
  cmaple::RealNumType default_blength, min_blength, max_blength,
      min_blength_mid, min_blength_sensitivity, half_max_blength,
      half_min_blength_mid, double_min_blength;

  /**
   Program parameters
   */
  // std::optional<cmaple::Params> params;
  std::unique_ptr<cmaple::Params> params = nullptr;

  /**
   Alignment
   */
  Alignment* aln = nullptr;

  /**
   Evolutionary model
   */
  ModelBase* model = nullptr;

  /**
   cumulative rates
   */
  cmaple::RealNumType* cumulative_rate = nullptr;

  /**
   cumulative bases
   */
  std::vector<std::vector<cmaple::PositionType>> cumulative_base;

  /**
   Vector of phylonodes
   */
  std::vector<PhyloNode> nodes;
    
    /**
     Number of existing nodes before sample placement
     */
    NumSeqsType num_exiting_nodes;

  /**
   Vector of likelihood contributions of internal nodes
   */
  std::vector<NodeLh> node_lhs;
    
    /**
     Vector of the annotations of nodes
     */
    std::vector<std::string> annotations;
    
    /**
     Vector of the SPRTA scores
     */
    std::vector<RealNumType> sprta_scores;
    
    /**
     Vector of root supports
     */
    std::vector<RealNumType> root_supports;
    
    /**
     Vector of alternative branches (when computing SPRTA)
     */
    std::vector<std::vector<cmaple::AltBranch>> sprta_alt_branches;
    
    /**
     The inverse of sprta_alt_branches
     highlight which nodes could be placed (with probability above threshold) on the branch above the current node
     */
    std::vector<std::vector<cmaple::AltBranch>> sprta_support_list;
    
    /**
     Vector of number of descendants of nodes
     */
    std::vector<NumSeqsType> num_descendants;
    
    /**
     Vector of internal node names
     */
    std::vector<NumSeqsType> internal_names;

  /**
   (vector) Index of root in the vector of phylonodes
   */
  cmaple::NumSeqsType root_vector_index;

  /**
   A backup of sequence names attached to the current tree, in cases that users
   re-read the alignment NHANLT: @TODO - avoid this redundant vector (sequence
   names are store in alignment and here)
   */
  std::vector<std::string> seq_names;

  /**
   a vector denote whether a sequence in the alignment is added to the tree or
   not
   */
  std::vector<bool> sequence_added;
    
  /**
   TRUE if branch support (i.e. aLRT-SH) computed
  */
  bool aLRT_SH_computed = false;


  /*!
   * Apply some minor changes (collapsing zero-branch leaves into less-info
   * sequences, re-estimating model parameters) to make the processes of
   * outputting then re-inputting a tree result in a consistent tree
   * @throw std::logic\_error if unexpected values/behaviors found during the
   * operations
   */
  void makeTreeInOutConsistent();

  /**
   * Parse type of tree search from a string
   * @param[in] tree_search_type Tree search type in string
   * @return a TreeSearchType
   */
  static TreeSearchType parseTreeSearchType(
      const std::string& tree_search_type);

  /**
   * Get tree search type in string
   * @param[in] tree_search_type a type of tree search
   * @return a  type of tree search from a string
   */
  static std::string getTreeSearchStr(const TreeSearchType tree_search_type);

  /**
   * Parse tree type from a string
   * @param tree_type_str a tree type in string
   * @return a TreeType
   */
  static TreeType parseTreeType(const std::string& tree_type_str);
    
    /*! \brief Export the phylogenetic tree  to a string in NEXUS format.
     * @param[in] tree_type The type of the output tree (optional): BIN_TREE
     * (bifurcating tree), MUL_TREE (multifurcating tree)
     * @param[in] show_branch_supports TRUE to output the branch supports (aLRT-SH
     * values)
     * @return A tree string in NEXUS format
     * @throw std::invalid\_argument if any of the following situations occur.
     * - tree\_type is unknown
     */
    std::string exportNexus(const TreeType tree_type = BIN_TREE,
                             const bool show_branch_supports = true);
    
    /**
     Export a TSV file that contains useful information from SPRTA
     */
    std::string exportTSV();
    
  /*! \endcond */

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
  typedef void (Tree::*DoInferencePtrType)(const int, const TreeSearchType,
                                                  const bool, std::ostream&);
  DoInferencePtrType doInferencePtr;

  /**
      Pointer  to doPlacement method
   */
  typedef void (Tree::*DoPlacementPtrType)(const int, std::ostream&);
  DoPlacementPtrType doPlacementPtr;

  /**
      Pointer  to applySPR method
   */
  typedef void (Tree::*ApplySPRPtrType)(const int, const TreeSearchType,
                                               const bool, std::ostream&);
  ApplySPRPtrType applySPRPtr;

  /**
      Pointer  to optimizeBranch method
   */
  typedef void (Tree::*OptimizeBranchPtrType)(std::ostream&);
  OptimizeBranchPtrType optimizeBranchPtr;

  /**
      Pointer  to computeLh method
   */
  typedef RealNumType (Tree::*computeLhPtrType)();
  computeLhPtrType computeLhPtr;

  /**
      Pointer  to computeBranchSupport method
   */
  typedef void (Tree::*computeBranchSupportPtrType)(const int,
                                                           const int,
                                                           const double,
                                                           const bool, std::ostream&);
  computeBranchSupportPtrType computeBranchSupportPtr;

  typedef void (Tree::*MakeTreeInOutConsistentPtrType)();
  MakeTreeInOutConsistentPtrType makeTreeInOutConsistentPtr;

  /*! Template of loadTree()
   @param[in] tree_stream A stream of the input tree
   @param[in] fixed_blengths TRUE to keep the input branch lengths unchanged
   (optional)
   */
  template <const cmaple::StateType num_states>
  void loadTreeTemplate(std::istream& tree_stream, const bool fixed_blengths);

  /*! Template of changeAln()
   */
  template <const cmaple::StateType num_states>
  void changeAlnTemplate(Alignment* aln);

  /*! Template of changeModel()
   */
  template <const cmaple::StateType num_states>
  void changeModelTemplate(Model* model);

  /*! Template of doInference()
   */
  template <const cmaple::StateType num_states>
  void doInferenceTemplate(const int num_threads, const TreeSearchType tree_search_type,
                                  const bool shallow_tree_search, std::ostream& out_stream);

  /*! Template of doPlacement()
   */
  template <const cmaple::StateType num_states>
  void doPlacementTemplate(const int num_threads, std::ostream& out_stream);

  /*! Template of applySPR()
   */
  template <const cmaple::StateType num_states>
  void applySPRTemplate(const int num_threads, const TreeSearchType tree_search_type,
                               const bool shallow_tree_search, std::ostream& out_stream);

  /*! Template of optimizeBranch()
   */
  template <const cmaple::StateType num_states>
  void optimizeBranchTemplate(std::ostream& out_stream);

  /*! Template of computeLh()
   */
  template <const cmaple::StateType num_states>
  RealNumType computeLhTemplate();

  /*! Template of computeBranchSupport()
   */
  template <const cmaple::StateType num_states>
  void computeBranchSupportTemplate(const int num_threads,
                                           const int num_replicates,
                                           const double epsilon,
                                           const bool allow_replacing_ML_tree,
                                           std::ostream& out_stream);

  /*! Template of makeTreeInOutConsistent()
   */
  template <const cmaple::StateType num_states>
  void makeTreeInOutConsistentTemplate();

  /*! Setup function pointers
   @throw std::invalid\_argument If the sequence type is unsupported (neither
   DNA (for nucleotide data) nor AA (for protein data))
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
   * - the sequence type is unsupported (neither DNA (for nucleotide data) nor
   * AA (for protein data))
   * - the alignment is empty
   * - model is unknown/unsupported
   *
   * @throw std::logic\_error if the reference genome is empty
   */
  void initTree(Alignment* aln,
                Model* model,
                std::unique_ptr<cmaple::Params>&& params);

  /**
   Compute cumulative rate of the ref genome
   @throw std::logic\_error if the reference genome is empty
   */
  void computeCumulativeRate();

  /*! Optimize the tree topology
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void optimizeTreeTopology(const int num_threads,
                            const TreeSearchType tree_search_type,
                            bool short_range_search = false);

  /**
   Traverse the intial tree from root to re-calculate all non-lower likelihoods
   regarding the latest/final estimated model parameters
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void refreshAllNonLowerLhs();
    
  /*! Expand data vectors after tree expansion
  */
  void expandVectorsAfterTreeExpansion();

  /**
   Try to improve a subtree rooted at node with SPR moves
   @return total improvement
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  cmaple::RealNumType improveSubTree(const cmaple::Index index,
                                     PhyloNode& node,
                                     const TreeSearchType tree_search_type,
                                     bool short_range_search,
                                     std::vector<std::pair<cmaple::Index, double>>& SPR_found_vec,
                                     bool SPR_search_only = false);

  /**
   Calculate derivative starting from coefficients.
   @return derivative
   */
  cmaple::RealNumType calculateDerivative(
      const std::vector<cmaple::RealNumType>& coefficient_vec,
      const cmaple::RealNumType delta_t);

  /**
   Examine placing a sample at a mid-branch point
   */
  template <const cmaple::StateType num_states>
  void examineSamplePlacementMidBranch(
      cmaple::Index& selected_node_index,
      const std::unique_ptr<SeqRegions>& mid_branch_lh,
      cmaple::RealNumType& best_lh_diff,
      bool& is_mid_branch,
      cmaple::RealNumType& lh_diff_mid_branch,
      TraversingNode& current_extended_node,
      const std::unique_ptr<SeqRegions>& sample_regions);

  /**
   Examine placing a sample as a descendant of an existing node
   */
  template <const cmaple::StateType num_states>
  void examineSamplePlacementAtNode(
      cmaple::Index& selected_node_index,
      const std::unique_ptr<SeqRegions>& total_lh,
      cmaple::RealNumType& best_lh_diff,
      bool& is_mid_branch,
      cmaple::RealNumType& lh_diff_at_node,
      cmaple::RealNumType& lh_diff_mid_branch,
      cmaple::RealNumType& best_up_lh_diff,
      cmaple::RealNumType& best_down_lh_diff,
      cmaple::Index& best_child_index,
      TraversingNode& current_extended_node,
      const std::unique_ptr<SeqRegions>& sample_regions);

  /**
   Traverse downwards polytomy for more fine-grained placement
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void finetuneSamplePlacementAtNode(
      const PhyloNode& selected_node,
      cmaple::RealNumType& best_down_lh_diff,
      cmaple::Index& best_child_index,
      const std::unique_ptr<SeqRegions>& sample_regions);
    
/**
 Check if the new placement is sufficiently different from the original one
*/
template <const cmaple::StateType num_states>
bool isDiffFromOrigPlacement(
    const cmaple::Index ori_parent_index,
    cmaple::Index& new_placement_index,
    const cmaple::RealNumType best_mid_top_blength,
    const cmaple::RealNumType best_mid_bottom_blength,
    bool& is_root_considered);

  /**
   Add start nodes for seeking a placement for a subtree
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void addStartingNodes(const cmaple::Index& node_index,
                        PhyloNode& node,
                        const cmaple::Index& other_child_node_index,
                        const cmaple::RealNumType best_lh_diff,
                        std::stack<std::unique_ptr<UpdatingNode>>& node_stack);
    
    /**
     Add start nodes for root assessment
     @throw std::logic\_error if unexpected values/behaviors found during the
     operations
     */
    template <const cmaple::StateType num_states>
    void addStartingRootCandidate(const NumSeqsType& node_vec_id,
                                  std::stack<std::unique_ptr<RootCandidate>>& node_stack);
    
    /**
     Add nodes for root assessment
     @throw std::logic\_error if unexpected values/behaviors found during the
     operations
     */
    template <const cmaple::StateType num_states>
    void addChildrenAsRootCandidate(const std::unique_ptr<SeqRegions>& incoming_regions_ref,
                          const cmaple::RealNumType branch_length,
                          const cmaple::RealNumType lh_deducted,
                          const cmaple::RealNumType lh_diff,
                          const short int failure_count,
                          PhyloNode& parent_node,
                          std::stack<std::unique_ptr<RootCandidate>>& node_stack);
    
    /**
     Reroot the current tree
     @throw std::logic\_error if unexpected values/behaviors found during the
     operations
     */
    template <const cmaple::StateType num_states>
    void reroot(const NumSeqsType& new_root_vec_id);
    
    /**
     Transfer annotations due to rerooting
     @throw std::logic\_error if unexpected values/behaviors found during the
     operations
     */
    void transferAnnotations(const NumSeqsType& new_root_vec_id);

  /**
   Examine placing a subtree at a mid-branch point
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  bool examineSubtreePlacementMidBranch(
      cmaple::Index& best_node_index,
      PhyloNode& current_node,
      cmaple::RealNumType& best_lh_diff,
      cmaple::RealNumType& best_lh_diff_before_bl_opt,
      bool& is_mid_branch,
      cmaple::RealNumType& lh_diff_at_node,
      cmaple::RealNumType& lh_diff_mid_branch,
      cmaple::RealNumType& best_up_lh_diff,
      cmaple::RealNumType& best_down_lh_diff,
      std::unique_ptr<UpdatingNode>& updating_node,
      const std::unique_ptr<SeqRegions>& subtree_regions,
      const cmaple::RealNumType threshold_prob,
      const cmaple::RealNumType removed_blength,
      const cmaple::Index top_node_index,
      const cmaple::Index ori_parent_index,
      std::unique_ptr<SeqRegions>& bottom_regions,
      RealNumType& opt_appending_blength,
      RealNumType& opt_mid_top_blength,
      RealNumType& opt_mid_bottom_blength,
      std::vector<AltBranch>& alt_branches,
      bool& is_root_considered);

  /**
   Examine placing a subtree as a descendant of an existing node
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  bool examineSubTreePlacementAtNode(
      cmaple::Index& best_node_index,
      PhyloNode& current_node,
      cmaple::RealNumType& best_lh_diff,
      bool& is_mid_branch,
      cmaple::RealNumType& lh_diff_at_node,
      cmaple::RealNumType& lh_diff_mid_branch,
      cmaple::RealNumType& best_up_lh_diff,
      cmaple::RealNumType& best_down_lh_diff,
      std::unique_ptr<UpdatingNode>& updating_node,
      const std::unique_ptr<SeqRegions>& subtree_regions,
      const cmaple::RealNumType threshold_prob,
      const cmaple::RealNumType removed_blength,
      const cmaple::Index top_node_index);

  /**
   Add a child node for downwards traversal when seeking a new subtree placement
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void addChildSeekSubtreePlacement(
      const cmaple::Index child_1_index,
      PhyloNode& child_1,
      PhyloNode& child_2,
      const cmaple::RealNumType& lh_diff_at_node,
      const std::unique_ptr<UpdatingNode>& updating_node,
      std::stack<std::unique_ptr<UpdatingNode>>& node_stack,
      const cmaple::RealNumType threshold_prob);

  /**
   Add neighbor nodes (parent/sibling) for traversal when seeking a new subtree
   placement
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  bool addNeighborsSeekSubtreePlacement(
      PhyloNode& current_node,
      const cmaple::Index other_child_index,
      std::unique_ptr<SeqRegions>&& bottom_regions,
      const cmaple::RealNumType& lh_diff_at_node,
      const std::unique_ptr<UpdatingNode>& updating_node,
      std::stack<std::unique_ptr<UpdatingNode>>& node_stack,
      const cmaple::RealNumType threshold_prob);

  /**
   Check whether we can obtain a higher likelihood with a shorter length for an
   existing branch
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states,
            cmaple::RealNumType (Tree::*calculatePlacementCost)(
                const std::unique_ptr<SeqRegions>&,
                const std::unique_ptr<SeqRegions>&,
                const cmaple::RealNumType)>
  bool tryShorterBranch(
      const cmaple::RealNumType current_blength,
      std::unique_ptr<SeqRegions>& best_child_regions,
      const std::unique_ptr<SeqRegions>& sample,
      const std::unique_ptr<SeqRegions>& upper_left_right_regions,
      const std::unique_ptr<SeqRegions>& lower_regions,
      cmaple::RealNumType& best_split_lh,
      cmaple::RealNumType& best_branch_length_split,
      const cmaple::RealNumType new_branch_length,
      const bool try_first_branch);

  /**
   Check whether we can obtain a higher likelihood with a shorter length at root
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void tryShorterBranchAtRoot(const std::unique_ptr<SeqRegions>& sample,
                              const std::unique_ptr<SeqRegions>& lower_regions,
                              std::unique_ptr<SeqRegions>& best_parent_regions,
                              cmaple::RealNumType& best_root_blength,
                              cmaple::RealNumType& best_parent_lh,
                              const cmaple::RealNumType fixed_blength);

  /**
   Check whether we can obtain a higher likelihood with a shorter length for the
   new branch at root
   */
  template <const cmaple::StateType num_states>
  bool tryShorterNewBranchAtRoot(
      const std::unique_ptr<SeqRegions>& sample,
      const std::unique_ptr<SeqRegions>& lower_regions,
      std::unique_ptr<SeqRegions>& best_parent_regions,
      cmaple::RealNumType& best_length,
      cmaple::RealNumType& best_parent_lh,
      const cmaple::RealNumType fixed_blength);

  /**
   Check whether we can obtain a higher likelihood with a longer length for the
   new branch at root
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  bool tryLongerNewBranchAtRoot(
      const std::unique_ptr<SeqRegions>& sample,
      const std::unique_ptr<SeqRegions>& lower_regions,
      std::unique_ptr<SeqRegions>& best_parent_regions,
      cmaple::RealNumType& best_length,
      cmaple::RealNumType& best_parent_lh,
      const cmaple::RealNumType fixed_blength);

  /**
   Estimate the length for a new branch at root
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void estimateLengthNewBranchAtRoot(
      const std::unique_ptr<SeqRegions>& sample,
      const std::unique_ptr<SeqRegions>& lower_regions,
      std::unique_ptr<SeqRegions>& best_parent_regions,
      cmaple::RealNumType& best_length,
      cmaple::RealNumType& best_parent_lh,
      const cmaple::RealNumType fixed_blength,
      const cmaple::RealNumType short_blength_thresh,
      const bool optional_check);

  /**
   Check whether we can obtain a higher likelihood with a shorter length for the
   new branch
   */
  template <cmaple::RealNumType (Tree::*calculatePlacementCost)(
      const std::unique_ptr<SeqRegions>&,
      const std::unique_ptr<SeqRegions>&,
      const cmaple::RealNumType)>
  bool tryShorterNewBranch(
      const std::unique_ptr<SeqRegions>& best_child_regions,
      const std::unique_ptr<SeqRegions>& sample,
      cmaple::RealNumType& best_blength,
      cmaple::RealNumType& new_branch_lh,
      const cmaple::RealNumType short_blength_thresh);

  /**
   Check whether we can obtain a higher likelihood with a longer length for the
   new branch
   */
  template <cmaple::RealNumType (Tree::*calculatePlacementCost)(
      const std::unique_ptr<SeqRegions>&,
      const std::unique_ptr<SeqRegions>&,
      const cmaple::RealNumType)>
  void tryLongerNewBranch(const std::unique_ptr<SeqRegions>& best_child_regions,
                          const std::unique_ptr<SeqRegions>& sample,
                          cmaple::RealNumType& best_blength,
                          cmaple::RealNumType& new_branch_lh,
                          const cmaple::RealNumType long_blength_thresh);

  /**
   Estimate the length for a new branch
   */
  template <cmaple::RealNumType (Tree::*calculatePlacementCost)(
      const std::unique_ptr<SeqRegions>&,
      const std::unique_ptr<SeqRegions>&,
      const cmaple::RealNumType)>
  void estimateLengthNewBranch(
      const cmaple::RealNumType best_split_lh,
      const std::unique_ptr<SeqRegions>& best_child_regions,
      const std::unique_ptr<SeqRegions>& sample,
      cmaple::RealNumType& best_blength,
      const cmaple::RealNumType long_blength_thresh,
      const cmaple::RealNumType short_blength_thresh,
      const bool optional_check);

  /**
   Connect a new sample to a branch
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void connectNewSample2Branch(
      std::unique_ptr<SeqRegions>& sample,
      const cmaple::NumSeqsType seq_name_index,
      const cmaple::Index sibling_node_index,
      PhyloNode& sibling_node,
      const cmaple::RealNumType top_distance,
      const cmaple::RealNumType down_distance,
      const cmaple::RealNumType best_blength,
      std::unique_ptr<SeqRegions>& best_child_regions,
      const std::unique_ptr<SeqRegions>& upper_left_right_regions);

  /**
   Connect a new sample to root
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void connectNewSample2Root(std::unique_ptr<SeqRegions>& sample,
                             const cmaple::NumSeqsType seq_name_index,
                             const cmaple::Index sibling_node_index,
                             PhyloNode& sibling_node,
                             const cmaple::RealNumType best_root_blength,
                             const cmaple::RealNumType best_length2,
                             std::unique_ptr<SeqRegions>& best_parent_regions);

  /**
   Place a subtree as a descendant of a node
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void placeSubTreeAtNode(const cmaple::Index selected_node_index,
                          const cmaple::Index subtree_index,
                          PhyloNode& subtree,
                          const std::unique_ptr<SeqRegions>& subtree_regions,
                          const cmaple::RealNumType new_branch_length,
                          const cmaple::RealNumType new_lh);

  /**
   Place a subtree at a mid-branch point
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void placeSubTreeMidBranch(const cmaple::Index selected_node_index,
                             const cmaple::Index subtree_index,
                             PhyloNode& subtree,
                             const std::unique_ptr<SeqRegions>& subtree_regions,
                             const cmaple::RealNumType new_branch_length,
                             const cmaple::RealNumType opt_appending_blength,
                             const cmaple::RealNumType opt_mid_top_blength,
                             const cmaple::RealNumType opt_mid_bottom_blength,
                             const cmaple::RealNumType new_lh);

  /**
   Connect a subtree to a branch
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <
      const cmaple::StateType num_states,
      void (Tree::*updateRegionsSubTree)(PhyloNode&,
                                         PhyloNode&,
                                         PhyloNode&,
                                         std::unique_ptr<SeqRegions>&&,
                                         const std::unique_ptr<SeqRegions>&,
                                         const std::unique_ptr<SeqRegions>&,
                                         const std::unique_ptr<SeqRegions>&,
                                         cmaple::RealNumType&)>
  void connectSubTree2Branch(
      const std::unique_ptr<SeqRegions>& subtree_regions,
      const std::unique_ptr<SeqRegions>& lower_regions,
      const cmaple::Index subtree_index,
      PhyloNode& subtree,
      const cmaple::Index sibling_node_index,
      PhyloNode& sibling_node,
      const cmaple::RealNumType top_distance,
      const cmaple::RealNumType down_distance,
      cmaple::RealNumType& best_blength,
      std::unique_ptr<SeqRegions>&& best_child_regions,
      const std::unique_ptr<SeqRegions>& upper_left_right_regions);

  /**
   Connect a subtree to root
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void connectSubTree2Root(const cmaple::Index subtree_index,
                           PhyloNode& subtree,
                           const std::unique_ptr<SeqRegions>& subtree_regions,
                           const std::unique_ptr<SeqRegions>& lower_regions,
                           const cmaple::Index sibling_node_index,
                           PhyloNode& sibling_node,
                           const cmaple::RealNumType best_root_blength,
                           const cmaple::RealNumType best_length2,
                           std::unique_ptr<SeqRegions>&& best_parent_regions);

  /**
   Update next_node_1->partial_lh and new_internal_node->partial_lh after
   placing a subtree in common cases (e.g., at a mid-branch point, under a node)

   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void updateRegionsPlaceSubTree(
      PhyloNode& subtree,
      PhyloNode& sibling_node,
      PhyloNode& internal,
      std::unique_ptr<SeqRegions>&& best_child_regions,
      const std::unique_ptr<SeqRegions>& subtree_regions,
      const std::unique_ptr<SeqRegions>& upper_left_right_regions,
      const std::unique_ptr<SeqRegions>& lower_regions,
      cmaple::RealNumType& best_blength);

  /**
   Update next_node_1->partial_lh and new_internal_node->partial_lh after
   placing a subtree in other cases (i.e., above a node)

   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void updateRegionsPlaceSubTreeAbove(
      PhyloNode& subtree,
      PhyloNode& sibling_node,
      PhyloNode& internal,
      std::unique_ptr<SeqRegions>&& best_child_regions,
      const std::unique_ptr<SeqRegions>& subtree_regions,
      const std::unique_ptr<SeqRegions>& upper_left_right_regions,
      const std::unique_ptr<SeqRegions>& lower_regions,
      cmaple::RealNumType& best_blength);

  /**
   Handle polytomy when placing a subtree
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void handlePolytomyPlaceSubTree(
      const cmaple::Index selected_node_index,
      PhyloNode& selected_node,
      const std::unique_ptr<SeqRegions>& subtree_regions,
      const cmaple::RealNumType new_branch_length,
      cmaple::RealNumType& best_down_lh_diff,
      cmaple::Index& best_child_index,
      cmaple::RealNumType& best_child_blength_split,
      std::unique_ptr<SeqRegions>& best_child_regions);

  /**
   Update likelihood at mid-branch point
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void updateMidBranchLh(
      const cmaple::Index node_index,
      PhyloNode& node,
      const std::unique_ptr<SeqRegions>& parent_upper_regions,
      std::stack<cmaple::Index>& node_stack,
      bool& update_blength);

  /**
   Compute Upper Left/Right regions at a node, updating the top branch length if
   neccessary
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  std::unique_ptr<SeqRegions> computeUpperLeftRightRegions(
      const cmaple::Index node_index,
      PhyloNode& node,
      const cmaple::MiniIndex next_node_mini,
      const std::unique_ptr<SeqRegions>& parent_upper_regions,
      std::stack<cmaple::Index>& node_stack,
      bool& update_blength);

  /**
   Update the PartialLh (seqregions) at a node if the new one is different from
   the current one
   */
  bool updateNewPartialIfDifferent(
      PhyloNode& node,
      const cmaple::MiniIndex next_node_mini,
      std::unique_ptr<SeqRegions>& upper_left_right_regions,
      std::stack<cmaple::Index>& node_stack,
      const cmaple::PositionType seq_length);

  /**
   Handle cases when the new seqregions is null/empty: (1) update the branch
   length; or (2) return an error message
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void inline handleNullNewRegions(const cmaple::Index index,
                                   PhyloNode& node,
                                   const bool do_update_zeroblength,
                                   std::stack<cmaple::Index>& node_stack,
                                   bool& update_blength,
                                   const std::string err_msg) {
    if (do_update_zeroblength) {
      updateZeroBlength<num_states>(index, node, node_stack);
      update_blength = true;
    } else
      throw std::logic_error(err_msg);
  }

  /**
   Update partial_lh comming from the parent node

   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void updatePartialLhFromParent(
      const cmaple::Index index,
      PhyloNode& node,
      std::stack<cmaple::Index>& node_stack,
      const std::unique_ptr<SeqRegions>& parent_upper_regions,
      const cmaple::PositionType seq_length);

  /**
   Update partial_lh comming from the children

   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void updatePartialLhFromChildren(
      const cmaple::Index index,
      PhyloNode& node,
      std::stack<cmaple::Index>& node_stack,
      const std::unique_ptr<SeqRegions>& parent_upper_regions,
      const bool is_non_root,
      const cmaple::PositionType seq_length);

  /**
   Compute the mid-branch region for a node/branch
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  inline void computeMidBranchRegions(
      PhyloNode& node,
      std::unique_ptr<SeqRegions>& regions_2_update,
      const SeqRegions& parent_upper_lr_lh) {
    std::unique_ptr<SeqRegions>& lower_lh = node.getPartialLh(cmaple::TOP);
    cmaple::RealNumType half_branch_length = node.getUpperLength() * 0.5;
    parent_upper_lr_lh.mergeUpperLower<num_states>(
        regions_2_update, half_branch_length, *lower_lh, half_branch_length,
        aln, model, params->threshold_prob);
  }

  /**
   Refresh all non-lowerlhs traversing from a parent node
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void refreshNonLowerLhsFromParent(cmaple::Index& node_index,
                                    cmaple::Index& last_node_index);

  /**
   Refresh upper left/right regions
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void refreshUpperLR(const cmaple::Index node_index,
                      PhyloNode& node,
                      const cmaple::Index neighbor_index,
                      std::unique_ptr<SeqRegions>& replaced_regions,
                      const std::unique_ptr<SeqRegions>& parent_upper_lr_lh);

  /**
   Calculate coefficients when merging R with O to estimate a branch length
   */
  template <const cmaple::StateType num_states>
  void estimateBlength_R_O(const SeqRegion& seq1_region,
                           const SeqRegion& seq2_region,
                           const cmaple::RealNumType total_blength,
                           const cmaple::PositionType end_pos,
                           cmaple::RealNumType& coefficient,
                           std::vector<cmaple::RealNumType>& coefficient_vec);

  /**
   Calculate coefficients when merging R with ACGT to estimate a branch length
   */
  void estimateBlength_R_ACGT(
      const SeqRegion& seq1_region,
      const cmaple::StateType seq2_state,
      const cmaple::RealNumType total_blength,
      const cmaple::PositionType end_pos,
      std::vector<cmaple::RealNumType>& coefficient_vec);

  /**
   Calculate coefficients when merging O with X(i.e., O, R, ACGT) to estimate a
   branch length
   */
  template <const cmaple::StateType num_states>
  void estimateBlength_O_X(const SeqRegion& seq1_region,
                           const SeqRegion& seq2_region,
                           const cmaple::RealNumType total_blength,
                           const cmaple::PositionType end_pos,
                           cmaple::RealNumType& coefficient,
                           std::vector<cmaple::RealNumType>& coefficient_vec);

  /**
   Calculate coefficients when merging ACGT with O to estimate a branch length
   */
  template <const cmaple::StateType num_states>
  void estimateBlength_ACGT_O(
      const SeqRegion& seq1_region,
      const SeqRegion& seq2_region,
      const cmaple::RealNumType total_blength,
      cmaple::RealNumType& coefficient,
      std::vector<cmaple::RealNumType>& coefficient_vec);

  /**
   Calculate coefficients when merging ACGT with R/ACGT to estimate a branch
   length
   */
  void estimateBlength_ACGT_RACGT(
      const SeqRegion& seq1_region,
      const SeqRegion& seq2_region,
      const cmaple::RealNumType total_blength,
      const cmaple::PositionType end_pos,
      std::vector<cmaple::RealNumType>& coefficient_vec);

  /**
   Estimate a branch length from coefficients
   */
  cmaple::RealNumType estimateBlengthFromCoeffs(
      cmaple::RealNumType& coefficient,
      const std::vector<cmaple::RealNumType> coefficient_vec);

  /**
   Handle branch length changed when improve a subtree
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void handleBlengthChanged(PhyloNode& node,
                            const cmaple::Index node_index,
                            const cmaple::RealNumType best_blength);

  /**
   Optimize a branch length before seeking an SPR move for a subtree
   */
  template <const cmaple::StateType num_states>
  void optimizeBlengthBeforeSeekingSPR(
      PhyloNode& node,
      cmaple::RealNumType& best_blength,
      cmaple::RealNumType& best_lh,
      bool& blength_changed,
      const std::unique_ptr<SeqRegions>& parent_upper_lr_lh,
      const std::unique_ptr<SeqRegions>& lower_lh);

  /**
   Check and apply SPR move
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void checkAndApplySPR(const cmaple::RealNumType best_lh_diff,
                        const cmaple::RealNumType best_blength,
                        const cmaple::RealNumType opt_appending_blength,
                        const cmaple::RealNumType opt_mid_top_blength,
                        const cmaple::RealNumType opt_mid_bottom_blength,
                        const cmaple::RealNumType best_lh,
                        const cmaple::Index node_index,
                        PhyloNode& node,
                        const cmaple::Index best_node_index,
                        const cmaple::Index parent_node_index,
                        const bool is_mid_node,
                        cmaple::RealNumType& total_improvement,
                        bool& topology_updated,
                        std::vector<std::pair<cmaple::Index, double>>& SPR_found_vec,
                        bool SPR_search_only);

  /**
   Create a new internal phylonode
   */
  inline void createAnInternalNode() { nodes.emplace_back(InternalNode()); }

  /**
   Create a new leaf phylonode
   */
  inline void createALeafNode(const cmaple::NumSeqsType new_seq_name_index) {
    nodes.emplace_back(LeafNode(
        new_seq_name_index));  //(PhyloNode(std::move(LeafNode(new_seq_name_index))));
  }

  /**
   Get partial_lh at a node by its index
   */
  std::unique_ptr<SeqRegions>& getPartialLhAtNode(const cmaple::Index index);

  /**
   Calculate the likelihood of an NNI neighbor
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  bool calculateNNILh(std::stack<cmaple::Index>& node_stack_aLRT,
                      cmaple::RealNumType& lh_diff,
                      PhyloNode& current_node,
                      PhyloNode& child_1,
                      PhyloNode& child_2,
                      PhyloNode& sibling,
                      PhyloNode& parent,
                      const cmaple::Index parent_index,
                      cmaple::RealNumType& lh_at_root,
                      const bool allow_replacing_ML_tree);

  /**
   Calculate the likelihood of an NNI neighbor on the branch connecting to root
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  bool calculateNNILhRoot(std::stack<cmaple::Index>& node_stack_aLRT,
                          cmaple::RealNumType& lh_diff,
                          std::unique_ptr<SeqRegions>& parent_new_lower_lh,
                          const cmaple::RealNumType& child_2_new_blength,
                          PhyloNode& current_node,
                          PhyloNode& child_1,
                          PhyloNode& child_2,
                          PhyloNode& sibling,
                          PhyloNode& parent,
                          const cmaple::Index parent_index,
                          cmaple::RealNumType& lh_at_root,
                          const bool allow_replacing_ML_tree);

  /**
   Calculate the likelihood of an NNI neighbor on the branch connecting to a
   non-root node
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  bool calculateNNILhNonRoot(std::stack<cmaple::Index>& node_stack_aLRT,
                             cmaple::RealNumType& lh_diff,
                             std::unique_ptr<SeqRegions>& parent_new_lower_lh,
                             const cmaple::RealNumType& child_2_new_blength,
                             PhyloNode& current_node,
                             PhyloNode& child_1,
                             PhyloNode& child_2,
                             PhyloNode& sibling,
                             PhyloNode& parent,
                             const cmaple::Index parent_index,
                             cmaple::RealNumType& lh_at_root,
                             const bool allow_replacing_ML_tree);

  /**
   Replace the current ML Tree by an NNI neighbor on a branch connecting to root
   @throw std::logic_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void replaceMLTreebyNNIRoot(std::stack<cmaple::Index>& node_stack_aLRT,
                              cmaple::RealNumType& lh_diff,
                              PhyloNode& current_node,
                              PhyloNode& child_1,
                              PhyloNode& child_2,
                              PhyloNode& sibling,
                              PhyloNode& parent,
                              cmaple::RealNumType& lh_at_root,
                              const cmaple::RealNumType child_1_best_blength,
                              const cmaple::RealNumType child_2_best_blength,
                              const cmaple::RealNumType sibling_best_blength,
                              const cmaple::RealNumType parent_best_blength);

  /**
   Replace the current ML Tree by an NNI neighbor on a branch connecting to a
   non-root node
   @throw std::logic_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void replaceMLTreebyNNINonRoot(
      std::stack<cmaple::Index>& node_stack_aLRT,
      cmaple::RealNumType& lh_diff,
      PhyloNode& current_node,
      PhyloNode& child_1,
      PhyloNode& child_2,
      PhyloNode& sibling,
      PhyloNode& parent,
      cmaple::RealNumType& lh_at_root,
      const cmaple::RealNumType child_1_best_blength,
      const cmaple::RealNumType child_2_best_blength,
      const cmaple::RealNumType sibling_best_blength,
      const cmaple::RealNumType parent_best_blength,
      const cmaple::RealNumType new_parent_best_blength);

  /**
   Traverse downward to update the upper_left/right_region until the changes is
   insignificant

   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void updateUpperLR(std::stack<cmaple::Index>& node_stack,
                     std::stack<cmaple::Index>& node_stack_aLRT);

  /**
   Add grand-children of a node to node_stack_aLRT to recompute their aLRT after
   replacing the ML tree with an NNI neighbor
   */
  void recompute_aLRT_GrandChildren(PhyloNode& parent,
                                    std::stack<cmaple::Index>& node_stack_aLRT);

  /**
   Calculate aLRT for each internal branches
   @param[in] allow_replacing_ML_tree TRUE to allow replacing the ML tree by a
   higher likelihood tree found when computing branch supports
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void calculate_aRLT(const bool allow_replacing_ML_tree);

  /**
   Perform a DFS to calculate the Site-lh-contribution
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  cmaple::RealNumType calculateSiteLhs(
      std::vector<cmaple::RealNumType>& site_lh_contributions,
      std::vector<cmaple::RealNumType>& site_lh_root);

  /**
   Calculate aLRT-SH for each internal branches
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void calculate_aRLT_SH(
      std::vector<cmaple::RealNumType>& site_lh_contributions,
      std::vector<cmaple::RealNumType>& site_lh_root,
      const cmaple::RealNumType& LT1);

  /**
   Count aLRT-SH for an internal branch
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  cmaple::PositionType count_aRLT_SH_branch(
      std::vector<cmaple::RealNumType>& site_lh_contributions,
      std::vector<cmaple::RealNumType>& site_lh_root,
      PhyloNode& node,
      const cmaple::RealNumType& LT1);

  /**
   Calculate the site-lh differences  between an NNI neighbor on the branch
   connecting to root and the ML tree
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void calSiteLhDiffRoot(std::vector<cmaple::RealNumType>& site_lh_diff,
                         std::vector<cmaple::RealNumType>& site_lh_root_diff,
                         const std::vector<cmaple::RealNumType>& site_lh_root,
                         std::unique_ptr<SeqRegions>& parent_new_lower_lh,
                         const cmaple::RealNumType& child_2_new_blength,
                         PhyloNode& current_node,
                         PhyloNode& child_1,
                         PhyloNode& child_2,
                         PhyloNode& sibling,
                         PhyloNode& parent,
                         const cmaple::Index parent_index);

  /**
   Calculate the site-lh differences  between an NNI neighbor on the branch
   connecting to a non-root node and the ML tree
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void calSiteLhDiffNonRoot(
      std::vector<cmaple::RealNumType>& site_lh_diff,
      std::vector<cmaple::RealNumType>& site_lh_root_diff,
      const std::vector<cmaple::RealNumType>& site_lh_root,
      std::unique_ptr<SeqRegions>& parent_new_lower_lh,
      const cmaple::RealNumType& child_2_new_blength,
      PhyloNode& current_node,
      PhyloNode& child_1,
      PhyloNode& child_2,
      PhyloNode& sibling,
      PhyloNode& parent,
      const cmaple::Index parent_index);

  /**
   Calculate the site-lh differences  between an NNI neighbor and the ML tree
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void calSiteLhDiff(std::vector<cmaple::RealNumType>& site_lh_diff,
                     std::vector<cmaple::RealNumType>& site_lh_root_diff,
                     const std::vector<cmaple::RealNumType>& site_lh_root,
                     PhyloNode& current_node,
                     PhyloNode& child_1,
                     PhyloNode& child_2,
                     PhyloNode& sibling,
                     PhyloNode& parent,
                     const cmaple::Index parent_index);

  /**
   Read the next character from the treefile
   */
  const char readNextChar(std::istream& in,
                          cmaple::PositionType& in_line,
                          cmaple::PositionType& in_column,
                          std::string& in_comment,
                          const char& current_ch = 0) const;

  /**
   Read string from tree file to create new nodes
   @throw std::logic\_error if any of the following situations occur.
   - any taxa in the tree is not found in the  alignment
   - unexpected values/behaviors found during the operations
   */
  cmaple::NumSeqsType parseFile(
      std::istream& infile,
      char& ch,
      cmaple::RealNumType& branch_len,
      cmaple::PositionType& in_line,
      cmaple::PositionType& in_column,
      std::string& in_comment,
      const std::map<std::string, cmaple::NumSeqsType>& map_seqname_index,
      bool& missing_blengths);

  /**
   Collapse leaves with zero-branch-lengths into a vector of less-info-seqs of a
   leaf
   */
  void collapseAllZeroLeave();

  /**
   Collapse one zero-branch-length leaf into its sibling's vector of
   less-info-seqs
   */
  void collapseOneZeroLeaf(PhyloNode& node,
                           cmaple::Index& node_index,
                           PhyloNode& neighbor_1,
                           const cmaple::Index neighbor_1_index,
                           PhyloNode& neighbor_2);

  /**
   Update the pesudocount of the model based on the sequence of a leaf
   */
  template <const cmaple::StateType num_states>
  void updatePesudoCountModel(PhyloNode& node,
                              const cmaple::Index node_index,
                              const cmaple::Index parent_index);
    
    /**
     Compute the number of descendants of a node
     */
    void computeNumDescendantsOfNode(PhyloNode& node,
                                const cmaple::Index node_index,
                                const cmaple::Index parent_index);
    
    /**
     Compute the number of descendants of all nodes
     */
    void computeNumDescendantsTree();


  /**
   Expand the tree by placing one less-info-seq
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void expandTreeByOneLessInfoSeq(PhyloNode& node,
                                  const cmaple::Index node_index,
                                  const cmaple::Index parent_index);

  /**
   Carefully update blength of a node when replacing the ML tree by an NNI
   neighbor Expand the new tree by adding one less-info -seq of the current node
   (if neccessary) to make sure we compute aLRT for all non-zero internal
   branches
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void updateBlengthReplaceMLTree(std::stack<cmaple::Index>& node_stack_aLRT,
                                  cmaple::RealNumType& lh_diff,
                                  PhyloNode& node,
                                  const cmaple::Index node_index,
                                  const cmaple::RealNumType best_blength);

  /**
   Expand the new tree by adding one less-info -seq of the current node after
   replacing the ML tree by an NNI neighbor to make sure we compute aLRT for all
   non-zero internal branches
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void addLessInfoSeqReplacingMLTree(std::stack<cmaple::Index>& node_stack_aLRT,
                                     cmaple::RealNumType& lh_diff,
                                     PhyloNode& node,
                                     const cmaple::Index node_index,
                                     const cmaple::Index parent_index);

  /**
   Export Node in Newick format
   @throw std::invalid\_argument if show\_branch\_supports = true but branch
   support values have yet been computed
   */
  std::string exportNodeString(const bool is_newick_format,
                               const bool binary,
                               const cmaple::NumSeqsType node_vec_index,
                               const bool print_internal_id,
                               const bool show_branch_supports);
    
    /**
     Export string of an alternative branch (for SPRTA)
     */
    std::string exportStringAltBranch(const AltBranch& alt_branch);

  /**
   Read an input tree from a stream
   @return TRUE if the tree contains any branch without a length
   @throw std::invalid\_argument if the tree in an incorrect format
   @throw std::logic\_error if any of the following situations occur.
   - any taxa in the tree is not found in the alignment
   - unexpected values/behaviors found during the operations

   @throw std::bad\_alloc if failing to allocate memory to store the tree
   */
  bool readTree(std::istream& tree_stream, PositionType& in_line);
    
    /**
     Read an input tree (in Nexus format) from a stream
     @return TRUE if the tree contains any branch without a length
     @throw std::invalid\_argument if the tree in an incorrect format
     @throw std::logic\_error if any of the following situations occur.
     - any taxa in the tree is not found in the alignment
     - unexpected values/behaviors found during the operations

     @throw std::bad\_alloc if failing to allocate memory to store the tree
     */
    bool readNexusTree(std::istream& tree_stream, PositionType& in_line);

  /**
   Check if the current tree is complete (i.e., containing all sequences from
   the alignment)
   @return TRUE if the tree is complete
   */
  bool isComplete();

  /**
   Update model according to alignment data
   @throw std::logic\_error if the reference genome is empty
   */
  void updateModelByAln();

  /**
   Update model params and partial likelihoods after loading a tree
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void updateModelLhAfterLoading();

  /**
   Initialize a mapping between sequence names and their index in the alignment
   */
  std::map<std::string, NumSeqsType> initMapSeqNameIndex();

  /**
   * Re-mark the sequences in the alignment, which already existed in the
   * current tree
   * @throw std::logic\_error if any taxa in the current tree is not found in
   * the new alignment
   */
  void remarkExistingSeqs();

  /**
   * Find and mark a sequence existed in the tree
   * @param[in] seq_name Name of the sequence
   * @param[in] map_name_index a mapping between sequence names and its index in
   * the alignment
   * @return the new index of the corresponding sequence in the new alignment
   * @throw std::logic\_error if the taxon named seq\_name is not found the
   * alignment
   */
  NumSeqsType markAnExistingSeq(
      const std::string& seq_name,
      const std::map<std::string, NumSeqsType>& map_name_index);

  /**
   * Mark all sequences (in the alignment) as not yet added to the current tree
   */
  void resetSeqAdded();

  /**
   Attach alignment and model
   @throw std::invalid\_argument If the sequence type is unsupported (neither
   DNA (for nucleotide data) nor AA (for protein data))
   @throw std::logic\_error if the reference genome is empty
   */
  void attachAlnModel(Alignment* aln, ModelBase* model);

  /**
   Export tree std::string in Newick format
   */
  std::string exportNewick(const bool binary, const bool print_internal_id,
                           const bool show_branch_supports);
    
    /**
     Export tree std::string in NEXUS format
     */
    std::string exportNexus(const bool binary,
                             const bool show_branch_supports);
    
    /**
     Traverse the tree to export TSV content
     */
    std::string exportTsvContent();

  /**
   Increase the length of a 0-length branch (connecting this node to its parent)
   to resolve the inconsistency when updating regions in updatePartialLh()
   */
  template <const cmaple::StateType num_states>
  void updateZeroBlength(const cmaple::Index index,
                         PhyloNode& node,
                         std::stack<cmaple::Index>& node_stack);

  /**
   Iteratively update partial_lh starting from the nodes in node_stack

   @param node_stack stack of nodes;
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void updatePartialLh(std::stack<cmaple::Index>& node_stack);

  /**
   Seek a position for a sample placement starting at the start_node

   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void seekSamplePlacement(const cmaple::Index start_node_index,
                           const cmaple::NumSeqsType seq_name_index,
                           const std::unique_ptr<SeqRegions>& sample_regions,
                           cmaple::Index& selected_node_index,
                           cmaple::RealNumType& best_lh_diff,
                           bool& is_mid_branch,
                           cmaple::RealNumType& best_up_lh_diff,
                           cmaple::RealNumType& best_down_lh_diff,
                           cmaple::Index& best_child_index);

  /**
   Seek a position for placing a subtree/sample starting at the start_node

   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void seekSubTreePlacement(
      cmaple::Index& best_node_index,
      cmaple::RealNumType& best_lh_diff,
      bool& is_mid_branch,
      cmaple::RealNumType& best_up_lh_diff,
      cmaple::RealNumType& best_down_lh_diff,
      cmaple::Index& best_child_index,
      const bool short_range_search,
      const cmaple::Index child_node_index,
      cmaple::RealNumType& removed_blength,
      cmaple::RealNumType& opt_appending_blength,
      cmaple::RealNumType& opt_mid_top_blength,
      cmaple::RealNumType& opt_mid_bottom_blength);
    
    /**
     Seek the best root position
     @throw std::logic\_error if unexpected values/behaviors found during the
     operations
     */
    template <const cmaple::StateType num_states>
    NumSeqsType seekBestRoot();
    
    /**
     Compute the root supports
     @throw std::logic\_error if unexpected values/behaviors found during the
     operations
     */
    void computeRootSupports(const NumSeqsType& best_node_vec_index,
                             const RealNumType& best_lh_diff,
                             std::vector<AltBranch>& alt_roots);

  /**
   Place a new sample at a mid-branch point
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void placeNewSampleMidBranch(const cmaple::Index& selected_node_index,
                               std::unique_ptr<SeqRegions>& sample,
                               const cmaple::NumSeqsType seq_name_index,
                               const cmaple::RealNumType best_lh_diff);

  /**
   Place a new sample as a descendant of a node
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void placeNewSampleAtNode(const cmaple::Index selected_node_index,
                            std::unique_ptr<SeqRegions>& sample,
                            const cmaple::NumSeqsType seq_name_index,
                            const cmaple::RealNumType best_lh_diff,
                            const cmaple::RealNumType best_up_lh_diff,
                            const cmaple::RealNumType best_down_lh_diff,
                            const cmaple::Index best_child_index);

  /**
   Apply a single SPR move
   pruning a subtree then regrafting it to a new position
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void applyOneSPR(const cmaple::Index subtree_index,
                   PhyloNode& subtree,
                   const cmaple::Index best_node_index,
                   const bool is_mid_branch,
                   const cmaple::RealNumType branch_length,
                   const cmaple::RealNumType opt_appending_blength,
                   const cmaple::RealNumType opt_mid_top_blength,
                   const cmaple::RealNumType opt_mid_bottom_blength,
                   const cmaple::RealNumType best_lh_diff);

  /**
   Traverse the intial tree from root to re-calculate all likelihoods regarding
   the latest/final estimated model parameters
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void refreshAllLhs(bool avoid_using_upper_lr_lhs = false);

  /**
   Reset the SPR flags
   @param update_outdated TRUE to update outdated
   @param n_outdated the new value of outdated
   */
  void resetSPRFlags(const bool update_outdated,
                     const bool n_outdated);

  /**
   Try to improve the entire tree with SPR moves
   @return total improvement
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  cmaple::RealNumType improveEntireTree(const TreeSearchType tree_search_type,
                                        bool short_range_search);
    
    /**
     The parallel version of improveEntireTree() - Try to improve the entire tree with SPR moves using multithreading
     @return total improvement
     @throw std::logic\_error if unexpected values/behaviors found during the
     operations
     */
    template <const cmaple::StateType num_states>
    cmaple::RealNumType improveEntireTreeParallel(const TreeSearchType tree_search_type,
                                          bool short_range_search);

  /**
   Try to optimize branch lengths of the tree by one round of tree traversal
   @return num of improvements
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  cmaple::PositionType optimizeBranchIter();

  /**
   Estimate the length of a branch using the derivative of the likelihood cost
   function wrt the branch length
   */
  template <const cmaple::StateType num_states>
  cmaple::RealNumType estimateBranchLength(
      const std::unique_ptr<SeqRegions>& parent_regions,
      const std::unique_ptr<SeqRegions>& child_regions);

  /**
   Estimate the length of a branch and check whether the new branch is different
   from the current one
   */
  template <const cmaple::StateType num_states>
  cmaple::RealNumType estimateBranchLengthWithCheck(
      const std::unique_ptr<SeqRegions>& upper_lr_regions,
      const std::unique_ptr<SeqRegions>& lower_regions,
      const cmaple::RealNumType current_blength);

  /**
   Calculate the placement cost of a sample
   @param child_regions: vector of regions of the new sample
   */
  template <const cmaple::StateType num_states>
  cmaple::RealNumType calculateSamplePlacementCost(
      const std::unique_ptr<SeqRegions>& parent_regions,
      const std::unique_ptr<SeqRegions>& child_regions,
      const cmaple::RealNumType blength);

  /**
   Calculate the placement cost of a subtree
   @param child_regions: vector of regions of the new sample
   */
  template <const cmaple::StateType num_states>
  cmaple::RealNumType calculateSubTreePlacementCost(
      const std::unique_ptr<SeqRegions>& parent_regions,
      const std::unique_ptr<SeqRegions>& child_regions,
      const cmaple::RealNumType blength);

  /**
   Update lower lh of a node
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void updateLowerLh(cmaple::RealNumType& total_lh,
                     std::unique_ptr<SeqRegions>& new_lower_lh,
                     PhyloNode& node,
                     const std::unique_ptr<SeqRegions>& lower_lh_1,
                     const std::unique_ptr<SeqRegions>& lower_lh_2,
                     const cmaple::Index neighbor_1_index,
                     PhyloNode& neighbor_1,
                     const cmaple::Index neighbor_2_index,
                     PhyloNode& neighbor_2,
                     const cmaple::PositionType& seq_length);

  /**
   Update lower lh of a node but avoid using UpperLeft/Right lhs to update
   zero-blength This function is called after reading a tree from an input file,
   thus, UpperLeft/Right lhs have not yet been computed
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void updateLowerLhAvoidUsingUpperLRLh(
      cmaple::RealNumType& total_lh,
      std::unique_ptr<SeqRegions>& new_lower_lh,
      PhyloNode& node,
      const std::unique_ptr<SeqRegions>& lower_lh_1,
      const std::unique_ptr<SeqRegions>& lower_lh_2,
      const cmaple::Index neighbor_1_index,
      PhyloNode& neighbor_1,
      const cmaple::Index neighbor_2_index,
      PhyloNode& neighbor_2,
      const cmaple::PositionType& seq_length);

  /**
   compute the likelihood contribution of (the upper branch of) a node
   @throw std::logic\_error if unexpected values/behaviors found during the
   operations
   */
  template <const cmaple::StateType num_states>
  void computeLhContribution(cmaple::RealNumType& total_lh,
                             std::unique_ptr<SeqRegions>& new_lower_lh,
                             PhyloNode& node,
                             const std::unique_ptr<SeqRegions>& lower_lh_1,
                             const std::unique_ptr<SeqRegions>& lower_lh_2,
                             const cmaple::Index neighbor_1_index,
                             PhyloNode& neighbor_1,
                             const cmaple::Index neighbor_2_index,
                             PhyloNode& neighbor_2,
                             const cmaple::PositionType& seq_length);

  /**
   Employ Depth First Search to do a task at internal nodes
   */
  template <void (Tree::*task)(cmaple::RealNumType&,
                               std::unique_ptr<SeqRegions>&,
                               PhyloNode&,
                               const std::unique_ptr<SeqRegions>&,
                               const std::unique_ptr<SeqRegions>&,
                               const cmaple::Index,
                               PhyloNode&,
                               const cmaple::Index,
                               PhyloNode&,
                               const cmaple::PositionType&)>
  cmaple::RealNumType performDFS();
    
    /**
     Employ Depth First Search to do a task at internal nodes
     The template is different from performDFS
     */
    template <void (Tree::*task)(PhyloNode&,
            const cmaple::Index, const cmaple::Index)>
    void performDFSv2();

  /**
   Update model parameters from an alignment and a tree
   @throw std::logic\_error if the reference genome is empty
   */
  template <const cmaple::StateType num_states>
  void updateModelParams();

  /**
   Employ Depth First Search to do a task at leaves
   */
  template <
      void (Tree::*task)(PhyloNode&, const cmaple::Index, const cmaple::Index)>
  void performDFSAtLeave();
   
  /**
   Check if we should traverse further down
  */
  bool keepTraversing(const RealNumType& best_lh_diff,
                        const RealNumType& lh_diff_at_node,
                        const bool& strict_stop_seeking_placement_subtree,
                        const short int& failure_count,
                        const int& failure_limit_subtree,
                        const RealNumType& thresh_log_lh_subtree,
                              const bool able2traverse = true);
    
    /**
        Generate internal names
    */
    void genIntNames();

  // NHANLT: Debug aLRT
  // void log_current(std::stack<cmaple::Index>& node_stack_aLRT);
};

/*!
 *  \addtogroup cmaple
 *  @{
 */

/** \brief Customized << operator to output the tree string in a (bifurcating)
 * NEWICK format to a stream
 */
std::ostream& operator<<(std::ostream& out_stream, cmaple::Tree& tree);

/** \brief Customized >> operator to read a tree from a stream
 */
std::istream& operator>>(std::istream& in_stream, cmaple::Tree& tree);

/*! @} End of Doxygen Groups*/

/*! \cond PRIVATE */
template <const StateType num_states>
void cmaple::Tree::refreshAllLhs(bool avoid_using_upper_lr_lhs) {
  assert(aln);
  assert(model);
  assert(cumulative_rate);
    
  // 1. update all the lower lhs along the tree
  if (avoid_using_upper_lr_lhs) {
    performDFS<&cmaple::Tree::updateLowerLhAvoidUsingUpperLRLh<num_states>>();
  } else {
    performDFS<&cmaple::Tree::updateLowerLh<num_states>>();
  }

  // 2. update all the non-lower lhs along the tree
  refreshAllNonLowerLhs<num_states>();
}

template <const StateType num_states>
RealNumType cmaple::Tree::improveEntireTree(const TreeSearchType tree_search_type,
                                            bool short_range_search) {
  assert(aln);
  assert(model);
  assert(cumulative_rate);
  assert(nodes.size() > 0);
    
  // start from the root
  std::stack<Index> node_stack;
  node_stack.push(Index(root_vector_index, TOP));

  // dummy variables
  RealNumType total_improvement = 0;
  PositionType num_nodes = 0;
  PositionType count_node_1K = 0;

  // traverse downward the tree
  while (!node_stack.empty()) {
    // pick the top node from the stack
    Index index = node_stack.top();
    node_stack.pop();
    PhyloNode& node = nodes[index.getVectorIndex()];
    // MiniIndex mini_index = index.getMiniIndex();

    // add all children of the current nodes to the stack for further traversing
    // later
    /*for (Index neighbor_index:node.getNeighborIndexes(mini_index))
        node_stack.push(neighbor_index);*/
    assert(index.getMiniIndex() == TOP);
    if (node.isInternal()) {
      node_stack.push(node.getNeighborIndex(RIGHT));
      node_stack.push(node.getNeighborIndex(LEFT));
    }

    // only process outdated node to avoid traversing the same part of the tree
    // multiple times
    if (node.isOutdated() && node.getSPRCount() <= 5) {
      node.setOutdated(false);

      // if checkEachSPR:
      //   root=node
      //  while root.up!=None:
      //      root=root.up
      //  #print("Pre-SPR tree: "+createBinaryNewick(root))
      //  oldTreeLK=calculateTreeLikelihood(root,mutMatrix,checkCorrectness=True)
      //  #print("Pre-SPR tree likelihood: "+str(oldTreeLK))
      // reCalculateAllGenomeLists(root,mutMatrix, checkExistingAreCorrect=True)

      // do SPR moves to improve the tree
      std::vector<std::pair<cmaple::Index, double>> SPR_found_vec;
      RealNumType improvement =
          improveSubTree<num_states>(index, node,
                                     tree_search_type, short_range_search, SPR_found_vec);

      // if checkEachSPR:
      //          #print(" apparent improvement "+str(improvement))
      //        root=node
      //       while root.up!=None:
      //         root=root.up
      //           #print("Post-SPR tree: "+createBinaryNewick(root))
      //         newTreeLK=calculateTreeLikelihood(root,mutMatrix)
      //       reCalculateAllGenomeLists(root,mutMatrix,
      //       checkExistingAreCorrect=True)
      //     #print("Post-SPR tree likelihood: "+str(newTreeLK))
      //   if newTreeLK-oldTreeLK < improvement-1.0:
      //              print("In startTopologyUpdates, LK score of improvement
      //              "+str(newTreeLK)+" - "+str(oldTreeLK)+" =
      //              "+str(newTreeLK-oldTreeLK)+" is less than //what is
      //              supposd to be "+str(improvement))
      //        exit()
      //

      // update total_improvement
      total_improvement += improvement;

      // NHANLT: LOGS FOR DEBUGGING
      /*if (params->debug && improvement > 0)
          std::cout << num_nodes << ": " << std::setprecision(20) <<
         total_improvement << std::endl;*/

      // Show log every 1000 nodes
      ++num_nodes;
      if (cmaple::verbose_mode >= cmaple::VB_MED && num_nodes - count_node_1K >= 1000
          && tree_search_type != FAST_TREE_SEARCH) {
        std::cout << "Processed topology for " << convertIntToString(num_nodes)
             << " nodes." << std::endl;
        count_node_1K = num_nodes;
      }
    }
  }

  return total_improvement;
}

template <const StateType num_states>
RealNumType cmaple::Tree::improveEntireTreeParallel(const TreeSearchType tree_search_type,
                                            bool short_range_search) {
  assert(aln);
  assert(model);
  assert(cumulative_rate);
  assert(nodes.size() > 0);
    
  // generate a vector/list of (outdated) nodes to search for SPR moves.
  std::vector<Index> outdated_nodes;
  outdated_nodes.reserve(nodes.size());
    
  // start from the root
  std::stack<Index> node_stack;
  node_stack.push(Index(root_vector_index, TOP));

  // traverse downward the tree
  while (!node_stack.empty()) {
      // pick the top node from the stack
      Index index = node_stack.top();
      node_stack.pop();
      PhyloNode& node = nodes[index.getVectorIndex()];

      // add all children of the current nodes to the stack for further traversing
      // later
      assert(index.getMiniIndex() == TOP);
      if (node.isInternal()) {
        node_stack.push(node.getNeighborIndex(RIGHT));
        node_stack.push(node.getNeighborIndex(LEFT));
      }

      // only process outdated node to avoid traversing the same part of the tree
      // multiple times
      if (node.isOutdated() && node.getSPRCount() <= 5) {
          outdated_nodes.emplace_back(index);
      }
    }
  
    if (cmaple::verbose_mode >= cmaple::VB_DEBUG)
    {
        std::cout << "The number of all nodes: " << nodes.size() << std::endl;
        std::cout << "The number of outdated nodes: " << outdated_nodes.size() << std::endl;
    }
    
    // don't allow blengths change during parallel SPR search
    const bool bk_fixed_blengths = fixed_blengths;
    fixed_blengths = true;
    
    std::vector<std::pair<cmaple::Index, double>> all_SPR_found_vec;
    #pragma omp parallel
    {
        std::vector<std::pair<cmaple::Index, double>> SPR_found_vec;
        SPR_found_vec.reserve(outdated_nodes.size());
        int thread_id = 0;
        #ifdef _OPENMP
        thread_id = omp_get_thread_num();
        #endif
        
        #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < outdated_nodes.size(); ++i)
        {
            cmaple::Index& index = outdated_nodes[i];
            PhyloNode& node = nodes[index.getVectorIndex()];
            node.setOutdated(false);

            // do SPR moves to improve the tree
            improveSubTree<num_states>(index, node, tree_search_type,
                                           short_range_search, SPR_found_vec, true);
        }
        
        // Merge results safely at the end
        #pragma omp critical
        {
            all_SPR_found_vec.insert(all_SPR_found_vec.end(), SPR_found_vec.begin(), SPR_found_vec.end());
            if (cmaple::verbose_mode >= cmaple::VB_DEBUG && tree_search_type != FAST_TREE_SEARCH)
            {
                std::cout << "Thread " << thread_id << " found " << SPR_found_vec.size() << " SPR moves." << std::endl;
            }
        }
    }
    
    if (cmaple::verbose_mode >= cmaple::VB_DEBUG && tree_search_type != FAST_TREE_SEARCH)
    {
        std::cout << "All threads found " << all_SPR_found_vec.size() << " SPR moves." << std::endl;
    }
    
    // sort all SPRs found in a descending order of the lh improvement
    std::sort(all_SPR_found_vec.begin(), all_SPR_found_vec.end(),
              [](const auto &a, const auto &b) {
                  return a.second > b.second;
              });
    
    // restore the fixed_blengths
    fixed_blengths = bk_fixed_blengths;
    
    // dummy variables
    RealNumType total_improvement = 0;
    
    // sequentially search and apply SPRs on nodes found in the above step
    if (tree_search_type != FAST_TREE_SEARCH)
    {
        for (const auto &item : all_SPR_found_vec) {
            const cmaple::Index& index = item.first;
            PhyloNode& node = nodes[index.getVectorIndex()];
            
            // search and apply SPR on that node
            std::vector<std::pair<cmaple::Index, double>> SPR_found_vec;
            RealNumType improvement =
            improveSubTree<num_states>(index, node,
                                       tree_search_type, short_range_search, SPR_found_vec);
            
            // update total_improvement
            total_improvement += improvement;
        }
    }
    

  return total_improvement;
}

template <const StateType num_states>
void cmaple::Tree::updateModelParams() {
  assert(aln);
  assert(model);
    
    // reset the pseudo count
    model->initMutationMat();
    
  // perform a DFS -> at each node, update the pesudoCount of the model based on
  // the sequence of that node
  performDFSv2<&cmaple::Tree::updatePesudoCountModel<num_states>>();

  // update model params based on the pseudo count
  if (model->updateMutationMatEmpirical()) {
    computeCumulativeRate();
  }
}

template <const StateType num_states>
void cmaple::Tree::seekSamplePlacement(
    const Index start_node_index,
    const NumSeqsType seq_name_index,
    const std::unique_ptr<SeqRegions>& sample_regions,
    Index& selected_node_index,
    RealNumType& best_lh_diff,
    bool& is_mid_branch,
    RealNumType& best_up_lh_diff,
    RealNumType& best_down_lh_diff,
    Index& best_child_index) {
  assert(sample_regions && sample_regions->size() > 0);
  assert(seq_name_index >= 0);
  assert(aln);
  assert(model);
  assert(cumulative_rate);
    
  // init variables
  // output variables
  selected_node_index = start_node_index;
  // dummy variables
  const bool collapse_only_ident_seqs = params->compute_SPRTA &&
    params->compute_SPRTA_zero_length_branches;
  RealNumType lh_diff_mid_branch = 0;
  RealNumType lh_diff_at_node = 0;
  PositionType seq_length = static_cast<PositionType>(aln->ref_seq.size());
  // stack of nodes to examine positions
  std::stack<TraversingNode> extended_node_stack;
  extended_node_stack.push(TraversingNode(start_node_index, 0, MIN_NEGATIVE));

  // recursively examine positions for placing the new sample
  while (!extended_node_stack.empty()) {
    TraversingNode current_extended_node = std::move(extended_node_stack.top());
    extended_node_stack.pop();
    const NumSeqsType current_node_vec =
        current_extended_node.getIndex().getVectorIndex();
    PhyloNode& current_node = nodes[current_node_vec];
    const bool& is_internal = current_node.isInternal();

    // NHANLT: debug
    // if (current_node->next && ((current_node->next->neighbor &&
    // current_node->next->neighbor->seq_name == "25")
    //                        || (current_node->next->next->neighbor &&
    //                        current_node->next->next->neighbor->seq_name ==
    //                        "25")))
    // cout << "fdsfsd";

    // if the current node is a leaf AND the new sample/sequence is strictly
    // less informative than the current node
    // -> add the new sequence into the list of minor sequences of the current
    // node + stop seeking the placement
    if ((!is_internal) &&
        (current_node.getPartialLh(TOP)->compareWithSample(
             *sample_regions, seq_length, aln,
            collapse_only_ident_seqs) == 1)) {
#ifdef _OPENMP
#pragma omp critical
#endif
      current_node.addLessInfoSeqs(seq_name_index);
      selected_node_index = Index();
      return;
    }

    const RealNumType current_node_blength = current_node.getUpperLength();

    // 1. try first placing as a descendant of the mid-branch point of the
    // branch above the current node
    if (root_vector_index != current_node_vec && current_node_blength > 0) {
      examineSamplePlacementMidBranch<num_states>(
          selected_node_index, current_node.getMidBranchLh(), best_lh_diff,
          is_mid_branch, lh_diff_mid_branch, current_extended_node,
          sample_regions);
    }
    // otherwise, don't consider mid-branch point
    else {
      lh_diff_mid_branch = MIN_NEGATIVE;
    }

    // 2. try to place as descendant of the current node (this is skipped if the
    // node has top branch length 0 and so is part of a polytomy).
    if (root_vector_index == current_node_vec || current_node_blength > 0) {
      examineSamplePlacementAtNode<num_states>(
          selected_node_index, current_node.getTotalLh(), best_lh_diff,
          is_mid_branch, lh_diff_at_node, lh_diff_mid_branch, best_up_lh_diff,
          best_down_lh_diff, best_child_index, current_extended_node,
          sample_regions);
    } else {
      lh_diff_at_node = current_extended_node.getLhDiff();
    }

    // keep trying to place at children nodes, unless the number of attempts has
    // reaches the failure limit
    const short int failure_count = current_extended_node.getFailureCount();
    if ((params->strict_stop_seeking_placement_sample &&
         failure_count <= params->failure_limit_sample &&
         lh_diff_at_node > (best_lh_diff - params->thresh_log_lh_sample)) ||
        (!params->strict_stop_seeking_placement_sample &&
         (failure_count <= params->failure_limit_sample ||
          lh_diff_at_node > (best_lh_diff - params->thresh_log_lh_sample)))) {
      /*for (Index neighbor_index:current_node.getNeighborIndexes(TOP))
          extended_node_stack.push(TraversingNode(neighbor_index,
         current_extended_node.getFailureCount(), lh_diff_at_node));*/
      if (is_internal) {
        extended_node_stack.push(
            TraversingNode(current_node.getNeighborIndex(RIGHT), failure_count,
                           lh_diff_at_node));
        extended_node_stack.push(
            TraversingNode(current_node.getNeighborIndex(LEFT), failure_count,
                           lh_diff_at_node));
      }
    }
  }

  // exploration of the tree is finished, and we are left with the node found so
  // far with the best appending likelihood cost. Now we explore placement just
  // below this node for more fine-grained placement within its descendant
  // branches.
  best_down_lh_diff = MIN_NEGATIVE;
  best_child_index = Index();

  // if best position so far is the descendant of a node -> explore further at
  // its children
  if (!is_mid_branch) {
    finetuneSamplePlacementAtNode<num_states>(
        nodes[selected_node_index.getVectorIndex()], best_down_lh_diff,
        best_child_index, sample_regions);
  }
}

template <const StateType num_states>
void cmaple::Tree::placeNewSampleAtNode(const Index selected_node_index,
                                        std::unique_ptr<SeqRegions>& sample,
                                        const NumSeqsType seq_name_index,
                                        const RealNumType best_lh_diff,
                                        const RealNumType best_up_lh_diff,
                                        const RealNumType best_down_lh_diff,
                                        const Index best_child_index) {
  // dummy variables
  RealNumType best_child_lh = MIN_NEGATIVE;
  RealNumType best_child_blength_split = 0;
  RealNumType best_parent_lh;
  RealNumType best_parent_blength_split = 0;
  RealNumType best_root_blength = -1;
  const RealNumType threshold_prob = params->threshold_prob;
  std::unique_ptr<SeqRegions> best_parent_regions = nullptr;
  std::unique_ptr<SeqRegions> best_child_regions = nullptr;

  assert(selected_node_index.getMiniIndex() == TOP);
  assert(sample && sample->size() > 0);
  assert(seq_name_index >= 0);
  assert(aln);
  assert(model);
  assert(cumulative_rate);
    
  const NumSeqsType selected_node_vec_index =
      selected_node_index.getVectorIndex();
  PhyloNode& selected_node = nodes[selected_node_vec_index];

  // place the new sample as a descendant of an existing node
  if (best_child_index.getMiniIndex() != UNDEFINED) {
    PhyloNode& best_child = nodes[best_child_index.getVectorIndex()];
    best_child_lh = best_down_lh_diff;
    assert(best_child_index.getMiniIndex() == TOP);
    best_child_blength_split =
        0.5 * best_child.getUpperLength();  // best_child->length;
    const std::unique_ptr<SeqRegions>& upper_left_right_regions =
        getPartialLhAtNode(best_child.getNeighborIndex(
            TOP));  // best_child->neighbor->getPartialLhAtNode(aln,
                    // model, threshold_prob);
    const std::unique_ptr<SeqRegions>& lower_regions = best_child.getPartialLh(
        TOP);  // ->getPartialLhAtNode(aln, model, threshold_prob);
    // best_child_regions = new SeqRegions(best_child->mid_branch_lh);
    // SeqRegions best_child_mid_clone =
    // SeqRegions(best_child.getMidBranchLh());
    best_child_regions =
        nullptr;  // cmaple::make_unique<SeqRegions>(std::move(best_child_mid_clone));

    // try a shorter split
    // tryShorterBranch<&cmaple::Tree::calculateSamplePlacementCost<num_states>>(best_child->length,
    // best_child_regions, sample, upper_left_right_regions, lower_regions,
    // best_child_lh, best_child_blength_split, default_blength, true);
    tryShorterBranch<num_states,
                     &cmaple::Tree::calculateSamplePlacementCost<num_states>>(
        best_child.getUpperLength(), best_child_regions, sample,
        upper_left_right_regions, lower_regions, best_child_lh,
        best_child_blength_split, default_blength, true);

    // Delay cloning SeqRegions
    if (!best_child_regions) {
      best_child_regions =
          cmaple::make_unique<SeqRegions>(SeqRegions(best_child.getMidBranchLh()));
    }
  }

  // if node is root, try to place as sibling of the current root.
  RealNumType old_root_lh = MIN_NEGATIVE;
  if (root_vector_index == selected_node_vec_index) {
    /*old_root_lh = selected_node->getPartialLhAtNode(aln, model,
    threshold_prob)->computeAbsoluteLhAtRoot(num_states, model); SeqRegions*
    lower_regions = selected_node->getPartialLhAtNode(aln, model,
    threshold_prob);*/
    const std::unique_ptr<SeqRegions>& lower_regions =
        selected_node.getPartialLh(TOP);
    old_root_lh = lower_regions->computeAbsoluteLhAtRoot<num_states>(
        model, cumulative_base);

    // merge 2 lower vector into one
    RealNumType new_root_lh = lower_regions->mergeTwoLowers<num_states>(
        best_parent_regions, default_blength, *sample, default_blength, aln,
        model, cumulative_rate, threshold_prob, true);

    new_root_lh += best_parent_regions->computeAbsoluteLhAtRoot<num_states>(
        model, cumulative_base);
    best_parent_lh = new_root_lh;

    // try shorter branch lengths
    best_root_blength = default_blength;
    tryShorterBranchAtRoot<num_states>(sample, lower_regions,
                                       best_parent_regions, best_root_blength,
                                       best_parent_lh, default_blength);

    // update best_parent_lh (taking into account old_root_lh)
    best_parent_lh -= old_root_lh;
  }
  // selected_node is not root
  else {
    best_parent_lh = best_up_lh_diff;
    best_parent_blength_split =
        0.5 * selected_node.getUpperLength();  // selected_node->length;
    /*SeqRegions* upper_left_right_regions =
    selected_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
    SeqRegions* lower_regions = selected_node->getPartialLhAtNode(aln, model,
    threshold_prob); best_parent_regions = new
    SeqRegions(selected_node->mid_branch_lh);*/
    const std::unique_ptr<SeqRegions>& upper_left_right_regions =
        getPartialLhAtNode(selected_node.getNeighborIndex(TOP));
    const std::unique_ptr<SeqRegions>& lower_regions =
        selected_node.getPartialLh(TOP);
    // SeqRegions seq_regions_clone =
    // SeqRegions(selected_node.getMidBranchLh());
    best_parent_regions =
        nullptr;  // cmaple::make_unique<SeqRegions>(std::move(seq_regions_clone));

    // try a shorter split
    // tryShorterBranch<&cmaple::Tree::calculateSamplePlacementCost<num_states>>(selected_node->length,
    // best_parent_regions, sample, upper_left_right_regions, lower_regions,
    // best_parent_lh, best_parent_blength_split, default_blength, false);
    tryShorterBranch<num_states,
                     &cmaple::Tree::calculateSamplePlacementCost<num_states>>(
        selected_node.getUpperLength(), best_parent_regions, sample,
        upper_left_right_regions, lower_regions, best_parent_lh,
        best_parent_blength_split, default_blength, false);

    // Delay cloning SeqRegions
    if (!best_parent_regions) {
      best_parent_regions = cmaple::make_unique<SeqRegions>(
          SeqRegions(selected_node.getMidBranchLh()));
    }
  }

  // if the best placement is below the selected_node => add an internal node
  // below the selected_node
  if (best_child_lh >= best_parent_lh && best_child_lh >= best_lh_diff) {
    assert(best_child_index.getMiniIndex() == TOP);
    PhyloNode& best_child = nodes[best_child_index.getVectorIndex()];
    const std::unique_ptr<SeqRegions>& upper_left_right_regions =
        getPartialLhAtNode(best_child.getNeighborIndex(
            TOP));  // best_child->neighbor->getPartialLhAtNode(aln,
                    // model, threshold_prob);

    // Estimate the length for the new branch
    RealNumType best_length = default_blength;
    estimateLengthNewBranch<
        &cmaple::Tree::calculateSamplePlacementCost<num_states>>(
        best_child_lh, best_child_regions, sample, best_length, max_blength,
        min_blength, false);

    // create new internal node and append child to it
    // connectNewSample2Branch(sample, seq_name_index, best_child,
    // best_child_blength_split, best_child->length - best_child_blength_split,
    // best_length, best_child_regions, upper_left_right_regions);
    connectNewSample2Branch<num_states>(
        sample, seq_name_index, best_child_index, best_child,
        best_child_blength_split,
        best_child.getUpperLength() - best_child_blength_split, best_length,
        best_child_regions, upper_left_right_regions);
  }
  // otherwise, add new parent to the selected_node
  else {
    // new parent is actually part of a polytomy since best placement is exactly
    // at the node
    if (best_lh_diff >= best_parent_lh) {
      best_root_blength = -1;
      best_parent_blength_split = -1;
      best_parent_lh = best_lh_diff;
      /*if (best_parent_regions) delete best_parent_regions;
      best_parent_regions = NULL;*/
      best_parent_regions = nullptr;

      if (root_vector_index == selected_node_vec_index) {
        // selected_node->getPartialLhAtNode(aln, model,
        // threshold_prob)->mergeTwoLowers<num_states>(best_parent_regions,
        // -1, *sample, default_blength, aln, model, threshold_prob);
        selected_node.getPartialLh(TOP)->mergeTwoLowers<num_states>(
            best_parent_regions, -1, *sample, default_blength, aln, model,
            cumulative_rate, threshold_prob);
      } else {
        // best_parent_regions = new SeqRegions(selected_node->total_lh);
        SeqRegions seq_regions_clone = SeqRegions(selected_node.getTotalLh());
        best_parent_regions =
            cmaple::make_unique<SeqRegions>(std::move(seq_regions_clone));
      }
    }

    // add parent to the root
    if (root_vector_index == selected_node_vec_index) {
      // now try different lengths for right branch
      best_parent_lh += old_root_lh;
      RealNumType best_length2 = default_blength;
      const std::unique_ptr<SeqRegions>& lower_regions =
          selected_node.getPartialLh(
              TOP);  // ->getPartialLhAtNode(aln, model, threshold_prob);

      estimateLengthNewBranchAtRoot<num_states>(
          sample, lower_regions, best_parent_regions, best_length2,
          best_parent_lh, best_root_blength, min_blength, false);

      // update best_parent_lh (taking into account old_root_lh)
      best_parent_lh -= old_root_lh;

      // add new sample to a new root
      connectNewSample2Root<num_states>(
          sample, seq_name_index, selected_node_index, selected_node,
          best_root_blength, best_length2, best_parent_regions);
    }
    // add parent to non-root node
    else {
      const std::unique_ptr<SeqRegions>& upper_left_right_regions =
          getPartialLhAtNode(selected_node.getNeighborIndex(
              TOP));  // selected_node->neighbor->getPartialLhAtNode(aln,
                      // model, threshold_prob);

      // now try different lengths for the new branch
      RealNumType best_length = default_blength;
      estimateLengthNewBranch<
          &cmaple::Tree::calculateSamplePlacementCost<num_states>>(
          best_parent_lh, best_parent_regions, sample, best_length, max_blength,
          min_blength, false);

      // now create new internal node and append child to it
      RealNumType down_distance = best_parent_blength_split;
      RealNumType top_distance =
          selected_node.getUpperLength() -
          down_distance;  // selected_node->length - down_distance;
      if (best_parent_blength_split < 0) {
        down_distance = -1;
        top_distance =
            selected_node.getUpperLength();  // selected_node->length;

        /*if (selected_node->total_lh) delete selected_node->total_lh;
        selected_node->total_lh = NULL;

        if (selected_node->mid_branch_lh) delete selected_node->mid_branch_lh;
        selected_node->mid_branch_lh = NULL;*/
        selected_node.setTotalLh(nullptr);
        selected_node.setMidBranchLh(nullptr);

        // node.furtherMidNodes=None
      }
      connectNewSample2Branch<num_states>(
          sample, seq_name_index, selected_node_index, selected_node,
          top_distance, down_distance, best_length, best_parent_regions,
          upper_left_right_regions);
    }
  }

  // delete best_parent_regions and best_child_regions
  /*if (best_parent_regions)
      delete best_parent_regions;
  if (best_child_regions)
      delete best_child_regions;*/
}

template <const StateType num_states>
void cmaple::Tree::placeNewSampleMidBranch(const Index& selected_node_index,
                                           std::unique_ptr<SeqRegions>& sample,
                                           const NumSeqsType seq_name_index,
                                           const RealNumType best_lh_diff) {
  // dummy variables
  // const RealNumType threshold_prob = params->threshold_prob;
  std::unique_ptr<SeqRegions> best_child_regions = nullptr;

  // selected_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
  // const MiniIndex seleted_node_mini_index =
  // selected_node_index.getMiniIndex();
  assert(selected_node_index.getMiniIndex() == TOP);
  assert(sample && sample->size() > 0);
  assert(seq_name_index >= 0);
  assert(aln);
  assert(model);
  assert(cumulative_rate);
  
// Use new method to optimize blengths
    // optmize the branch lengths
    RealNumType best_appending_blength;
    RealNumType best_mid_top_blength;
    RealNumType best_mid_bottom_blength;
    const RealNumType threshold_prob = params->threshold_prob;
    PhyloNode& selected_node = nodes[selected_node_index.getVectorIndex()];
    
    // optimize the appending branch
    best_appending_blength =
        estimateBranchLength<num_states>(selected_node.getMidBranchLh(), sample);
        
    // optimize the mid_top blength
    const std::unique_ptr<SeqRegions>& upper_lr_regions =
        getPartialLhAtNode(selected_node.getNeighborIndex(TOP));
    const std::unique_ptr<SeqRegions>& lower_regions =
        selected_node.getPartialLh(TOP);
    const RealNumType mid_branch_length =
        selected_node.getUpperLength() * 0.5;
    std::unique_ptr<SeqRegions> two_lower_regions = nullptr;
    lower_regions->mergeTwoLowers<num_states>(two_lower_regions, mid_branch_length,
        *sample, best_appending_blength,
        aln, model, cumulative_rate, threshold_prob);
    best_mid_top_blength =
        estimateBranchLength<num_states>(upper_lr_regions, two_lower_regions);
        
    // optimize the mid_bottom blength
    std::unique_ptr<SeqRegions> tmp_upper_lr_regions = nullptr;
    upper_lr_regions->mergeUpperLower<num_states>(
        tmp_upper_lr_regions, best_mid_top_blength, *sample,
        best_appending_blength, aln, model, threshold_prob);
    best_mid_bottom_blength =
        estimateBranchLength<num_states>(tmp_upper_lr_regions, lower_regions);
        
    // compute new mid-branch regions
    std::unique_ptr<SeqRegions> new_mid_branch_regions = nullptr;
    upper_lr_regions->mergeUpperLower<num_states>(
        new_mid_branch_regions, best_mid_top_blength, *lower_regions,
        best_mid_bottom_blength, aln, model, threshold_prob);
    
    // create new internal node and append child to it
    connectNewSample2Branch<num_states>(
        sample, seq_name_index, selected_node_index, selected_node,
        best_mid_top_blength,
        best_mid_bottom_blength, best_appending_blength,
        new_mid_branch_regions, upper_lr_regions);
        
    
    /** Old method to optimize blengths -> comment out
  PhyloNode& selected_node = nodes[selected_node_index.getVectorIndex()];
  const std::unique_ptr<SeqRegions>& upper_left_right_regions =
      getPartialLhAtNode(selected_node.getNeighborIndex(TOP));
  RealNumType best_split_lh = best_lh_diff;
  const RealNumType selected_node_blength = selected_node.getUpperLength();
  RealNumType best_branch_length_split = 0.5 * selected_node_blength;
  best_child_regions =
      nullptr;  // cmaple::make_unique<SeqRegions>(SeqRegions(selected_node.getMidBranchLh()));
  // selected_node->getPartialLhAtNode(aln, model, threshold_prob);
  const std::unique_ptr<SeqRegions>& lower_regions =
      selected_node.getPartialLh(TOP);

  // try different positions on the existing branch
  bool found_new_split =
      tryShorterBranch<num_states,
                       &cmaple::Tree::calculateSamplePlacementCost<num_states>>(
          selected_node_blength, best_child_regions, sample,
          upper_left_right_regions, lower_regions, best_split_lh,
          best_branch_length_split, default_blength, true);

  if (!found_new_split) {
    // try on the second branch
    found_new_split = tryShorterBranch<
        num_states, &cmaple::Tree::calculateSamplePlacementCost<num_states>>(
        selected_node_blength, best_child_regions, sample,
        upper_left_right_regions, lower_regions, best_split_lh,
        best_branch_length_split, default_blength, false);

    if (found_new_split) {
      best_branch_length_split =
          selected_node_blength - best_branch_length_split;
    }
  }

  // Delay cloning SeqRegions
  if (!best_child_regions) {
    best_child_regions = cmaple::make_unique<SeqRegions>(
        SeqRegions(selected_node.getMidBranchLh()));
  }

  // now try different lengths for the new branch
  RealNumType best_blength = default_blength;
  estimateLengthNewBranch<
      &cmaple::Tree::calculateSamplePlacementCost<num_states>>(
      best_split_lh, best_child_regions, sample, best_blength, max_blength,
      min_blength, false);

  // create new internal node and append child to it
  connectNewSample2Branch<num_states>(
      sample, seq_name_index, selected_node_index, selected_node,
      best_branch_length_split,
      selected_node_blength - best_branch_length_split, best_blength,
      best_child_regions, upper_left_right_regions);
     */
}
/*! \endcond */
}  // namespace cmaple
