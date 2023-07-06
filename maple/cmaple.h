#include "../tree/treebase.h"
#include "../utils/timeutil.h"
#include "../alignment/seqregions.h"
#include "../tree/tree.h"

#pragma once

/** Class CMaple specifies APIs */
class CMaple {
    
    /**
        Status of the CMaple instance
     */
    enum CMapleStatus {
        NEW_INSTANCE, INFERENCE_DONE, BRANCH_SUPPORT_DONE
    };
    
public:
    
    /////////////////////////////////////////////////////////////////////////////////////////////
    /// KEY APIs
    //////////////////////////////////////////////////////////////////////////////////////////////
    // The same for set Alignment
    
    /*! \brief CMaple constructor
     *
     * Custom CMaple constructor
     * @param[in] aln_filename Name of an alignment file
     * @param[in] format Alignment format (optional): "", "MAPLE", "PHYLIP", "FASTA"
     * @param[in] seqtype Data type of sequences (optional): "", "DNA", "AA"
     */
    CMaple(const std::string& aln_filename = "", const std::string& format = "", const std::string& seqtype = "");
    
    /*! \brief CMaple constructor
     *
     * Custom CMaple constructor
     * @param[in] aln_stream Stream of an input alignment
     * @param[in] format Alignment format (optional): "", "MAPLE", "PHYLIP", "FASTA"
     * @param[in] seqtype Data type of sequences (optional): "", "DNA", "AA"
     */
    CMaple(const std::istream& aln_stream, const std::string& format = "", const std::string& seqtype = "");
    
    /*! \brief Specify an input alignment file
     *
     * Specify an input alignment in MAPLE, PHYLIP, or FASTA format
     * @param[in] aln_filename Name of an alignment file
     * @param[in] format Alignment format (optional): "", "MAPLE", "PHYLIP", "FASTA"
     * @param[in] seqtype Data type of sequences (optional): "", "DNA", "AA"
     * @return a status code: 0 - success, non-zero - failed with errors
     */
    int setAlignment(const std::string& aln_filename, const std::string& format = "", const std::string& seqtype = "");
    
    /*! \brief Specify an input alignment file
     *
     * Specify an input alignment in MAPLE, PHYLIP, or FASTA format
     * @param[in] aln_stream Stream of an input alignment
     * @param[in] format Alignment format (optional): "", "MAPLE", "PHYLIP", "FASTA"
     * @param[in] seqtype Data type of sequences (optional): "", "DNA", "AA"
     * @return a status code: 0 - success, non-zero - failed with errors
     */
    int setAlignment(const std::istream& aln_stream, const std::string& format = "", const std::string& seqtype = "");

    /*! \brief Set the subtitution model
     *
     * Set the substitution model by its name
     * @param[in] model_name Name of a substitution model
     * @return a status code: 0 - success, non-zero - failed with errors
     */
    int setModel(const std::string& model_name);
    
    /*! \brief Specify an input tree file in NEWICK format
     *
     * Specify a strictly bifurcating tree with may or may not contain all taxa from a file in NEWICK format
     * @param[in] tree_filename Name of a tree file
     * @return a status code: 0 - success, non-zero - failed with errors
     */
    int setTree(const std::string& tree_filename);
    
    /*! \brief Set the tree search type
     *
     * Set the tree search type:
     * - "FAST": no tree search (placement only).
     * - "NORMAL": only consider pruning branches at newly-added nodes when seeking SPR moves.
     * - "SLOW": consider all nodes when seeking SPR moves.
     *
     * @param[in] tree_search_type a type of tree search
     * @param[in] shallow_tree_search TRUE ton enable a shallow tree search before a deeper tree search
     * @return a status code: 0 - success, non-zero - failed with errors
     */
    int setTreeSearchType(const std::string& tree_search_type = "NORMAL", const bool shallow_tree_search = false);
    
    /*! \brief Set the prefix for the output files
     *
     * Set the prefix for the output filenames
     * @param[in] prefix A prefix for the output filenames
     * @return a status code: 0 - success, non-zero - failed with errors
     */
    int setPrefix(const std::string& prefix = "");
    
    /*! \brief Set the minimum branch length
     *
     * Set the minimum length of branches
     * @param[in] min_blength a positive minimum branch length
     * @return a status code: 0 - success, non-zero - failed with errors
     */
    int setMinBlength(const double min_blength);
    
    /*! \brief Set the approximation threshold
     *
     * Set the approximation threshold, which is used for many approximations in CMaple algorithm. Default: 1e-8
     * @param[in] thresh_prob a positive threshold
     * @return a status code: 0 - success, non-zero - failed with errors
     */
    int setThreshProb(const double thresh_prob);
    
    /*! \brief Allow CMaple to overwrite existing output files
     *
     * Allow CMaple to overwrite existing output files
     * @param[in] enable TRUE to allow overwrite the output files
     * @return a status code: 0 - success, non-zero - failed with errors
     */
    int overwriteOutputs(const bool enable);
    
    /*! \brief Specify a seed number for random generators
     *
     * Specify a seed number for random generators. Default: the clock of the PC. Be careful! To make the CMaple reproducible, users should specify the seed number.
     * @param[in] seed a seed number
     * @return a status code: 0 - success, non-zero - failed with errors
     */
    int setRandomSeed(const int seed);
    
    /*! \brief Run the inference
     *
     * Run the inference to infer a phylogenetic tree
     * @param[in] force_rerun TRUE to force rerunning the inference
     * @param[in] tree_type the type of the output tree (optional): "BIN" - bifurcating tree or "MUL" - multifurcating tree
     * @return a status code: 0 - success, non-zero - failed with errors
     */
    int inferTree(const bool force_rerun = false, const std::string& tree_type = "BIN");
    
    /*! \brief Compute branch supports
     *
     * Compute branch supports (aLRT-SH)
     * @param[in] force_rerun TRUE to force rerunning the inference
     * @param[in] num_threads number of threads (optional)
     * @param[in] num_replicates a positive number of replicates (optional)
     * @param[in] epsilon a positive epsilon (optional), // NHANLT- TODO: explain epsilon
     * @return a status code: 0 - success, non-zero - failed with errors
     */
    int computeBranchSupports(const bool force_rerun = false, const int num_threads = 1, const int num_replicates = 1000, const double epsilon = 0.05);
    
    /*! \brief Extract a MAPLE file from an alignment
     *
     * Convert an alignment in FASTA/PHYLIP format into MAPLE format. Note: the reference sequence is user-specified (via setRef) or automatically generated by CMaple
     * @param[in] aln_filename Name of an alignment file in FASTA or PHYLIP format
     * @param[in] output_filename Name of the output MAPLE file
     * @return a status code: 0 - success, non-zero - failed with errors
     */
    //int extractMaple(const std::string& aln_filename, const std::string& output_filename);
    
    /*! \brief Extract a FASTA alignment from a MAPLE file
     *
     * Convert an alignment in MAPLE into FASTA format.
     * @param[in] aln_filename Name of an alignment file in MAPLE format
     * @param[in] output_filename Name of the output FASTA file
     * @return a status code: 0 - success, non-zero - failed with errors
     */
    int extractFASTA(const std::string& aln_filename, const std::string& output_filename);
    
    /*! \brief Export the tree string in NEWICK format.
     *
     * Export the phylogenetic tree in a NEWICK string.
     * @param[in] tree_type the type of the output tree (optional): "BIN" - bifurcating tree or "MUL" - multifurcating tree
     * @param[in] show_branch_supports TRUE to output branch supports (aLRT-SH)
     * @return a tree string in NEWICK format
     */
    std::string getTreeString(const std::string& tree_type = "BIN", const bool show_branch_supports = false);
    
    /*! \brief Export the substitution model and its parameters in a dictionary
     *
     * Export the substitution model and its parameters in a dictionary
     * @return a dictionary (std::map<std::string, std::string>) with the keys and values as in the following.
     * key: "model_name", value: the name of the model in string
     * key: "model_freqs", value: the state frequencies in string
     * key: "model_rates", value: the mutation matrix in string
     */
    std::map<std::string, std::string> getModelParams();
    
    /////////////////////////////////////////////////////////////////////////////////////////////
    /// ADDITIONAL APIs
    //////////////////////////////////////////////////////////////////////////////////////////////
    
    /*! \brief Get the version of CMaple (APIs)
     *
     * Get the version of CMaple (APIs)
     * @return CMAPLE APIs version
     */
    std::string getVersion();
    
    /*! \brief Get the citation string
     *
     * Get the citation strings for CMaple/Maple manuscripts
     * @return citation strings
     */
    std::string getCitations();
    
    
    
    /*! \brief Get an instance of CMaple settings
     *
     * Get CMaple settings (in an instance of class Params). Users can use that Params instance to change other (minor) settings of CMaple (which is not yet supported via the above APIs).
     * @return an instance of cmaple::Params
     */
    cmaple::Params& getSettings();
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////
    /// PRIVATE FUNCTIONS
    //////////////////////////////////////////////////////////////////////////////////////////////
private:
    /** The phylogenetic tree */
    cmaple::TreeBase tree;
    
    /** Status of this CMaple instance */
    CMapleStatus status;
    
    /**
        Pointer  to doInference method
     */
    typedef void (CMaple::*DoInferencePtrType)();
    DoInferencePtrType doInferencePtr;
    
    /**
        Pointer  to postInference method
     */
    typedef void (CMaple::*PostInferencePtrType)();
    PostInferencePtrType postInferencePtr;
    
    /*! Template of doInference()
     */
    template <const cmaple::StateType  num_states>
    void doInferenceTemplate();
    
    /*! Template of postInference()
     */
    template <const cmaple::StateType  num_states>
    void postInferenceTemplate();
    
    /*! Setup function pointers
     */
    void setupFuncPtrs(const cmaple::StateType  num_states);
    
    /*! Build an Initial Tree
     */
    template <const cmaple::StateType  num_states>
    void buildInitialTree();
    
    /*! Optimize the current tree
     */
    template <const cmaple::StateType  num_states>
    void optimizeTree();
    
    /*! Optimize the tree topology
     */
    template <const cmaple::StateType  num_states>
    void optimizeTreeTopology(bool short_range_search = false);
    
    /*! Optimize the branch lengths of the current tree
     */
    template <const cmaple::StateType  num_states>
    void optimizeBranchLengthsOfTree();

    /*! \brief Export output files
     *
     * @param[in] filename name of the output file(s)
     * @param[in] show_branch_support TRUE to output the branch supports (aLRT-SH)
     */
    void exportOutput(const std::string &filename, const bool show_branch_support = false);
    
    /*! Calculate branch supports
     */
    template <const cmaple::StateType  num_states>
    void calculateBranchSupports();
    
    /*! Load an input tree
     */
    template <const cmaple::StateType  num_states>
    void loadInputTree();
    
    /*! Load input data
     */
    void loadInput();
    
    /*! Prepare for the inference
     */
    void preInference();
    
    /*! Do the inference
     */
    void doInference();
    
    /*! Complete the inference
     */
    void postInference();
    
    /*! \brief Reset CMaple status
     * Delete all related objects, except params
     */
    void resetStatus();
};
