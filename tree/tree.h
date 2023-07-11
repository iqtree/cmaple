#pragma once
#include "treebase.h"
#include "../alignment/alignment.h"
#include "../model/model.h"

namespace cmaple
{
    /** The tree structure */
    class Tree {
    public:
        /*! \brief Tree constructor
         *
         * Tree constructor from a (bifurcating or multifurcating) tree, which may or may not contain all taxa from a file in NEWICK format
         * @param[in] aln an alignment
         * @param[in] model a substitution model
         * @param[in] tree_stream A stream of the input tree
         * @param[in] fixed_blengths TRUE to keep the input branch lengths unchanged (optional)
         */
        Tree(Alignment& aln, Model& model, std::istream& tree_stream, const bool fixed_blengths = false);
        
        /*! \brief Tree constructor
         *
         * Tree constructor from a (bifurcating or multifurcating) tree, which may or may not contain all taxa from a file in NEWICK format
         * @param[in] aln an alignment
         * @param[in] model a substitution model
         * @param[in] tree_filename Name of a tree file (optinal)
         * @param[in] fixed_blengths TRUE to keep the input branch lengths unchanged (optional)
         */
        Tree(Alignment& aln, Model& model, const std::string& tree_filename = "", const bool fixed_blengths = false);
        
        /*! \brief Tree destructor
         *
         * Tree destructor
         */
        ~Tree();
        
        /*! \brief Infer a tree from an alignment using a substitution model
         *
         * Infer a phylogenetic tree from an alignment using a substitution model.
         * - If users didn't supply an input tree or supplied an incomplete tree (which doesn't contain all the taxa in the alignment) when initializing the tree (by Tree() constructor), this function:
         * + does placement (i.e., adding missing taxa from the alignment to the tree)
         * + applies a NORMAL tree search (which does SPR moves only on newly-added nodes)
         * + optimizes all branch lengths
         * - If users already supplied a complete tree, this function:
         * + by default, does no tree search but optimizes all branch lengths.
         * + If users want to keep the branch lengths fixed, they should set fixed_blengths = true when initializing the tree (by Tree() constructor);
         * - If users want to use the input tree as a starting tree (then does SPR moves and optimizes branch lengths), they should set tree_search_type = MORE_ACCURATE
         *
         * @param[in] tree_search_type one of the following tree search:
         * - "FAST": no tree search (placement only).
         * - "NORMAL": only consider pruning branches at newly-added nodes when seeking SPR moves.
         * - "MORE_ACCURATE": consider all nodes when seeking SPR moves.
         * @param[in] shallow_tree_search TRUE ton enable a shallow tree search before a deeper tree search
         * @Return a string contains all messages redirected from std::cout
         */
        std::string infer(const std::string& tree_search_type = "NORMAL", const bool shallow_tree_search = false);
        
        /*! \brief Compute the likelihood of the current tree
         *
         * Compute the likelihood of the current tree
         * @return the likelihood of the tree
         */
        RealNumType computeLh();
        
        /*! \brief Compute branch supports of the current tree
         *
         * Compute branch supports (aLRT-SH) of the current tree
         * @param[in] num_threads number of threads (optional)
         * @param[in] num_replicates a positive number of replicates (optional)
         * @param[in] epsilon a positive epsilon, which is used to avoid rounding effects, when the best and second best NNI trees have nearly identical site log-likelihood values (see Guindon et al., 2010) (optional)
         * @param[in] allow_replacing_ML_tree TRUE to allow replacing the ML tree by a higher likelihood tree found when computing branch supports (optional)
         * @Return a string contains all messages redirected from std::cout
         */
        std::string computeBranchSupports(const int num_threads = 1, const int num_replicates = 1000, const double epsilon = 0.1, const bool allow_replacing_ML_tree = true);
        
        /*! \brief Export the tree string in NEWICK format.
         *
         * Export the phylogenetic tree in a NEWICK string.
         * @param[in] tree_type the type of the output tree (optional): "BIN" - bifurcating tree or "MUL" - multifurcating tree
         * @param[in] show_branch_supports TRUE to output branch supports (aLRT-SH)
         * @return a tree string in NEWICK format
         */
        std::string exportString(const std::string& tree_type = "BIN", const bool show_branch_supports = false);
        
        /*! \brief Get an instance of Params, which stores all parameter settings
         *
         * Get an instance of Params, which stores all parameter settings. Users can use that Params instance to change other (minor) settings of CMaple (which is not yet supported via the above APIs).
         * @return an instance of cmaple::Params
         */
        cmaple::Params& getParams();
        
    private:
        /**
         A tree base instance
         */
        TreeBase* tree_base;
        
        /*! \brief Initialize a tree base instance
         *
         * Initialize tree base instance
         * @param[in] aln an alignment
         * @param[in] model a substitution model
         */
        void initTree(Alignment& aln, Model& model);
    };
}
