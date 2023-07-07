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
         * @param[in] final_blength_optimization TRUE to optimize all the branch lengths at the end of the inference
         */
        Tree(Alignment& aln, Model& model, std::istream& tree_stream, const bool final_blength_optimization = true);
        
        /*! \brief Tree constructor
         *
         * Tree constructor from a (bifurcating or multifurcating) tree, which may or may not contain all taxa from a file in NEWICK format
         * @param[in] aln an alignment
         * @param[in] model a substitution model
         * @param[in] tree_filename Name of a tree file (optinal)
         * @param[in] final_blength_optimization TRUE to optimize all the branch lengths at the end of the inference
         */
        Tree(Alignment& aln, Model& model, const std::string& tree_filename = "", const bool final_blength_optimization = true);
        
        /*! \brief Tree destructor
         *
         * Tree destructor
         */
        ~Tree();
        
        /*! \brief Infer a tree from an alignment using a substitution model
         *
         * Infer a phylogenetic tree from an alignment using a substitution model.
         * - If users didn't supply an input tree or supplied an incomplete tree (which doesn't contain all taxa in the alignment), this function will:
         * + do placement (add missing taxa from the alignment to the tree)
         * + apply a NORMAL tree search (which do SPR moves on newly-added nodes)
         * + optimize all branch lengths
         * - If users already supplied a complete tree, this function optimize all branch lengths (by default)
         */
        void infer();
        
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
         * @return a status code: 0 - success, non-zero - failed with errors
         */
        int computeBranchSupports(const int num_threads = 1, const int num_replicates = 1000, const double epsilon = 0.1);
        
        /*! \brief Export the tree string in NEWICK format.
         *
         * Export the phylogenetic tree in a NEWICK string.
         * @param[in] tree_type the type of the output tree (optional): "BIN" - bifurcating tree or "MUL" - multifurcating tree
         * @param[in] show_branch_supports TRUE to output branch supports (aLRT-SH)
         * @return a tree string in NEWICK format
         */
        std::string exportString(const std::string& tree_type = "BIN", const bool show_branch_supports = false);
        
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
         * @param[in] final_blength_optimization TRUE to optimize all the branch lengths at the end of the inference
         */
        void initTree(Alignment& aln, Model& model, const bool final_blength_optimization);
    };
}
