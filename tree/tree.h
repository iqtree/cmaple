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
         * Tree constructor from a strictly bifurcating tree with may or may not contain all taxa from a file in NEWICK format
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
         * If users already inputted a tree:
         * - If it is an incomplete tree (which doesn't contain all taxa in the alignment), this function will:
         * + do placement (add missing taxa from the alignment to the tree)
         * + apply a NORMAL tree search (which do SPR moves on newly-added nodes)
         * + optimize all branch lengths
         * - If it is a complete tree, this function optimize all branch lengths (by default)
         */
        void infer();
        
        /*! \brief Compute the likelihood of the current tree
         *
         * Compute the likelihood of the current tree
         * @return the likelihood of the tree
         */
        RealNumType computeLh();
        
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
    };
}
