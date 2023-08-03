#pragma once
#include "tree.h"
#include "../alignment/alignment.h"
#include "../model/model.h"

namespace cmaple
{
    /** Class represents a phylogenetic tree */
    class TreeWrapper {
    public:
        /*! \brief Constructor from a stream of a (bifurcating or multifurcating) tree (with/without branch lengths in NEWICK format), which may or may not contain all taxa in the alignment
         * @param[in] aln An alignment
         * @param[in] model A substitution model
         * @param[in] tree_stream A stream of an input tree
         * @param[in] fixed_blengths TRUE to keep the input branch lengths unchanged (optional)
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
        TreeWrapper(Alignment* aln, Model* model, std::istream& tree_stream, const bool fixed_blengths = false);
        
        /*! \brief Constructor from an optional (bifurcating or multifurcating) tree (with/without branch lengths in NEWICK format), which may or may not contain all taxa in the alignment
         * @param[in] aln An alignment
         * @param[in] model A substitution model
         * @param[in] tree_filename Name of a tree file (optinal)
         * @param[in] fixed_blengths TRUE to keep the input branch lengths unchanged (optional)
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
        TreeWrapper(Alignment* aln, Model* model, const std::string& tree_filename = "", const bool fixed_blengths = false);
        
        /*! \brief Destructor
         */
        ~TreeWrapper();
        
        /*! \brief Load a tree from a stream of a (bifurcating or multifurcating) tree (with/without branch lengths) in NEWICK format, which may or may not contain all taxa in the alignment
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
        
        /*! \brief Load a tree from a (bifurcating or multifurcating) tree (with/without branch lengths) in NEWICK format, which may or may not contain all taxa in the alignment
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
        
        /*! \brief Infer a phylogenetic tree
         * - If users didn't supply an input tree or supplied an incomplete tree (which doesn't contain all the taxa in the alignment) when initializing the tree (by Tree() constructor), this function:
         *  + performs placement (i.e., adding missing taxa from the alignment to the tree)
         *  + applies a NORMAL tree search (which does SPR moves only on newly-added nodes)
         *  + optimizes all branch lengths
         * - If users already supplied a complete tree, this function:
         *  + by default, does neither placment nor tree search, but it optimizes all branch lengths.
         *  + If users want to keep the branch lengths fixed, they should set fixed_blengths = true when initializing the tree (by Tree() constructor);
         *  + If users want to use the input tree as a starting tree (then performs SPR moves and optimizes branch lengths), they should set tree_search_type = MORE_ACCURATE
         *
         * @param[in] tree_search_type One of the following tree search:
         * <br><em>FAST_TREE_SEARCH</em>: no tree search (placement only).
         * <br><em>NORMAL_TREE_SEARCH</em>: only consider pruning branches at newly-added nodes when seeking SPR moves.
         * <br><em>MORE_ACCURATE_TREE_SEARCH</em>: consider all nodes when seeking SPR moves.
         * @param[in] shallow_tree_search TRUE ton enable a shallow tree search before a deeper tree search
         * @return a string contains all messages redirected from std::cout (for information and debugging purpuses only). To output the tree in NEWICK format, one could call exportString() later
         * @throw std::invalid\_argument if tree\_search\_type is unknown
         * @throw std::logic\_error if any of the following situations occur.
         * - the attached substitution model is unknown/unsupported
         * - unexpected values/behaviors found during the operations
         */
        std::string infer(const TreeSearchType tree_search_type = NORMAL_TREE_SEARCH, const bool shallow_tree_search = false);
        
        /*! \brief Compute the log likelihood of the current tree, which may or may not contain all taxa in the alignment
         * @return The log likelihood of the current tree
         * @throw std::logic\_error if any of the following situations occur.
         * - the tree is empty
         * - unexpected values/behaviors found during the operations
         */
        RealNumType computeLh();
        
        /*! \brief Compute branch supports ([aLRT-SH](https://academic.oup.com/sysbio/article/59/3/307/1702850)) of the current tree, which may or may not contain all taxa in the alignment
         * @param[in] num_threads The number of threads (optional)
         * @param[in] num_replicates A positive number of replicates (optional)
         * @param[in] epsilon A positive epsilon (optional), which is used to avoid rounding effects, when the best and second best NNI trees have nearly identical site log-likelihood values (see [Guindon et al., 2010](https://academic.oup.com/sysbio/article/59/3/307/1702850))
         * @param[in] allow_replacing_ML_tree TRUE to allow replacing the ML tree by a higher likelihood tree found when computing branch supports (optional)
         * @return A string contains all messages redirected from std::cout (for information and debugging purpuses only). To output the branch supports values, one could call exportString("BIN", true) later
         * @throw std::invalid\_argument if any of the following situations occur.
         * - num_threads < 0 or num_threads > the number of CPU cores
         * - num_replicates <= 0
         * - epsilon < 0
         *
         * @throw std::logic\_error if any of the following situations occur.
         * - the tree is empty
         * - unexpected values/behaviors found during the operations
         */
        std::string computeBranchSupports(const int num_threads = 1, const int num_replicates = 1000, const double epsilon = 0.1, const bool allow_replacing_ML_tree = true);
        
        /*! \brief Export the phylogenetic tree  to a string in NEWICK format.
         * @param[in] tree_type The type of the output tree (optional): BIN_TREE (bifurcating tree), MUL_TREE (multifurcating tree)
         * @param[in] show_branch_supports TRUE to output the branch supports (aLRT-SH values)
         * @return A tree string in NEWICK format
         * @throw std::invalid\_argument if any of the following situations occur.
         * - n\_tree\_type is unknown
         * - show\_branch\_supports = true but branch support values have yet been computed
         */
        std::string exportString(const TreeType tree_type = BIN_TREE, const bool show_branch_supports = false) const;
        
        /*! \brief Get an instance of cmaple::Params, which stores all parameter settings. Users can use that params instance to change [**other minor settings**](classcmaple_1_1_params.html) of CMaple (which are not yet supported via the APIs)
         * @return An instance of cmaple::Params
         */
        cmaple::Params& getParams();
        
    private:
        /**
         A tree base instance
         */
        Tree* tree_base;
        
        /*! \brief Initialize tree base instance
         * @param[in] aln An alignment
         * @param[in] model A substitution model
         * @throw std::invalid\_argument If any of the following situation occurs.
         * - the sequence type is unsupported (neither DNA (for nucleotide data) nor AA (for protein data))
         * - the alignment is empty
         * - model is unknown/unsupported
         *
         * @throw std::logic\_error if the reference genome is empty
         */
        void initTree(Alignment* aln, Model* model);
    };

    /** \brief Customized << operator to output the tree string in a (bifurcating) NEWICK format
     */
    std::ostream& operator<<(std::ostream& out_stream, const cmaple::TreeWrapper& tree);

    /** \brief Customized >> operator to read the tree from a stream
     */
    std::istream& operator>>(std::istream& in_stream, cmaple::TreeWrapper& tree);
}
