#include "tree.h"

using namespace std;
using namespace cmaple;

cmaple::Tree::Tree(Alignment& aln, Model& model, const std::string& tree_filename, const bool final_blength_optimization):tree_base(nullptr)
{
    // Init the tree base
    tree_base = new TreeBase();
    
    // TODO: preserve this option
    // set skip_final_blength_optimization
    tree_base->params->optimize_blength = final_blength_optimization;
    
    // Attach alignment and model
    tree_base->attachAlnModel(aln.aln_base, model.model_base);
    
    // load the input tree (if users specify it)
    if (tree_filename.length())
        tree_base->loadTree(tree_filename);
}

cmaple::Tree::~Tree()
{
    if (tree_base)
    {
        delete tree_base;
        tree_base = nullptr;
    }
}

void cmaple::Tree::infer()
{
    ASSERT(tree_base);
    tree_base->doInference();
}

RealNumType cmaple::Tree::computeLh()
{
    ASSERT(tree_base);
    return tree_base->calculateLh();
}

int cmaple::Tree::computeBranchSupports(const int num_threads, const int num_replicates, const double epsilon)
{
    ASSERT(tree_base);
    
    // validate inputs
    if (num_threads < 0)
    {
        std::cout << "Number of threads cannot be negative!" << std::endl;
        return CODE_ERROR_1;
    }
    if (num_replicates <= 0)
    {
        std::cout << "Number of replicates must be positive!" << std::endl;
        return CODE_ERROR_1;
    }
    if (epsilon < 0)
    {
        std::cout << "Epsilon cannot be negative!" << std::endl;
        return CODE_ERROR_1;
    }
    
    // Only compute the branch supports
    tree_base->calculateBranchSupport(num_threads, num_replicates, epsilon);
    
    return CODE_SUCCESS;
}

std::string cmaple::Tree::exportString(const std::string& tree_type, const bool show_branch_supports)
{
    return tree_base->exportTreeString(tree_type, show_branch_supports);
}
