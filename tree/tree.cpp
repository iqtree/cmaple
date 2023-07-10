#include "tree.h"

using namespace std;
using namespace cmaple;

void cmaple::Tree::initTree(Alignment& aln, Model& model, const bool final_blength_optimization)
{
    // Init the tree base
    tree_base = new TreeBase();
    
    // TODO: preserve this option
    // set skip_final_blength_optimization
    tree_base->params->optimize_blength = final_blength_optimization;
    
    // Attach alignment and model
    tree_base->attachAlnModel(aln.aln_base, model.model_base);
}

cmaple::Tree::Tree(Alignment& aln, Model& model, std::istream& tree_stream, const bool final_blength_optimization):tree_base(nullptr)
{
    // Initialize a tree
    initTree(aln, model, final_blength_optimization);
    
    // load the input tree from a stream
    tree_base->loadTree(tree_stream);
}

cmaple::Tree::Tree(Alignment& aln, Model& model, const std::string& tree_filename, const bool final_blength_optimization):tree_base(nullptr)
{
    // Initialize a tree
    initTree(aln, model, final_blength_optimization);
    
    // load the input tree from a tree file (if users specify it)
    if (tree_filename.length())
    {
        // open the treefile
        std::ifstream tree_stream;
        try {
            tree_stream.exceptions(ios::failbit | ios::badbit);
            tree_stream.open(tree_filename);
            tree_base->loadTree(tree_stream);
            tree_stream.close();
        } catch (ios::failure) {
            outError(ERR_READ_INPUT, tree_filename);
        }
    }
}

cmaple::Tree::~Tree()
{
    if (tree_base)
    {
        delete tree_base;
        tree_base = nullptr;
    }
}

std::string cmaple::Tree::infer()
{
    ASSERT(tree_base);
    
    // Redirect the original src_cout to the target_cout
    streambuf* src_cout = cout.rdbuf();
    ostringstream target_cout;
    cout.rdbuf(target_cout.rdbuf());
    
    tree_base->doInference();
    
    // Restore the source cout
    cout.rdbuf(src_cout);

    // Will output our Hello World! from above.
    return target_cout.str();
}

RealNumType cmaple::Tree::computeLh()
{
    ASSERT(tree_base);
    return tree_base->calculateLh();
}

std::string cmaple::Tree::computeBranchSupports(const int num_threads, const int num_replicates, const double epsilon)
{
    ASSERT(tree_base);
    
    // Redirect the original src_cout to the target_cout
    streambuf* src_cout = cout.rdbuf();
    ostringstream target_cout;
    cout.rdbuf(target_cout.rdbuf());
    
    // validate inputs
    if (num_threads < 0)
        outError("Number of threads cannot be negative!");
    if (num_replicates <= 0)
        outError("Number of replicates must be positive!");
    if (epsilon < 0)
        outError("Epsilon cannot be negative!");
    
    // Only compute the branch supports
    tree_base->calculateBranchSupport(num_threads, num_replicates, epsilon);
    
     // Restore the source cout
     cout.rdbuf(src_cout);

     // Will output our Hello World! from above.
     return target_cout.str();
}

std::string cmaple::Tree::exportString(const std::string& tree_type, const bool show_branch_supports)
{
    return tree_base->exportTreeString(tree_type, show_branch_supports);
}
