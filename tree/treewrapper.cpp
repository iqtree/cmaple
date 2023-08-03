#include "treewrapper.h"

using namespace std;
using namespace cmaple;

void cmaple::TreeWrapper::initTree(Alignment* aln, Model* model)
{    
    // Validate input aln
    if (!aln || !aln->data.size())
        throw std::invalid_argument("Alignment is empty. Please call read(...) first!");
    
    // Validate input model
    if (!model || !model->model_base)
        throw std::invalid_argument("Model is null!");
    
    // Init the tree base
    tree_base = new Tree();
    
    // Attach alignment and model
    tree_base->attachAlnModel(aln, model->model_base);
}

cmaple::TreeWrapper::TreeWrapper(Alignment* aln, Model* model, std::istream& tree_stream, const bool fixed_blengths):tree_base(nullptr)
{
    // Initialize a tree
    initTree(aln, model);
    
    // load the input tree from a stream
    load(tree_stream, fixed_blengths);
}

cmaple::TreeWrapper::TreeWrapper(Alignment* aln, Model* model, const std::string& tree_filename, const bool fixed_blengths):tree_base(nullptr)
{
    // Initialize a tree
    initTree(aln, model);
    
    // load the input tree from a tree file (if users specify it)
    if (tree_filename.length())
        load(tree_filename, fixed_blengths);
    // Unable to keep blengths fixed if users don't input a tree
    else if (fixed_blengths && cmaple::verbose_mode > cmaple::VB_QUIET)
        outWarning("Disable the option to keep the branch lengths fixed because users didn't supply an input tree.");
}

cmaple::TreeWrapper::~TreeWrapper()
{
    if (tree_base)
    {
        delete tree_base;
        tree_base = nullptr;
    }
}

void cmaple::TreeWrapper::load(std::istream& tree_stream, const bool fixed_blengths)
{
    // Make sure tree_base is not null
    ASSERT(tree_base && "tree_base is null");
    
    // load the input tree from a stream
    tree_base->loadTree(tree_stream, fixed_blengths);
}

void cmaple::TreeWrapper::load(const std::string& tree_filename, const bool fixed_blengths)
{
    // Validate input
    if (!tree_filename.length())
        throw std::invalid_argument("The tree file name is empty");
    
    // open the treefile
    std::ifstream tree_stream;
    try {
        tree_stream.exceptions(ios::failbit | ios::badbit);
        tree_stream.open(tree_filename);
        load(tree_stream, fixed_blengths);
        tree_stream.close();
    } catch (ios::failure const &e) {
        std::string error_msg(ERR_READ_INPUT);
        throw ios::failure(error_msg + tree_filename);
    }
}

void cmaple::TreeWrapper::changeAln(Alignment* n_aln)
{
    ASSERT(tree_base);
    
    // Validate input aln
    if (!n_aln || !n_aln->data.size())
        throw std::invalid_argument("Alignment is empty. Please call read(...) first!");
    
    // change the aln_base
    tree_base->changeAln(n_aln);
}

void cmaple::TreeWrapper::changeModel(Model* n_model)
{
    ASSERT(tree_base);
    
    // Validate input model
    if (!n_model || !n_model->model_base)
        throw std::invalid_argument("Model is null!");
    
    // change the model_base
    tree_base->changeModel(n_model->model_base);
}

std::string cmaple::TreeWrapper::infer(const TreeSearchType tree_search_type, const bool shallow_tree_search)
{
    ASSERT(tree_base);
    
    // Validate tree_search_type
    if (tree_search_type == UNKNOWN_TREE_SEARCH)
        throw std::invalid_argument("Unknown tree search type. Please check and try again!");
    
    // Redirect the original src_cout to the target_cout
    streambuf* src_cout = cout.rdbuf();
    ostringstream target_cout;
    cout.rdbuf(target_cout.rdbuf());
    
    tree_base->doInference(tree_search_type, shallow_tree_search);
    
    // Restore the source cout
    cout.rdbuf(src_cout);

    // Will output our Hello World! from above.
    return target_cout.str();
}

RealNumType cmaple::TreeWrapper::computeLh()
{
    ASSERT(tree_base);
    return tree_base->calculateLh();
}

std::string cmaple::TreeWrapper::computeBranchSupports(const int num_threads, const int num_replicates, const double epsilon, const bool allow_replacing_ML_tree)
{
    ASSERT(tree_base);
        
    // Redirect the original src_cout to the target_cout
    streambuf* src_cout = cout.rdbuf();
    ostringstream target_cout;
    cout.rdbuf(target_cout.rdbuf());
    
    // record the start time
    auto start = getRealTime();
    if (cmaple::verbose_mode >= cmaple::VB_MED)
        cout << "Calculating branch supports" << endl;
    
    // validate inputs
    if (num_threads < 0)
        throw std::invalid_argument("Number of threads cannot be negative!");
    if (num_replicates <= 0)
        throw std::invalid_argument("Number of replicates must be positive!");
    if (epsilon < 0)
        throw std::invalid_argument("Epsilon cannot be negative!");
    
    // Only compute the branch supports
    tree_base->calculateBranchSupport(num_threads, num_replicates, epsilon, allow_replacing_ML_tree);
    
    // show the runtime for calculating branch supports
    auto end = getRealTime();
    if (cmaple::verbose_mode >= cmaple::VB_MAX)
        cout << " - Time spent on calculating branch supports: " << std::setprecision(3) << end - start << endl;
    
     // Restore the source cout
     cout.rdbuf(src_cout);

     // return the redirected messages
     return target_cout.str();
}

std::string cmaple::TreeWrapper::exportString(const TreeType tree_type, const bool show_branch_supports) const
{
    ASSERT(tree_base);
    
    // validate tree_type
    if (tree_type == UNKNOWN_TREE)
        throw std::invalid_argument("Unknown tree type. Please use BIN_TREE or MUL_TREE!");
    
    return tree_base->exportTreeString(tree_type, show_branch_supports);
}

cmaple::Params& cmaple::TreeWrapper::getParams()
{
    ASSERT(tree_base && tree_base->params);
    return *tree_base->params;
}

std::ostream& cmaple::operator<<(std::ostream& out_stream, const cmaple::TreeWrapper& tree)
{
    out_stream << tree.exportString();
    return out_stream;
}

std::istream& cmaple::operator>>(std::istream& in_stream, cmaple::TreeWrapper& tree)
{
    // go back to the beginning og the stream
    in_stream.clear();
    in_stream.seekg(0);
    
    // read the stream
    tree.load(in_stream);
    return in_stream;
}
