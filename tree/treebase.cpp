#include "treebase.h"

#include <cassert>
#include <utils/matrix.h>

using namespace std;
using namespace cmaple;

// explicit instantiation of templates
template void cmaple::TreeBase::refreshAllLhs<4>(bool);
template void cmaple::TreeBase::refreshAllLhs<20>(bool);

template double cmaple::TreeBase::improveEntireTree<4>(bool);
template double cmaple::TreeBase::improveEntireTree<20>(bool);

template void cmaple::TreeBase::updateModelParams<4>();
template void cmaple::TreeBase::updateModelParams<20>();

template void cmaple::TreeBase::seekSamplePlacement<4>(const Index, const NumSeqsType, const std::unique_ptr<SeqRegions>&, Index&, RealNumType&, bool&, RealNumType&, RealNumType&, Index&);
template void cmaple::TreeBase::seekSamplePlacement<20>(const Index, const NumSeqsType, const std::unique_ptr<SeqRegions>&, Index&, RealNumType&, bool&, RealNumType&, RealNumType&, Index&);

template void cmaple::TreeBase::placeNewSampleAtNode<4>(const Index, std::unique_ptr<SeqRegions>&, const NumSeqsType, const RealNumType, const RealNumType, const RealNumType, const Index);
template void cmaple::TreeBase::placeNewSampleAtNode<20>(const Index, std::unique_ptr<SeqRegions>&, const NumSeqsType, const RealNumType, const RealNumType, const RealNumType, const Index);

template int cmaple::TreeBase::optimizeBranchLengths<4>();
template int cmaple::TreeBase::optimizeBranchLengths<20>();

template void cmaple::TreeBase::placeNewSampleMidBranch<4>(const Index&, std::unique_ptr<SeqRegions>&, const NumSeqsType, const RealNumType);
template void cmaple::TreeBase::placeNewSampleMidBranch<20>(const Index&, std::unique_ptr<SeqRegions>&, const NumSeqsType, const RealNumType);

void cmaple::TreeBase::setupFuncPtrs()
{
    ASSERT(aln);
    
    switch (aln->num_states) {
        case 4:
            loadTreePtr = &TreeBase::loadTreeTemplate<4>;
            changeAlnPtr = &TreeBase::changeAlnTemplate<4>;
            changeModelPtr = &TreeBase::changeModelTemplate<4>;
            doInferencePtr = &TreeBase::doInferenceTemplate<4>;
            calculateLhPtr = &TreeBase::calculateLhTemplate<4>;
            calculateBranchSupportPtr = &TreeBase::calculateBranchSupportTemplate<4>;
            break;
        case 20:
            loadTreePtr = &TreeBase::loadTreeTemplate<20>;
            changeAlnPtr = &TreeBase::changeAlnTemplate<20>;
            changeModelPtr = &TreeBase::changeModelTemplate<20>;
            doInferencePtr = &TreeBase::doInferenceTemplate<20>;
            calculateLhPtr = &TreeBase::calculateLhTemplate<20>;
            calculateBranchSupportPtr = &TreeBase::calculateBranchSupportTemplate<20>;
            break;
            
        default:
            outError("Sorry! currently we only support DNA and Protein data!");
            break;
    }
}

void cmaple::TreeBase::setupBlengthThresh()
{
    default_blength = 1.0 / aln->ref_seq.size();
    min_blength = params->fixed_min_blength == -1 ? params->min_blength_factor * default_blength : params->fixed_min_blength;
    max_blength = params->max_blength_factor * default_blength;
    min_blength_mid = params->min_blength_mid_factor * default_blength;
    min_blength_sensitivity = min_blength * 1e-5;
    half_min_blength_mid = min_blength_mid * 0.5;
    half_max_blength = max_blength * 0.5;
    double_min_blength = min_blength + min_blength;
    
    // compute thresholds for approximations
    params->threshold_prob2 = params->threshold_prob * params->threshold_prob;
}

void cmaple::TreeBase::resetSeqAdded()
{
    ASSERT(aln);
    const NumSeqsType num_sequences = aln->data.size();
    sequence_added.resize(num_sequences);
    for (auto i = 0; i < num_sequences; ++i)
        sequence_added[i] = false;
}

void cmaple::TreeBase::attachAlnModel(AlignmentBase* n_aln, ModelBase* n_model)
{
    // Validate inputs
    ASSERT(n_aln && n_model);
    // Sequence type (~num_states) of model must match that of alignment
    ASSERT(n_aln->num_states == n_model->num_states_ && "Sequence type (~num_states) of model must match that of alignment");
    
    // Attach aln
    aln = n_aln;
    // record the current tree in the list of trees that the alignment is attached to
    n_aln->attached_trees.insert(this);
    
    // Attach model
    model = n_model;
    
    // reserve spaces for nodes
    const NumSeqsType num_seqs = aln->data.size();
    nodes.reserve(num_seqs + num_seqs);
    
    // Initialize sequence_added -> all sequences has yet added to the tree
    resetSeqAdded();
    
    // init params & thresholds
    setupBlengthThresh();
    
    // extract a backup vector of sequence names
    seq_names.resize(num_seqs);
    for (auto i = 0; i < num_seqs; ++i)
        seq_names[i] = aln->data[i].seq_name;
    
    // update model according to the data from the alignment
    updateModelByAln();
    
    // setup function pointers
    setupFuncPtrs();
}

void cmaple::TreeBase::updateModelByAln()
{
    // extract related info (freqs, log_freqs) of the ref sequence
    model->extractRefInfo(aln);
    
    // update the mutation matrix
    model->updateMutationMatEmpirical(aln);
    
    // compute cumulative rates of the ref sequence
    model->computeCumulativeRate(aln);
}

void cmaple::TreeBase::loadTree(std::istream& tree_stream, const bool fixed_blengths)
{
    (this->*loadTreePtr)(tree_stream, fixed_blengths);
}

template <const StateType num_states>
void cmaple::TreeBase::loadTreeTemplate(std::istream& tree_stream, const bool n_fixed_blengths)
{
    // reset variables in tree
    ASSERT(aln);
    const NumSeqsType num_seqs = aln->data.size();
    // reset nodes
    nodes.clear();
    nodes.reserve(num_seqs + num_seqs);
    // reset node_lhs
    node_lhs.clear();
    node_lhs.reserve(num_seqs);
    node_lhs.emplace_back(0);
    // reset sequence_added
    resetSeqAdded();
    
    // read tree from the input treefile
    bool missing_blength = readTree(tree_stream);
    
    // make sure users can only keep the blengths fixed if they input a complete tree with branch lengths
    if (n_fixed_blengths && (!isComplete() || missing_blength))
    {
        if (cmaple::verbose_mode > cmaple::VB_QUIET)
            outWarning("Disable the option to keep the branch lengths fixed because the input tree is incomplete (i.e., not containing all taxa from the alignment) or contains missing branch length(s).");
        fixed_blengths = false;
    }
    else
        fixed_blengths = n_fixed_blengths;
    
    // update model params & partial lhs along tree after loading the tree from a file
    updateModelLhAfterLoading<num_states>();
    
    // If the tree contains any missing blengths -> re-estimate the blengths
    if (missing_blength)
    {
        if (cmaple::verbose_mode >= cmaple::VB_MED)
            std::cout << "The input tree contains missing branch lengths. Re-estimate all branch lengths." << std::endl;
        
        // optimize blengths
        optimizeBranchLengths<num_states>();
        
        // refresh all likelihoods
        refreshAllLhs<num_states>();
    }
    
    // set outdated = false at all nodes to avoid considering SPR moves at those nodes
    resetSPRFlags(false, true, false);
}

template <const cmaple::StateType  num_states>
void cmaple::TreeBase::updateModelLhAfterLoading()
{
    // do nothing on an empty tree
    if (!nodes.size()) return;
    
    // calculate all lower, upper left/right likelihoods
    refreshAllLhs<num_states>(true);
    
    // update model params
    /*cout << " - Model params before updating: " << endl;
     tree.showModelParams();*/
    updateModelParams<num_states>();
    /*cout << " - Model params after updating: " << endl;
     tree.showModelParams();*/
    
    // refresh all lower after updating model params
    performDFS<&TreeBase::updateLowerLh<num_states>>();
}

void cmaple::TreeBase::changeAln(AlignmentBase* aln)
{
    (this->*changeAlnPtr)(aln);
}

template <const cmaple::StateType  num_states>
void cmaple::TreeBase::changeAlnTemplate(AlignmentBase* n_aln)
{
    // make sure both alignments have the same seqtype
    if (aln->getSeqType() != n_aln->getSeqType())
        outError("Sorry! we have yet supported changing to a new alignment with a sequence type different from the current one.");
    
    // show information
    if (cmaple::verbose_mode >= cmaple::VB_MED)
        std::cout << "Changing the alignment" << std::endl;
    
    // change the alignment
    aln = n_aln;
    // record the current tree in the list of trees that the alignment is attached to
    n_aln->attached_trees.insert(this);
    
    // re-mark taxa in the new alignment, which already existed in the current tree
    remarkExistingSeqs();
    
    // extract a backup vector of sequence names
    seq_names.resize(aln->data.size());
    for (auto i = 0; i < aln->data.size(); ++i)
        seq_names[i] = aln->data[i].seq_name;
    
    // update model according to the data in the new alignment
    updateModelByAln();
    
    // make sure users can only keep the blengths fixed if they input a complete tree with branch lengths
    if (fixed_blengths && !isComplete())
    {
        if (cmaple::verbose_mode > cmaple::VB_QUIET)
            outWarning("Disable the option to keep the branch lengths fixed because the input tree is incomplete (i.e., not containing all taxa from the alignment).");
        fixed_blengths = false;
    }
    
    // update Model params and partial lhs along the current tree (if any)
    updateModelLhAfterLoading<num_states>();
}

void cmaple::TreeBase::changeModel(ModelBase* model)
{
    (this->*changeModelPtr)(model);
}

template <const cmaple::StateType  num_states>
void cmaple::TreeBase::changeModelTemplate(ModelBase* n_model)
{
    ASSERT(n_model && model && aln);
    
    // Make sure we use the updated alignment (in case users re-read the alignment from a new file after attaching the alignment to the tree)
    if (aln->attached_trees.find(this) == aln->attached_trees.end())
        changeAln(aln);
    
    // make sure both models have the same seqtype
    if (model->num_states_ != n_model->num_states_)
        outError("Sorry! we have yet supported changing to a new model with a sequence type different from the current one.");
    
    // do nothing if new model is the same as the current one
    if (model == n_model) return;
    
    // show information
    if (cmaple::verbose_mode >= cmaple::VB_MED)
        std::cout << "Changing the model" << std::endl;
    
    // change the model
    model = n_model;
    
    // update model according to the data in the alignment
    updateModelByAln();
    
    // update Model params and partial lhs along the current tree (if any)
    updateModelLhAfterLoading<num_states>();
}

void cmaple::TreeBase::doInference(const TreeSearchType tree_search_type, const bool shallow_tree_search)
{
    (this->*doInferencePtr)(tree_search_type, shallow_tree_search);
}

template <const StateType num_states>
void cmaple::TreeBase::doInferenceTemplate(const TreeSearchType tree_search_type, const bool shallow_tree_search)
{
    // validate input
    ASSERT(aln->ref_seq.size() > 0 && "Reference sequence is not found!");
    ASSERT(aln->data.size() >= 3 && "The number of input sequences must be at least 3! Please check and try again!");
    
    // Validate tree search type
    if (tree_search_type == UNKNOWN_TREE_SEARCH)
        outError("Unknown tree search type!");
    
    // Make sure we use the updated alignment (in case users re-read the alignment from a new file after attaching the alignment to the tree)
    if (aln->attached_trees.find(this) == aln->attached_trees.end())
        changeAln(aln);
    
    //  Check whether we infer a phologeny from an input tree
    const bool from_input_tree = nodes.size() > 0;
        
    // 1. Build an initial tree
    buildInitialTree<num_states>(from_input_tree);
    
    // 2. Optimize the tree with SPR if there is any new nodes added to the tree
    optimizeTree<num_states>(from_input_tree, tree_search_type, shallow_tree_search);
}

template <const StateType num_states>
void cmaple::TreeBase::buildInitialTree(const bool from_input_tree)
{
    // show information
    if (cmaple::verbose_mode >= cmaple::VB_MED)
        std::cout << "Performing placement" << std::endl;
    
    // record the start time
    auto start = getRealTime();
    
    // dummy variables
    const PositionType seq_length = aln->ref_seq.size();
    const PositionType num_seqs = aln->data.size();
    PositionType num_new_sequences = num_seqs;
    Sequence* sequence = &aln->data.front();
    // make sure we allocate enough space to store all nodes
    if (nodes.capacity() < num_seqs + num_seqs)
        nodes.reserve(num_seqs + num_seqs);
    NumSeqsType i = 0;
    
    // if users don't input a tree -> create the root from the first sequence
    if (!from_input_tree)
    {
        // place the root node
        root_vector_index = 0;
        nodes.emplace_back(LeafNode(0));
        PhyloNode& root = nodes[0];
        root.setPartialLh(TOP, std::move(sequence->getLowerLhVector(seq_length, num_states, aln->getSeqType())));
        root.getPartialLh(TOP)->computeTotalLhAtRoot<num_states>(root.getTotalLh(), model);
        root.setUpperLength(0);
        
        // move to the next sequence in the alignment
        sequence_added[i] = true;
        ++sequence;
        ++i;
    }
    
    // iteratively place other samples (sequences)
    for (; i < num_seqs; ++i, ++sequence)
    {
        // don't add sequence that was already added in the input tree
        if (from_input_tree && sequence_added[i])
        {
            --num_new_sequences;
            continue;
        }
        // otherwise, mark the current sequence as added
        else
            sequence_added[i] = true;
        
        // get the lower likelihood vector of the current sequence
        std::unique_ptr<SeqRegions> lower_regions = sequence->getLowerLhVector(seq_length, num_states, aln->getSeqType());
        
        // update the mutation matrix from empirical number of mutations observed from the recent sequences
        if (i % params->mutation_update_period == 0)
            model->updateMutationMatEmpirical(aln);
        
        // NHANLT: debug
        //if ((*sequence)->seq_name == "39")
        //    cout << "debug" <<endl;
        
        // seek a position for new sample placement
        Index selected_node_index;
        RealNumType best_lh_diff = MIN_NEGATIVE;
        bool is_mid_branch = false;
        RealNumType best_up_lh_diff = MIN_NEGATIVE;
        RealNumType best_down_lh_diff = MIN_NEGATIVE;
        Index best_child_index;
        seekSamplePlacement<num_states>(Index(root_vector_index, TOP), i, lower_regions, selected_node_index, best_lh_diff, is_mid_branch, best_up_lh_diff, best_down_lh_diff, best_child_index);
        
        // if new sample is not less informative than existing nodes (~selected_node != NULL) -> place the new sample in the existing tree
        if (selected_node_index.getMiniIndex() != UNDEFINED)
        {
            // place new sample as a descendant of a mid-branch point
            if (is_mid_branch)
                placeNewSampleMidBranch<num_states>(selected_node_index, lower_regions, i, best_lh_diff);
            // otherwise, best lk so far is for appending directly to existing node
            else
                placeNewSampleAtNode<num_states>(selected_node_index, lower_regions, i, best_lh_diff, best_up_lh_diff, best_down_lh_diff, best_child_index);
        }
        
        // NHANLT: debug
        //cout << "Added node " << (*sequence)->seq_name << endl;
        //cout << (*sequence)->seq_name << endl;
        //cout << tree.exportTreeString() << endl;
        
        //if ((*sequence)->seq_name == "2219")
        //{
            //cout << tree.exportTreeString() << endl;
            //string output_file(tree.params->output_prefix);
            //exportOutput(output_file + "_init.treefile");
            //exit(0);
        //}
    }
    
    // flag denotes whether there is any new nodes added
    // show the number of new sequences added to the tree
    if (num_new_sequences > 0)
    {
        if (cmaple::verbose_mode > cmaple::VB_QUIET)
            std::cout << num_new_sequences << " sequences have been added to the tree." << std::endl;
        
        // traverse the intial tree from root to re-calculate all likelihoods regarding the latest/final estimated model parameters
        refreshAllLhs<num_states>();
    }
    else if (cmaple::verbose_mode > cmaple::VB_QUIET)
        std::cout << "All sequences were presented in the input tree. No new sequence has been added!" << std::endl;
    
    // show the runtime for building an initial tree
    auto end = getRealTime();
    if (cmaple::verbose_mode >= cmaple::VB_MAX)
        cout << " - Time spent on building an initial tree: " << std::setprecision(3) << end - start << endl;
}

template <const StateType num_states>
void cmaple::TreeBase::optimizeTree(const bool from_input_tree, const TreeSearchType tree_search_type, const bool shallow_tree_search)
{
    // show information
    if (tree_search_type == FAST_TREE_SEARCH && cmaple::verbose_mode >= cmaple::VB_MED)
        std::cout << "No tree search is invoked." << std::endl;
    // tree.params->debug = true;
    // string output_file(params->output_prefix);
    // exportOutput(output_file + "_init.treefile");
    
    // run a shallow (short range) search for tree topology improvement (if neccessary)
    // Don't apply shallow tree search if users chose no tree search
    if (shallow_tree_search && tree_search_type != FAST_TREE_SEARCH)
    {
        // don't apply shallow tree search if users inputted a tree and wanted to run a partial tree search later (because outdated flags may become incorrect due to shallow tree search)
        if (from_input_tree && tree_search_type == NORMAL_TREE_SEARCH)
        {
            if (cmaple::verbose_mode > cmaple::VB_QUIET)
                outWarning("Disable shallow tree search because it's not supported if users input a tree and want to apply a normal tree search.");
        }
        // otherwise, apply a shallow tree search
        else
        {
            if (cmaple::verbose_mode >= cmaple::VB_MED)
                std::cout << "Applying a shallow tree search" << std::endl;
            
            // apply short-range SPR search
            optimizeTreeTopology<num_states>(true);
            // exportOutput(output_file + "_short_search.treefile");
            
            // reset the SPR flags so that we can start a deeper SPR search later
            resetSPRFlags(false, true, true);
        }
    }
    
    // output log-likelihood of the tree
    if (cmaple::verbose_mode >= cmaple::VB_DEBUG)
        std::cout << std::setprecision(10) << "Tree log likelihood (before a deeper tree search): " << calculateLh() << std::endl;
    
    // run a normal search for tree topology improvement
    if (tree_search_type != FAST_TREE_SEARCH)
    {
        if (cmaple::verbose_mode >= cmaple::VB_MED)
        {
            std::string tree_search_str = getTreeSearchStr(tree_search_type);
            transform(tree_search_str.begin(), tree_search_str.end(), tree_search_str.begin(), ::tolower);
            std::cout << "Applying a " + tree_search_str + " tree search" << std::endl;
        }
        
        optimizeTreeTopology<num_states>();
        // exportOutput(output_file + "_topo.treefile");
    }
    
    // traverse the tree from root to re-calculate all likelihoods after optimizing the tree topology
    refreshAllLhs<num_states>();
    
    // output log-likelihood of the tree
    if (cmaple::verbose_mode >= cmaple::VB_DEBUG)
        std::cout << std::setprecision(10) << "Tree log likelihood (after the deeper tree search (if any)): " << calculateLh() << std::endl;
    
    // do further optimization on branch lengths (if needed)
    if (!fixed_blengths)
        optimizeBranchLengthsOfTree<num_states>();
    
    // Make writing and re-reading tree consistently
    // colapse zero branch lengths in leave
    collapseAllZeroLeave();
    // 0. Traverse tree using DFS, at each non-zero-branch leaf -> expand the tree by adding one less-info-seq
    performDFSAtLeave<&cmaple::TreeBase::expandTreeByOneLessInfoSeq<num_states>>();
    
    // NhanLT: update the model params
    if (aln->getSeqType() == SEQ_DNA)
    {
        model->initMutationMat();
        updateModelParams<num_states>();
    }
    
    // traverse the tree from root to re-calculate all lower likelihoods after optimizing branch lengths
    performDFS<&TreeBase::updateLowerLh<num_states>>();
    
    // output log-likelihood of the tree
    if (cmaple::verbose_mode >= cmaple::VB_DEBUG)
        std::cout << std::setprecision(10) << "Tree log likelihood (at the end of the inference): " << calculateLh() << std::endl;
}

template <const StateType num_states>
void cmaple::TreeBase::optimizeTreeTopology(bool short_range_search)
{
    // record the start time
    auto start = getRealTime();
    int num_tree_improvement = short_range_search ? 1 : params->num_tree_improvement;
    
    for (int i = 0; i < num_tree_improvement; ++i)
    {
        // first, set all nodes outdated
        // no need to do so anymore as new nodes were already marked as outdated
        // resetSPRFlags(false, true, true);
        
        // traverse the tree from root to try improvements on the entire tree
        RealNumType improvement = improveEntireTree<num_states>(short_range_search);
        
        // stop trying if the improvement is so small
        if (improvement < params->thresh_entire_tree_improvement)
        {
            if (cmaple::verbose_mode >= cmaple::VB_DEBUG)
                cout << "Small improvement, stopping topological search." << endl;
            break;
        }
        
        // run improvements only on the nodes that have been affected by some changes in the last round, and so on
        for (int j = 0; j < 20; ++j)
        {
            // forget SPR_applied flag to allow new SPR moves
            resetSPRFlags(false, false, true);
            
            improvement = improveEntireTree<num_states>(short_range_search);
            if (cmaple::verbose_mode >= cmaple::VB_DEBUG)
                cout << "Tree was improved by " + convertDoubleToString(improvement) + " at subround " + convertIntToString(j + 1) << endl;
            
            // stop trying if the improvement is so small
            if (improvement < params->thresh_entire_tree_improvement)
                break;
        }
            
    }
    
    // show the runtime for optimize the tree
    auto end = getRealTime();
    if (cmaple::verbose_mode >= cmaple::VB_MAX)
    {
        cout << " - Time spent on";
        cout << (short_range_search ? " a shallow search for" : "");
        cout << " optimizing the tree topology: " << std::setprecision(3) << end - start << endl;
    }
}

template <const StateType num_states>
void cmaple::TreeBase::optimizeBranchLengthsOfTree()
{
    // record the start time
    auto start = getRealTime();
    
    if (cmaple::verbose_mode >= cmaple::VB_MED)
        cout << "Optimizing branch lengths" << endl;
    
    // first, set all nodes outdated
    resetSPRFlags(false, true, true);
   
    // traverse the tree from root to optimize branch lengths
    PositionType num_improvement = optimizeBranchLengths<num_states>();
   
    // run improvements only on the nodes that have been affected by some changes in the last round, and so on
    for (int j = 0; j < 20; ++j)
    {
        // stop trying if the improvement is so small
        if (num_improvement < params->thresh_entire_tree_improvement)
        //if (num_improvement == 0)
            break;
        
        // traverse the tree from root to optimize branch lengths
        num_improvement = optimizeBranchLengths<num_states>();
    }

    // show the runtime for optimize the branch lengths
    auto end = getRealTime();
    if (cmaple::verbose_mode >= cmaple::VB_MAX)
        cout << " - Time spent on optimizing the branch lengths: " << std::setprecision(3) << end - start << endl;
}

RealNumType cmaple::TreeBase::calculateLh()
{
    return (this->*calculateLhPtr)();
}

template <const StateType num_states>
RealNumType cmaple::TreeBase::calculateLhTemplate()
{
    // Make sure the tree is not empty
    if (!nodes.size())
    {
        outError("Tree is empty. Please call infer() to infer a tree from the alignment before computing the likelihood!", false);
        return 0;
    }
    
    // Make sure we use the updated alignment (in case users re-read the alignment from a new file after attaching the alignment to the tree)
    if (aln->attached_trees.find(this) == aln->attached_trees.end())
        changeAln(aln);
    
    // traverse the tree from root to re-calculate all likelihoods after optimizing the tree topology
    refreshAllLhs<num_states>();
    
    // initialize the total_lh by the likelihood from root
    RealNumType total_lh = nodes[root_vector_index].getPartialLh(TOP)->computeAbsoluteLhAtRoot<num_states>(model);
    
    // perform a DFS to add likelihood contributions from each internal nodes
    total_lh += performDFS<&cmaple::TreeBase::computeLhContribution<num_states>>();
   
    // return total_lh
    return total_lh;
}

void cmaple::TreeBase::calculateBranchSupport(const int num_threads, const int num_replicates, const double epsilon, const bool allow_replacing_ML_tree)
{
    return (this->*calculateBranchSupportPtr)(num_threads, num_replicates, epsilon, allow_replacing_ML_tree);
}

template <const StateType num_states>
void cmaple::TreeBase::calculateBranchSupportTemplate(const int num_threads, const int num_replicates, const double epsilon, const bool allow_replacing_ML_tree)
{
    // Make sure the tree is not empty
    if (!nodes.size())
    {
        outError("Tree is empty. Please call infer() to infer a tree from the alignment first!", false);
        return;
    }
    
    // set num_threads
    ASSERT(num_threads >= 0 && "Number of threads must be non-negative!");
    setNumThreads(num_threads);
    
    // set num_replicates
    ASSERT(num_replicates > 0 && "Number of replicates must be positive!");
    params->aLRT_SH_replicates = num_replicates;
    
    // set epsilon
    ASSERT(epsilon >= 0 && "Epsilon must be non-negative!");
    params->aLRT_SH_half_epsilon = epsilon * 0.5;
    
    // Make sure we use the updated alignment (in case users re-read the alignment from a new file after attaching the alignment to the tree)
    if (aln->attached_trees.find(this) == aln->attached_trees.end())
        changeAln(aln);
    
    // refresh all upper left/right lhs before calculating aLRT-SH
    // refreshAllNonLowerLhs<num_states>();
    refreshAllLhs<num_states>();
    
    // 0. Traverse tree using DFS, at each leaf -> expand the tree by adding one less-info-seq to make sure all we compute the aLRT of all internal branches
    performDFSAtLeave<&cmaple::TreeBase::expandTreeByOneLessInfoSeq<num_states>>();
    
    // 1. calculate aLRT for each internal branches, replacing the ML tree if a higher ML NNI neighbor was found
    calculate_aRLT<num_states>(allow_replacing_ML_tree);
    
    // 2. calculate the site lh contributions
    std::vector<RealNumType> site_lh_contributions, site_lh_at_root;
    RealNumType total_lh = calculateSiteLhs<num_states>(site_lh_contributions, site_lh_at_root);
    
    // validate the result
    /*RealNumType tmp_total_lh = 0;
    for (RealNumType site_lh : site_lh_contributions)
        tmp_total_lh += site_lh;
    for (RealNumType site_lh : site_lh_at_root)
        tmp_total_lh += site_lh;
    ASSERT(fabs(tmp_total_lh - total_lh) < 1e-6);*/
    
    // calculate aLRT-SH for all internal branches
    calculate_aRLT_SH<num_states>(site_lh_contributions, site_lh_at_root, total_lh);
    
    // refresh all non-lower likelihoods
    refreshAllNonLowerLhs<num_states>();
}

const string cmaple::TreeBase::exportNodeString(const bool binary, const NumSeqsType node_vec_index, const bool show_branch_supports)
{
    PhyloNode& node = nodes[node_vec_index];
    string output = "(";
    
    // if it's a leaf
    if (!node.isInternal())
        return node.exportString(binary, seq_names, show_branch_supports);
    // if it's an internal node
    else
    {
        /*bool add_comma = false;
         for (Index neighbor_index:node.getNeighborIndexes(TOP))
         {
         if (!add_comma)
         add_comma = true;
         else
         output += ",";
         output += exportNodeString(binary, neighbor_index.getVectorIndex());
         }*/
        output += exportNodeString(binary, node.getNeighborIndex(RIGHT).getVectorIndex(), show_branch_supports);
        output += ",";
        output += exportNodeString(binary, node.getNeighborIndex(LEFT).getVectorIndex(), show_branch_supports);
    }

    string branch_support = "";
    if (show_branch_supports)
    {
        // Make sure Branch supports have been computed
        if (!node.getNodelhIndex())
            outError("Branch supports is not available. Please compute them first!");
        
        branch_support = convertDoubleToString(node_lhs[node.getNodelhIndex()].get_aLRT_SH());
    }
    string length = node.getUpperLength() < 0 ? "0" : convertDoubleToString(node.getUpperLength(), 20);
    output += ")" + branch_support + ":" + length;
    
    return output;
}

const std::string cmaple::TreeBase::exportTreeString(const TreeType tree_type, const bool show_branch_supports)
{
    // output the tree according to its type
    switch (tree_type)
    {
        case BIN_TREE:
            return exportTreeString(true, show_branch_supports);
            break;
        case MUL_TREE:
            return exportTreeString(false, show_branch_supports);
            break;
        default:
            outError("Unknown tree type. Please use BIN_TREE or MUL_TREE");
            break;
    }
}

const std::string cmaple::TreeBase::exportTreeString(const bool binary, const bool show_branch_supports)
{
    // make sure tree is not empty
    if (nodes.size() < 3)
        return "";
    
    return exportNodeString(binary, root_vector_index, show_branch_supports) + ";" ;
}

template <const StateType num_states>
void cmaple::TreeBase::updateMidBranchLh(const Index node_index, PhyloNode& node, const std::unique_ptr<SeqRegions>& parent_upper_regions, stack<Index> &node_stack, bool &update_blength)
{
    // update vector of regions at mid-branch point
    std::unique_ptr<SeqRegions> mid_branch_regions = nullptr;
    computeMidBranchRegions<num_states>(node, mid_branch_regions, *parent_upper_regions);
    
    if (!mid_branch_regions)
        handleNullNewRegions<num_states>(node_index, node, (node.getUpperLength() <= 1e-100), node_stack, update_blength, "inside updatePartialLh(), from parent: should not have happened since node->length > 0");
    // update likelihood at the mid-branch point
    else
        node.setMidBranchLh(std::move(mid_branch_regions));
}

template <const StateType num_states>
std::unique_ptr<SeqRegions> cmaple::TreeBase::computeUpperLeftRightRegions(const Index node_index, PhyloNode& node, const MiniIndex next_node_mini, const std::unique_ptr<SeqRegions>& parent_upper_regions, std::stack<Index> &node_stack, bool &update_blength)
{
    std::unique_ptr<SeqRegions> upper_left_right_regions = nullptr;
    const std::unique_ptr<SeqRegions>& lower_regions = getPartialLhAtNode(node.getNeighborIndex(next_node_mini));  // next_node->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob);
    // parent_upper_regions->mergeUpperLower<num_states>(upper_left_right_regions, node.getUpperLength(), *lower_regions, next_node->length, aln, model, params->threshold_prob);
    parent_upper_regions->mergeUpperLower<num_states>(upper_left_right_regions, node.getUpperLength(), *lower_regions, node.getCorrespondingLength(next_node_mini, nodes), aln, model, params->threshold_prob);
    
    // handle cases when new regions is null/empty
    if (!upper_left_right_regions || upper_left_right_regions->size() == 0)
        handleNullNewRegions<num_states>(node_index, node, (node.getUpperLength() <= 0 && node.getCorrespondingLength(next_node_mini, nodes) <= 0), node_stack, update_blength, "Strange: None vector from non-zero distances in updatePartialLh() from parent direction.");
    
    return upper_left_right_regions;
}

bool cmaple::TreeBase::updateNewPartialIfDifferent(PhyloNode& node, const MiniIndex next_node_mini, std::unique_ptr<SeqRegions>& upper_left_right_regions, std::stack<Index> &node_stack, const PositionType seq_length)
{
    // ASSERT(params.has_value());
    ASSERT(params);
    if (!node.getPartialLh(next_node_mini) || node.getPartialLh(next_node_mini)->areDiffFrom(*upper_left_right_regions, seq_length, aln->num_states, *params)) //(next_node->getPartialLhAtNode(aln, model, params->threshold_prob)->areDiffFrom(*upper_left_right_regions, seq_length, aln->num_states, &params.value()))
    {
        /*replacePartialLH(next_node->partial_lh, upper_left_right_regions);
        node_stack.push(next_node->neighbor);*/
        node.setPartialLh(next_node_mini, std::move(upper_left_right_regions));
        node_stack.push(node.getNeighborIndex(next_node_mini));
        return true;
    }
    
    return false;
}

template <const StateType num_states>
void cmaple::TreeBase::updatePartialLhFromParent(const Index index, PhyloNode& node, stack<Index> &node_stack, const std::unique_ptr<SeqRegions>& parent_upper_regions, const PositionType seq_length)
{
    bool update_blength = false;
    
    // if necessary, update the total probabilities at the mid node.
    if (node.getUpperLength() > 0)
    {
        // update vector of regions at mid-branch point
        updateMidBranchLh<num_states>(index, node, parent_upper_regions, node_stack, update_blength);
        
        // if necessary, update the total probability vector.
        if (!update_blength)
        {
            // node->computeTotalLhAtNode(aln, model, params->threshold_prob, node == root);
            node.computeTotalLhAtNode<num_states>(node.getTotalLh(), nodes[node.getNeighborIndex(TOP).getVectorIndex()], aln, model, params->threshold_prob, root_vector_index == index.getVectorIndex());
            
            if (!node.getTotalLh() || node.getTotalLh()->size() == 0)
                outError("inside updatePartialLh(), from parent 2: should not have happened since node->length > 0");
        }
    }
    
    // at valid internal node, update upLeft and upRight, and if necessary add children to node_stack.
    if (node.isInternal() && !update_blength)//(node->next && !update_blength)
    {
        /*Node* next_node_1 = node->next;
        Node* next_node_2 = next_node_1->next;*/
        
        // compute new upper left/right for next_node_1
        std::unique_ptr<SeqRegions> upper_left_right_regions_1 = computeUpperLeftRightRegions<num_states>(index, node, RIGHT, parent_upper_regions, node_stack, update_blength);//computeUpperLeftRightRegions(next_node_1, node, parent_upper_regions, node_stack, update_blength);
        std::unique_ptr<SeqRegions> upper_left_right_regions_2 = nullptr;
        
        // compute new upper left/right for next_node_1
        if (!update_blength)
            upper_left_right_regions_2 = computeUpperLeftRightRegions<num_states>(index, node, LEFT, parent_upper_regions, node_stack, update_blength);
        
        if (!update_blength)
        {
            // update new partiallh for next_node_1
            updateNewPartialIfDifferent(node, RIGHT, upper_left_right_regions_2, node_stack, seq_length);
            
            // update new partiallh for next_node_2
            updateNewPartialIfDifferent(node, LEFT, upper_left_right_regions_1, node_stack, seq_length);
        }
    }
}

template <const StateType num_states>
void cmaple::TreeBase::updatePartialLhFromChildren(const Index index, PhyloNode& node, std::stack<Index> &node_stack, const std::unique_ptr<SeqRegions>& parent_upper_regions, const bool is_non_root, const PositionType seq_length)
{
    bool update_blength = false;
    /*Node* top_node = NULL;
    Node* other_next_node = NULL;
    Node* next_node = NULL;
    FOR_NEXT(node, next_node)
    {
        if (next_node->is_top)
            top_node = next_node;
        else
            other_next_node = next_node;
    }
    
    ASSERT(top_node && other_next_node);*/
    
    const NumSeqsType node_vec_index = index.getVectorIndex();
    const Index top_node_index = Index(node_vec_index, TOP);
    const MiniIndex node_mini = index.getMiniIndex();
    const MiniIndex other_next_node_mini = index.getFlipMiniIndex();
    
    const RealNumType this_node_distance = node.getCorrespondingLength(node_mini, nodes); //node->length;
    const RealNumType other_next_node_distance = node.getCorrespondingLength(other_next_node_mini, nodes); //other_next_node->length;
    const Index neighbor_index = node.getNeighborIndex(node_mini);
    PhyloNode& neighbor = nodes[neighbor_index.getVectorIndex()];
    const std::unique_ptr<SeqRegions>& this_node_lower_regions = neighbor.getPartialLh(neighbor_index.getMiniIndex()); // node->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob);
    
    // update lower likelihoods
    std::unique_ptr<SeqRegions> merged_two_lower_regions = nullptr;
    std::unique_ptr<SeqRegions> old_lower_regions = nullptr;
    // other_next_node->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob)->mergeTwoLowers<num_states>(merged_two_lower_regions, other_next_node_distance, *this_node_lower_regions, this_node_distance, aln, model, params->threshold_prob);
    getPartialLhAtNode(node.getNeighborIndex(other_next_node_mini))->mergeTwoLowers<num_states>(merged_two_lower_regions, other_next_node_distance, *this_node_lower_regions, this_node_distance, aln, model, params->threshold_prob);
    
    if (!merged_two_lower_regions || merged_two_lower_regions->size() == 0)
    {
        //handleNullNewRegions(node->neighbor, (this_node_distance <= 0 && other_next_node_distance <= 0), node_stack, update_blength, "Strange: None vector from non-zero distances in updatePartialLh() from child direction.");
        handleNullNewRegions<num_states>(neighbor_index, neighbor, (this_node_distance <= 0 && other_next_node_distance <= 0), node_stack, update_blength, "Strange: None vector from non-zero distances in updatePartialLh() from child direction.");
    }
    else
    {
        /*replacePartialLH(old_lower_regions, top_node->partial_lh);
        top_node->partial_lh = merged_two_lower_regions;
        merged_two_lower_regions = NULL;*/
        old_lower_regions = std::move(node.getPartialLh(TOP));
        node.setPartialLh(TOP, std::move(merged_two_lower_regions));
    }

    // update total likelihood
    if (!update_blength)
    {
        if (node.getUpperLength() > 0 || root_vector_index == node_vec_index) //(top_node->length > 0 || top_node == root)
        {
            // SeqRegions* new_total_lh_regions = top_node->computeTotalLhAtNode(aln, model, params->threshold_prob, top_node == root, false);
            std::unique_ptr<SeqRegions> new_total_lh_regions = nullptr;
            node.computeTotalLhAtNode<num_states>(new_total_lh_regions, nodes[node.getNeighborIndex(TOP).getVectorIndex()], aln, model, params->threshold_prob, root_vector_index == node_vec_index);
            
            if (!new_total_lh_regions)
            {
                // handleNullNewRegions(top_node, (node.getUpperLength() <= 0), node_stack, update_blength, "Strange: None vector from non-zero distances in updatePartialLh() from child direction while doing overall likelihood.");
                handleNullNewRegions<num_states>(top_node_index, node, (node.getUpperLength() <= 0), node_stack, update_blength, "Strange: None vector from non-zero distances in updatePartialLh() from child direction while doing overall likelihood.");
            }
            else
            {
                // replacePartialLH(top_node->total_lh, new_total_lh_regions);
                node.setTotalLh(std::move(new_total_lh_regions));
            }
        }
    }
    
    // update total mid-branches likelihood
    if (!update_blength && node.getUpperLength() > 0 && is_non_root) //(!update_blength && top_node->length > 0 && is_non_root)
        //updateMidBranchLh(top_node, parent_upper_regions, node_stack, update_blength);
        updateMidBranchLh<num_states>(top_node_index, node, parent_upper_regions, node_stack, update_blength);
    
    if (!update_blength)
    {
        // update likelihoods at parent node
        // ASSERT(params.has_value());
        ASSERT(params);
        if (node.getPartialLh(TOP)->areDiffFrom(*old_lower_regions, seq_length, num_states, *params) && root_vector_index != node_vec_index) //(top_node->getPartialLhAtNode(aln, model, params->threshold_prob)->areDiffFrom(*old_lower_regions, seq_length, aln->num_states, &params.value()) && root != top_node)
            //node_stack.push(top_node->neighbor);
            node_stack.push(node.getNeighborIndex(TOP));

        // update likelihoods at sibling node
        std::unique_ptr<SeqRegions> new_upper_regions = nullptr;
        if (is_non_root)
            parent_upper_regions->mergeUpperLower<num_states>(new_upper_regions, node.getUpperLength(), *this_node_lower_regions, this_node_distance, aln, model, params->threshold_prob);
        else
            //new_upper_regions = node->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob)->computeTotalLhAtRoot(aln->num_states, model, this_node_distance);
             getPartialLhAtNode(neighbor_index)->computeTotalLhAtRoot<num_states>(new_upper_regions, model, this_node_distance);
        
        if (!new_upper_regions || new_upper_regions->size() == 0)
        {
            //handleNullNewRegions(top_node, (top_node->length <= 0 && this_node_distance <= 0), node_stack, update_blength, "Strange: None vector from non-zero distances in updatePartialLh() from child direction, new_upper_regions.");
            handleNullNewRegions<num_states>(top_node_index, node, (node.getUpperLength() <= 0 && this_node_distance <= 0), node_stack, update_blength, "Strange: None vector from non-zero distances in updatePartialLh() from child direction, new_upper_regions.");
        }
        // update partiallh of other_next_node if the new one is different from the current one
        else
            // updateNewPartialIfDifferent(other_next_node, new_upper_regions, node_stack, seq_length);
            updateNewPartialIfDifferent(node, other_next_node_mini, new_upper_regions, node_stack, seq_length);
            
        // delete new_upper_regions
       // if (new_upper_regions) delete new_upper_regions;
    }
            
    // delete old_lower_regions
    // if (old_lower_regions) delete old_lower_regions;
}

template <const StateType num_states>
void cmaple::TreeBase::updatePartialLh(stack<Index> &node_stack)
{
    const PositionType seq_length = aln->ref_seq.size();
    
    while (!node_stack.empty())
    {
        Index node_index = node_stack.top();
        node_stack.pop();
        PhyloNode& node = nodes[node_index.getVectorIndex()];
        
        // NHANLT: debug
        //if (node->next && (node->neighbor->seq_name == "25" || (node->next->neighbor && node->next->neighbor->seq_name == "25") || //(node->next->next->neighbor && node->next->next->neighbor->seq_name == "25")))
        //if (node->seq_name == "25")
         //   cout << "dsdas";
        
        node.setOutdated(true);
        
        std::unique_ptr<SeqRegions> null_seqregions_ptr = nullptr;
        bool is_non_root = root_vector_index != node_index.getVectorIndex();
        /*if (is_non_root)
        {
            //parent_upper_regions = node->getTopNode()->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob);
            SeqRegions parent_upper_regions_clone = SeqRegions(getPartialLhAtNode(node.getNeighborIndex(TOP)));
            parent_upper_regions = std::make_unique<SeqRegions>(std::move(parent_upper_regions_clone));
        }*/
        const std::unique_ptr<SeqRegions>& parent_upper_regions = is_non_root ? getPartialLhAtNode(node.getNeighborIndex(TOP)) :  null_seqregions_ptr;
            
        // change in likelihoods is coming from parent node
        if (node_index.getMiniIndex() == TOP)
        {
            ASSERT(is_non_root);
            updatePartialLhFromParent<num_states>(node_index, node, node_stack, parent_upper_regions, seq_length);
        }
        // otherwise, change in likelihoods is coming from a child.
        else
        {
            updatePartialLhFromChildren<num_states>(node_index, node, node_stack, parent_upper_regions, is_non_root, seq_length);
        }
    }
}

template <const StateType num_states>
void cmaple::TreeBase::examineSamplePlacementMidBranch(Index& selected_node_index, const std::unique_ptr<SeqRegions>& mid_branch_lh, RealNumType &best_lh_diff, bool& is_mid_branch, RealNumType& lh_diff_mid_branch, TraversingNode& current_extended_node, const std::unique_ptr<SeqRegions>& sample_regions)
{
    // compute the placement cost
    lh_diff_mid_branch = calculateSamplePlacementCost<num_states>(nodes[current_extended_node.getIndex().getVectorIndex()].getMidBranchLh(), sample_regions, default_blength);
    
    // record the best_lh_diff if lh_diff_mid_branch is greater than the best_lh_diff ever
    if (lh_diff_mid_branch > best_lh_diff)
    {
        best_lh_diff = lh_diff_mid_branch;
        selected_node_index = current_extended_node.getIndex();
        current_extended_node.setFailureCount(0);
        is_mid_branch = true;
    }
}

template <const StateType num_states>
void cmaple::TreeBase::examineSamplePlacementAtNode(Index& selected_node_index, const std::unique_ptr<SeqRegions>& total_lh, RealNumType &best_lh_diff, bool& is_mid_branch, RealNumType& lh_diff_at_node, RealNumType& lh_diff_mid_branch, RealNumType &best_up_lh_diff, RealNumType &best_down_lh_diff, Index& best_child_index, TraversingNode& current_extended_node, const std::unique_ptr<SeqRegions>& sample_regions)
{
    // compute the placement cost
    lh_diff_at_node = calculateSamplePlacementCost<num_states>(total_lh, sample_regions, default_blength);
    
    // record the best_lh_diff if lh_diff_at_node is greater than the best_lh_diff ever
    if (lh_diff_at_node > best_lh_diff)
    {
        best_lh_diff = lh_diff_at_node;
        selected_node_index = current_extended_node.getIndex();
        current_extended_node.setFailureCount(0);
        is_mid_branch = false;
        best_up_lh_diff = lh_diff_mid_branch;
    }
    else if (lh_diff_mid_branch >= (best_lh_diff - params->threshold_prob))
    {
        best_up_lh_diff = current_extended_node.getLhDiff();
        best_down_lh_diff = lh_diff_at_node;
        best_child_index = current_extended_node.getIndex();
    }
    // placement at current node is considered failed if placement likelihood is not improved by a certain margin compared to best placement so far for the nodes above it.
    else if (lh_diff_at_node < (current_extended_node.getLhDiff() - params->thresh_log_lh_failure))
        current_extended_node.increaseFailureCount();
}

template <const StateType num_states>
void cmaple::TreeBase::finetuneSamplePlacementAtNode(const PhyloNode& selected_node, RealNumType &best_down_lh_diff, Index& best_child_index, const std::unique_ptr<SeqRegions>& sample_regions)
{
    // current node might be part of a polytomy (represented by 0 branch lengths) so we want to explore all the children of the current node to find out if the best placement is actually in any of the branches below the current node.
    // Node* neighbor_node;
    stack<Index> node_stack;
    /*for (Index neighbor_index:nodes[selected_node_index.getVectorIndex()].getNeighborIndexes(selected_node_index.getMiniIndex()))
        node_stack.push(neighbor_index);*/
    // ASSERT(selected_node_index.getMiniIndex() == TOP);
    if (selected_node.isInternal())
    {
        node_stack.push(selected_node.getNeighborIndex(RIGHT));
        node_stack.push(selected_node.getNeighborIndex(LEFT));
    }
    
    
    while (!node_stack.empty())
    {
        Index node_index = node_stack.top();
        // MiniIndex node_mini_index = node_index.getMiniIndex();
        node_stack.pop();
        ASSERT(node_index.getMiniIndex() == TOP);
        PhyloNode& node = nodes[node_index.getVectorIndex()];
        // const RealNumType current_blength = node.getCorrespondingLength(node_mini_index, nodes);
        const RealNumType current_blength = node.getUpperLength();

        if (current_blength <= 0)
        {
            /*for (Index neighbor_index:node.getNeighborIndexes(node_index.getMiniIndex()))
                node_stack.push(neighbor_index);*/
            if (node.isInternal())
            {
                node_stack.push(node.getNeighborIndex(RIGHT));
                node_stack.push(node.getNeighborIndex(LEFT));
            }
        }
        else
        {
            // now try to place on the current branch below the best node, at an height above the mid-branch.
            RealNumType new_blength = current_blength * 0.5;
            RealNumType new_best_lh_mid_branch = MIN_NEGATIVE;
            // node->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob);
            // const std::unique_ptr<SeqRegions>& upper_lr_regions = getPartialLhAtNode(node.getNeighborIndex(node_mini_index));
            const std::unique_ptr<SeqRegions>& upper_lr_regions = getPartialLhAtNode(node.getNeighborIndex(TOP));
            //SeqRegions* lower_regions = node->getPartialLhAtNode(aln, model, params->threshold_prob);
            // const std::unique_ptr<SeqRegions>& lower_regions = node.getPartialLh(node_mini_index);
            const std::unique_ptr<SeqRegions>& lower_regions = node.getPartialLh(TOP);
            RealNumType new_lh_mid_branch = calculateSamplePlacementCost<num_states>(node.getMidBranchLh(), sample_regions, default_blength);
            std::unique_ptr<SeqRegions> mid_branch_regions = nullptr;

            // try to place new sample along the upper half of the current branch
            while (true)
            {
                // record new_best_lh_mid_branch
                if (new_lh_mid_branch > new_best_lh_mid_branch)
                    new_best_lh_mid_branch = new_lh_mid_branch;
                // otherwise, stop trying along the current branch
                else
                    break;
                
                // stop trying if reaching the minimum branch length
                if (new_blength <= min_blength_mid)
                    break;
             
                // try at different position along the current branch
                new_blength *= 0.5;

                // get new mid_branch_regions based on the new_blength
                upper_lr_regions->mergeUpperLower<num_states>(mid_branch_regions, new_blength, *lower_regions, current_blength - new_blength, aln, model, params->threshold_prob);
                
                // compute the placement cost
                new_lh_mid_branch = calculateSamplePlacementCost<num_states>(mid_branch_regions, sample_regions, default_blength);
            }
            
            // record new best_down_lh_diff
            if (new_best_lh_mid_branch > best_down_lh_diff)
            {
                best_down_lh_diff = new_best_lh_mid_branch;
                best_child_index = node_index;
            }
        }
    }
}

template <const StateType num_states>
void cmaple::TreeBase::seekSamplePlacement(const Index start_node_index, const NumSeqsType seq_name_index, const std::unique_ptr<SeqRegions>& sample_regions, Index& selected_node_index, RealNumType &best_lh_diff, bool &is_mid_branch, RealNumType &best_up_lh_diff, RealNumType &best_down_lh_diff, Index& best_child_index)
{
    // init variables
    // output variables
    selected_node_index = start_node_index;
    // dummy variables
    RealNumType lh_diff_mid_branch = 0;
    RealNumType lh_diff_at_node = 0;
    // stack of nodes to examine positions
    stack<TraversingNode> extended_node_stack;
    extended_node_stack.push(TraversingNode(start_node_index, 0, MIN_NEGATIVE));
    
    // recursively examine positions for placing the new sample
    while (!extended_node_stack.empty())
    {
        TraversingNode current_extended_node = std::move(extended_node_stack.top());
        extended_node_stack.pop();
        const NumSeqsType current_node_vec = current_extended_node.getIndex().getVectorIndex();
        PhyloNode& current_node = nodes[current_node_vec];
        const bool& is_internal = current_node.isInternal();
        
        // NHANLT: debug
        //if (current_node->next && ((current_node->next->neighbor && current_node->next->neighbor->seq_name == "25")
          //                        || (current_node->next->next->neighbor && current_node->next->next->neighbor->seq_name == "25")))
            //cout << "fdsfsd";
    
        // if the current node is a leaf AND the new sample/sequence is strictly less informative than the current node
        // -> add the new sequence into the list of minor sequences of the current node + stop seeking the placement
        if ((!is_internal) && (current_node.getPartialLh(TOP)->compareWithSample(*sample_regions, aln->ref_seq.size(), num_states) == 1))
        {
            current_node.addLessInfoSeqs(seq_name_index);
            selected_node_index = Index();
            return;
        }
        
        const RealNumType current_node_blength = current_node.getUpperLength();
        
        // 1. try first placing as a descendant of the mid-branch point of the branch above the current node
        if (root_vector_index != current_node_vec && current_node_blength > 0)
        {
            examineSamplePlacementMidBranch<num_states>(selected_node_index, current_node.getMidBranchLh(), best_lh_diff, is_mid_branch, lh_diff_mid_branch, current_extended_node, sample_regions);
        }
        // otherwise, don't consider mid-branch point
        else
            lh_diff_mid_branch = MIN_NEGATIVE;

        // 2. try to place as descendant of the current node (this is skipped if the node has top branch length 0 and so is part of a polytomy).
        if (root_vector_index == current_node_vec || current_node_blength > 0)
        {
            examineSamplePlacementAtNode<num_states>(selected_node_index, current_node.getTotalLh(), best_lh_diff, is_mid_branch, lh_diff_at_node, lh_diff_mid_branch, best_up_lh_diff, best_down_lh_diff, best_child_index, current_extended_node, sample_regions);
        }
        else
            lh_diff_at_node = current_extended_node.getLhDiff();
        
        // keep trying to place at children nodes, unless the number of attempts has reaches the failure limit
        const short int failure_count = current_extended_node.getFailureCount();
        if ((params->strict_stop_seeking_placement_sample
             && failure_count < params->failure_limit_sample
             && lh_diff_at_node > (best_lh_diff - params->thresh_log_lh_sample))
            || (!params->strict_stop_seeking_placement_sample
                && (failure_count < params->failure_limit_sample
                    || lh_diff_at_node > (best_lh_diff - params->thresh_log_lh_sample))))
        {
            /*for (Index neighbor_index:current_node.getNeighborIndexes(TOP))
                extended_node_stack.push(TraversingNode(neighbor_index, current_extended_node.getFailureCount(), lh_diff_at_node));*/
            if (is_internal)
            {
                extended_node_stack.push(TraversingNode(current_node.getNeighborIndex(RIGHT), failure_count, lh_diff_at_node));
                extended_node_stack.push(TraversingNode(current_node.getNeighborIndex(LEFT), failure_count, lh_diff_at_node));
            }
        }
    }

    // exploration of the tree is finished, and we are left with the node found so far with the best appending likelihood cost. Now we explore placement just below this node for more fine-grained placement within its descendant branches.
    best_down_lh_diff = MIN_NEGATIVE;
    best_child_index = Index();;
    
    // if best position so far is the descendant of a node -> explore further at its children
    if (!is_mid_branch)
    {
        finetuneSamplePlacementAtNode<num_states>(nodes[selected_node_index.getVectorIndex()], best_down_lh_diff, best_child_index, sample_regions);
    }
}

template <const StateType num_states>
void cmaple::TreeBase::addStartingNodes(const Index& node_index, PhyloNode& node, const Index& other_child_node_index, const RealNumType best_lh_diff, std::stack<std::unique_ptr<UpdatingNode>> &node_stack)
{
    // dummy variables
    const NumSeqsType vec_index = node_index.getVectorIndex();
    PhyloNode& other_child_node = nodes[other_child_node_index.getVectorIndex()];
    
    // node is not the root
    if (root_vector_index != vec_index)
    {
        const Index parent_index = node.getNeighborIndex(TOP);
        std::unique_ptr<SeqRegions>& parent_upper_lr_regions = getPartialLhAtNode(parent_index); // node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
        std::unique_ptr<SeqRegions>& other_child_node_regions = other_child_node.getPartialLh(TOP); // other_child_node->getPartialLhAtNode(aln, model, threshold_prob);
       
        // add nodes (sibling and parent of the current node) into node_stack which we will need to traverse to update their regions due to the removal of the sub tree
        RealNumType branch_length = other_child_node.getUpperLength(); // other_child_node->length;
        if (node.getUpperLength() > 0)
            branch_length = branch_length > 0 ? branch_length + node.getUpperLength() : node.getUpperLength();
        
        /*node_stack.push(new UpdatingNode(node->neighbor, other_child_node_regions, branch_length, true, best_lh_diff, 0, false));
        node_stack.push(new UpdatingNode(other_child_node, parent_upper_lr_regions, branch_length, true, best_lh_diff, 0, false));*/
        std::unique_ptr<SeqRegions> null_seqregions_ptr1 = nullptr;
        node_stack.push(std::make_unique<UpdatingNode>(UpdatingNode(parent_index, std::move(null_seqregions_ptr1), other_child_node_regions, branch_length, true, best_lh_diff, 0)));
        std::unique_ptr<SeqRegions> null_seqregions_ptr2 = nullptr;
        node_stack.push(std::make_unique<UpdatingNode>(UpdatingNode(other_child_node_index, std::move(null_seqregions_ptr2),  parent_upper_lr_regions, branch_length, true, best_lh_diff, 0)));
    }
    // node is root
    else
    {
        // there is only one sample outside of the subtree doesn't need to be considered
        if (other_child_node.isInternal()) // other_child_node->next)
        {
            // add nodes (children of the sibling of the current node) into node_stack which we will need to traverse to update their regions due to the removal of the sub tree
            const Index grand_child_1_index = other_child_node.getNeighborIndex(RIGHT);
            const Index grand_child_2_index = other_child_node.getNeighborIndex(LEFT);
            PhyloNode& grand_child_1 = nodes[grand_child_1_index.getVectorIndex()]; // other_child_node->next->neighbor;
            PhyloNode& grand_child_2 = nodes[grand_child_2_index.getVectorIndex()]; // other_child_node->next->next->neighbor;
            
            // always unique_ptr<SeqRegions>& -> always automatically delete
            // SeqRegions* up_lr_regions_1 = grand_child_2->computeTotalLhAtNode(aln, model, threshold_prob, true, false, grand_child_2->length);
            // std::unique_ptr<SeqRegions> up_lr_regions_1 = grand_child_2.computeTotalLhAtNode(other_child_node, aln, model, threshold_prob, true, grand_child_2.getUpperLength());
            std::unique_ptr<SeqRegions> up_lr_regions_1 = nullptr;
            grand_child_2.getPartialLh(TOP)->computeTotalLhAtRoot<num_states>(up_lr_regions_1, model, grand_child_2.getUpperLength());
            
            // node_stack.push(new UpdatingNode(grand_child_1, up_lr_regions_1, grand_child_1->length, true, best_lh_diff, 0, true));
            std::unique_ptr<SeqRegions> null_seqregions_ptr1 = nullptr;
            node_stack.push(std::make_unique<UpdatingNode>(UpdatingNode(grand_child_1_index, std::move(up_lr_regions_1), null_seqregions_ptr1, grand_child_1.getUpperLength(), true, best_lh_diff, 0)));
            
            
            // SeqRegions* up_lr_regions_2 = grand_child_1->computeTotalLhAtNode(aln, model, threshold_prob, true, false, grand_child_1->length);
            std::unique_ptr<SeqRegions> up_lr_regions_2 = nullptr;
            grand_child_1.getPartialLh(TOP)->computeTotalLhAtRoot<num_states>(up_lr_regions_2, model, grand_child_1.getUpperLength());
            
            // node_stack.push(new UpdatingNode(grand_child_2, up_lr_regions_2, grand_child_2->length, true, best_lh_diff, 0, true));
            std::unique_ptr<SeqRegions> null_seqregions_ptr2 = nullptr;
            node_stack.push(std::make_unique<UpdatingNode>(UpdatingNode(grand_child_2_index, std::move(up_lr_regions_2), null_seqregions_ptr2, grand_child_2.getUpperLength(), true, best_lh_diff, 0)));
        }
    }
}

// NOTE: top_node != null <=> case when crawling up from child to parent
// otherwise, top_node == null <=> case we are moving from a parent to a child
template <const StateType num_states>
bool cmaple::TreeBase::examineSubtreePlacementMidBranch(Index& best_node_index, PhyloNode& current_node, RealNumType& best_lh_diff, bool& is_mid_branch, RealNumType& lh_diff_at_node, RealNumType& lh_diff_mid_branch, RealNumType& best_up_lh_diff, RealNumType& best_down_lh_diff, std::unique_ptr<UpdatingNode>& updating_node, const std::unique_ptr<SeqRegions>& subtree_regions, const RealNumType threshold_prob, const RealNumType removed_blength, const Index top_node_index, std::unique_ptr<SeqRegions>& bottom_regions)
{
    const bool top_node_exists = (top_node_index.getMiniIndex() != UNDEFINED);
    const Index updating_node_index = updating_node->getIndex();
    const Index at_node_index = top_node_exists ? top_node_index : updating_node_index;
    const NumSeqsType at_node_vec = at_node_index.getVectorIndex();
    PhyloNode& at_node = top_node_exists ? nodes[at_node_vec] : current_node;
    
    std::unique_ptr<SeqRegions> new_mid_branch_regions = nullptr;
    // get or recompute the lh regions at the mid-branch position
    if (updating_node->needUpdate())
    {
        // recompute mid_branch_regions in case when crawling up from child to parent
        if (top_node_exists)
        {
            /*Node* other_child = updating_node->node->getOtherNextNode()->neighbor;
            SeqRegions* other_child_lower_regions = other_child->getPartialLhAtNode(aln, model, threshold_prob);
            other_child_lower_regions->mergeTwoLowers<num_states>(bottom_regions, other_child->length, *updating_node->incoming_regions, updating_node->branch_length, aln, model, threshold_prob);*/
            const Index other_child_index = current_node.getNeighborIndex(updating_node_index.getFlipMiniIndex());
            PhyloNode& other_child = nodes[other_child_index.getVectorIndex()]; // updating_node->node->getOtherNextNode()->neighbor;
            const std::unique_ptr<SeqRegions>& other_child_lower_regions = other_child.getPartialLh(TOP); // other_child->getPartialLhAtNode(aln, model, threshold_prob);
            // other_child_lower_regions->mergeTwoLowers<num_states>(bottom_regions, other_child->length, *updating_node->incoming_regions, updating_node->branch_length, aln, model, threshold_prob);
            other_child_lower_regions->mergeTwoLowers<num_states>(bottom_regions, other_child.getUpperLength(), *updating_node->getIncomingRegions(), updating_node->getBranchLength(), aln, model, threshold_prob);
                                
            // skip if bottom_regions is null (inconsistent)
            if (!bottom_regions)
            {
                // delete updating_node
                // delete updating_node;

                return false; // continue;
            }
           
            // compute new mid-branch regions
            /*SeqRegions* upper_lr_regions = top_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
            RealNumType mid_branch_length = top_node->length * 0.5;
            upper_lr_regions->mergeUpperLower<num_states>(mid_branch_regions, mid_branch_length, *bottom_regions, mid_branch_length, aln, model, threshold_prob);*/
            
            const std::unique_ptr<SeqRegions>& upper_lr_regions = getPartialLhAtNode(at_node.getNeighborIndex(TOP)); // top_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
            const RealNumType mid_branch_length = at_node.getUpperLength() * 0.5;
            upper_lr_regions->mergeUpperLower<num_states>(new_mid_branch_regions, mid_branch_length, *bottom_regions, mid_branch_length, aln, model, threshold_prob);
        }
        // recompute mid_branch_regions in case we are moving from a parent to a child
        else
        {
            /*SeqRegions* lower_regions = updating_node->node->getPartialLhAtNode(aln, model, threshold_prob);
            RealNumType mid_branch_length = updating_node->branch_length * 0.5;
            updating_node->incoming_regions->mergeUpperLower<num_states>(new_mid_branch_regions, mid_branch_length, *lower_regions, mid_branch_length, aln, model, threshold_prob);*/
            const std::unique_ptr<SeqRegions>& lower_regions = current_node.getPartialLh(TOP); // getPartialLhAtNode(updating_node_index);
            const RealNumType mid_branch_length = updating_node->getBranchLength() * 0.5;
            updating_node->getIncomingRegions()->mergeUpperLower<num_states>(new_mid_branch_regions, mid_branch_length, *lower_regions, mid_branch_length, aln, model, threshold_prob);
        }
    }
    
    std::unique_ptr<SeqRegions>& mid_branch_regions = updating_node->needUpdate() ? new_mid_branch_regions : at_node.getMidBranchLh();
    
    // skip if mid_branch_regions is null (branch length == 0)
    if (!mid_branch_regions)
    {
        // delete bottom_regions if it's existed
        // if (bottom_regions) delete bottom_regions;
                            
        // delete updating_node
        // delete updating_node;
        
        return false; // continue;
    }
    
    // compute the placement cost
    // if (search_subtree_placement)
    lh_diff_mid_branch = calculateSubTreePlacementCost<num_states>(mid_branch_regions, subtree_regions, removed_blength);
    //else
       // lh_diff_mid_branch = calculateSamplePlacementCost( mid_branch_regions, subtree_regions, removed_blength);
    
    if (top_node_exists && best_node_index.getVectorIndex() == at_node_vec)// top_node && best_node == top_node) // only update in case when crawling up from child to parent
        best_up_lh_diff = lh_diff_mid_branch;
    
    // if this position is better than the best position found so far -> record it
    if (lh_diff_mid_branch > best_lh_diff)
    {
        best_node_index = at_node_index;
        best_lh_diff = lh_diff_mid_branch;
        is_mid_branch = true;
        updating_node->setFailureCount(0);
        if (top_node_exists) best_down_lh_diff = lh_diff_at_node; // only update in case when crawling up from child to parent
    }
    else if (top_node_exists && lh_diff_at_node >= (best_lh_diff - threshold_prob))
        best_up_lh_diff = lh_diff_mid_branch;
        
    // delete mid_branch_regions
    // if (updating_node->need_updating) delete mid_branch_regions;
    
    // no error
    return true;
}

// NOTE: top_node != null <=> case when crawling up from child to parent
// otherwise, top_node == null <=> case we are moving from a parent to a child
template <const StateType num_states>
bool cmaple::TreeBase::examineSubTreePlacementAtNode(Index& best_node_index, PhyloNode& current_node, RealNumType &best_lh_diff, bool& is_mid_branch, RealNumType& lh_diff_at_node, RealNumType& lh_diff_mid_branch, RealNumType &best_up_lh_diff, RealNumType &best_down_lh_diff, std::unique_ptr<UpdatingNode>& updating_node, const std::unique_ptr<SeqRegions>& subtree_regions, const RealNumType threshold_prob, const RealNumType removed_blength, const Index top_node_index)
{
    const PositionType seq_length = aln->ref_seq.size();
    
    // Node* at_node = top_node ? top_node: updating_node->node;
    const bool top_node_exits = top_node_index.getMiniIndex() != UNDEFINED;
    const Index updating_node_index = updating_node->getIndex();
    const Index at_node_index = top_node_exits ? top_node_index : updating_node_index;
    PhyloNode& at_node = top_node_exits ? nodes[at_node_index.getVectorIndex()] : current_node;
    
    std::unique_ptr<SeqRegions> new_at_node_regions = nullptr;
    const bool need_updating = updating_node->needUpdate();
    if (updating_node->needUpdate())
    {
        // get or recompute the lh regions at the current node position
        const std::unique_ptr<SeqRegions>& updating_node_partial = getPartialLhAtNode(updating_node_index); // updating_node->node->getPartialLhAtNode(aln, model, threshold_prob);
        if (top_node_exits)
        {
            updating_node_partial->mergeUpperLower<num_states>(new_at_node_regions, -1, *updating_node->getIncomingRegions(), updating_node->getBranchLength(), aln, model, threshold_prob);
        }
        else
        {
            updating_node->getIncomingRegions()->mergeUpperLower<num_states>(new_at_node_regions, updating_node->getBranchLength(), *updating_node_partial, -1, aln, model, threshold_prob);
        }
        
        // skip if at_node_regions is null (branch length == 0)
        if (!new_at_node_regions)
        {
            // delete updating_node
            // delete updating_node;
            
            // continue;
            return false;
        }
        
        // stop updating if the difference between the new and old regions is insignificant
        // ASSERT(params.has_value());
        ASSERT(params);
        if  (!new_at_node_regions->areDiffFrom(*at_node.getTotalLh(), seq_length, num_states, *params))
            updating_node->setUpdate(false);
    }
    // else
    const std::unique_ptr<SeqRegions>& at_node_regions = need_updating ? new_at_node_regions : at_node.getTotalLh();
    
    //if (search_subtree_placement)
    lh_diff_at_node = calculateSubTreePlacementCost<num_states>(at_node_regions, subtree_regions, removed_blength);
    //else
        //lh_diff_at_node = calculateSamplePlacementCost(at_node_regions, subtree_regions, removed_blength);
    
    // if this position is better than the best position found so far -> record it
    if (lh_diff_at_node > best_lh_diff)
    {
        best_node_index = at_node_index;
        best_lh_diff = lh_diff_at_node;
        is_mid_branch = false;
        updating_node->setFailureCount(0);
        if (!top_node_exits) best_up_lh_diff = lh_diff_mid_branch; // only update in case we are moving from a parent to a child
    }
    else if (!top_node_exits && lh_diff_mid_branch >= (best_lh_diff - threshold_prob)) // only update in case we are moving from a parent to a child
    {
        best_up_lh_diff = updating_node->getLhDiff();
        best_down_lh_diff = lh_diff_at_node;
    }
    // placement at current node is considered failed if placement likelihood is not improved by a certain margin compared to best placement so far for the nodes above it.
    else if (lh_diff_at_node < (updating_node->getLhDiff() - params->thresh_log_lh_failure))
        updating_node->increaseFailureCount();
        
    // delete at_node_regions
    // if (delete_at_node_regions) delete at_node_regions;
    
    // no error
    return true;
}

bool keepTraversing(const RealNumType& best_lh_diff, const RealNumType& lh_diff_at_node, const bool& strict_stop_seeking_placement_subtree, const std::unique_ptr<UpdatingNode>& updating_node, const int& failure_limit_subtree, const RealNumType& thresh_log_lh_subtree, const bool able2traverse = true)
{
    //if (search_subtree_placement)
    //{
    if (strict_stop_seeking_placement_subtree)
    {
        if (updating_node->getFailureCount() <= failure_limit_subtree && lh_diff_at_node > (best_lh_diff - thresh_log_lh_subtree) && able2traverse)
            return true;
    }
    else
    {
        if ((updating_node->getFailureCount() <= failure_limit_subtree || lh_diff_at_node > (best_lh_diff - thresh_log_lh_subtree))
            && able2traverse)
            return true;
    }
    //}
    //else
    //{
      //  if (params->strict_stop_seeking_placement_sample)
        //{
          //  if (updating_node->failure_count <= params->failure_limit_sample && lh_diff_at_node > (best_lh_diff - params->thresh_log_lh_sample)
            //    && updating_node->node->next)
              //      return true;
       // }
       // else
        //{
          //  if ((updating_node->failure_count <= params->failure_limit_sample || lh_diff_at_node > (best_lh_diff - params->thresh_log_lh_sample))
            //    && updating_node->node->next)
              //      return true;
        //}
    //}
    
    // default
    return false;
}

template <const StateType num_states>
void cmaple::TreeBase::addChildSeekSubtreePlacement(const Index child_1_index, const Index child_2_index, PhyloNode& child_1, PhyloNode& child_2, const RealNumType& lh_diff_at_node, const std::unique_ptr<UpdatingNode>& updating_node, std::stack<std::unique_ptr<UpdatingNode>>& node_stack, const RealNumType threshold_prob)
{
    // get or recompute the upper left/right regions of the children node
    if (updating_node->needUpdate())
    {
        std::unique_ptr<SeqRegions> upper_lr_regions = nullptr;
        const std::unique_ptr<SeqRegions>& lower_regions = child_2.getPartialLh(TOP); // ->getPartialLhAtNode(aln, model, threshold_prob);
        
        updating_node->getIncomingRegions()->mergeUpperLower<num_states>(upper_lr_regions, updating_node->getBranchLength(), *lower_regions, child_2.getUpperLength(), aln, model, threshold_prob);
        
        // traverse to this child's subtree
        if (upper_lr_regions)
        {
            // node_stack.push(new UpdatingNode(child_1, upper_lr_regions, child_1->length, updating_node->need_updating, lh_diff_at_node, updating_node->failure_count, updating_node->need_updating));
            
            std::unique_ptr<SeqRegions> null_seqregions_ptr = nullptr;
            node_stack.push(std::make_unique<UpdatingNode>(UpdatingNode(child_1_index, std::move(upper_lr_regions), null_seqregions_ptr, child_1.getUpperLength(), updating_node->needUpdate(), lh_diff_at_node, updating_node->getFailureCount())));
        }
    }
    else
    {
        std::unique_ptr<SeqRegions>& upper_lr_regions = getPartialLhAtNode(child_1.getNeighborIndex(TOP));
        const std::unique_ptr<SeqRegions> null_seqregions_ptr = nullptr;
        /*if (child_1->neighbor->partial_lh)
            upper_lr_regions = child_1->neighbor->getPartialLhAtNode(aln, model, threshold_prob);*/
        
        // traverse to this child's subtree
        if (upper_lr_regions)
        {
            // node_stack.push(new UpdatingNode(child_1, upper_lr_regions, child_1->length, updating_node->need_updating, lh_diff_at_node, updating_node->failure_count, updating_node->need_updating));
            std::unique_ptr<SeqRegions> null_seqregions_ptr = nullptr;
            node_stack.push(std::make_unique<UpdatingNode>(UpdatingNode(child_1_index, std::move(null_seqregions_ptr), upper_lr_regions, child_1.getUpperLength(), updating_node->needUpdate(), lh_diff_at_node, updating_node->getFailureCount())));
        }
    }
}

template <const StateType num_states>
bool cmaple::TreeBase::addNeighborsSeekSubtreePlacement(PhyloNode& current_node, const Index other_child_index, std::unique_ptr<SeqRegions>&& bottom_regions, const RealNumType& lh_diff_at_node, const std::unique_ptr<UpdatingNode>& updating_node, std::stack<std::unique_ptr<UpdatingNode>>& node_stack, const RealNumType threshold_prob)
{
    ASSERT(other_child_index.getMiniIndex() == TOP);
    PhyloNode& other_child = nodes[other_child_index.getVectorIndex()];
    const Index updating_node_index = updating_node->getIndex();
    const MiniIndex updating_node_mini = updating_node_index.getMiniIndex();
    
    // keep crawling up into parent and sibling node
    // case the node is not the root
    if (root_vector_index != updating_node_index.getVectorIndex())
    {
        // first pass the crawling down the other child (sibling)
        
        // get or recompute the upper left/right regions of the sibling node
        if (updating_node->needUpdate())
        {
            const std::unique_ptr<SeqRegions>& parent_upper_lr_regions = getPartialLhAtNode(current_node.getNeighborIndex(TOP)); // top_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
            std::unique_ptr<SeqRegions> upper_lr_regions = nullptr;
            
            parent_upper_lr_regions->mergeUpperLower<num_states>(upper_lr_regions, current_node.getUpperLength(), *updating_node->getIncomingRegions(), updating_node->getBranchLength(), aln, model, threshold_prob);
            
            if (!upper_lr_regions)
            {
                // delete bottom_regions if it's existed
                // if (bottom_regions) delete bottom_regions;
                
                // delete updating_node
                // delete updating_node;
                
                return false; // continue;
            }
            else
            {
                std::unique_ptr<SeqRegions> null_seqregions_ptr = nullptr;
                //node_stack.push(std::make_unique<UpdatingNode>(UpdatingNode(other_child_index, std::move(upper_lr_regions), null_seqregions_ptr, other_child->length, updating_node->need_updating_, lh_diff_at_node, updating_node->failure_count_)));
                node_stack.push(std::make_unique<UpdatingNode>(UpdatingNode(other_child_index, std::move(upper_lr_regions), null_seqregions_ptr, other_child.getUpperLength(), updating_node->needUpdate(), lh_diff_at_node, updating_node->getFailureCount())));
            }
        }
        else
        {
            std::unique_ptr<SeqRegions>& upper_lr_regions = current_node.getPartialLh(updating_node_mini);
            
            if (!upper_lr_regions)//updating_node->node->partial_lh)
            {
                // delete bottom_regions if it's existed
                // if (bottom_regions) delete bottom_regions;
                
                // delete updating_node
                // delete updating_node;
                
                return false; // continue;
            }
            else
            {
                /*upper_lr_regions = updating_node->node->getPartialLhAtNode(aln, model, threshold_prob);
                node_stack.push(new UpdatingNode(other_child, upper_lr_regions, other_child->length, updating_node->need_updating, lh_diff_at_node, updating_node->failure_count, updating_node->need_updating));*/
                std::unique_ptr<SeqRegions> null_seqregions_ptr = nullptr;
                node_stack.push(std::make_unique<UpdatingNode>(UpdatingNode(other_child_index, std::move(null_seqregions_ptr), upper_lr_regions, other_child.getUpperLength(), updating_node->needUpdate(), lh_diff_at_node, updating_node->getFailureCount())));
            }
        }

        // add sibling node to node_stack for traversing later; skip if upper_lr_regions is null (inconsistent)
        /*if (!upper_lr_regions)
        {
            // delete bottom_regions if it's existed
            // if (bottom_regions) delete bottom_regions;
            
            // delete updating_node
            // delete updating_node;
            
            return false; // continue;
        }
        else
            node_stack.push(new UpdatingNode(other_child, upper_lr_regions, other_child->length, updating_node->need_updating, lh_diff_at_node, updating_node->failure_count, updating_node->need_updating));*/
        
        // now pass the crawling up to the parent node
        // get or recompute the bottom regions (comming from 2 children) of the parent node
        if (updating_node->needUpdate())
        {
            if (!bottom_regions)
            {
                const std::unique_ptr<SeqRegions>& other_child_lower_regions = other_child.getPartialLh(TOP); // other_child->getPartialLhAtNode(aln, model, threshold_prob);
                // other_child_lower_regions->mergeTwoLowers<num_states>(bottom_regions, other_child->length, *updating_node->incoming_regions, updating_node->branch_length_, aln, model, threshold_prob);
                other_child_lower_regions->mergeTwoLowers<num_states>(bottom_regions, other_child.getUpperLength(), *updating_node->getIncomingRegions(), updating_node->getBranchLength(), aln, model, threshold_prob);
                
                // skip if bottom_regions is null (inconsistent)
                if (!bottom_regions)
                {
                    // delete updating_node
                    // delete updating_node;
                    
                    return false; // continue;
                }
            }
            
            // node_stack.push(std::make_unique<UpdatingNode>(UpdatingNode(top_node->neighbor, bottom_regions, top_node->length, updating_node->need_updating, lh_diff_at_node, updating_node->failure_count, updating_node->need_updating)));
            std::unique_ptr<SeqRegions> null_seqregions_ptr = nullptr;
            node_stack.push(std::make_unique<UpdatingNode>(UpdatingNode(current_node.getNeighborIndex(TOP), std::move(bottom_regions), null_seqregions_ptr, current_node.getUpperLength(), updating_node->needUpdate(), lh_diff_at_node, updating_node->getFailureCount())));
        }
        else
        {
            // if (bottom_regions) delete bottom_regions;
            bottom_regions = nullptr;
            std::unique_ptr<SeqRegions>& bottom_regions_ref = current_node.getPartialLh(TOP); // top_node->getPartialLhAtNode(aln, model, threshold_prob);
            
            std::unique_ptr<SeqRegions> null_seqregions_ptr = nullptr;
            // node_stack.push(std::make_unique<UpdatingNode>(UpdatingNode(top_node->neighbor, bottom_regions, top_node->length, updating_node->need_updating, lh_diff_at_node, updating_node->failure_count, updating_node->need_updating)));
            node_stack.push(std::make_unique<UpdatingNode>(UpdatingNode(current_node.getNeighborIndex(TOP), std::move(null_seqregions_ptr), bottom_regions_ref, current_node.getUpperLength(), updating_node->needUpdate(), lh_diff_at_node, updating_node->getFailureCount())));
        }
        
        // add the parent node to node_stack for traversing later
        /* node_stack.push(new UpdatingNode(top_node->neighbor, bottom_regions, top_node->length, updating_node->need_updating, lh_diff_at_node, updating_node->failure_count, updating_node->need_updating));*/
    }
    // now consider case of root node -> only need to care about the sibling node
    else
    {
        // get or recompute the upper left/right regions of the sibling node
        if (updating_node->needUpdate())
        {
            std::unique_ptr<SeqRegions> upper_lr_regions = nullptr;
            updating_node->getIncomingRegions()->computeTotalLhAtRoot<num_states>(upper_lr_regions, model, updating_node->getBranchLength());
            
            std::unique_ptr<SeqRegions> null_seqregions_ptr = nullptr;
            node_stack.push(std::make_unique<UpdatingNode>(UpdatingNode(other_child_index, std::move(upper_lr_regions), null_seqregions_ptr, other_child.getUpperLength(), updating_node->needUpdate(), lh_diff_at_node, updating_node->getFailureCount())));
        }
        else
        {
            std::unique_ptr<SeqRegions>& upper_lr_regions = current_node.getPartialLh(updating_node_mini); // updating_node->node->getPartialLhAtNode(aln, model, threshold_prob);
            
            std::unique_ptr<SeqRegions> null_seqregions_ptr = nullptr;
            node_stack.push(std::make_unique<UpdatingNode>(UpdatingNode(other_child_index, std::move(null_seqregions_ptr), upper_lr_regions, other_child.getUpperLength(), updating_node->needUpdate(), lh_diff_at_node, updating_node->getFailureCount())));
        }
        
        // add the sibling node to node_stack for traversing later
        // node_stack.push(new UpdatingNode(other_child, upper_lr_regions, other_child->length, updating_node->need_updating, lh_diff_at_node, updating_node->failure_count, updating_node->need_updating));
        
        // delete bottom_regions
        // if (bottom_regions) delete bottom_regions;
    }
    
    return true;
}

template <const StateType num_states>
void cmaple::TreeBase::seekSubTreePlacement(Index& best_node_index, RealNumType &best_lh_diff, bool &is_mid_branch, RealNumType &best_up_lh_diff, RealNumType &best_down_lh_diff, Index& best_child_index, const bool short_range_search, const Index child_node_index, RealNumType &removed_blength) //, bool search_subtree_placement, SeqRegions* sample_regions)
{
    // init variables
    PhyloNode& child_node = nodes[child_node_index.getVectorIndex()];
    const Index node_index = child_node.getNeighborIndex(TOP);
    const NumSeqsType vec_index = node_index.getVectorIndex();
    PhyloNode& node = nodes[vec_index]; // child_node->neighbor->getTopNode();
    const Index other_child_node_index = node.getNeighborIndex(node_index.getFlipMiniIndex()); // child_node->neighbor->getOtherNextNode()->neighbor;
    best_node_index = node_index;
    const std::unique_ptr<SeqRegions>& subtree_regions = child_node.getPartialLh(TOP); // child_node->getPartialLhAtNode(aln, model, threshold_prob); nullptr;
    // stack of nodes to examine positions
    stack<std::unique_ptr<UpdatingNode>> node_stack;
    // dummy variables
    const RealNumType threshold_prob = params->threshold_prob;
    RealNumType lh_diff_mid_branch = 0;
    RealNumType lh_diff_at_node = 0;
    // const std::unique_ptr<SeqRegions> null_seqregions_ptr = nullptr;
    // const std::unique_ptr<SeqRegions>& parent_upper_lr_regions = root_vector_index == vec_index ? null_seqregions_ptr : getPartialLhAtNode(node.getNeighborIndex(TOP));
    const PositionType seq_length = aln->ref_seq.size();
    
    // get/init approximation params
    bool strict_stop_seeking_placement_subtree = params->strict_stop_seeking_placement_subtree;
    int failure_limit_subtree = params->failure_limit_subtree;
    RealNumType thresh_log_lh_subtree = params->thresh_log_lh_subtree;
    
    // for short range topology search
    if (short_range_search)
    {
        strict_stop_seeking_placement_subtree = params->strict_stop_seeking_placement_subtree_short_search;
        failure_limit_subtree = params->failure_limit_subtree_short_search;
        thresh_log_lh_subtree = params->thresh_log_lh_subtree_short_search;
    }

    // search a placement for a subtree
    //if (search_subtree_placement)
    //{
    // get the lower regions of the child node
    // subtree_regions = child_node->getPartialLhAtNode(aln, model, threshold_prob);
    
    // add starting nodes to start seek placement for the subtree
    addStartingNodes<num_states>(node_index, node, other_child_node_index, best_lh_diff, node_stack);
    
    //}
     // search a placement for a new sample
    //else
    //{
        // get the regions of the input sample
      //  subtree_regions = sample_regions;
        //RealNumType down_lh = is_mid_branch ? best_down_lh_diff : best_lh_diff;
        
        // node is not the root
  //      if (node != root)
    //    {
      //      SeqRegions* lower_regions = new SeqRegions(node->getPartialLhAtNode(aln, model, threshold_prob), num_states);
        //    parent_upper_lr_regions = node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
            
            // add the parent node of the current node into node_stack for traversing to seek the placement for the new sample
    //        node_stack.push(new UpdatingNode(node->neighbor, lower_regions, node->length, false, down_lh, 0));
      //  }
        
        // add the children nodes of the current node into node_stack for traversing to seek the placement for the new sample
    //    Node* neighbor_node;
      //  FOR_NEIGHBOR(node, neighbor_node)
        //{
         //   SeqRegions* upper_lr_regions = new SeqRegions(neighbor_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob), num_states);
          //  node_stack.push(new UpdatingNode(neighbor_node, upper_lr_regions, neighbor_node->length, false, down_lh, 0));
        //}
     //}
    
    // examine each node in the node stack to seek the "best" placement
    while (!node_stack.empty())
    {
        // extract updating_node from stack
        std::unique_ptr<UpdatingNode> updating_node = std::move(node_stack.top());
        node_stack.pop();
        const Index current_node_index = updating_node->getIndex();
        const NumSeqsType current_node_vec = current_node_index.getVectorIndex();
        PhyloNode& current_node = nodes[current_node_vec];
        
        // consider the case we are moving from a parent to a child
        if (current_node_index.getMiniIndex() == TOP)
        {
            if (current_node.getUpperLength() > 0) // updating_node->node->length > 0)
            {
                //  try to append mid-branch
                // avoid adding to the old position where the subtree was just been removed from
                if (root_vector_index != current_node_vec && current_node.getNeighborIndex(TOP).getVectorIndex() != vec_index) // updating_node->node != root && updating_node->node->neighbor->getTopNode() != node)
                {
                    std::unique_ptr<SeqRegions> bottom_regions = nullptr;
                    if (!examineSubtreePlacementMidBranch<num_states>(best_node_index, current_node, best_lh_diff, is_mid_branch, lh_diff_at_node, lh_diff_mid_branch, best_up_lh_diff, best_down_lh_diff, updating_node, subtree_regions, threshold_prob, removed_blength, Index(), bottom_regions)) continue;
                }
                // set the placement cost at the mid-branch position the most negative value if branch length is zero -> we can't place the subtree on that branch
                else
                    lh_diff_mid_branch = MIN_NEGATIVE;
                    
                // now try appending exactly at node
                if(!examineSubTreePlacementAtNode<num_states>(best_node_index, current_node, best_lh_diff, is_mid_branch, lh_diff_at_node, lh_diff_mid_branch, best_up_lh_diff, best_down_lh_diff, updating_node, subtree_regions, threshold_prob, removed_blength, Index())) continue;
            }
            // set the placement cost at the current node position at the most negative value if branch length is zero -> we can't place the subtree on that branch
            else
                lh_diff_at_node = updating_node->getLhDiff();
            
            // keep crawling down into children nodes unless the stop criteria for the traversal are satisfied.
            // check the stop criteria
            // keep traversing further down to the children
            if (keepTraversing(best_lh_diff, lh_diff_at_node, strict_stop_seeking_placement_subtree, updating_node, failure_limit_subtree, thresh_log_lh_subtree, current_node.isInternal()))// updating_node->node->next))
            {
                /*Node* child_1 = updating_node->node->getOtherNextNode()->neighbor;
                Node* child_2 = child_1->neighbor->getOtherNextNode()->neighbor;*/
                const Index child_1_index = current_node.getNeighborIndex(RIGHT); // updating_node->node->getOtherNextNode()->neighbor;
                const Index child_2_index = current_node.getNeighborIndex(LEFT); // child_1->neighbor->getOtherNextNode()->neighbor;
                PhyloNode& child_1 = nodes[child_1_index.getVectorIndex()];
                PhyloNode& child_2 = nodes[child_2_index.getVectorIndex()];
                
                // add child_1 to node_stack
                addChildSeekSubtreePlacement<num_states>(child_1_index, child_2_index, child_1, child_2, lh_diff_at_node, updating_node, node_stack, threshold_prob);
                
                // add child_2 to node_stack
                addChildSeekSubtreePlacement<num_states>(child_2_index, child_1_index, child_2, child_1, lh_diff_at_node, updating_node, node_stack, threshold_prob);
            }
        }
        // case when crawling up from child to parent
        else
        {
            // Node* top_node = updating_node->node->getTopNode();
            const Index top_node_index = Index(current_node_vec, TOP);
            
            // append directly at the node
            if (current_node.getUpperLength() > 0 || root_vector_index == current_node_vec) // top_node->length > 0 || top_node == root)
            {
                if (!examineSubTreePlacementAtNode<num_states>(best_node_index, current_node, best_lh_diff, is_mid_branch, lh_diff_at_node, lh_diff_mid_branch, best_up_lh_diff, best_down_lh_diff, updating_node, subtree_regions, threshold_prob, removed_blength, top_node_index)) continue;
            }
            // if placement cost at new position gets worse -> restore to the old one
            else
                lh_diff_at_node = updating_node->getLhDiff();

            // try appending mid-branch
            const Index other_child_index = current_node.getNeighborIndex(current_node_index.getFlipMiniIndex()); // updating_node->node->getOtherNextNode()->neighbor;
            std::unique_ptr<SeqRegions> bottom_regions = nullptr;
            if (current_node.getUpperLength() > 0 && root_vector_index != current_node_vec) // top_node->length > 0 && top_node != root)
            {
                if (!examineSubtreePlacementMidBranch<num_states>(best_node_index, current_node, best_lh_diff, is_mid_branch, lh_diff_at_node, lh_diff_mid_branch, best_up_lh_diff, best_down_lh_diff, updating_node, subtree_regions, threshold_prob, removed_blength, top_node_index, bottom_regions)) continue;
            }
            // set the placement cost at the mid-branch position at the most negative value if branch length is zero -> we can't place the subtree on that branch
            // NHANLT: we actually don't need to do that since lh_diff_mid_branch will never be read
            // else
                // lh_diff_mid_branch = MIN_NEGATIVE;
            
            // check stop rule of the traversal process
            // keep traversing upwards
            if (keepTraversing(best_lh_diff, lh_diff_at_node, strict_stop_seeking_placement_subtree, updating_node, failure_limit_subtree, thresh_log_lh_subtree))
            {
                // if(!addNeighborsSeekSubtreePlacement(top_node_index, other_child_index, parent_upper_lr_regions, bottom_regions, lh_diff_at_node, updating_node, node_stack, threshold_prob)) continue;
                if(!addNeighborsSeekSubtreePlacement<num_states>(current_node, other_child_index, std::move(bottom_regions), lh_diff_at_node, updating_node, node_stack, threshold_prob)) continue;
            }
            /*else
            {
                // delete bottom_regions if it's existed
                if (bottom_regions) delete bottom_regions;
            }*/
        }
        
        // delete updating_node
        // delete updating_node;
    }
    
    // ############ KEEP this section DISABLE/COMMENTED OUT ############
    // exploration of the tree is finished, and we are left with the node found so far with the best appending likelihood cost.
    // Now we explore placement just below this node for more fine-grained placement within its descendant branches.
    /*if (!search_subtree_placement)
    {
        best_down_lh_diff = MIN_NEGATIVE;
        best_child = NULL;
        
        if (is_mid_branch)
        {
            // go upward until we reach the parent node of a polytomy
            Node* parent_node = best_node->neighbor->getTopNode();
            while (parent_node->length <= 0 && parent_node != root)
                parent_node = parent_node->neighbor->getTopNode();
            
            best_up_lh_diff = calculateSamplePlacementCost(parent_node->total_lh, subtree_regions, removed_blength);
            child_node = best_node;
        }
        else
        {
            // current node might be part of a polytomy (represented by 0 branch lengths) so we want to explore all the children of the current node to find out if the best placement is actually in any of the branches below the current node.
            Node* neighbor_node;
            stack<Node*> new_node_stack;
            FOR_NEIGHBOR(best_node, neighbor_node)
                new_node_stack.push(neighbor_node);
            
            while (!new_node_stack.empty())
            {
                Node* node = new_node_stack.top();
                new_node_stack.pop();
                
                if (node->length <= 0)
                {
                    FOR_NEIGHBOR(node, neighbor_node)
                        new_node_stack.push(neighbor_node);
                }
                // now try to place on the current branch below the best node, at an height above the mid-branch.
                else
                {
                    // now try to place on the current branch below the best node, at an height above the mid-branch.
                    RealNumType new_blength = node->length * 0.5;
                    RealNumType new_best_lh_mid_branch = MIN_NEGATIVE;
                    SeqRegions* upper_lr_regions = node->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob);
                    SeqRegions* lower_regions = node->getPartialLhAtNode(aln, model, params->threshold_prob);
                    SeqRegions* mid_branch_regions = new SeqRegions(node->mid_branch_lh, aln->num_states);

                    // try to place new sample along the upper half of the current branch
                    while (true)
                    {
                        // compute the placement cost
                        RealNumType new_lh_mid_branch = calculateSamplePlacementCost(mid_branch_regions, subtree_regions, removed_blength);
                        
                        // record new_best_lh_mid_branch
                        if (new_lh_mid_branch > new_best_lh_mid_branch)
                            new_best_lh_mid_branch = new_lh_mid_branch;
                        // otherwise, stop trying along the current branch
                        else
                            break;
                        
                        // stop trying if reaching the minimum branch length
                        if (new_blength <= min_blength_mid)
                            break;
     
                         // try at different position along the current branch
                         new_blength *= 0.5;
                     
                        // get new mid_branch_regions based on the new_blength
                        upper_lr_regions->mergeUpperLower<num_states>(mid_branch_regions, new_blength, lower_regions, node->length - new_blength, aln, model, params->threshold_prob);
                    }
                    
                    //RealNumType new_best_lh_mid_branch = calculateSamplePlacementCost(node->mid_branch_lh, sample_regions, default_blength);
                    
                    // record new best_down_lh_diff
                    if (new_best_lh_mid_branch > best_down_lh_diff)
                    {
                        best_down_lh_diff = new_best_lh_mid_branch;
                        best_child = node;
                    }
                }
            }
        }
    }*/
    // ############ KEEP this section DISABLE/COMMENTED OUT ############
}

template <const StateType num_states>
void cmaple::TreeBase::applySPR(const Index subtree_index, PhyloNode& subtree, const Index best_node_index, const bool is_mid_branch, const RealNumType branch_length, const RealNumType best_lh_diff)
{
    // record the SPR applied at this subtree
    subtree.setSPRApplied(true);
    // remove subtree from the tree
    const Index parent_index = subtree.getNeighborIndex(TOP);
    PhyloNode& parent_subtree = nodes[parent_index.getVectorIndex()]; // subtree->neighbor->getTopNode();
    const Index sibling_index = parent_subtree.getNeighborIndex(parent_index.getFlipMiniIndex());
    PhyloNode& sibling_subtree = nodes[sibling_index.getVectorIndex()]; // subtree->neighbor->getOtherNextNode()->neighbor;
    const RealNumType threshold_prob = params->threshold_prob;
    
    // connect grandparent to sibling
    const Index grandparent_index = parent_subtree.getNeighborIndex(TOP);
    if (root_vector_index != parent_index.getVectorIndex())
    {
        // parent_subtree->neighbor->neighbor = sibling_subtree;
        nodes[grandparent_index.getVectorIndex()].setNeighborIndex(grandparent_index.getMiniIndex(), sibling_index);
    }
    // sibling_subtree->neighbor = parent_subtree->neighbor;
    sibling_subtree.setNeighborIndex(TOP, grandparent_index);
    
    // update the length of the branch connecting grandparent to sibling
    if (sibling_subtree.getUpperLength() > 0)// sibling_subtree->length > 0)
    {
        if (parent_subtree.getUpperLength() > 0)
            sibling_subtree.setUpperLength(sibling_subtree.getUpperLength() + parent_subtree.getUpperLength()); // sibling_subtree->length += parent_subtree->length;
    }
    else
        sibling_subtree.setUpperLength(parent_subtree.getUpperLength()); // sibling_subtree->length = parent_subtree->length;
    
    // update likelihood lists after subtree removal
    // case when the sibling_subtree becomes the new root
    if (root_vector_index == parent_index.getVectorIndex()) //!sibling_subtree->neighbor)
    {
        // update root
        root_vector_index = sibling_index.getVectorIndex();
        
        // delete mid_branch_lh
        /*if (root->mid_branch_lh)
        {
            delete sibling_subtree->mid_branch_lh;
            sibling_subtree->mid_branch_lh = NULL;
        }*/
        sibling_subtree.setMidBranchLh(nullptr);
        
        // reset branch length (to 0) if sibling_subtree is root
        sibling_subtree.setUpperLength(0);

        // recompute the total lh regions at sibling
        //sibling_subtree->computeTotalLhAtNode(aln, model, threshold_prob, true);
        sibling_subtree.getPartialLh(TOP)->computeTotalLhAtRoot<num_states>(sibling_subtree.getTotalLh(), model);
        
        // traverse downwards (to childrens of the sibling) to update their lh regions
        if (sibling_subtree.isInternal()) //->next)
        {
            // update upper left/right regions
            /*Node* next_node_1 = sibling_subtree->next; // right
            Node* next_node_2 = next_node_1->next;*/ // left
            const Index neighbor_1_index = sibling_subtree.getNeighborIndex(RIGHT);
            PhyloNode& neighbor_1 = nodes[neighbor_1_index.getVectorIndex()];
            const Index neighbor_2_index = sibling_subtree.getNeighborIndex(LEFT);
            PhyloNode& neighbor_2 = nodes[neighbor_2_index.getVectorIndex()];
            
            // if (next_node_1->partial_lh) delete next_node_1->partial_lh;
            const std::unique_ptr<SeqRegions>& lower_regions_2 = neighbor_2.getPartialLh(TOP); // next_node_2->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
            // next_node_1->partial_lh = lower_reions->computeTotalLhAtRoot(num_states, model, next_node_2->length);
            lower_regions_2->computeTotalLhAtRoot<num_states>(sibling_subtree.getPartialLh(RIGHT), model, neighbor_2.getUpperLength());
            
            // if (next_node_2->partial_lh) delete next_node_2->partial_lh;
            /*lower_reions = next_node_1->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
            next_node_2->partial_lh = lower_reions->computeTotalLhAtRoot(num_states, model, next_node_1->length);*/
            const std::unique_ptr<SeqRegions>& lower_regions_1 = neighbor_1.getPartialLh(TOP);
            lower_regions_1->computeTotalLhAtRoot<num_states>(sibling_subtree.getPartialLh(LEFT), model, neighbor_1.getUpperLength());
            
            // add children to node_stack for further traversing and updating likelihood regions
            stack<Index> node_stack;
            node_stack.push(neighbor_1_index);
            node_stack.push(neighbor_2_index);
            updatePartialLh<num_states>(node_stack);
        }
    }
    // case when the sibling_subtree is non-root node
    else
    {
        // update branch length from the grandparent side
        // sibling_subtree->neighbor->length = sibling_subtree->length;
        
        // traverse to parent and sibling node to update their likelihood regions due to subtree remova
        stack<Index> node_stack;
        node_stack.push(sibling_index);
        node_stack.push(grandparent_index); // sibling_subtree->neighbor);
        updatePartialLh<num_states>(node_stack);
    }
    
    // replace the node and re-update the vector lists
    const std::unique_ptr<SeqRegions>& subtree_lower_regions = subtree.getPartialLh(TOP); // ->getPartialLhAtNode(aln, model, threshold_prob);
    // try to place the new sample as a descendant of a mid-branch point
    if (is_mid_branch && root_vector_index != best_node_index.getVectorIndex())
        placeSubTreeMidBranch<num_states>(best_node_index, subtree_index, subtree, subtree_lower_regions, branch_length, best_lh_diff);
    // otherwise, best lk so far is for appending directly to existing node
    else
        placeSubTreeAtNode<num_states>(best_node_index, subtree_index, subtree, subtree_lower_regions, branch_length, best_lh_diff);
}

template <const StateType num_states>
void cmaple::TreeBase::updateRegionsPlaceSubTree(PhyloNode& subtree, PhyloNode& sibling_node, PhyloNode& internal, std::unique_ptr<SeqRegions>&& best_child_regions, const std::unique_ptr<SeqRegions>& subtree_regions, const std::unique_ptr<SeqRegions>& upper_left_right_regions, const std::unique_ptr<SeqRegions>& lower_regions, RealNumType& best_blength)
{
    // update next_node_1->partial_lh
    // replacePartialLH(next_node_1->partial_lh, best_child_regions);
    internal.setPartialLh(LEFT, std::move(best_child_regions));

    // update new_internal_node->partial_lh
    // sibling_node.getPartialLhAtNode(aln, model, params->threshold_prob)->mergeTwoLowers<num_states>(new_internal_node->partial_lh, sibling_node->length, *subtree_regions, best_blength, aln, model, params->threshold_prob);
    sibling_node.getPartialLh(TOP)->mergeTwoLowers<num_states>(internal.getPartialLh(TOP), sibling_node.getUpperLength(), *subtree_regions, best_blength, aln, model, params->threshold_prob);
}

template <const StateType num_states>
void cmaple::TreeBase::updateRegionsPlaceSubTreeAbove(PhyloNode& subtree, PhyloNode& sibling_node, PhyloNode& internal, std::unique_ptr<SeqRegions>&& best_child_regions, const std::unique_ptr<SeqRegions>& subtree_regions, const std::unique_ptr<SeqRegions>& upper_left_right_regions, const std::unique_ptr<SeqRegions>& lower_regions, RealNumType& best_length)
{
    // sibling_node->getPartialLhAtNode(aln, model, params->threshold_prob)->mergeTwoLowers<num_states>(new_internal_node->partial_lh, sibling_node->length, *subtree_regions, best_length, aln, model, params->threshold_prob);
    sibling_node.getPartialLh(TOP)->mergeTwoLowers<num_states>(internal.getPartialLh(TOP), sibling_node.getUpperLength(), *subtree_regions, best_length, aln, model, params->threshold_prob);
    
    if (!internal.getPartialLh(TOP)) // new_internal_node->partial_lh)
    {
        if (cmaple::verbose_mode >= cmaple::VB_DEBUG)
            outWarning("Problem, non lower likelihood while placing subtree -> set best branch length to min length");
        best_length = min_blength;
        /*subtree->length = best_length;
        subtree->neighbor->length = best_length;*/
        subtree.setUpperLength(best_length);
        // lower_regions->mergeTwoLowers<num_states>(new_internal_node->partial_lh, sibling_node->length, *subtree_regions, best_length, aln, model, params->threshold_prob);
        lower_regions->mergeTwoLowers<num_states>(internal.getPartialLh(TOP), sibling_node.getUpperLength(), *subtree_regions, best_length, aln, model, params->threshold_prob);
    }
    // upper_left_right_regions->mergeUpperLower<num_states>(next_node_1->partial_lh, new_internal_node->length, *lower_regions, sibling_node->length, aln, model, params->threshold_prob);
    upper_left_right_regions->mergeUpperLower<num_states>(internal.getPartialLh(LEFT), internal.getUpperLength(), *lower_regions, sibling_node.getUpperLength(), aln, model, params->threshold_prob);
}

template<const StateType num_states, void (cmaple::TreeBase::*updateRegionsSubTree)(PhyloNode&, PhyloNode&, PhyloNode&, std::unique_ptr<SeqRegions>&&, const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, RealNumType&)>
void cmaple::TreeBase::connectSubTree2Branch(const std::unique_ptr<SeqRegions>& subtree_regions, const std::unique_ptr<SeqRegions>& lower_regions, const Index subtree_index, PhyloNode& subtree, const Index sibling_node_index, PhyloNode& sibling_node, const RealNumType top_distance, const RealNumType down_distance, RealNumType &best_blength, std::unique_ptr<SeqRegions>&& best_child_regions, const std::unique_ptr<SeqRegions>& upper_left_right_regions)
{
    const RealNumType threshold_prob = params->threshold_prob;
    ASSERT(sibling_node_index.getMiniIndex() == TOP);
    const NumSeqsType internal_vec = subtree.getNeighborIndex(TOP).getVectorIndex();
    PhyloNode& internal = nodes[internal_vec];
        
    // re-use internal nodes
    /*Node* next_node_1 = subtree->neighbor;
    Node* new_internal_node = next_node_1->getTopNode();
    Node* next_node_2 = next_node_1->getOtherNextNode();*/
    
    // NHANLT NOTES: UNNECESSARY
    // re-order next circle (not neccessary, just to make it consistent with Python code)
    /*new_internal_node->next = next_node_2;
    next_node_2->next = next_node_1;
    next_node_1->next = new_internal_node;*/
    
    // connect new_internal_node to the parent of the selected node
    /*new_internal_node->outdated = true;
    sibling_node->neighbor->neighbor = new_internal_node;
    new_internal_node->neighbor = sibling_node->neighbor;
    new_internal_node->length = top_distance;
    new_internal_node->neighbor->length = top_distance;*/
    
    internal.setOutdated(true);
    const Index parent_index = sibling_node.getNeighborIndex(TOP);
    const NumSeqsType parent_vec = parent_index.getVectorIndex();
    nodes[parent_vec].setNeighborIndex(parent_index.getMiniIndex(), Index(internal_vec, TOP));
    internal.setNeighborIndex(TOP, parent_index);
    internal.setUpperLength(top_distance);
    
    // connect the selected_node to new_internal_node (via next_node_2)
    /*sibling_node->neighbor = next_node_2;
    next_node_2->neighbor = sibling_node;
    sibling_node->length = down_distance;
    sibling_node->neighbor->length = down_distance;*/
    
    sibling_node.setNeighborIndex(TOP, Index(internal_vec, RIGHT));
    internal.setNeighborIndex(RIGHT, sibling_node_index);
    sibling_node.setUpperLength(down_distance);
    
    // subtree already connected to new_internal_node (via next_node_1)
    /*subtree->length = best_blength;
    subtree->neighbor->length = best_blength;*/
    
    subtree.setNeighborIndex(TOP, Index(internal_vec, LEFT));
    internal.setNeighborIndex(LEFT, subtree_index);
    subtree.setUpperLength(best_blength);
            
    // update all likelihood regions
    (this->*updateRegionsSubTree)(subtree, sibling_node, internal, std::move(best_child_regions), subtree_regions, upper_left_right_regions, lower_regions, best_blength);
    
    // upper_left_right_regions->mergeUpperLower<num_states>(next_node_2->partial_lh, new_internal_node->length, *subtree_regions, best_blength, aln, model, threshold_prob);
    upper_left_right_regions->mergeUpperLower<num_states>(internal.getPartialLh(RIGHT), internal.getUpperLength(), *subtree_regions, best_blength, aln, model, threshold_prob);
    
    RealNumType mid_branch_length = internal.getUpperLength() * 0.5; // new_internal_node->length * 0.5;
    // upper_left_right_regions->mergeUpperLower<num_states>(new_internal_node->mid_branch_lh, mid_branch_length, *new_internal_node->partial_lh, mid_branch_length, aln, model, threshold_prob);
    upper_left_right_regions->mergeUpperLower<num_states>(internal.getMidBranchLh(), mid_branch_length, *internal.getPartialLh(TOP), mid_branch_length, aln, model, threshold_prob);
    
    // new_internal_node->computeTotalLhAtNode(aln, model, threshold_prob, new_internal_node == root);
    internal.computeTotalLhAtNode<num_states>(internal.getTotalLh(), nodes[parent_vec], aln, model, threshold_prob, root_vector_index == internal_vec);
    
    if (!internal.getTotalLh()) //->total_lh || new_internal_node->total_lh->size() == 0)
        outError("Problem, None vector when re-placing sample, placing subtree at mid-branch point");
    
    //if distTop>=2*min_blengthForMidNode:
    // createFurtherMidNodes(newInternalNode,upper_left_right_regions)
    
    // NHANLT: LOGS FOR DEBUGGING
    /*if (params->debug)
    {
        cout << "2Branch " << (best_blength > 0 ? subtree.getTotalLh()->size():0)<< " " << (best_blength > 0 ? subtree.getMidBranchLh()->size():0)<< " " << subtree.getPartialLh(TOP)->size() << " " << internal.getTotalLh()->size() << " " << internal.getMidBranchLh()->size()<< " " << internal.getPartialLh(TOP)->size() << " " << internal.getPartialLh(LEFT)->size() << " " << internal.getPartialLh(RIGHT)->size() << std::endl;
         cout << std::setprecision(20) << internal.getUpperLength() << " " << internal.getCorrespondingLength(RIGHT, nodes) << " " << internal.getCorrespondingLength(LEFT, nodes) << std::endl;
    }*/

    // iteratively traverse the tree to update partials from the current node
    stack<Index> node_stack;
    node_stack.push(sibling_node_index);
    node_stack.push(subtree_index);
    node_stack.push(parent_index); // new_internal_node->neighbor);
    updatePartialLh<num_states>(node_stack);
}

template <const StateType num_states>
void cmaple::TreeBase::placeSubTreeMidBranch(const Index selected_node_index, const Index subtree_index, PhyloNode& subtree, const std::unique_ptr<SeqRegions>& subtree_regions, const RealNumType new_branch_length, const RealNumType new_lh)
{
    PhyloNode& selected_node = nodes[selected_node_index.getVectorIndex()];
    const RealNumType threshold_prob = params->threshold_prob;
    const std::unique_ptr<SeqRegions>& upper_left_right_regions = getPartialLhAtNode(selected_node.getNeighborIndex(TOP)); // selected_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
    // RealNumType best_split = 0.5;
    RealNumType best_blength_split = selected_node.getUpperLength() * 0.5;
    RealNumType best_split_lh = new_lh;
    // RealNumType new_split = 0.25;
    std::unique_ptr<SeqRegions> best_child_regions = nullptr; // std::make_unique<SeqRegions>(SeqRegions(selected_node.getMidBranchLh()));
    const std::unique_ptr<SeqRegions>& lower_regions = selected_node.getPartialLh(TOP);

    // try different positions on the existing branch
    bool found_new_split = tryShorterBranch<num_states, &cmaple::TreeBase::calculateSubTreePlacementCost<num_states>>(selected_node.getUpperLength(), best_child_regions, subtree_regions, upper_left_right_regions, lower_regions, best_split_lh, best_blength_split, new_branch_length, true);
    
    if (!found_new_split)
    {
        found_new_split = tryShorterBranch<num_states, &cmaple::TreeBase::calculateSubTreePlacementCost<num_states>>(selected_node.getUpperLength(), best_child_regions, subtree_regions, upper_left_right_regions, lower_regions, best_split_lh, best_blength_split, new_branch_length, false);
        
        if (found_new_split)
            best_blength_split = selected_node.getUpperLength() - best_blength_split;
    }
    
    // Delay cloning SeqRegions
    if (!best_child_regions)
        best_child_regions = std::make_unique<SeqRegions>(SeqRegions(selected_node.getMidBranchLh()));
    
    // now try different lengths for the new branch
    RealNumType best_blength = new_branch_length;
    estimateLengthNewBranch<&cmaple::TreeBase::calculateSubTreePlacementCost<num_states>>(best_split_lh, best_child_regions, subtree_regions, best_blength, max_blength, double_min_blength, (new_branch_length <= 0));
    
    // attach subtree to the branch above the selected node
    connectSubTree2Branch<num_states, &cmaple::TreeBase::updateRegionsPlaceSubTree<num_states>>(subtree_regions, nullptr, subtree_index, subtree, selected_node_index, selected_node, best_blength_split, selected_node.getUpperLength() - best_blength_split, best_blength, std::move(best_child_regions), upper_left_right_regions);
    
    // delete best_child_regions
    /*if (best_child_regions)
        delete best_child_regions;*/
}

template <const StateType num_states>
void cmaple::TreeBase::connectSubTree2Root(const Index subtree_index, PhyloNode& subtree, const std::unique_ptr<SeqRegions>& subtree_regions, const std::unique_ptr<SeqRegions>& lower_regions, const Index sibling_node_index, PhyloNode& sibling_node, const RealNumType best_root_blength, const RealNumType best_length2, std::unique_ptr<SeqRegions>&& best_parent_regions)
{
    ASSERT(sibling_node_index.getMiniIndex() == TOP);
    const NumSeqsType new_root_vec = subtree.getNeighborIndex(TOP).getVectorIndex();
    PhyloNode& new_root = nodes[new_root_vec];
    
    // re-use internal nodes
    /*Node* next_node_1 = subtree->neighbor;
    Node* new_root = next_node_1->getTopNode();
    Node* next_node_2 = next_node_1->getOtherNextNode();
    
    // NHANLT NOTES: UNNECESSARY
    // re-order next circle (not neccessary, just to make it consistent with Python code)
    new_root->next = next_node_2;
    next_node_2->next = next_node_1;
    next_node_1->next = new_root;*/
    
    // connect new_internal_node to the parent of the selected node
    /*new_root->outdated = true;
    new_root->neighbor = sibling_node->neighbor; // actually NULL since selected_node is root
    new_root->length = 0;*/
    
    new_root.setOutdated(true);
    new_root.setNeighborIndex(TOP, Index());
    new_root.setUpperLength(0);
    
    // connect the selected_node to new_internal_node (via next_node_2)
    /*sibling_node->neighbor = next_node_2;
    next_node_2->neighbor = sibling_node;
    sibling_node->length = best_root_blength;
    sibling_node->neighbor->length = best_root_blength;*/
    
    sibling_node.setNeighborIndex(TOP, Index(new_root_vec, RIGHT));
    new_root.setNeighborIndex(RIGHT, sibling_node_index);
    sibling_node.setUpperLength(best_root_blength);
    
    if (best_root_blength <= 0)
    {
        /*delete sibling_node->total_lh;
        sibling_node->total_lh = NULL;
        
        if (sibling_node->mid_branch_lh) delete sibling_node->mid_branch_lh;
        sibling_node->mid_branch_lh = NULL;*/
        
        sibling_node.setTotalLh(nullptr);
        sibling_node.setMidBranchLh(nullptr);
        
        //selected_node.furtherMidNodes=None
    }
    
    // subtree already connected to new_internal_node (via next_node_1)
    /*subtree->length = best_length2;
    subtree->neighbor->length = best_length2;*/
    
    subtree.setNeighborIndex(TOP, Index(new_root_vec, LEFT));
    new_root.setNeighborIndex(LEFT, subtree_index);
    subtree.setUpperLength(best_length2);
            
    // update all likelihood regions
    /*replacePartialLH(new_root->partial_lh, best_parent_regions);
    if (new_root->mid_branch_lh) delete new_root->mid_branch_lh;
    new_root->mid_branch_lh = NULL;*/
    
    new_root.setPartialLh(TOP, std::move(best_parent_regions));
    new_root.setMidBranchLh(nullptr);
    
    //new_root->computeTotalLhAtNode(aln, model, params->threshold_prob, true);
    new_root.getPartialLh(TOP)->computeTotalLhAtRoot<num_states>(new_root.getTotalLh(), model);

    /*if (next_node_1->partial_lh) delete next_node_1->partial_lh;
    next_node_1->partial_lh = lower_regions->computeTotalLhAtRoot(aln->num_states, model, best_root_blength);*/
    lower_regions->computeTotalLhAtRoot<num_states>(new_root.getPartialLh(LEFT), model, best_root_blength);
    
    /*if (next_node_2->partial_lh) delete next_node_2->partial_lh;
    next_node_2->partial_lh = subtree_regions->computeTotalLhAtRoot(aln->num_states, model, best_length2);*/
    subtree_regions->computeTotalLhAtRoot<num_states>(new_root.getPartialLh(RIGHT), model, best_length2);
    
    if (!new_root.getTotalLh() && cmaple::verbose_mode >= cmaple::VB_DEBUG) // ->total_lh || new_root->total_lh->size() == 0)
        outWarning("Problem, None vector when re-placing sample, position root");
    
    // NHANLT: LOGS FOR DEBUGGING
    /*if (params->debug)
    {
        cout << "2Root " << (best_length2 > 0 ? subtree.getTotalLh()->size():0)<< " " << (best_length2 > 0 ? subtree.getMidBranchLh()->size():0)<< " " << subtree.getPartialLh(TOP)->size()<< " " << new_root.getTotalLh()->size()<< " " << new_root.getPartialLh(TOP)->size() << " " << new_root.getPartialLh(LEFT)->size() << " " << new_root.getPartialLh(RIGHT)->size() << std::endl;
        cout << std::setprecision(20) << new_root.getCorrespondingLength(RIGHT, nodes) << " " << new_root.getCorrespondingLength(LEFT, nodes) << std::endl;
    }*/
    
    // update tree->root;
    root_vector_index = new_root_vec;
    
    // iteratively traverse the tree to update partials from the current node
    stack<Index> node_stack;
    node_stack.push(sibling_node_index);
    node_stack.push(subtree_index);
    updatePartialLh<num_states>(node_stack);
}

template <const StateType num_states>
void cmaple::TreeBase::handlePolytomyPlaceSubTree(const Index selected_node_index, PhyloNode& selected_node, const std::unique_ptr<SeqRegions>& subtree_regions, const RealNumType new_branch_length, RealNumType& best_down_lh_diff, Index& best_child_index, RealNumType& best_child_blength_split, std::unique_ptr<SeqRegions>& best_child_regions)
{
    // current node might be part of a polytomy (represented by 0 branch lengths) so we want to explore all the children of the current node to find out if the best placement is actually in any of the branches below the current node.
    const RealNumType threshold_prob = params->threshold_prob;
    stack<Index> new_node_stack;
    /*for (Index neighbor_index:nodes[selected_node_index.getVectorIndex()].getNeighborIndexes(TOP))
        new_node_stack.push(neighbor_index);*/
    if (selected_node.isInternal())
    {
        new_node_stack.push(selected_node.getNeighborIndex(RIGHT));
        new_node_stack.push(selected_node.getNeighborIndex(LEFT));
    }
    
    while (!new_node_stack.empty())
    {
        const Index node_index = new_node_stack.top();
        new_node_stack.pop();
        PhyloNode& node = nodes[node_index.getVectorIndex()];


        // add all nodes in polytomy
        if (node.getUpperLength() <= 0) //node->length <= 0)
        {
            /*FOR_NEIGHBOR(node, neighbor_node)
                new_node_stack.push(neighbor_node);*/
            /*for (Index neighbor_index:node.getNeighborIndexes(TOP))
                new_node_stack.push(neighbor_index);*/
            if (node.isInternal())
            {
                new_node_stack.push(node.getNeighborIndex(RIGHT));
                new_node_stack.push(node.getNeighborIndex(LEFT));
            }
        }
        else
        {
            // now try to place on the current branch below the best node, at an height above or equal to the mid-branch.
            std::unique_ptr<SeqRegions> mid_branch_regions = nullptr; // std::make_unique<SeqRegions>(SeqRegions(node.getMidBranchLh())); // new SeqRegions(node->mid_branch_lh);
            const std::unique_ptr<SeqRegions>& parent_upper_lr_regions = getPartialLhAtNode(node.getNeighborIndex(TOP)); // node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
            const std::unique_ptr<SeqRegions>& lower_regions = node.getPartialLh(TOP); // node->getPartialLhAtNode(aln, model, threshold_prob);
            RealNumType new_branch_length_split = 0.5 * node.getUpperLength(); // node->length;
            RealNumType tmp_lh_diff = calculateSubTreePlacementCost<num_states>(node.getMidBranchLh(), subtree_regions, new_branch_length);;
            
            while (true)
            {
                // if better placement found -> record it
                if (tmp_lh_diff > best_down_lh_diff)
                {
                    best_down_lh_diff = tmp_lh_diff;
                    best_child_index = node_index;
                    best_child_blength_split = new_branch_length_split;
                    new_branch_length_split *= 0.5;
                    
                    // replacePartialLH(best_child_regions, mid_branch_regions);
                    // Delay cloning SeqRegions
                    best_child_regions = mid_branch_regions ? std::move(mid_branch_regions) : std::make_unique<SeqRegions>(SeqRegions(node.getMidBranchLh()));
                    
                    if (new_branch_length_split <= half_min_blength_mid)
                        break;
                    
                    // compute mid_branch_regions
                    parent_upper_lr_regions->mergeUpperLower<num_states>(mid_branch_regions, new_branch_length_split, *lower_regions, node.getUpperLength() - new_branch_length_split, aln, model, threshold_prob);
                    
                    tmp_lh_diff = calculateSubTreePlacementCost<num_states>(mid_branch_regions, subtree_regions, new_branch_length);
                }
                else
                    break;
            }
            
            // delete mid_branch_regions
            // if (mid_branch_regions) delete mid_branch_regions;
        }
    }
}

template <const StateType num_states>
void cmaple::TreeBase::placeSubTreeAtNode(const Index selected_node_index, const Index subtree_index, PhyloNode& subtree, const std::unique_ptr<SeqRegions>& subtree_regions, const RealNumType new_branch_length, const RealNumType new_lh)
{
    // dummy variables
    const RealNumType threshold_prob = params->threshold_prob;
    RealNumType best_child_lh;
    RealNumType best_child_blength_split = -1;
    RealNumType best_parent_lh;
    RealNumType best_parent_blength_split = 0;
    std::unique_ptr<SeqRegions> best_parent_regions = nullptr;
    RealNumType best_root_blength = -1;
    std::unique_ptr<SeqRegions> best_child_regions = nullptr;
    RealNumType best_down_lh_diff = MIN_NEGATIVE;
    Index best_child_index;
    const NumSeqsType selected_node_vec = selected_node_index.getVectorIndex();
    PhyloNode& selected_node = nodes[selected_node_vec];
    
    // We first explore placement just below the best placement node for more fine-grained placement within its descendant branches (accounting for polytomies).
    handlePolytomyPlaceSubTree<num_states>(selected_node_index, selected_node, subtree_regions, new_branch_length, best_down_lh_diff, best_child_index, best_child_blength_split, best_child_regions);
    
    // place the new sample as a descendant of an existing node
    if (best_child_index.getMiniIndex() != UNDEFINED)
    {
        PhyloNode& best_child = nodes[best_child_index.getVectorIndex()];
        const std::unique_ptr<SeqRegions>& upper_left_right_regions = getPartialLhAtNode(best_child.getNeighborIndex(TOP)); // best_child->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
        const std::unique_ptr<SeqRegions>& lower_regions = best_child.getPartialLh(TOP); // ->getPartialLhAtNode(aln, model, threshold_prob);
        best_child_lh = best_down_lh_diff;
        best_child_blength_split = (best_child_blength_split == -1) ? (0.5 * best_child.getUpperLength()) : best_child_blength_split;
        
        // try with a shorter branch length
        tryShorterBranch<num_states, &cmaple::TreeBase::calculateSubTreePlacementCost<num_states>>(best_child.getUpperLength(), best_child_regions, subtree_regions, upper_left_right_regions, lower_regions, best_child_lh, best_child_blength_split, new_branch_length, true);
    }
    else
        best_child_lh = MIN_NEGATIVE;
    
    // if node is root, try to place as sibling of the current root.
    RealNumType old_root_lh = MIN_NEGATIVE;
    if (root_vector_index == selected_node_vec)
    {
        const std::unique_ptr<SeqRegions>& lower_regions = selected_node.getPartialLh(TOP); // ->getPartialLhAtNode(aln, model, threshold_prob);
        old_root_lh = lower_regions->computeAbsoluteLhAtRoot<num_states>(model);
        
        // merge 2 lower vector into one
        best_parent_lh = lower_regions->mergeTwoLowers<num_states>(best_parent_regions, default_blength, *subtree_regions, new_branch_length, aln, model, threshold_prob, true);
        
        best_parent_lh += best_parent_regions->computeAbsoluteLhAtRoot<num_states>(model);
        
        // Try shorter branch lengths at root
        best_root_blength = default_blength;
        tryShorterBranchAtRoot<num_states>(subtree_regions, lower_regions, best_parent_regions, best_root_blength, best_parent_lh, new_branch_length);
        
        // update best_parent_lh (taking into account old_root_lh)
        best_parent_lh -= old_root_lh;
    }
    // selected_node is not root
    // try to append just above node
    else
    {
        const std::unique_ptr<SeqRegions>& upper_left_right_regions = getPartialLhAtNode(selected_node.getNeighborIndex(TOP)); // selected_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
        const std::unique_ptr<SeqRegions>& lower_regions = selected_node.getPartialLh(TOP); // selected_node->getPartialLhAtNode(aln, model, threshold_prob);
        best_parent_regions = nullptr; // std::make_unique<SeqRegions>(SeqRegions(selected_node.getMidBranchLh()));
        best_parent_lh = calculateSubTreePlacementCost<num_states>(selected_node.getMidBranchLh(), subtree_regions, new_branch_length);
        best_parent_blength_split = 0.5 * selected_node.getUpperLength();
        
        // try with a shorter split
        tryShorterBranch<num_states, &cmaple::TreeBase::calculateSubTreePlacementCost<num_states>>(selected_node.getUpperLength(), best_parent_regions, subtree_regions, upper_left_right_regions, lower_regions, best_parent_lh, best_parent_blength_split, new_branch_length, false);
        
        // Delay cloning SeqRegions
        if (!best_parent_regions)
            best_parent_regions = std::make_unique<SeqRegions>(SeqRegions(selected_node.getMidBranchLh()));
    }
    
    // if the best placement is below the selected_node => add an internal node below the selected_node
    // now we have three likelihood costs,
    // best_child_lh is the likelihood score of appending below node;
    // best_parent_lh is the likelihood score of appending above node;
    // new_lh is the likelihood cost of appending exactly at node.
    if (best_child_lh >= best_parent_lh && best_child_lh >= new_lh)
    {
        ASSERT(best_child_index.getMiniIndex() != UNDEFINED);
        PhyloNode& best_child = nodes[best_child_index.getVectorIndex()];
        const std::unique_ptr<SeqRegions>& upper_left_right_regions = getPartialLhAtNode(best_child.getNeighborIndex(TOP)); // best_child->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
        
        // now try different lengths for the new branch
        RealNumType best_length = new_branch_length;
        estimateLengthNewBranch<&cmaple::TreeBase::calculateSubTreePlacementCost<num_states>>(best_child_lh, best_child_regions, subtree_regions, best_length, max_blength, double_min_blength, (new_branch_length <= 0));
        
        // attach subtree to the phylogenetic tree (below the selected_node ~ above the child node)
        connectSubTree2Branch<num_states, &cmaple::TreeBase::updateRegionsPlaceSubTree<num_states>>(subtree_regions, nullptr, subtree_index, subtree, best_child_index, best_child, best_child_blength_split, best_child.getUpperLength() - best_child_blength_split, best_length, std::move(best_child_regions), upper_left_right_regions);

    }
    // otherwise, add new parent to the selected_node
    else
    {
        const std::unique_ptr<SeqRegions>& lower_regions = selected_node.getPartialLh(TOP); // ->getPartialLhAtNode(aln, model, threshold_prob);
        
        // new parent is actually part of a polytomy since best placement is exactly at the node
        if (new_lh >= best_parent_lh)
        {
            best_root_blength = -1;
            best_parent_blength_split = -1;
            best_parent_lh = new_lh;
            /*if (best_parent_regions) delete best_parent_regions;
            best_parent_regions = NULL;*/
            best_parent_regions = nullptr;
            
            if (root_vector_index == selected_node_vec)
                lower_regions->mergeTwoLowers<num_states>(best_parent_regions, -1, *subtree_regions, new_branch_length, aln, model, threshold_prob);
            else
                best_parent_regions = std::make_unique<SeqRegions>(std::move(selected_node.getTotalLh()));
        }

        // add parent to the root
        if (root_vector_index == selected_node_vec)
        {
            // remove old_root_lh from best_parent_lh before estimating the new branch length
            best_parent_lh += old_root_lh;
            
            // estimate the new branch length
            RealNumType best_length2 = new_branch_length;
            estimateLengthNewBranchAtRoot<num_states>(subtree_regions, lower_regions, best_parent_regions, best_length2, best_parent_lh, best_root_blength, double_min_blength, new_branch_length <= 0);
            
            // update best_parent_lh (taking into account old_root_lh)
            best_parent_lh -= old_root_lh;
            
            // attach subtree to the phylogenetic tree (exactly at the seleted root node)
            connectSubTree2Root<num_states>(subtree_index, subtree, subtree_regions, lower_regions, selected_node_index, selected_node, best_root_blength, best_length2, std::move(best_parent_regions));
        }
        //add parent to non-root node (place subtree exactly at the selected non-root node)
        else
        {
            const std::unique_ptr<SeqRegions>& upper_left_right_regions = getPartialLhAtNode(selected_node.getNeighborIndex(TOP)); // selected_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
            
            // estimate the length for the new branch
            RealNumType best_length = new_branch_length;
            estimateLengthNewBranch<&cmaple::TreeBase::calculateSubTreePlacementCost<num_states>>(best_parent_lh, best_parent_regions, subtree_regions, best_length, new_branch_length * 10, double_min_blength, (new_branch_length <= 0));
            
            // attach subtree to the phylogenetic tree (exactly at the selected non-root node)
            RealNumType down_distance = best_parent_blength_split;
            RealNumType top_distance = selected_node.getUpperLength() - down_distance; // selected_node->length - down_distance;
            if (best_parent_blength_split <= 0)
            {
                down_distance = -1;
                top_distance = selected_node.getUpperLength(); // ->length;
                
                /*if (selected_node->total_lh) delete selected_node->total_lh;
                selected_node->total_lh = NULL;*/
                selected_node.setTotalLh(nullptr);
                
                /*if (selected_node->mid_branch_lh) delete selected_node->mid_branch_lh;
                selected_node->mid_branch_lh = NULL;*/
                selected_node.setMidBranchLh(nullptr);
            }
            
            connectSubTree2Branch<num_states, &cmaple::TreeBase::updateRegionsPlaceSubTreeAbove<num_states>>(subtree_regions, lower_regions, subtree_index, subtree, selected_node_index, selected_node, top_distance, down_distance, best_length, std::move(best_child_regions), upper_left_right_regions);
        }
    }
    
    // delete best_parent_regions and best_child_regions
    /*if (best_parent_regions)
        delete best_parent_regions;
    if (best_child_regions)
        delete best_child_regions;*/
}

template <const StateType num_states, RealNumType(cmaple::TreeBase::*calculatePlacementCost)(const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, const RealNumType)>
bool cmaple::TreeBase::tryShorterBranch(const RealNumType current_blength, std::unique_ptr<SeqRegions>& best_child_regions, const std::unique_ptr<SeqRegions>& sample, const std::unique_ptr<SeqRegions>& upper_left_right_regions, const std::unique_ptr<SeqRegions>& lower_regions, RealNumType &best_split_lh, RealNumType &best_branch_length_split, const RealNumType new_branch_length, const bool try_first_branch)
{
    std::unique_ptr<SeqRegions> new_parent_regions = nullptr;
    bool found_new_split = false;
    RealNumType new_branch_length_split = 0.5 * best_branch_length_split;
    
    while (new_branch_length_split > min_blength)
    {
        // try on the first or second branch
        if (try_first_branch)
            upper_left_right_regions->mergeUpperLower<num_states>(new_parent_regions, new_branch_length_split, *lower_regions,  current_blength - new_branch_length_split, aln, model, params->threshold_prob);
        else
            upper_left_right_regions->mergeUpperLower<num_states>(new_parent_regions, current_blength - new_branch_length_split, *lower_regions, new_branch_length_split, aln, model, params->threshold_prob);
        
        // calculate placement_cost
        RealNumType placement_cost = (this->*calculatePlacementCost)(new_parent_regions, sample, new_branch_length);
        
        if (placement_cost > best_split_lh)
        {
            best_split_lh = placement_cost;
            best_branch_length_split = new_branch_length_split;
            new_branch_length_split *= 0.5;
            found_new_split = true;
            
            best_child_regions = std::move(new_parent_regions);
        }
        else
            break;
    }
    
    return found_new_split;
}

template <RealNumType(cmaple::TreeBase::*calculatePlacementCost)(const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, const RealNumType)>
bool cmaple::TreeBase::tryShorterNewBranch(const std::unique_ptr<SeqRegions>& best_child_regions, const std::unique_ptr<SeqRegions>& sample, RealNumType &best_blength, RealNumType &new_branch_lh, const RealNumType short_blength_thresh)
{
    bool found_new_split = false;
    RealNumType new_blength = best_blength;
    RealNumType placement_cost;
    
    while (best_blength > short_blength_thresh)
    {
        new_blength *= 0.5;
        placement_cost = (this->*calculatePlacementCost)(best_child_regions, sample, new_blength);
        
        if (placement_cost > new_branch_lh)
        {
            new_branch_lh = placement_cost;
            best_blength = new_blength;
            found_new_split = true;
        }
        else
            break;
    }
    
    return found_new_split;
}

template <RealNumType(cmaple::TreeBase::*calculatePlacementCost)(const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, const RealNumType)>
void cmaple::TreeBase::tryLongerNewBranch(const std::unique_ptr<SeqRegions>& best_child_regions, const std::unique_ptr<SeqRegions>& sample, RealNumType &best_blength, RealNumType &new_branch_lh, const RealNumType long_blength_thresh)
{
    RealNumType new_blength = best_blength;
    RealNumType placement_cost;
    
    while (best_blength < long_blength_thresh)
    {
        new_blength += new_blength;
        placement_cost = (this->*calculatePlacementCost)(best_child_regions, sample, new_blength);
        if (placement_cost > new_branch_lh)
        {
            new_branch_lh = placement_cost;
            best_blength = new_blength;
        }
        else
            break;
    }
}

template <RealNumType(cmaple::TreeBase::*calculatePlacementCost)(const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, const RealNumType)>
void cmaple::TreeBase::estimateLengthNewBranch(const RealNumType best_split_lh, const std::unique_ptr<SeqRegions>& best_child_regions, const std::unique_ptr<SeqRegions>& sample, RealNumType &best_blength, const RealNumType long_blength_thresh, const RealNumType short_blength_thresh, const bool optional_check)
{
    RealNumType new_branch_lh = best_split_lh;
    
    // change zero branch length to min branch length
    if (optional_check)
    {
        best_blength = min_blength;
        new_branch_lh = (this->*calculatePlacementCost)(best_child_regions, sample, best_blength);
    }
    
    // try shorter lengths for the new branch
    bool found_new_blength = tryShorterNewBranch<calculatePlacementCost>(best_child_regions, sample, best_blength, new_branch_lh, min_blength);
    
    // try longer lengths for the new branch
    if (optional_check || !found_new_blength) // (best_blength > 0.7 * default_blength)
        tryLongerNewBranch<calculatePlacementCost>(best_child_regions, sample, best_blength, new_branch_lh, long_blength_thresh);
    
    // try zero-length for the new branch
    if (best_blength < short_blength_thresh)
    {
        RealNumType zero_branch_lh = (this->*calculatePlacementCost)(best_child_regions, sample, -1);
        if (zero_branch_lh > new_branch_lh)
            best_blength = -1;
    }
}

template <const StateType num_states>
void cmaple::TreeBase::connectNewSample2Branch(std::unique_ptr<SeqRegions>& sample, const NumSeqsType seq_name_index, const Index sibling_node_index, PhyloNode& sibling_node, const RealNumType top_distance, const RealNumType down_distance, const RealNumType best_blength, std::unique_ptr<SeqRegions>& best_child_regions, const std::unique_ptr<SeqRegions>& upper_left_right_regions)
{
    const RealNumType threshold_prob = params->threshold_prob;
    
    // create new internal node and append child to it
    /*Node* new_internal_node = new Node(true);
    Node* next_node_1 = new Node();
    Node* next_node_2 = new Node();
    Node* new_sample_node = new Node(seq_name);
    
    new_internal_node->next = next_node_2;
    next_node_2->next = next_node_1;
    next_node_1->next = new_internal_node;*/
    createAnInternalNode();
    createALeafNode(seq_name_index);
    const NumSeqsType leaf_vec_index = nodes.size() - 1;
    PhyloNode& leaf = nodes[leaf_vec_index];
    const NumSeqsType internal_vec_index = leaf_vec_index - 1;
    PhyloNode& internal = nodes[internal_vec_index];
    
    /*new_internal_node->neighbor = sibling_node->neighbor;
     sibling_node->neighbor->neighbor = new_internal_node;
     new_internal_node->length = top_distance;
     new_internal_node->neighbor->length = top_distance;*/
    
    Index parent_index = sibling_node.getNeighborIndex(TOP);
    PhyloNode& parent_node = nodes[parent_index.getVectorIndex()];
    internal.setNeighborIndex(TOP, parent_index);
    parent_node.setNeighborIndex(parent_index.getMiniIndex(), Index(internal_vec_index, TOP));
    internal.setUpperLength(top_distance);
    
    /*sibling_node->neighbor = next_node_2;
     next_node_2->neighbor = sibling_node;
     sibling_node->length = down_distance;
     sibling_node->neighbor->length = down_distance;*/
    
    internal.setNeighborIndex(RIGHT, sibling_node_index);
    sibling_node.setNeighborIndex(TOP, Index(internal_vec_index, RIGHT));
    sibling_node.setUpperLength(down_distance);
    
    /*new_sample_node->neighbor = next_node_1;
     next_node_1->neighbor = new_sample_node;
     new_sample_node->length = best_blength;
     new_sample_node->neighbor->length = best_blength;*/
    
    internal.setNeighborIndex(LEFT, Index(leaf_vec_index, TOP));
    leaf.setNeighborIndex(TOP, Index(internal_vec_index, LEFT));
    leaf.setUpperLength(best_blength);
    
    /*new_sample_node->partial_lh = sample;
    next_node_1->partial_lh = best_child_regions;
    best_child_regions = NULL;
    upper_left_right_regions->mergeUpperLower<num_states>(next_node_2->partial_lh, new_internal_node->length, *sample, best_blength, aln, model, threshold_prob);
    sibling_node->getPartialLhAtNode(aln, model, threshold_prob)->mergeTwoLowers<num_states>(new_internal_node->partial_lh, sibling_node->length, *sample, best_blength, aln, model, threshold_prob);
    RealNumType half_branch_length = new_internal_node->length * 0.5;
    upper_left_right_regions->mergeUpperLower<num_states>(new_internal_node->mid_branch_lh, half_branch_length, *new_internal_node->partial_lh, half_branch_length, aln, model, threshold_prob);*/
    
    leaf.setPartialLh(TOP, std::move(sample));
    SeqRegions& leaf_lower_regions = *leaf.getPartialLh(TOP);
    internal.setPartialLh(LEFT, std::move(best_child_regions));
    upper_left_right_regions->mergeUpperLower<num_states>(internal.getPartialLh(RIGHT), internal.getUpperLength(), leaf_lower_regions, best_blength, aln, model, threshold_prob);
    sibling_node.getPartialLh(TOP)->mergeTwoLowers<num_states>(internal.getPartialLh(TOP), sibling_node.getUpperLength(), leaf_lower_regions, best_blength, aln, model, threshold_prob);
    RealNumType half_branch_length = internal.getUpperLength() * 0.5;
    upper_left_right_regions->mergeUpperLower<num_states>(internal.getMidBranchLh(), half_branch_length, *(internal.getPartialLh(TOP)), half_branch_length, aln, model, threshold_prob);
    
    //new_internal_node->computeTotalLhAtNode(aln, model, threshold_prob, new_internal_node == root);
    internal.computeTotalLhAtNode<num_states>(internal.getTotalLh(), parent_node, aln, model, threshold_prob, root_vector_index == internal_vec_index);
    
    //if (!internal.getTotalLh() || internal.getTotalLh()->empty())
    if (!internal.getTotalLh())
        outError("Problem, None vector when placing sample, below node");
    
    if (best_blength > 0)
    {
        // new_sample_node->computeTotalLhAtNode(aln, model, threshold_prob, new_sample_node == root);
        leaf.computeTotalLhAtNode<num_states>(leaf.getTotalLh(), internal, aln, model, threshold_prob, root_vector_index == leaf_vec_index);
        
        /*RealNumType half_branch_length = new_sample_node->length * 0.5;
        next_node_1->getPartialLhAtNode(aln, model, threshold_prob)->mergeUpperLower<num_states>(new_sample_node->mid_branch_lh, half_branch_length, *sample, half_branch_length, aln, model, threshold_prob);*/
        RealNumType half_branch_length = leaf.getUpperLength() * 0.5;
        internal.getPartialLh(LEFT)->mergeUpperLower<num_states>(leaf.getMidBranchLh(), half_branch_length, leaf_lower_regions, half_branch_length, aln, model, threshold_prob);
    }
    
    // NHANLT: LOGS FOR DEBUGGING
    /*if (params->debug)
    {
        cout << "2Branch " << aln->data[seq_name_index].seq_name << " " << (best_blength > 0 ? leaf.getTotalLh()->size():0)<< " " << (best_blength > 0 ? leaf.getMidBranchLh()->size():0)<< " " << leaf.getPartialLh(TOP)->size() << " " << internal.getTotalLh()->size() << " " << internal.getMidBranchLh()->size()<< " " << internal.getPartialLh(TOP)->size() << " " << internal.getPartialLh(LEFT)->size() << " " << internal.getPartialLh(RIGHT)->size() << std::endl;
        cout << std::setprecision(20) << internal.getUpperLength() << " " << internal.getCorrespondingLength(RIGHT, nodes) << " " << internal.getCorrespondingLength(LEFT, nodes) << std::endl;
    }*/
    
    // update pseudo_count
    model->updatePesudoCount(aln, *internal.getPartialLh(LEFT), *leaf.getPartialLh(TOP));

    // iteratively traverse the tree to update partials from the current node
    stack<Index> node_stack;
    node_stack.push(sibling_node_index);
    node_stack.push(parent_index);
    updatePartialLh<num_states>(node_stack);
}

template <const StateType num_states>
void cmaple::TreeBase::placeNewSampleMidBranch(const Index& selected_node_index, std::unique_ptr<SeqRegions>& sample, const NumSeqsType seq_name_index, const RealNumType best_lh_diff)
{
    // dummy variables
    // const RealNumType threshold_prob = params->threshold_prob;
    std::unique_ptr<SeqRegions> best_child_regions = nullptr;
    
    // selected_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
    // const MiniIndex seleted_node_mini_index = selected_node_index.getMiniIndex();
    ASSERT(selected_node_index.getMiniIndex() == TOP);
    PhyloNode& selected_node = nodes[selected_node_index.getVectorIndex()];
    const std::unique_ptr<SeqRegions>& upper_left_right_regions = getPartialLhAtNode(selected_node.getNeighborIndex(TOP));
    RealNumType best_split_lh = best_lh_diff;
    const RealNumType selected_node_blength = selected_node.getUpperLength();
    RealNumType best_branch_length_split = 0.5 * selected_node_blength;
    best_child_regions = nullptr; // std::make_unique<SeqRegions>(SeqRegions(selected_node.getMidBranchLh()));
    //selected_node->getPartialLhAtNode(aln, model, threshold_prob);
    const std::unique_ptr<SeqRegions>& lower_regions = selected_node.getPartialLh(TOP);
    
    // try different positions on the existing branch
    bool found_new_split = tryShorterBranch<num_states, &cmaple::TreeBase::calculateSamplePlacementCost<num_states>>(selected_node_blength, best_child_regions, sample, upper_left_right_regions, lower_regions, best_split_lh, best_branch_length_split, default_blength, true);
    
    if (!found_new_split)
    {
        // try on the second branch
        found_new_split = tryShorterBranch<num_states, &cmaple::TreeBase::calculateSamplePlacementCost<num_states>>(selected_node_blength, best_child_regions, sample, upper_left_right_regions, lower_regions, best_split_lh, best_branch_length_split, default_blength, false);
        
        if (found_new_split)
            best_branch_length_split = selected_node_blength - best_branch_length_split;
    }
    
    // Delay cloning SeqRegions
    if (!best_child_regions)
        best_child_regions = std::make_unique<SeqRegions>(SeqRegions(selected_node.getMidBranchLh()));
    
    // now try different lengths for the new branch
    RealNumType best_blength = default_blength;
    estimateLengthNewBranch<&cmaple::TreeBase::calculateSamplePlacementCost<num_states>>(best_split_lh, best_child_regions, sample, best_blength, max_blength, min_blength, false);
    
    // create new internal node and append child to it
    connectNewSample2Branch<num_states>(sample, seq_name_index, selected_node_index, selected_node, best_branch_length_split, selected_node_blength - best_branch_length_split, best_blength, best_child_regions, upper_left_right_regions);
}

template <const StateType num_states>
void cmaple::TreeBase::tryShorterBranchAtRoot(const std::unique_ptr<SeqRegions>& sample, const std::unique_ptr<SeqRegions>& lower_regions, std::unique_ptr<SeqRegions>& best_parent_regions, RealNumType &best_root_blength, RealNumType &best_parent_lh, const RealNumType fixed_blength)
{
    std::unique_ptr<SeqRegions> merged_root_sample_regions = nullptr;
    RealNumType new_blength = 0.5 * best_root_blength;
    RealNumType new_root_lh;
    
    while (new_blength > min_blength)
    {
        // merge 2 lower vector into one
        new_root_lh = lower_regions->mergeTwoLowers<num_states>(merged_root_sample_regions, new_blength, *sample, fixed_blength, aln, model, params->threshold_prob, true);
        new_root_lh += merged_root_sample_regions->computeAbsoluteLhAtRoot<num_states>(model);
        
        if (new_root_lh > best_parent_lh)
        {
            best_parent_lh = new_root_lh;
            best_root_blength = new_blength;
            new_blength *= 0.5;
            
            // replacePartialLH(best_parent_regions, merged_root_sample_regions);
            best_parent_regions = std::move(merged_root_sample_regions);
        }
        else
            break;
    }
    
    // delete merged_root_sample_regions
    // if (merged_root_sample_regions) delete merged_root_sample_regions;
}

template <const StateType num_states>
bool cmaple::TreeBase::tryShorterNewBranchAtRoot(const std::unique_ptr<SeqRegions>& sample, const std::unique_ptr<SeqRegions>& lower_regions, std::unique_ptr<SeqRegions>&best_parent_regions, RealNumType &best_length, RealNumType &best_parent_lh, const RealNumType fixed_blength)
{
    std::unique_ptr<SeqRegions> new_root_lower_regions = nullptr;
    bool found_new_split = false;
    RealNumType new_blength = best_length;
    RealNumType new_root_lh = 0;
    
    while (best_length > min_blength)
    {
        new_blength *= 0.5;
        
        new_root_lh = lower_regions->mergeTwoLowers<num_states>(new_root_lower_regions, fixed_blength, *sample, new_blength, aln, model, params->threshold_prob, true);
        new_root_lh += new_root_lower_regions->computeAbsoluteLhAtRoot<num_states>(model);
        
        if (new_root_lh > best_parent_lh)
        {
            best_parent_lh = new_root_lh;
            best_length = new_blength;
            found_new_split = true;
            
            // replacePartialLH(best_parent_regions, new_root_lower_regions);
            best_parent_regions = std::move(new_root_lower_regions);
        }
        else
            break;
    }
    
    // delete new_root_lower_regions
    // if (new_root_lower_regions) delete new_root_lower_regions;
    
    return found_new_split;
}

template <const StateType num_states>
bool cmaple::TreeBase::tryLongerNewBranchAtRoot(const std::unique_ptr<SeqRegions>& sample, const std::unique_ptr<SeqRegions>& lower_regions, std::unique_ptr<SeqRegions>& best_parent_regions, RealNumType &best_length, RealNumType &best_parent_lh, const RealNumType fixed_blength)
{
    std::unique_ptr<SeqRegions> new_root_lower_regions = nullptr;
    bool found_new_split = false;
    RealNumType new_blength = best_length;
    RealNumType new_root_lh = 0;
    
    while (best_length < max_blength)
    {
        new_blength += new_blength;
        
        new_root_lh = lower_regions->mergeTwoLowers<num_states>(new_root_lower_regions, fixed_blength, *sample, new_blength, aln, model, params->threshold_prob, true);
        new_root_lh += new_root_lower_regions->computeAbsoluteLhAtRoot<num_states>(model);
        
        if (new_root_lh > best_parent_lh)
        {
            best_parent_lh = new_root_lh;
            best_length = new_blength;
            found_new_split = true;
            
            // replacePartialLH(best_parent_regions, new_root_lower_regions);
            best_parent_regions = std::move(new_root_lower_regions);
        }
        else
            break;
    }
    
    // delete new_root_lower_regions
    // if (new_root_lower_regions) delete new_root_lower_regions;
    
    return found_new_split;
}

template <const StateType num_states>
void cmaple::TreeBase::estimateLengthNewBranchAtRoot(const std::unique_ptr<SeqRegions>& sample, const std::unique_ptr<SeqRegions>& lower_regions, std::unique_ptr<SeqRegions>& best_parent_regions, RealNumType &best_length, RealNumType &best_parent_lh, const RealNumType fixed_blength, const RealNumType short_blength_thresh, const bool optional_check)
{
    if (optional_check)
    {
        std::unique_ptr<SeqRegions> new_root_lower_regions = nullptr;
        best_length = min_blength;
        best_parent_lh = lower_regions->mergeTwoLowers<num_states>(new_root_lower_regions, fixed_blength, *sample, best_length, aln, model, params->threshold_prob, true);
        
        best_parent_lh += new_root_lower_regions->computeAbsoluteLhAtRoot<num_states>(model);
        
        // replacePartialLH(best_parent_regions, new_root_lower_regions);
        best_parent_regions = std::move(new_root_lower_regions);
    }
    
    // try shorter lengths
    bool found_new_split = tryShorterNewBranchAtRoot<num_states>(sample, lower_regions, best_parent_regions, best_length, best_parent_lh, fixed_blength);
    
    // try longer lengths
    if (optional_check || !found_new_split)
        tryLongerNewBranchAtRoot<num_states>(sample, lower_regions, best_parent_regions, best_length, best_parent_lh, fixed_blength);
    
    
    // try with length zero
    if (best_length < short_blength_thresh)
    {
        std::unique_ptr<SeqRegions> new_root_lower_regions = nullptr;
        
        RealNumType new_root_lh = lower_regions->mergeTwoLowers<num_states>(new_root_lower_regions, fixed_blength, *sample, -1, aln, model, params->threshold_prob, true);
        new_root_lh += new_root_lower_regions->computeAbsoluteLhAtRoot<num_states>(model);
        
        if (new_root_lh > best_parent_lh)
        {
            best_length = -1;
            // replacePartialLH(best_parent_regions, new_root_lower_regions);
            best_parent_regions = std::move(new_root_lower_regions);
        }
        
        // delete new_root_lower_regions
        // if (new_root_lower_regions) delete new_root_lower_regions;
    }
}

template <const StateType num_states>
void cmaple::TreeBase::connectNewSample2Root(std::unique_ptr<SeqRegions>& sample, const NumSeqsType seq_name_index, const Index sibling_node_index, PhyloNode& sibling_node, const RealNumType best_root_blength, const RealNumType best_length2, std::unique_ptr<SeqRegions>& best_parent_regions)
{
    const RealNumType threshold_prob = params->threshold_prob;
    // const MiniIndex sibling_node_mini_index = sibling_node_index.getMiniIndex();
    
    /*Node* new_root = new Node(true);
    Node* next_node_1 = new Node();
    Node* next_node_2 = new Node();
    Node* new_sample_node = new Node(seq_name);
    
    new_root->next = next_node_2;
    next_node_2->next = next_node_1;
    next_node_1->next = new_root;*/
    
    createAnInternalNode();
    createALeafNode(seq_name_index);
    const NumSeqsType leaf_vec_index = nodes.size() - 1;
    PhyloNode& leaf = nodes[leaf_vec_index];
    const NumSeqsType new_root_vec_index = leaf_vec_index - 1;
    PhyloNode& new_root = nodes[new_root_vec_index];
    
    new_root.setNeighborIndex(TOP, Index());
    
    // attach the left child
    /*sibling_node->neighbor = next_node_2;
    next_node_2->neighbor = sibling_node;
    sibling_node->length = best_root_blength;
    sibling_node->neighbor->length = best_root_blength;*/
    
    new_root.setNeighborIndex(RIGHT, sibling_node_index);
    sibling_node.setNeighborIndex(TOP, Index(new_root_vec_index, RIGHT));
    sibling_node.setUpperLength(best_root_blength);
    
    if (best_root_blength <= 0)
    {
        /*if (sibling_node->total_lh) delete sibling_node->total_lh;
        sibling_node->total_lh = NULL;
        
        if (sibling_node->mid_branch_lh) delete sibling_node->mid_branch_lh;
        sibling_node->mid_branch_lh = NULL;*/
        sibling_node.setTotalLh(nullptr);
        sibling_node.setMidBranchLh(nullptr);
        //selected_node.furtherMidNodes=None
    }
    
    // attach the right child
    /*new_sample_node->neighbor = next_node_1;
    next_node_1->neighbor = new_sample_node;
    new_sample_node->length = best_length2;
    new_sample_node->neighbor->length = best_length2;*/
    
    new_root.setNeighborIndex(LEFT, Index(leaf_vec_index, TOP));
    leaf.setNeighborIndex(TOP, Index(new_root_vec_index, LEFT));
    leaf.setUpperLength(best_length2);
    
    /*new_root->partial_lh = best_parent_regions;
    best_parent_regions = NULL;*/
    new_root.setPartialLh(TOP, std::move(best_parent_regions));
    
    // new_root->total_lh = new_root->computeTotalLhAtNode(aln, model, threshold_prob, true);
    new_root.getPartialLh(TOP)->computeTotalLhAtRoot<num_states>(new_root.getTotalLh(), model);

    /*next_node_1->partial_lh = sibling_node->getPartialLhAtNode(aln, model, threshold_prob)->computeTotalLhAtRoot(aln->num_states, model, best_root_blength);
    next_node_2->partial_lh = sample->computeTotalLhAtRoot(aln->num_states, model, best_length2);*/
    sibling_node.getPartialLh(TOP)->computeTotalLhAtRoot<num_states>(new_root.getPartialLh(LEFT), model, best_root_blength);
    sample->computeTotalLhAtRoot<num_states>(new_root.getPartialLh(RIGHT), model, best_length2);
    
    // new_sample_node->partial_lh = sample;
    leaf.setPartialLh(TOP, std::move(sample));
    
    if (!new_root.getTotalLh()) //(!new_root.getTotalLh() || new_root->total_lh->size() == 0)
    {
        outError("Problem, None vector when placing sample, new root");
        //print(merged_root_sample_regions)
        //print(node.probVect)
        //print(sample)
        //print(best_length2)
        //print(best_root_blength)
    }
    
    if (best_length2 > 0)
    {
        //new_sample_node->computeTotalLhAtNode(aln, model, threshold_prob, new_sample_node == root);
        leaf.computeTotalLhAtNode<num_states>(leaf.getTotalLh(), new_root, aln, model, threshold_prob, root_vector_index == leaf_vec_index);
        
        RealNumType half_branch_length = leaf.getUpperLength() * 0.5; // new_sample_node->length * 0.5;
        // next_node_1->getPartialLhAtNode(aln, model, threshold_prob)->mergeUpperLower<num_states>(new_sample_node->mid_branch_lh, half_branch_length, *sample, half_branch_length, aln, model, threshold_prob);
        new_root.getPartialLh(LEFT)->mergeUpperLower<num_states>(leaf.getMidBranchLh(), half_branch_length, *leaf.getPartialLh(TOP), half_branch_length, aln, model, threshold_prob);
        
        //if best_length2>=2*min_blengthForMidNode:
          //  createFurtherMidNodes(new_root.children[1],new_root.probVectUpLeft)
    }
    
    // NHANLT: LOGS FOR DEBUGGING
    /*if (params->debug)
    {
        cout << "2Root " << aln->data[seq_name_index].seq_name << " " << (best_length2 > 0 ? leaf.getTotalLh()->size():0)<< " " << (best_length2 > 0 ? leaf.getMidBranchLh()->size():0)<< " " << leaf.getPartialLh(TOP)->size()<< " " << new_root.getTotalLh()->size()<< " " << new_root.getPartialLh(TOP)->size() << " " << new_root.getPartialLh(LEFT)->size() << " " << new_root.getPartialLh(RIGHT)->size() << std::endl;
        cout << std::setprecision(20) << new_root.getCorrespondingLength(RIGHT, nodes) << " " << new_root.getCorrespondingLength(LEFT, nodes) << std::endl;
    }*/
    
    // update tree->root;
    root_vector_index = new_root_vec_index;
    
    // iteratively traverse the tree to update partials from the current node
    stack<Index> node_stack;
    node_stack.push(sibling_node_index);
    updatePartialLh<num_states>(node_stack);
}

template <const StateType num_states>
void cmaple::TreeBase::placeNewSampleAtNode(const Index selected_node_index, std::unique_ptr<SeqRegions>& sample, const NumSeqsType seq_name_index, const RealNumType best_lh_diff, const RealNumType best_up_lh_diff, const RealNumType best_down_lh_diff, const Index best_child_index)
{
    // dummy variables
    RealNumType best_child_lh = MIN_NEGATIVE;
    RealNumType best_child_blength_split = 0;
    RealNumType best_parent_lh;
    RealNumType best_parent_blength_split = 0;
    RealNumType best_root_blength = -1;
    const RealNumType threshold_prob = params->threshold_prob;
    std::unique_ptr<SeqRegions> best_parent_regions = nullptr;
    std::unique_ptr<SeqRegions> best_child_regions = nullptr;
    
    ASSERT(selected_node_index.getMiniIndex() == TOP);
    const NumSeqsType selected_node_vec_index = selected_node_index.getVectorIndex();
    PhyloNode& selected_node = nodes[selected_node_vec_index];
    
    // place the new sample as a descendant of an existing node
    if (best_child_index.getMiniIndex() != UNDEFINED)
    {
        PhyloNode& best_child = nodes[best_child_index.getVectorIndex()];
        best_child_lh = best_down_lh_diff;
        ASSERT(best_child_index.getMiniIndex() == TOP);
        best_child_blength_split = 0.5 * best_child.getUpperLength(); // best_child->length;
        const std::unique_ptr<SeqRegions>& upper_left_right_regions = getPartialLhAtNode(best_child.getNeighborIndex(TOP)); // best_child->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
        const std::unique_ptr<SeqRegions>& lower_regions = best_child.getPartialLh(TOP); // ->getPartialLhAtNode(aln, model, threshold_prob);
        // best_child_regions = new SeqRegions(best_child->mid_branch_lh);
        // SeqRegions best_child_mid_clone = SeqRegions(best_child.getMidBranchLh());
        best_child_regions = nullptr; // std::make_unique<SeqRegions>(std::move(best_child_mid_clone));
        
        // try a shorter split
        // tryShorterBranch<&cmaple::TreeBase::calculateSamplePlacementCost<num_states>>(best_child->length, best_child_regions, sample, upper_left_right_regions, lower_regions, best_child_lh, best_child_blength_split, default_blength, true);
        tryShorterBranch<num_states, &cmaple::TreeBase::calculateSamplePlacementCost<num_states>>(best_child.getUpperLength(), best_child_regions, sample, upper_left_right_regions, lower_regions, best_child_lh, best_child_blength_split, default_blength, true);
        
        // Delay cloning SeqRegions
        if (!best_child_regions)
            best_child_regions = std::make_unique<SeqRegions>(SeqRegions(best_child.getMidBranchLh()));
    }
    
    // if node is root, try to place as sibling of the current root.
    RealNumType old_root_lh = MIN_NEGATIVE;
    if (root_vector_index == selected_node_vec_index)
    {
        /*old_root_lh = selected_node->getPartialLhAtNode(aln, model, threshold_prob)->computeAbsoluteLhAtRoot(num_states, model);
        SeqRegions* lower_regions = selected_node->getPartialLhAtNode(aln, model, threshold_prob);*/
        const std::unique_ptr<SeqRegions>& lower_regions = selected_node.getPartialLh(TOP);
        old_root_lh = lower_regions->computeAbsoluteLhAtRoot<num_states>(model);
        
        // merge 2 lower vector into one
        RealNumType new_root_lh = lower_regions->mergeTwoLowers<num_states>(best_parent_regions, default_blength, *sample, default_blength, aln, model, threshold_prob, true);
        
        new_root_lh += best_parent_regions->computeAbsoluteLhAtRoot<num_states>(model);
        best_parent_lh = new_root_lh;
        
        // try shorter branch lengths
        best_root_blength = default_blength;
        tryShorterBranchAtRoot<num_states>(sample, lower_regions, best_parent_regions, best_root_blength, best_parent_lh, default_blength);
        
        // update best_parent_lh (taking into account old_root_lh)
        best_parent_lh -= old_root_lh;
    }
    // selected_node is not root
    else
    {
        best_parent_lh = best_up_lh_diff;
        best_parent_blength_split = 0.5 * selected_node.getUpperLength(); // selected_node->length;
        /*SeqRegions* upper_left_right_regions = selected_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
        SeqRegions* lower_regions = selected_node->getPartialLhAtNode(aln, model, threshold_prob);
        best_parent_regions = new SeqRegions(selected_node->mid_branch_lh);*/
        const std::unique_ptr<SeqRegions>& upper_left_right_regions = getPartialLhAtNode(selected_node.getNeighborIndex(TOP));
        const std::unique_ptr<SeqRegions>& lower_regions = selected_node.getPartialLh(TOP);
        // SeqRegions seq_regions_clone = SeqRegions(selected_node.getMidBranchLh());
        best_parent_regions = nullptr; // std::make_unique<SeqRegions>(std::move(seq_regions_clone));
        
        // try a shorter split
        // tryShorterBranch<&cmaple::TreeBase::calculateSamplePlacementCost<num_states>>(selected_node->length, best_parent_regions, sample, upper_left_right_regions, lower_regions, best_parent_lh, best_parent_blength_split, default_blength, false);
        tryShorterBranch<num_states, &cmaple::TreeBase::calculateSamplePlacementCost<num_states>>(selected_node.getUpperLength(), best_parent_regions, sample, upper_left_right_regions, lower_regions, best_parent_lh, best_parent_blength_split, default_blength, false);
        
        // Delay cloning SeqRegions
        if (!best_parent_regions)
            best_parent_regions = std::make_unique<SeqRegions>(SeqRegions(selected_node.getMidBranchLh()));
    }
    
    // if the best placement is below the selected_node => add an internal node below the selected_node
    if (best_child_lh >= best_parent_lh && best_child_lh >= best_lh_diff)
    {
        ASSERT(best_child_index.getMiniIndex() == TOP);
        PhyloNode& best_child = nodes[best_child_index.getVectorIndex()];
        const std::unique_ptr<SeqRegions>& upper_left_right_regions = getPartialLhAtNode(best_child.getNeighborIndex(TOP)); // best_child->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
        
        // Estimate the length for the new branch
        RealNumType best_length = default_blength;
        estimateLengthNewBranch<&cmaple::TreeBase::calculateSamplePlacementCost<num_states>>(best_child_lh, best_child_regions, sample, best_length, max_blength, min_blength, false);
        
        // create new internal node and append child to it
        // connectNewSample2Branch(sample, seq_name_index, best_child, best_child_blength_split, best_child->length - best_child_blength_split, best_length, best_child_regions, upper_left_right_regions);
        connectNewSample2Branch<num_states>(sample, seq_name_index, best_child_index, best_child, best_child_blength_split, best_child.getUpperLength() - best_child_blength_split, best_length, best_child_regions, upper_left_right_regions);
    }
    // otherwise, add new parent to the selected_node
    else
    {
        // new parent is actually part of a polytomy since best placement is exactly at the node
        if (best_lh_diff >= best_parent_lh)
        {
            best_root_blength = -1;
            best_parent_blength_split = -1;
            best_parent_lh = best_lh_diff;
            /*if (best_parent_regions) delete best_parent_regions;
            best_parent_regions = NULL;*/
            best_parent_regions = nullptr;
            
            if (root_vector_index == selected_node_vec_index)
                // selected_node->getPartialLhAtNode(aln, model, threshold_prob)->mergeTwoLowers<num_states>(best_parent_regions, -1, *sample, default_blength, aln, model, threshold_prob);
                selected_node.getPartialLh(TOP)->mergeTwoLowers<num_states>(best_parent_regions, -1, *sample, default_blength, aln, model, threshold_prob);
            else
            {
                // best_parent_regions = new SeqRegions(selected_node->total_lh);
                SeqRegions seq_regions_clone = SeqRegions(selected_node.getTotalLh());
                best_parent_regions = std::make_unique<SeqRegions>(std::move(seq_regions_clone));
            }
        }

        // add parent to the root
        if (root_vector_index == selected_node_vec_index)
        {
            // now try different lengths for right branch
            best_parent_lh += old_root_lh;
            RealNumType best_length2 = default_blength;
            const std::unique_ptr<SeqRegions>& lower_regions = selected_node.getPartialLh(TOP); // ->getPartialLhAtNode(aln, model, threshold_prob);
            
            estimateLengthNewBranchAtRoot<num_states>(sample, lower_regions, best_parent_regions, best_length2, best_parent_lh, best_root_blength, min_blength, false);
            
            // update best_parent_lh (taking into account old_root_lh)
            best_parent_lh -= old_root_lh;
            
            // add new sample to a new root
            connectNewSample2Root<num_states>(sample, seq_name_index, selected_node_index, selected_node, best_root_blength, best_length2, best_parent_regions);
        }
        //add parent to non-root node
        else
        {
            const std::unique_ptr<SeqRegions>& upper_left_right_regions = getPartialLhAtNode(selected_node.getNeighborIndex(TOP)); // selected_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
            
            // now try different lengths for the new branch
            RealNumType new_branch_length_lh = best_parent_lh;
            RealNumType best_length = default_blength;
            estimateLengthNewBranch<&cmaple::TreeBase::calculateSamplePlacementCost<num_states>>(best_parent_lh, best_parent_regions, sample, best_length, max_blength, min_blength, false);
            
            // now create new internal node and append child to it
            RealNumType down_distance = best_parent_blength_split;
            RealNumType top_distance = selected_node.getUpperLength() - down_distance; // selected_node->length - down_distance;
            if (best_parent_blength_split < 0)
            {
                down_distance = -1;
                top_distance = selected_node.getUpperLength(); // selected_node->length;
                
                /*if (selected_node->total_lh) delete selected_node->total_lh;
                selected_node->total_lh = NULL;
                
                if (selected_node->mid_branch_lh) delete selected_node->mid_branch_lh;
                selected_node->mid_branch_lh = NULL;*/
                selected_node.setTotalLh(nullptr);
                selected_node.setMidBranchLh(nullptr);
                
                // node.furtherMidNodes=None
            }
            connectNewSample2Branch<num_states>(sample, seq_name_index, selected_node_index, selected_node, top_distance, down_distance, best_length, best_parent_regions, upper_left_right_regions);
        }
    }
    
    // delete best_parent_regions and best_child_regions
    /*if (best_parent_regions)
        delete best_parent_regions;
    if (best_child_regions)
        delete best_child_regions;*/
}

template <const StateType num_states>
void cmaple::TreeBase::refreshAllLhs(bool avoid_using_upper_lr_lhs)
{
    // 1. update all the lower lhs along the tree
    if (avoid_using_upper_lr_lhs)
        performDFS<&cmaple::TreeBase::updateLowerLhAvoidUsingUpperLRLh<num_states>>();
    else
        performDFS<&cmaple::TreeBase::updateLowerLh<num_states>>();
    
    // 2. update all the non-lower lhs along the tree
    refreshAllNonLowerLhs<num_states>();
}

template <const StateType num_states>
void cmaple::TreeBase::refreshUpperLR(const Index node_index, PhyloNode& node, const Index neighbor_index, std::unique_ptr<SeqRegions>& replaced_regions, const SeqRegions& parent_upper_lr_lh)
{
    // recalculate the upper left/right lh of the current node
    std::unique_ptr<SeqRegions> new_upper_lr_lh = nullptr;
    PhyloNode& neighbor = nodes[neighbor_index.getVectorIndex()];
    const std::unique_ptr<SeqRegions>& lower_lh = neighbor.getPartialLh(TOP); // next_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
    // parent_upper_lr_lh.mergeUpperLower<num_states>(new_upper_lr_lh, node->length, *lower_lh, next_node->length, aln, model, threshold_prob);
    parent_upper_lr_lh.mergeUpperLower<num_states>(new_upper_lr_lh, node.getUpperLength(), *lower_lh, neighbor.getUpperLength(), aln, model, params->threshold_prob);
    
    // if the upper left/right lh is null -> try to increase the branch length
    if (!new_upper_lr_lh)
    {
        if (neighbor.getUpperLength() <= 0) // next_node->length <= 0)
        {
            stack<Index> node_stack;
            // updateZeroBlength<num_states>(next_node->neighbor, node_stack, threshold_prob);
            updateZeroBlength<num_states>(neighbor_index, neighbor, node_stack);
            updatePartialLh<num_states>(node_stack);
        }
        else if (node.getUpperLength() <= 0) // node->length <= 0)
        {
            stack<Index> node_stack;
            // updateZeroBlength<num_states>(node, node_stack, threshold_prob);
            updateZeroBlength<num_states>(node_index, node, node_stack);
            updatePartialLh<num_states>(node_stack);
        }
        else
            outError("Strange, inconsistent upper left/right lh creation in refreshAllNonLowerLhs()");
    }
    // otherwise, everything is good -> update upper left/right lh of the current node
    else
        // replacePartialLH(replaced_regions, new_upper_lr_lh);
        replaced_regions = std::move(new_upper_lr_lh);
    
    // delete new_upper_lr_lh
    // if (new_upper_lr_lh) delete new_upper_lr_lh;
}

template <const StateType num_states>
void cmaple::TreeBase::refreshNonLowerLhsFromParent(Index& node_index, Index& last_node_index)
{
    PhyloNode& node = nodes[node_index.getVectorIndex()];
    const RealNumType threshold_prob = params->threshold_prob;
    const Index parent_index = node.getNeighborIndex(TOP);
    PhyloNode& parent_node = nodes[parent_index.getVectorIndex()];
    const std::unique_ptr<SeqRegions>& parent_upper_lr_lh = parent_node.getPartialLh(parent_index.getMiniIndex()); // node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
    
    // update the total lh, total lh at the mid-branch point of the current node
    if (node.getUpperLength() > 0) // node->length > 0)
    {
        // update the total lh
        // node->computeTotalLhAtNode(aln, model, threshold_prob, node == root);
        node.computeTotalLhAtNode<num_states>(node.getTotalLh(), parent_node, aln, model, threshold_prob, root_vector_index == node_index.getVectorIndex());
        
        if (!node.getTotalLh())
            outError("Strange, inconsistent total lh creation in refreshAllNonLowerLhs()");

        // update mid_branch_lh
        computeMidBranchRegions<num_states>(node, node.getMidBranchLh(), *parent_upper_lr_lh);
        
        // NHANLT: LOGS FOR DEBUGGING
        /*if (params->debug)
        {
            std::cout << "TotalLh " << node.getTotalLh()->size() << std::endl;
            std::cout << "MidbranchLh " << node.getMidBranchLh()->size() << std::endl;
        }*/
    }
    
    // if the current node is an internal node (~having children) -> update its upper left/right lh then traverse downward to update non-lower lhs of other nodes
    if (node.isInternal()) // !node->isLeave())
    {
        /*Node* next_node_1 = node->next;
        Node* next_node_2 = next_node_1->next;*/
        const Index neighbor_1_index = node.getNeighborIndex(RIGHT);
        const Index neighbor_2_index = node.getNeighborIndex(LEFT);
        
        // recalculate the FIRST upper left/right lh of the current node
        // refreshUpperLR(node, next_node_2, next_node_1->partial_lh, *parent_upper_lr_lh);
        refreshUpperLR<num_states>(node_index, node, neighbor_2_index, node.getPartialLh(RIGHT), *parent_upper_lr_lh);
        
        // recalculate the SECOND upper left/right lh of the current node
        // refreshUpperLR(node, next_node_1, next_node_2->partial_lh, *parent_upper_lr_lh);
        refreshUpperLR<num_states>(node_index, node, neighbor_1_index, node.getPartialLh(LEFT), *parent_upper_lr_lh);
        
        // NHANLT: LOGS FOR DEBUGGING
        /*if (params->debug)
        {
            std::cout << "next_node_1->partial_lh " << node.getPartialLh(RIGHT)->size() << std::endl;
            std::cout << "next_node_2->partial_lh " << node.getPartialLh(LEFT)->size() << std::endl;
        }*/
        
        // keep traversing downward to its firt child
        node_index = neighbor_1_index;
    }
    // if the current node is a leaf -> traverse upward to its parent
    else
    {
        last_node_index = node_index;
        node_index = node.getNeighborIndex(TOP); // node->neighbor;
    }
}

template <const StateType num_states>
void cmaple::TreeBase::refreshAllNonLowerLhs()
{
    // start from the root
    // update the total lh at root
    //node->computeTotalLhAtNode(aln, model, params->threshold_prob, true);
    PhyloNode& root = nodes[root_vector_index];
    root.getPartialLh(TOP)->computeTotalLhAtRoot<num_states>(root.getTotalLh(), model);
    
    // if the root has children -> update its upper left/right lh then traverse downward to update non-lower lhs of other nodes
    if (root.isInternal()) // !node->isLeave())
    {
        // update upper left/right lh of the root
        /*Node* next_node_1 = node->next;
        Node* next_node_2 = next_node_1->next;*/
        Index neighbor_1_index = root.getNeighborIndex(RIGHT);
        PhyloNode& neighbor_1 = nodes[neighbor_1_index.getVectorIndex()];
        PhyloNode& neighbor_2 = nodes[root.getNeighborIndex(LEFT).getVectorIndex()];
        
        /*delete next_node_1->partial_lh;
        next_node_1->partial_lh = next_node_2->neighbor->getPartialLhAtNode(aln, model, threshold_prob)->computeTotalLhAtRoot(num_states, model, next_node_2->length);*/
        neighbor_2.getPartialLh(TOP)->computeTotalLhAtRoot<num_states>(root.getPartialLh(RIGHT), model, neighbor_2.getUpperLength());
        
        /*delete next_node_2->partial_lh;
        next_node_2->partial_lh = next_node_1->neighbor->getPartialLhAtNode(aln, model, threshold_prob)->computeTotalLhAtRoot(num_states, model, next_node_1->length);*/
        neighbor_1.getPartialLh(TOP)->computeTotalLhAtRoot<num_states>(root.getPartialLh(LEFT), model, neighbor_1.getUpperLength());
        
        // NHANLT: LOGS FOR DEBUGGING
        /*if (params->debug)
        {
            std::cout << "root.setPartialLh " << root.getPartialLh(RIGHT)->size() << std::endl;
            std::cout << "root.setPartialLh " << root.getPartialLh(LEFT)->size() << std::endl;
        }*/
        
        // traverse the tree downward and update the non-lower genome lists for all other nodes of the tree.
        /*Node* last_node = NULL;
        node = next_node_1->neighbor;*/
        Index last_node_index;
        Index node_index = neighbor_1_index;
        while (node_index.getMiniIndex() != UNDEFINED)
        {
            // we reach a top node by a downward traversing
            if (node_index.getMiniIndex() == TOP) // node->is_top)
                refreshNonLowerLhsFromParent<num_states>(node_index, last_node_index);
            // we reach the current node by an upward traversing from its children
            else
            {
                /*Node* top_node = node->getTopNode();
                Node* next_node_1 = top_node->next;
                Node* next_node_2 = next_node_1->next;*/
                NumSeqsType node_vec = node_index.getVectorIndex();
                PhyloNode& node = nodes[node_vec];
                Index neighbor_1_index = node.getNeighborIndex(RIGHT);
                Index neighbor_2_index = node.getNeighborIndex(LEFT);
                
                // if we reach the current node by an upward traversing from its first children -> traversing downward to its second children
                if (last_node_index == neighbor_1_index) // next_node_1->neighbor)
                    node_index = neighbor_2_index; // node = next_node_2->neighbor;
                // otherwise, all children of the current node are updated -> update the lower lh of the current node
                else
                {
                    /*last_node_index = top_node;
                    node = top_node->neighbor;*/
                    last_node_index = Index(node_vec, TOP);
                    node_index = node.getNeighborIndex(TOP);
                }
            }
        }
    }
}

template <const StateType num_states>
RealNumType cmaple::TreeBase::improveEntireTree(bool short_range_search)
{
    // start from the root
    stack<Index> node_stack;
    node_stack.push(Index(root_vector_index, TOP));
    
    // dummy variables
    RealNumType total_improvement = 0;
    PositionType num_nodes = 0;
    
    // traverse downward the tree
    while (!node_stack.empty())
    {
        // pick the top node from the stack
        Index index = node_stack.top();
        node_stack.pop();
        PhyloNode& node = nodes[index.getVectorIndex()];
        // MiniIndex mini_index = index.getMiniIndex();

        // add all children of the current nodes to the stack for further traversing later
        /*for (Index neighbor_index:node.getNeighborIndexes(mini_index))
            node_stack.push(neighbor_index);*/
        ASSERT(index.getMiniIndex() == TOP);
        if (node.isInternal())
        {
            node_stack.push(node.getNeighborIndex(RIGHT));
            node_stack.push(node.getNeighborIndex(LEFT));
        }
        
        // only process outdated node to avoid traversing the same part of the tree multiple times
        if (node.isOutdated() && !node.isSPRApplied())
        {
            node.setOutdated(false);
            
            //if checkEachSPR:
              //  root=node
               // while root.up!=None:
               //     root=root.up
               // #print("Pre-SPR tree: "+createBinaryNewick(root))
               // oldTreeLK=calculateTreeLikelihood(root,mutMatrix,checkCorrectness=True)
               // #print("Pre-SPR tree likelihood: "+str(oldTreeLK))
                //reCalculateAllGenomeLists(root,mutMatrix, checkExistingAreCorrect=True)
            
            // do SPR moves to improve the tree
            RealNumType improvement = improveSubTree<num_states>(index, node, short_range_search);
            
            //if checkEachSPR:
      //          #print(" apparent improvement "+str(improvement))
        //        root=node
         //       while root.up!=None:
           //         root=root.up
     //           #print("Post-SPR tree: "+createBinaryNewick(root))
       //         newTreeLK=calculateTreeLikelihood(root,mutMatrix)
         //       reCalculateAllGenomeLists(root,mutMatrix, checkExistingAreCorrect=True)
           //     #print("Post-SPR tree likelihood: "+str(newTreeLK))
             //   if newTreeLK-oldTreeLK < improvement-1.0:
      //              print("In startTopologyUpdates, LK score of improvement "+str(newTreeLK)+" - "+str(oldTreeLK)+" = "+str(newTreeLK-oldTreeLK)+" is less than //what is supposd to be "+str(improvement))
            //        exit()
             //
            
            // update total_improvement
            total_improvement += improvement;
            
            // NHANLT: LOGS FOR DEBUGGING
            /*if (params->debug && improvement > 0)
                std::cout << num_nodes << ": " << std::setprecision(20) << total_improvement << std::endl;*/
            
            // Show log every 1000 nodes
            num_nodes += 1;
            if (cmaple::verbose_mode >= cmaple::VB_DEBUG && num_nodes % 1000 == 0)
                cout << "Processed topology for " << convertIntToString(num_nodes) << " nodes." << endl;
        }
    }
    
    return total_improvement;
}

template <const StateType num_states>
PositionType cmaple::TreeBase::optimizeBranchLengths()
{
    // start from the root's children
    stack<Index> node_stack;
    PhyloNode& root = nodes[root_vector_index];
    if (!root.isInternal()) // !root || !root->next)
        return 0;
    /*Node* neighbor_node = NULL;
    FOR_NEIGHBOR(root, neighbor_node)
        node_stack.push(neighbor_node);*/
    node_stack.push(root.getNeighborIndex(RIGHT));
    node_stack.push(root.getNeighborIndex(LEFT));
    
    // dummy variables
    PositionType num_improvement = 0;
    const RealNumType threshold_prob = params->threshold_prob;
    
    // traverse downward the tree
    while (!node_stack.empty())
    {
        // pick the top node from the stack
        Index node_index = node_stack.top();
        node_stack.pop();
        PhyloNode& node = nodes[node_index.getVectorIndex()];
        
        const Index parent_index = node.getNeighborIndex(TOP);
        const std::unique_ptr<SeqRegions>& upper_lr_regions = getPartialLhAtNode(parent_index); // node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
        const std::unique_ptr<SeqRegions>& lower_regions = node.getPartialLh(TOP); // node->getPartialLhAtNode(aln, model, threshold_prob);
        
        // add all children of the current nodes to the stack for further traversing later
        /*neighbor_node = NULL;
        FOR_NEIGHBOR(node, neighbor_node)
            node_stack.push(neighbor_node);*/
        /*for (Index neighbor_index:node.getNeighborIndexes(TOP))
            node_stack.push(neighbor_index);*/
        if (node.isInternal())
        {
            node_stack.push(node.getNeighborIndex(RIGHT));
            node_stack.push(node.getNeighborIndex(LEFT));
        }
        
        // only process outdated node to avoid traversing the same part of the tree multiple times
        if (node.isOutdated())
        {
            // estimate the branch length
            RealNumType best_length = estimateBranchLength<num_states>(upper_lr_regions, lower_regions);
            
            if (best_length > 0 || node.getUpperLength() > 0)
            {
                RealNumType diff_thresh = 0.01 * best_length;
                if (best_length <= 0 || node.getUpperLength() <= 0 || (node.getUpperLength() > (best_length + diff_thresh)) || (node.getUpperLength() < (best_length - diff_thresh)))
                {
                    node.setUpperLength(best_length);
                    // node->neighbor->length = node->length;
                    ++num_improvement;
                    
                    // update partial likelihood regions
                    stack<Index> new_node_stack;
                    new_node_stack.push(node_index);
                    new_node_stack.push(parent_index);
                    updatePartialLh<num_states>(new_node_stack);
                }
            }
        }
    }
    
    return num_improvement;
}

template <const StateType num_states>
void cmaple::TreeBase::estimateBlength_R_O(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType total_blength, const PositionType end_pos, RealNumType &coefficient, std::vector<RealNumType> &coefficient_vec)
{
    const StateType seq1_state = aln->ref_seq[end_pos];
    RealNumType* mutation_mat_row = model->mutation_mat + model->row_index[seq1_state];
    RealNumType coeff0 = seq2_region.getLH(seq1_state);
    RealNumType coeff1 = 0;

    if (seq1_region.plength_observation2root >= 0)
    {
      coeff0 *= model->root_freqs[seq1_state];

      RealNumType* transposed_mut_mat_row = model->transposed_mut_mat + model->row_index[seq1_state];

      updateCoeffs<num_states>(model->root_freqs, transposed_mut_mat_row, &(*seq2_region.likelihood)[0], mutation_mat_row, seq1_region.plength_observation2node, coeff0, coeff1);

      coeff1 *= model->root_freqs[seq1_state];
    }
    else
    {
      // NHANLT NOTES:
      // x = seq1_state
      // l = log(1 + q_xx * t + sum(q_xy * t)
      // l' = [q_xx + sum(q_xy)]/[1 + q_xx * t + sum(q_xy * t)]
      // coeff1 = numerator = q_xx + sum(q_xy)
        coeff1 += dotProduct<num_states>(&(*seq2_region.likelihood)[0], mutation_mat_row);
    }

    // NHANLT NOTES:
    // l = log(1 + q_xx * t + sum(q_xy * t)
    // l' = [q_xx + sum(q_xy)]/[1 + q_xx * t + sum(q_xy * t)]
    // coeff0 = denominator = 1 + q_xx * t + sum(q_xy * t)
    if (total_blength > 0)
      coeff0 += coeff1 * total_blength;

    // NHANLT NOTES:
    // l' = [q_xx + sum(q_xy)]/[1 + q_xx * t + sum(q_xy * t)] = coeff1 / coeff0
    if (coeff1 < 0)
      coefficient += coeff1 / coeff0;
    else
      coefficient_vec.push_back(coeff0 / coeff1);
}

void cmaple::TreeBase::estimateBlength_R_ACGT(const SeqRegion& seq1_region, const StateType seq2_state, const RealNumType total_blength, const PositionType end_pos, std::vector<RealNumType> &coefficient_vec)
{
    if (seq1_region.plength_observation2root >= 0)
    {
        StateType seq1_state = aln->ref_seq[end_pos];
        
        RealNumType coeff1 = model->root_freqs[seq1_state] * model->mutation_mat[model->row_index[seq1_state] + seq2_state];
        RealNumType coeff0 = model->root_freqs[seq2_state] * model->mutation_mat[model->row_index[seq2_state] + seq1_state] * seq1_region.plength_observation2node;
        
        if (total_blength > 0)
            coeff0 += coeff1 * total_blength;
        
        coefficient_vec.push_back(coeff0 / coeff1);
    }
    // NHANLT: add else here, otherwise, coefficient_vec.push_back(total_blength > 0 ? total_blength : 0) is called even when (seq1_region->plength_observation2root >= 0)
    else
        // NHANLT NOTES:
        // l = log(q_xy * t)
        // l' = q_xy / (q_xy * t) = 1 / t
        coefficient_vec.push_back(total_blength > 0 ? total_blength : 0);
}

template <const StateType num_states>
void cmaple::TreeBase::estimateBlength_O_X(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType total_blength, const PositionType end_pos, RealNumType &coefficient, std::vector<RealNumType> &coefficient_vec)
{
    RealNumType coeff0 = 0;
    RealNumType coeff1 = 0;
    
    // 3.1. e1.type = O and e2.type = O
    if (seq2_region.type == TYPE_O)
    {
        RealNumType* mutation_mat_row = model->mutation_mat;
        
        // NHANLT NOTES:
        // l = log(sum_x(1 + q_xx * t + sum_y(q_xy * t)))
        // l' = [sum_x(q_xx + sum_y(q_xy))]/[sum_x(1 + q_xx * t + sum_y(q_xy * t))]
        // coeff1 = numerator = sum_x(q_xx + sum_y(q_xy))
        // coeff0 = denominator = sum_x(1 + q_xx * t + sum_y(q_xy * t))
        for (StateType i = 0; i < num_states; ++i, mutation_mat_row += num_states)
        {
            RealNumType seq1_lh_i = seq1_region.getLH(i);
            coeff0 += seq1_lh_i * seq2_region.getLH(i);
            
            for (StateType j = 0; j < num_states; ++j)
                coeff1 += seq1_lh_i * seq2_region.getLH(j) * mutation_mat_row[j];
        }
    }
    // 3.2. e1.type = O and e2.type = R or A/C/G/T
    else
    {
        StateType seq2_state = seq2_region.type;
        if (seq2_state == TYPE_R)
            seq2_state = aln->ref_seq[end_pos];
        
        coeff0 = seq1_region.getLH(seq2_state);

        // NHANLT NOTES:
        // y = seq2_state
        // l = log(1 + q_yy * t + sum_x(q_xy * t)
        // l' = [q_yy + sum_x(q_xy))]/[1 + q_xx * t + sum_y(q_xy * t)]
        // coeff1 = numerator = q_yy + sum_x(q_xy))
        // coeff0 = denominator = 1 + q_xx * t + sum_y(q_xy * t)
        RealNumType* transposed_mut_mat_row = model->transposed_mut_mat + model->row_index[seq2_state];
        coeff1 += dotProduct<num_states>(&(*seq1_region.likelihood)[0], transposed_mut_mat_row);
    }
    
    if (total_blength > 0)
        coeff0 += coeff1 * total_blength;
    
    // NHANLT NOTES:
    // l' = coeff1 / coeff0
    if (coeff1 < 0)
        coefficient += coeff1 / coeff0;
    else
        coefficient_vec.push_back(coeff0 / coeff1);
}

template <const StateType num_states>
void cmaple::TreeBase::estimateBlength_ACGT_O(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType total_blength, RealNumType &coefficient, std::vector<RealNumType> &coefficient_vec)
{
    StateType seq1_state = seq1_region.type;
    RealNumType coeff0 = seq2_region.getLH(seq1_state);
    RealNumType coeff1 = 0;
    
    RealNumType* mutation_mat_row = model->mutation_mat + model->row_index[seq1_state];
    
    if (seq1_region.plength_observation2root >= 0)
    {
        coeff0 *= model->root_freqs[seq1_state];

        RealNumType* transposed_mut_mat_row = model->transposed_mut_mat + model->row_index[seq1_state];
                            
        updateCoeffs<num_states>(model->root_freqs, transposed_mut_mat_row, &(*seq2_region.likelihood)[0], mutation_mat_row, seq1_region.plength_observation2node, coeff0, coeff1);
        
        coeff1 *= model->root_freqs[seq1_state];
    }
    else
    {
        // NHANLT NOTES:
        // x = seq1_state
        // l = log(1 + q_xx * t + sum(q_xy * t)
        // l' = [q_xx + sum(q_xy)]/[1 + q_xx * t + sum(q_xy * t)]
        // coeff1 = numerator = q_xx + sum(q_xy)
        coeff1 += dotProduct<num_states>(&(*seq2_region.likelihood)[0], mutation_mat_row);
    }
    
    // NHANLT NOTES:
    // l = log(1 + q_xx * t + sum(q_xy * t)
    // l' = [q_xx + sum(q_xy)]/[1 + q_xx * t + sum(q_xy * t)]
    // coeff0 = denominator = 1 + q_xx * t + sum(q_xy * t)
    if (total_blength > 0)
        coeff0 += coeff1 * total_blength;
    
    // NHANLT NOTES:
    // l' = [q_xx + sum(q_xy)]/[1 + q_xx * t + sum(q_xy * t)] = coeff1 / coeff0;
    if (coeff1 < 0)
        coefficient += coeff1 / coeff0;
    else
        coefficient_vec.push_back(coeff0 / coeff1);
}

void cmaple::TreeBase::estimateBlength_ACGT_RACGT(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType total_blength, const PositionType end_pos, std::vector<RealNumType> &coefficient_vec)
{
    RealNumType coeff0 = 0;
    StateType seq1_state = seq1_region.type;
    StateType seq2_state = seq2_region.type;
    if (seq2_state == TYPE_R)
        seq2_state = aln->ref_seq[end_pos];
    
    if (seq1_region.plength_observation2root >= 0)
    {
        coeff0 = model->root_freqs[seq2_state] * model->mutation_mat[model->row_index[seq2_state] + seq1_state] * seq1_region.plength_observation2node;
        RealNumType coeff1 = model->root_freqs[seq1_state] * model->mutation_mat[model->row_index[seq1_state] + seq2_state];
        
        if (total_blength > 0)
            coeff0 += coeff1 * total_blength;
        
        coeff0 /= coeff1;
    }
    // NHANLT NOTES:
    // l = log(q_xy * t)
    // l' = q_xy / (q_xy * t) = 1 / t
    else if (total_blength > 0)
        coeff0 = total_blength;
    
    coefficient_vec.push_back(coeff0);
}

RealNumType cmaple::TreeBase::estimateBlengthFromCoeffs(RealNumType &coefficient, const std::vector<RealNumType> coefficient_vec)
{
    coefficient = -coefficient;
    PositionType num_coefficients = coefficient_vec.size();
    if (num_coefficients == 0)
        return -1;

    // Get min and max coefficients
    RealNumType min_coefficient = coefficient_vec[0];
    RealNumType max_coefficient = coefficient_vec[0];
    for (PositionType i = 1; i < num_coefficients; ++i)
    {
        RealNumType coefficient_i = coefficient_vec[i];
        if (coefficient_i < min_coefficient)
            min_coefficient = coefficient_i;
        if (coefficient_i > max_coefficient)
            max_coefficient = coefficient_i;
    }
    
    RealNumType num_coefficients_over_coefficient = num_coefficients / coefficient;
    RealNumType tDown = num_coefficients_over_coefficient - min_coefficient;
    if (tDown <= 0)
        return 0;
    RealNumType derivative_tDown = calculateDerivative(coefficient_vec, tDown);
    
    RealNumType tUp = num_coefficients_over_coefficient - max_coefficient;
    if (tUp < 0)
    {
        if (min_coefficient > 0)
            tUp = 0;
        else
            tUp = min_blength_sensitivity;
    }
    RealNumType derivative_tUp = calculateDerivative(coefficient_vec, tUp);
    
    if ((derivative_tDown > coefficient + min_blength_sensitivity) || (derivative_tUp < coefficient - min_blength_sensitivity))
        if ((derivative_tUp < coefficient - min_blength_sensitivity) && (tUp == 0))
            return 0;
    
    while (tDown - tUp > min_blength_sensitivity)
    {
        RealNumType tMiddle = (tUp + tDown) * 0.5;
        RealNumType derivative_tMiddle = calculateDerivative(coefficient_vec, tMiddle);
        
        if (derivative_tMiddle > coefficient)
            tUp = tMiddle;
        else
            tDown = tMiddle;
    }
    
    return tUp;
}

using DoubleState = uint16_t;
static constexpr DoubleState RR = (DoubleState(TYPE_R) << 8) | TYPE_R;
static constexpr DoubleState RO = (DoubleState(TYPE_R) << 8) | TYPE_O;
static constexpr DoubleState OO = (DoubleState(TYPE_O) << 8) | TYPE_O;

template <const StateType num_states>
RealNumType cmaple::TreeBase::estimateBranchLength(const std::unique_ptr<SeqRegions>& parent_regions, const std::unique_ptr<SeqRegions>& child_regions)
{
    // init dummy variables
    RealNumType coefficient = 0;
    vector<RealNumType> coefficient_vec;
    PositionType pos = 0;
    const SeqRegions& seq1_regions = *parent_regions;
    const SeqRegions& seq2_regions = *child_regions;
    size_t iseq1 = 0;
    size_t iseq2 = 0;
    const PositionType seq_length = aln->ref_seq.size();
    const RealNumType* const &cumulative_rate = model->cumulative_rate;
    
    // avoid reallocations
    coefficient_vec.reserve(parent_regions->countSharedSegments(seq2_regions, seq_length)); // avoid realloc of vector data

    while (pos < seq_length)
    {
        PositionType end_pos;
        RealNumType total_blength;
        
        // get the next shared segment in the two sequences
        SeqRegions::getNextSharedSegment(pos, seq1_regions, seq2_regions, iseq1, iseq2, end_pos);
        const auto* seq1_region = &seq1_regions[iseq1];
        const auto* seq2_region = &seq2_regions[iseq2];
        
        // 1. e1.type = N || e2.type = N
        if ((seq2_region->type == TYPE_N) + (seq1_region->type == TYPE_N))
        {
            pos = end_pos + 1;
            continue;
        }
        
        // e1.type != N && e2.type != N
        const DoubleState s1s2 = (DoubleState(seq1_region->type) << 8) | seq2_region->type;

        // total_blength will be here the total length from the root or from the upper node, down to the down node.
        if (seq1_region->plength_observation2root >= 0)
            total_blength = seq1_region->plength_observation2root;
        else if (seq1_region->plength_observation2node >= 0)
            total_blength = seq1_region->plength_observation2node;
        else
            total_blength = 0;
            
        if (seq2_region->plength_observation2node >= 0)
            total_blength = total_blength + seq2_region->plength_observation2node;
            //total_blength = (total_blength > 0 ? total_blength : 0) + seq2_region->plength_observation2node;
        
        // 2. e1.type = R
        // 2.1.e1.type = R and e2.type = R
        if (s1s2 == RR)
        {
            // NHANLT NOTES:
            // coefficient = derivative of log likelihood function wrt t
            // l = log(1 + q_xx * t) ~ q_xx * t
            // => l' = q_xx
            //if (seq2_region->type == TYPE_R)
              coefficient += cumulative_rate[end_pos + 1] - cumulative_rate[pos];
        }
        // 2.2. e1.type = R and e2.type = O
        else if (s1s2 == RO)
        {
            estimateBlength_R_O<num_states>(*seq1_region, *seq2_region, total_blength, end_pos, coefficient, coefficient_vec);
        }
        // 2.3. e1.type = R and e2.type = A/C/G/T
        else if (seq1_region->type == TYPE_R)
        {
            estimateBlength_R_ACGT(*seq1_region, seq2_region->type, total_blength, end_pos, coefficient_vec);
        }
        // 3. e1.type = O
        else if (seq1_region->type == TYPE_O)
        {
            estimateBlength_O_X<num_states>(*seq1_region, *seq2_region, total_blength, end_pos, coefficient, coefficient_vec);
        }
        // 4. e1.type = A/C/G/T
        // 4.1. e1.type =  e2.type
        // NHANLT NOTES:
        // coefficient = derivative of log likelihood function wrt t
        // l = log(1 + q_xx * t) ~ q_xx * t
        // => l' = q_xx
        else if (seq1_region->type == seq2_region->type)
            coefficient += model->diagonal_mut_mat[seq1_region->type];
        // e1.type = A/C/G/T and e2.type = O/A/C/G/T
        // 4.2. e1.type = A/C/G/T and e2.type = O
        else if (seq2_region->type == TYPE_O)
        {
            estimateBlength_ACGT_O<num_states>(*seq1_region, *seq2_region, total_blength, coefficient, coefficient_vec);
        }
        // 4.3. e1.type = A/C/G/T and e2.type = R or A/C/G/T
        else
        {
            estimateBlength_ACGT_RACGT(*seq1_region, *seq2_region, total_blength, end_pos, coefficient_vec);
        }
        
        // update pos
        pos = end_pos + 1;
    }
    
    // now optimized branch length based on coefficients
    return estimateBlengthFromCoeffs(coefficient, coefficient_vec);
}

RealNumType cmaple::TreeBase::calculateDerivative(const vector<RealNumType> &coefficient_vec, const RealNumType delta_t)
{
    RealNumType result = 0;
    
    for (RealNumType coefficient : coefficient_vec)
        result += 1.0 / (coefficient + delta_t);
        
    return result;
}

template <const StateType num_states>
void cmaple::TreeBase::handleBlengthChanged(PhyloNode& node, const Index node_index, const RealNumType best_blength)
{
    /*node->length = best_blength;
    node->neighbor->length = node->length;*/
    node.setUpperLength(best_blength);
    
    stack<Index> node_stack;
    node_stack.push(node_index);
    node_stack.push(node.getNeighborIndex(TOP));
    // node_stack.push(node->neighbor);
    updatePartialLh<num_states>(node_stack);
}

template <const StateType num_states>
void cmaple::TreeBase::optimizeBlengthBeforeSeekingSPR(PhyloNode& node, RealNumType &best_blength, RealNumType &best_lh, bool &blength_changed, const std::unique_ptr<SeqRegions>& parent_upper_lr_lh, const std::unique_ptr<SeqRegions>& lower_lh)
{
    RealNumType original_lh = best_lh;
    
    // try different branch lengths for the current node placement (just in case branch length can be improved, in which case it counts both as tree improvment and better strategy to find a new placement).
    if (node.getUpperLength() <= 0)
    {
        best_blength = min_blength;
        best_lh = calculateSubTreePlacementCost<num_states>(parent_upper_lr_lh, lower_lh, best_blength);
    }
    
    // cache best_blength
    const RealNumType cached_blength = best_blength;

    // try shorter branch lengths
    bool found_new_blength = tryShorterNewBranch<&cmaple::TreeBase::calculateSubTreePlacementCost<num_states>>(parent_upper_lr_lh, lower_lh, best_blength, best_lh, double_min_blength);
    
    // try longer branch lengths
    if (!found_new_blength)
        tryLongerNewBranch<&cmaple::TreeBase::calculateSubTreePlacementCost<num_states>>(parent_upper_lr_lh, lower_lh, best_blength, best_lh, half_max_blength);
    
    // update blength_changed
    if (cached_blength != best_blength)
        blength_changed = true;
    
    if (node.getUpperLength() <= 0 && original_lh > best_lh)
        best_lh = original_lh;
}

template <const StateType num_states>
void cmaple::TreeBase::checkAndApplySPR(const RealNumType best_lh_diff, const RealNumType best_blength, const RealNumType best_lh, const Index node_index, PhyloNode& node, const Index best_node_index, const Index parent_node_index, const bool is_mid_node, RealNumType& total_improvement, bool& topology_updated)
{
    const NumSeqsType parent_node_vec = parent_node_index.getVectorIndex();
    const NumSeqsType best_node_vec = best_node_index.getVectorIndex();
    PhyloNode& parent_node = nodes[parent_node_vec];
    if (best_node_vec == parent_node_vec)
    {
        if (cmaple::verbose_mode >= cmaple::VB_DEBUG)
            std::cout << "Strange, re-placement is at same node" << std::endl;
    }
    else if ((best_node_vec == parent_node.getNeighborIndex(LEFT).getVectorIndex() || best_node_vec == parent_node.getNeighborIndex(RIGHT).getVectorIndex()) && is_mid_node) // ((best_node == parent_node->next->neighbor || best_node == parent_node->next->next->neighbor) && is_mid_node)
    {
        if (cmaple::verbose_mode >= cmaple::VB_DEBUG)
            std::cout << "Re-placement is above sibling node" << std::endl;
    }
    else [[likely]]
    {
        // reach the top of a multifurcation, which is the only place in a multifurcatio where placement is allowed.
        NumSeqsType top_polytomy_vec = best_node_vec;

        while (nodes[top_polytomy_vec].getUpperLength() <= 0 && root_vector_index != top_polytomy_vec)// top_polytomy->length <= 0 && top_polytomy!= root)
        {
            // top_polytomy = top_polytomy->neighbor->getTopNode();
            top_polytomy_vec = nodes[top_polytomy_vec].getNeighborIndex(TOP).getVectorIndex();
        }
        
        if (top_polytomy_vec != best_node_vec && cmaple::verbose_mode >= cmaple::VB_DEBUG)
            std::cout << "Strange, placement node not at top of polytomy" << std::endl;
        
        // reach the top of the multifurcation of the parent
        /*Node* parent_top_polytomy = parent_node;
        while (parent_top_polytomy->length <= 0 && parent_top_polytomy != root)
            parent_top_polytomy = parent_top_polytomy->neighbor->getTopNode();*/
        NumSeqsType parent_top_polytomy_vec = parent_node_vec;

        while (nodes[parent_top_polytomy_vec].getUpperLength() <= 0 && root_vector_index != parent_top_polytomy_vec)
            parent_top_polytomy_vec = nodes[parent_top_polytomy_vec].getNeighborIndex(TOP).getVectorIndex();
        
        if (!(parent_top_polytomy_vec == top_polytomy_vec && !is_mid_node))
        {
            total_improvement = best_lh_diff - best_lh;
            
            if (verbose_mode == VB_DEBUG)
                cout << "In improveSubTree() found SPR move with improvement " << total_improvement << endl;
            
            // apply an SPR move
            applySPR<num_states>(node_index, node, best_node_index, is_mid_node, best_blength, best_lh_diff);
            
            topology_updated = true;
        }
    }
}

template <const StateType num_states>
RealNumType cmaple::TreeBase::improveSubTree(const Index node_index, PhyloNode& node, bool short_range_search)
{
    // dummy variables
    ASSERT(node_index.getMiniIndex() == TOP);
    const NumSeqsType vec_index = node_index.getVectorIndex();
    const RealNumType threshold_prob = params->threshold_prob;
    const RealNumType thresh_placement_cost = short_range_search ? params->thresh_placement_cost_short_search : params->thresh_placement_cost;
    RealNumType total_improvement = 0;
    bool blength_changed = false; // true if a branch length has been changed
    
    // we avoid the root node since it cannot be re-placed with SPR moves
    if (root_vector_index != vec_index)
    {
        // evaluate current placement
        const std::unique_ptr<SeqRegions>& parent_upper_lr_lh = getPartialLhAtNode(node.getNeighborIndex(TOP)); // node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
        const std::unique_ptr<SeqRegions>& lower_lh = node.getPartialLh(TOP); // node->getPartialLhAtNode(aln, model, threshold_prob);
        RealNumType best_blength = node.getUpperLength(); // node->length;
        RealNumType best_lh = calculateSubTreePlacementCost<num_states>(parent_upper_lr_lh, lower_lh, best_blength);
        
        // optimize branch length
        if (best_lh < thresh_placement_cost)
            optimizeBlengthBeforeSeekingSPR<num_states>(node, best_blength, best_lh, blength_changed, parent_upper_lr_lh, lower_lh);
           
        // find new placement
        if (best_lh < thresh_placement_cost)
        {
            // now find the best place on the tree where to re-attach the subtree rooted at "node" but to do that we need to consider new vector probabilities after removing the node that we want to replace this is done using findBestParentTopology().
            bool topology_updated = false;
            const Index parent_index = node.getNeighborIndex(TOP);
            // PhyloNode parent_node = nodes[parent_index.getVectorIndex()]; // node->neighbor->getTopNode();
            Index best_node_index;
            RealNumType best_lh_diff = best_lh;
            bool is_mid_node = false;
            RealNumType best_up_lh_diff = MIN_NEGATIVE;
            RealNumType best_down_lh_diff = MIN_NEGATIVE;
            Index best_child_index;
            
            // seek a new placement for the subtree
            seekSubTreePlacement<num_states>(best_node_index, best_lh_diff, is_mid_node, best_up_lh_diff, best_down_lh_diff, best_child_index, short_range_search, node_index, best_blength); // , true, NULL);
            
            // validate the new placement cost
            if (best_lh_diff > params->threshold_prob2)
                outError("Strange, lh cost is positive");
            else if (best_lh_diff < -1e50)
                outError("Likelihood cost is very heavy, this might mean that the reference used is not the same used to generate the input MAPLE file");
            
            if (best_lh_diff + thresh_placement_cost > best_lh)
            {
                // check and apply SPR move
                checkAndApplySPR<num_states>(best_lh_diff, best_blength, best_lh, node_index, node, best_node_index, parent_index, is_mid_node, total_improvement, topology_updated);
                
                if (!topology_updated && blength_changed)
                    handleBlengthChanged<num_states>(node, node_index, best_blength);
            }
            else if (blength_changed)
                handleBlengthChanged<num_states>(node, node_index, best_blength);
        }
        else if (blength_changed)
            handleBlengthChanged<num_states>(node, node_index, best_blength);
    }
    
                        
    return total_improvement;
}

void calculateSubtreeCost_R_R(const SeqRegion& seq1_region, const RealNumType* const &cumulative_rate, RealNumType& total_blength, const PositionType pos, const PositionType end_pos, RealNumType& lh_cost)
{
    if (seq1_region.plength_observation2root >= 0)
    {
        ASSERT(total_blength >= 0);
        total_blength += seq1_region.plength_observation2node;
    }

    // NHANLT NOTE:
    // approximation log(1+x)~x
    if (total_blength > 0)
      lh_cost += total_blength * (cumulative_rate[end_pos + 1] - cumulative_rate[pos]);
}

template <const StateType num_states>
void calculateSubtreeCost_R_O(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType total_blength, const StateType seq1_state, RealNumType& total_factor, const ModelBase* model)
{
    RealNumType tot = 0;
    
    if (seq1_region.plength_observation2root >= 0)
    {
        RealNumType* transposed_mut_mat_row = model->transposed_mut_mat + model->row_index[seq1_state];
        RealNumType* mutation_mat_row = model->mutation_mat;
                            
        for (StateType i = 0; i < num_states; ++i, mutation_mat_row += num_states)
        {
            // NHANLT NOTE: UNSURE
            // tot2: likelihood that we can observe seq1_state elvoving from i at root (account for the fact that the observation might have occurred on the other side of the phylogeny with respect to the root)
            // tot2 = root_freqs[seq1_state] * (1 + mut[seq1_state,seq1_state] * plength_observation2node) + root_freqs[i] * mut[i,seq1_state] * plength_observation2node
            RealNumType tot2 = model->root_freqs[i] * transposed_mut_mat_row[i] * seq1_region.plength_observation2node + (seq1_state == i ? model->root_freqs[i] : 0);
            
            // NHANLT NOTE:
            // tot3: likelihood of i evolves to j
            // tot3 = (1 + mut[i,i] * total_blength) * lh(seq2,i) + mut[i,j] * total_blength * lh(seq2,j)
            RealNumType tot3 = total_blength > 0 ? (total_blength * dotProduct<num_states>(mutation_mat_row, &((*seq2_region.likelihood)[0]))) : 0;
            
            // NHANLT NOTE:
            // tot = tot2 * tot3
            tot += tot2 * (seq2_region.getLH(i) + tot3);
        }
        
        // NHANLT NOTE: UNCLEAR
        // why we need to divide tot by root_freqs[seq1_state]
        tot *= model->inverse_root_freqs[seq1_state];
    }
    else
    {
        // NHANLT NOTE:
        // (1 + mut[seq1_state,seq1_state] * total_blength) * lh(seq2,seq1_state) + mut[seq1_state,j] * total_blength * lh(seq2,j)
        if (total_blength > 0)
        {
            const RealNumType* mutation_mat_row = model->mutation_mat + model->row_index[seq1_state];
            tot += dotProduct<num_states>(mutation_mat_row, &((*seq2_region.likelihood)[0]));
            tot *= total_blength;
        }
        tot += seq2_region.getLH(seq1_state);
    }
    total_factor *= tot;
}

bool calculateSubtreeCost_R_ACGT(const SeqRegion& seq1_region, const RealNumType total_blength, const StateType seq1_state, const StateType seq2_state, RealNumType& total_factor, const ModelBase* model)
{
    if (seq1_region.plength_observation2root >= 0)
    {
        if (total_blength > 0)
        {
            // NHANLT NOTE: UNSURE
            // seq1_state_evolves_seq2_state = (1) the likelihood that seq1_state evolves to seq2_state * (2) the likelihood that seq1_state unchanges from the observing position
            // (1) = model->mutation_mat[model->row_index[seq1_state] + seq2_state] * total_blength
            // (2) = (1.0 + model->diagonal_mut_mat[seq1_state] * seq1_region.plength_observation2node)
            RealNumType seq1_state_evolves_seq2_state = model->mutation_mat[model->row_index[seq1_state] + seq2_state] * total_blength * (1.0 + model->diagonal_mut_mat[seq1_state] * seq1_region.plength_observation2node);
            
            // NHANLT NOTE: UNCLEAR
            // consider the inverse process of the above
            // seq2_state_evolves_seq1_state = (1) the likelihood that seq2_state evolves to seq1_state * (2) the likelihood that seq2_state unchanges from the observing position
            // (1) = root_freqs[seq2_state] / root_freqs[seq1_state] * mutation_mat[model->row_index[seq2_state] + seq1_state] * seq1_region.plength_observation2node
            // (2) = (1.0 + model->diagonal_mut_mat[seq2_state] * total_blength)
            RealNumType seq2_state_evolves_seq1_state = model->freqi_freqj_qij[model->row_index[seq2_state] + seq1_state] * seq1_region.plength_observation2node * (1.0 + model->diagonal_mut_mat[seq2_state] * total_blength);
            
            total_factor *= seq1_state_evolves_seq2_state + seq2_state_evolves_seq1_state;
        }
        // NHANLT NOTE:
        // the same as above but total_blength = 0 then we simplify the formula to save the runtime (avoid multiplying with 0)
        else
            total_factor *= model->freqi_freqj_qij[model->row_index[seq2_state] + seq1_state] * seq1_region.plength_observation2node;
    }
    // NHANLT NOTE:
    // add the likelihood that seq1_state evoles to seq2_state = mut[seq1_state,seq2_state] * total_blength
    else if (total_blength > 0)
        total_factor *= model->mutation_mat[model->row_index[seq1_state] + seq2_state] * total_blength;
    else
        return false; // return MIN_NEGATIVE;

    // no error
    return true;
}

template <const StateType num_states>
void calculateSubtreeCost_O_O(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType total_blength, RealNumType& total_factor, const ModelBase* model)
{
    if (total_blength > 0)
    {
      total_factor *= matrixEvolve<num_states>(&((*seq1_region.likelihood)[0]), &((*seq2_region.likelihood)[0]), model->mutation_mat, total_blength);
    }
    // NHANLT NOTE:
    // the same as above but total_blength = 0 then we simplify the formula to save the runtime (avoid multiplying with 0)
    else
    {
      total_factor *= dotProduct<num_states>(&((*seq1_region.likelihood)[0]), &((*seq2_region.likelihood)[0]));
    }
}

template <const StateType num_states>
void calculateSubtreeCost_O_RACGT(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType total_blength, const PositionType end_pos, RealNumType& total_factor, const AlignmentBase* aln, const ModelBase* model)
{
    StateType seq2_state = seq2_region.type;
    if (seq2_state == TYPE_R)
        seq2_state = aln->ref_seq[end_pos];
    
    if (total_blength > 0)
    {
        // NHANLT NOTE:
        // tot2: likelihood of i evolves to seq2_state
        // tot2 = (1 + mut[seq2_state,seq2_state] * total_blength) * lh(seq1,seq2_state) + lh(seq1,i) * mut[i,seq2_state] * total_blength
        RealNumType* transposed_mut_mat_row = model->transposed_mut_mat + model->row_index[seq2_state];
        RealNumType tot2 = dotProduct<num_states>(&((*seq1_region.likelihood)[0]), transposed_mut_mat_row);
        total_factor *= seq1_region.getLH(seq2_state) + total_blength * tot2;
    }
    // NHANLT NOTE:
    // the same as above but total_blength = 0 then we simplify the formula to save the runtime (avoid multiplying with 0)
    else
        total_factor *= seq1_region.getLH(seq2_state);
}

void calculateSubtreeCost_identicalACGT(const SeqRegion& seq1_region, RealNumType& total_blength, RealNumType& lh_cost, const ModelBase* model)
{
    if (seq1_region.plength_observation2root >= 0)
        total_blength += seq1_region.plength_observation2node;
    
    // NHANLT NOTE:
    // the likelihood that seq1_state unchanges
    if (total_blength > 0)
        lh_cost += model->diagonal_mut_mat[seq1_region.type] * total_blength;
}

template <const StateType num_states>
void calculateSubtreeCost_ACGT_O(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType total_blength, RealNumType& total_factor, const ModelBase* model)
{
    StateType seq1_state = seq1_region.type;
    if (seq1_region.plength_observation2root >= 0)
    {
        RealNumType* transposed_mut_mat_row = model->transposed_mut_mat + model->row_index[seq1_state];
        RealNumType* mutation_mat_row = model->mutation_mat;
        RealNumType tot = matrixEvolveRoot<num_states>(&((*seq2_region.likelihood)[0]), seq1_state,
          model->root_freqs, transposed_mut_mat_row, mutation_mat_row, total_blength, seq1_region.plength_observation2node);
        // NHANLT NOTE: UNCLEAR
        // why we need to divide tot by root_freqs[seq1_state]
        total_factor *= (tot * model->inverse_root_freqs[seq1_state]);
    }
    else
    {
        RealNumType* mutation_mat_row = model->mutation_mat + model->row_index[seq1_state];
        
        // NHANLT NOTE:
        // tot = the likelihood of seq1_state evolving to j
        // (1 + mut[seq1_state,seq1_state] * total_blength) * lh(seq2,seq1_state) + mut[seq1_state,j] * total_blength * lh(seq2,j)
        RealNumType tot = dotProduct<num_states>(mutation_mat_row, &((*seq2_region.likelihood)[0]));
        tot *= total_blength;
        tot += seq2_region.getLH(seq1_state);
        total_factor *= tot;
    }
}

bool calculateSubtreeCost_ACGT_RACGT(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType total_blength, const PositionType end_pos, RealNumType& total_factor, const AlignmentBase* aln, const ModelBase* model)
{
    StateType seq1_state = seq1_region.type;
    StateType seq2_state = seq2_region.type;
    if (seq2_state == TYPE_R)
        seq2_state = aln->ref_seq[end_pos];
    
    if (seq1_region.plength_observation2root >= 0)
    {
        if (total_blength > 0)
        {
            // NHANLT NOTE: UNSURE
            // seq1_state_evolves_seq2_state = (1) the likelihood that seq1_state evolves to seq2_state * (2) the likelihood that seq1_state unchanges from the observing position
            // (1) = model->mutation_mat[model->row_index[seq1_state] + seq2_state] * total_blength
            // (2) = (1.0 + model->diagonal_mut_mat[seq1_state] * seq1_region.plength_observation2node)
            RealNumType seq1_state_evolves_seq2_state = model->mutation_mat[model->row_index[seq1_state] + seq2_state] * total_blength * (1.0 + model->diagonal_mut_mat[seq1_state] * seq1_region.plength_observation2node);
            
            // NHANLT NOTE: UNCLEAR
            // consider the inverse process of the above
            // seq2_state_evolves_seq1_state = (1) the likelihood that seq2_state evolves to seq1_state * (2) the likelihood that seq2_state unchanges from the observing position
            // (1) = root_freqs[seq2_state] / root_freqs[seq1_state] * mutation_mat[model->row_index[seq2_state] + seq1_state] * seq1_region.plength_observation2node
            // (2) = (1.0 + model->diagonal_mut_mat[seq2_state] * total_blength)
            RealNumType seq2_state_evolves_seq1_state = model->freqi_freqj_qij[model->row_index[seq2_state] + seq1_state] * seq1_region.plength_observation2node * (1.0 + model->diagonal_mut_mat[seq2_state] * total_blength);
            
            total_factor *= seq1_state_evolves_seq2_state + seq2_state_evolves_seq1_state;
        }
        // NHANLT NOTE:
        // the same as above but total_blength = 0 then we simplify the formula to save the runtime (avoid multiplying with 0)
        else
            total_factor *= model->freqi_freqj_qij[model->row_index[seq2_state] + seq1_state] * seq1_region.plength_observation2node;
    }
    // NHANLT NOTE:
    // add the likelihood that seq1_state evoles to seq2_state = mut[seq1_state,seq2_state] * total_blength
    else if (total_blength > 0)
        total_factor *= model->mutation_mat[model->row_index[seq1_state] + seq2_state] * total_blength;
    else
        return false; // return MIN_NEGATIVE;

    // no error
    return true;
}

// this implementation derives from appendProbNode
template <const StateType num_states>
RealNumType cmaple::TreeBase::calculateSubTreePlacementCost(const std::unique_ptr<SeqRegions>& parent_regions, const std::unique_ptr<SeqRegions>& child_regions, const RealNumType blength)
{
    // NHANLT BUG FIXED -> not sure it's the best way to due with cases where parent_regions is null
    if (!parent_regions)
        return MIN_NEGATIVE;
    
    // 55% of runtime
    // init dummy variables
    RealNumType lh_cost = 0;
    PositionType pos = 0;
    RealNumType total_factor = 1;
    const SeqRegions& seq1_regions = *parent_regions;
    const SeqRegions& seq2_regions = *child_regions;
    size_t iseq1 = 0;
    size_t iseq2 = 0;
    const PositionType seq_length = aln->ref_seq.size();
    
    while (pos < seq_length)
    {
        PositionType end_pos;
        RealNumType total_blength;
        
        // get the next shared segment in the two sequences
        SeqRegions::getNextSharedSegment(pos, seq1_regions, seq2_regions, iseq1, iseq2, end_pos);
        const auto* const seq1_region = &seq1_regions[iseq1];
        const auto* const seq2_region = &seq2_regions[iseq2];

        // 1. e1.type = N || e2.type = N
        if ((seq2_region->type == TYPE_N) + (seq1_region->type == TYPE_N))
        {
            pos = end_pos + 1;
            continue;
        }
       
        // e1.type != N && e2.type != N
        const DoubleState s1s2 = (DoubleState(seq1_region->type) << 8) | seq2_region->type;
        
        // total_blength will be here the total length from the root or from the upper node, down to the down node.
        if (seq1_region->plength_observation2root >= 0)
            total_blength = seq1_region->plength_observation2root + (blength >= 0 ? blength : 0);
        else if (seq1_region->plength_observation2node >= 0)
            total_blength = seq1_region->plength_observation2node + (blength >= 0 ? blength : 0);
        else
            total_blength = blength;
            
        if (seq2_region->plength_observation2node >= 0)
            total_blength = (total_blength > 0 ? total_blength : 0) + seq2_region->plength_observation2node;
        
        
        //assert(total_blength >= 0); // can be -1 ..

        // 2.1. e1.type = R and e2.type = R
        if (s1s2 == RR) [[likely]]
        {
            calculateSubtreeCost_R_R(*seq1_region, model->cumulative_rate, total_blength, pos, end_pos, lh_cost);
        }
        // 2.2. e1.type = R and e2.type = O
        else if (s1s2 == RO)
        {
            calculateSubtreeCost_R_O<num_states>(*seq1_region, *seq2_region, total_blength, aln->ref_seq[end_pos], total_factor, model);
        }
        // 2.3. e1.type = R and e2.type = A/C/G/T
        else if (seq1_region->type == TYPE_R)
        {
            if (!calculateSubtreeCost_R_ACGT(*seq1_region, total_blength, aln->ref_seq[end_pos], seq2_region->type, total_factor, model)) return MIN_NEGATIVE;
        }
        // 3. e1.type = O
        // 3.1. e1.type = O and e2.type = O
        else if (s1s2 == OO)
        {
            calculateSubtreeCost_O_O<num_states>(*seq1_region, *seq2_region, total_blength, total_factor, model);
        }
        // 3.2. e1.type = O and e2.type = R or A/C/G/T
        else if (seq1_region->type == TYPE_O)
        {
            calculateSubtreeCost_O_RACGT<num_states>(*seq1_region, *seq2_region, total_blength, end_pos, total_factor, aln, model);
        }
        // 4. e1.type = A/C/G/T
        // 4.1. e1.type =  e2.type
        else if (seq1_region->type == seq2_region->type)
        {
            calculateSubtreeCost_identicalACGT(*seq1_region, total_blength, lh_cost, model);
        }
        // e1.type = A/C/G/T and e2.type = O/A/C/G/T
        // 4.2. e1.type = A/C/G/T and e2.type = O
        else if (seq2_region->type == TYPE_O)
        {
            calculateSubtreeCost_ACGT_O<num_states>(*seq1_region, *seq2_region, total_blength, total_factor, model);
        }
        // 4.3. e1.type = A/C/G/T and e2.type = R or A/C/G/T
        else
        {
            if (!calculateSubtreeCost_ACGT_RACGT(*seq1_region, *seq2_region, total_blength, end_pos, total_factor, aln, model)) return MIN_NEGATIVE;
        }
         
        // avoid underflow on total_factor
        // approximately update lh_cost and total_factor
        if (total_factor <= MIN_CARRY_OVER)
        {
            if (total_factor < MIN_POSITIVE)
                return MIN_NEGATIVE;
            
            //lh_cost += log(total_factor);
            //total_factor = 1.0;
            total_factor *= MAX_POSITIVE;
            lh_cost -= LOG_MAX_POSITIVE;
        }
        
        // update pos
        pos = end_pos + 1;
    }
    
    return lh_cost + log(total_factor);
}

void calculateSampleCost_R_R(const SeqRegion& seq1_region, const RealNumType* const &cumulative_rate, const RealNumType blength, const PositionType pos, const PositionType end_pos, RealNumType& lh_cost)
{
    if (seq1_region.plength_observation2node < 0 && seq1_region.plength_observation2root < 0)
        lh_cost += blength * (cumulative_rate[end_pos + 1] - cumulative_rate[pos]);
    else
    {
        RealNumType total_blength = blength + seq1_region.plength_observation2node;
        if (seq1_region.plength_observation2root < 0)
            lh_cost += total_blength * (cumulative_rate[end_pos + 1] - cumulative_rate[pos]);
        else
            // here contribution from root frequency gets added and subtracted so it's ignored
            lh_cost += (total_blength + seq1_region.plength_observation2root) * (cumulative_rate[end_pos + 1] - cumulative_rate[pos]);
    }
}

template <const StateType num_states>
void calculateSampleCost_R_O(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType blength, const StateType seq1_state, RealNumType& lh_cost, RealNumType& total_factor, const ModelBase* model)
{
    if (seq1_region.plength_observation2root >= 0)
    {
        RealNumType total_blength = seq1_region.plength_observation2root + blength;
        
        if (seq2_region.getLH(seq1_state) > 0.1)
        {
            total_blength += seq1_region.plength_observation2node;
            
            // here contribution from root frequency can also be also ignored
            lh_cost += model->diagonal_mut_mat[seq1_state] * total_blength;
        }
        else
        {
            RealNumType tot = 0;
            RealNumType* freq_j_transposed_ij_row = model->freq_j_transposed_ij + model->row_index[seq1_state];
            RealNumType* mutation_mat_row = model->mutation_mat;
                                
            for (StateType i = 0; i < num_states; ++i, mutation_mat_row += num_states)
            {
                RealNumType tot2 = freq_j_transposed_ij_row[i] * seq1_region.plength_observation2node + ((seq1_state == i) ? model->root_freqs[i] : 0);
                RealNumType tot3 = ((seq2_region.getLH(i) > 0.1) ? 1 : 0) + sumMutationByLh<num_states>(&(*seq2_region.likelihood)[0], mutation_mat_row);
                
                tot += tot2 * tot3 * total_blength;
            }
            
            total_factor *= tot * model->inverse_root_freqs[seq1_state];
        }
    }
    else
    {
        if (seq2_region.getLH(seq1_state) > 0.1)
        {
            if (seq1_region.plength_observation2node >= 0)
                lh_cost += model->diagonal_mut_mat[seq1_state] * (blength + seq1_region.plength_observation2node);
            else
                lh_cost += model->diagonal_mut_mat[seq1_state] * blength;
        }
        else
        {
            RealNumType tot = 0;
            RealNumType* mutation_mat_row = model->mutation_mat + model->row_index[seq1_state];
            
            tot += sumMutationByLh<num_states>(&(*seq2_region.likelihood)[0], mutation_mat_row);
            
            if (seq1_region.plength_observation2node >= 0)
                total_factor *= tot * (blength + seq1_region.plength_observation2node);
            else
                total_factor *= tot * blength;
        }
    }
}

void calculateSampleCost_R_ACGT(const SeqRegion& seq1_region, const RealNumType blength, const StateType seq1_state, const StateType seq2_state, RealNumType& total_factor, const ModelBase* model)
{
    if (seq1_region.plength_observation2root >= 0)
    {
        // TODO: can cache model->mutation_mat[model->row_index[seq1_state] * model->diagonal_mut_mat[seq1_state]
        // TODO: can cache  model->freqi_freqj_qij[model->row_index[seq2_state] + seq1_state] * model->diagonal_mut_mat[seq2_state]
        RealNumType seq1_state_evolves_seq2_state = model->mutation_mat[model->row_index[seq1_state] + seq2_state] * blength * (1.0 + model->diagonal_mut_mat[seq1_state] * seq1_region.plength_observation2node);
        
        RealNumType seq2_state_evolves_seq1_state = model->freqi_freqj_qij[model->row_index[seq2_state] + seq1_state] * seq1_region.plength_observation2node * (1.0 + model->diagonal_mut_mat[seq2_state] * (blength + seq1_region.plength_observation2root));
                                                                                                                                                                                    
        total_factor *= seq1_state_evolves_seq2_state + seq2_state_evolves_seq1_state;
    }
    else
    {
        total_factor *= model->mutation_mat[model->row_index[seq1_state] + seq2_state] * (blength + (seq1_region.plength_observation2node < 0 ? 0 : seq1_region.plength_observation2node));
    }
}

template <const StateType num_states>
void calculateSampleCost_O_O(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType blength, RealNumType& total_factor, const ModelBase* model)
{
    RealNumType blength13 = blength;
    if (seq1_region.plength_observation2node >= 0)
    {
        blength13 = seq1_region.plength_observation2node;
        if (blength > 0)
            blength13 += blength;
    }
    
    RealNumType tot = 0;
    
    RealNumType* mutation_mat_row = model->mutation_mat;
                        
    for (StateType i = 0; i < num_states; ++i, mutation_mat_row += num_states)
    {
        RealNumType tot2 = blength13 * sumMutationByLh<num_states>(&(*seq2_region.likelihood)[0], mutation_mat_row);
        
        tot += (tot2 + (seq2_region.getLH(i) > 0.1 ? 1 : 0)) * seq1_region.getLH(i);
    }
        
    total_factor *= tot;
}

template <const StateType num_states>
void calculateSampleCost_O_RACGT(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType blength, const PositionType end_pos, RealNumType& total_factor, const AlignmentBase* aln, const ModelBase* model)
{
    RealNumType blength13 = blength;
    if (seq1_region.plength_observation2node >= 0)
    {
        blength13 = seq1_region.plength_observation2node;
        if (blength > 0)
            blength13 += blength;
    }
    
    StateType seq2_state = seq2_region.type;
    if (seq2_state == TYPE_R)
        seq2_state = aln->ref_seq[end_pos];
    
    RealNumType *transposed_mut_mat_row = model->transposed_mut_mat + model->row_index[seq2_state];
    RealNumType tot2 = dotProduct<num_states>(transposed_mut_mat_row, &((*seq1_region.likelihood)[0]));
    total_factor *= seq1_region.getLH(seq2_state) + blength13 * tot2;
}

void calculateSampleCost_identicalACGT(const SeqRegion& seq1_region, const RealNumType blength, RealNumType& lh_cost, const ModelBase* model)
{
    RealNumType total_blength = blength;
    total_blength += (seq1_region.plength_observation2node < 0 ? 0 : seq1_region.plength_observation2node);
    total_blength += (seq1_region.plength_observation2root < 0 ? 0 : seq1_region.plength_observation2root);

    lh_cost += model->diagonal_mut_mat[seq1_region.type] * total_blength;
}

template <const StateType num_states>
void calculateSampleCost_ACGT_O(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType blength, RealNumType& lh_cost, RealNumType& total_factor, const ModelBase* model)
{
    StateType seq1_state = seq1_region.type;
    RealNumType tot = 0.0;
    
    if (seq1_region.plength_observation2root >= 0)
    {
        RealNumType blength15 = blength + seq1_region.plength_observation2root;
        
        if (seq2_region.getLH(seq1_state) > 0.1)
            lh_cost += model->diagonal_mut_mat[seq1_state] * (blength15 + seq1_region.plength_observation2node);
        else
        {
            RealNumType* freq_j_transposed_ij_row = model->freq_j_transposed_ij + model->row_index[seq1_state];
            RealNumType* mutation_mat_row = model->mutation_mat;
                                
            for (StateType i = 0; i < num_states; ++i, mutation_mat_row += num_states)
            {
                RealNumType tot2 = freq_j_transposed_ij_row[i] * seq1_region.plength_observation2node + ((seq1_state == i) ? model->root_freqs[i] : 0);
                    
                RealNumType tot3 = sumMutationByLh<num_states>(&(*seq2_region.likelihood)[0], mutation_mat_row);
                
                tot += tot2 * blength15 * tot3 + (seq2_region.getLH(i) > 0.1 ? tot2 : 0);
            }
            
            total_factor *= (tot * model->inverse_root_freqs[seq1_state]);
        }
    }
    else
    {
        RealNumType tmp_blength = blength + (seq1_region.plength_observation2node < 0 ? 0 : seq1_region.plength_observation2node);
        if (seq2_region.getLH(seq1_state) > 0.1)
            lh_cost += model->diagonal_mut_mat[seq1_state] * tmp_blength;
        else
        {
            RealNumType* mutation_mat_row = model->mutation_mat + model->row_index[seq1_state];
            tot += sumMutationByLh<num_states>(&(*seq2_region.likelihood)[0], mutation_mat_row);
            
            total_factor *= tot * tmp_blength;
        }
    }
}

void calculateSampleCost_ACGT_RACGT(const SeqRegion& seq1_region, const SeqRegion& seq2_region, const RealNumType blength, const PositionType end_pos, RealNumType& total_factor, const AlignmentBase* aln, const ModelBase* model)
{
    StateType seq1_state = seq1_region.type;
    StateType seq2_state = seq2_region.type;
    if (seq2_state == TYPE_R)
        seq2_state = aln->ref_seq[end_pos];
    
    if (seq1_region.plength_observation2root >= 0)
    {
        // here we ignore contribution of non-parsimonious mutational histories
        RealNumType seq1_state_evoloves_seq2_state = model->mutation_mat[model->row_index[seq1_state] + seq2_state] * (blength + seq1_region.plength_observation2root) * (1.0 + model->diagonal_mut_mat[seq1_state] * seq1_region.plength_observation2node);
        
        RealNumType seq2_state_evolves_seq1_state = model->freqi_freqj_qij[model->row_index[seq2_state] + seq1_state] * seq1_region.plength_observation2node * (1.0 + model->diagonal_mut_mat[seq2_state] * (blength + seq1_region.plength_observation2root));
        
        total_factor *= (seq1_state_evoloves_seq2_state + seq2_state_evolves_seq1_state);
    }
    else
    {
        RealNumType tmp_blength = ((seq1_region.plength_observation2node < 0) ? blength : blength + seq1_region.plength_observation2node);
        
        total_factor *= model->mutation_mat[model->row_index[seq1_state] + seq2_state] * tmp_blength;
    }
}

// this implementation derives from appendProb
template <const StateType num_states>
RealNumType cmaple::TreeBase::calculateSamplePlacementCost(const std::unique_ptr<SeqRegions>& parent_regions, const std::unique_ptr<SeqRegions>& child_regions, const RealNumType input_blength)
{
    // NHANLT BUG FIXED -> not sure it's the best way to due with cases where parent_regions is null
    if (!parent_regions)
        return MIN_NEGATIVE;
    
    // 10% of total runtime
    // init dummy variables
    RealNumType lh_cost = 0;
    PositionType pos = 0;
    RealNumType total_factor = 1;
    const SeqRegions& seq1_regions = *parent_regions;
    const SeqRegions& seq2_regions = *child_regions;
    size_t iseq1 = 0;
    size_t iseq2 = 0;
    RealNumType blength = input_blength;
    if (blength < 0) blength = 0;
    const PositionType seq_length = aln->ref_seq.size();
    
    while (pos < seq_length)
    {
        PositionType end_pos;
        RealNumType total_blength = blength;
        
        // get the next shared segment in the two sequences
        SeqRegions::getNextSharedSegment(pos, seq1_regions, seq2_regions, iseq1, iseq2, end_pos);
        const auto* seq1_region = &seq1_regions[iseq1];
        const auto* seq2_region = &seq2_regions[iseq2]; 
        
        // 1. e1.type = N || e2.type = N
        if ((seq2_region->type == TYPE_N) + (seq1_region->type == TYPE_N))
        {
            pos = end_pos + 1;
            continue;
        }
        
        // e1.type != N && e2.type != N
        // A,C,G,T
        // R -> same as the reference
        // N -> gaps
        // O -> vector of 4 probability to observe A, C, G, T
        const DoubleState s1s2 = (DoubleState(seq1_region->type) << 8) | seq2_region->type;
        
        // 2. e1.type = R
        // 2.1. e1.type = R and e2.type = R
        if (s1s2 == RR) [[likely]]
        {
            calculateSampleCost_R_R(*seq1_region, model->cumulative_rate, blength, pos, end_pos, lh_cost);
        }
        // 2.2. e1.type = R and e2.type = O
        else if (s1s2 == RO)
        {
            calculateSampleCost_R_O<num_states>(*seq1_region, *seq2_region, blength, aln->ref_seq[end_pos], lh_cost, total_factor, model);
        }
        // 2.3. e1.type = R and e2.type = A/C/G/T
        else if (seq1_region->type == TYPE_R)
        {
            calculateSampleCost_R_ACGT(*seq1_region, blength, aln->ref_seq[end_pos], seq2_region->type, total_factor, model);
        }
        // 3. e1.type = O
        // 3.1. e1.type = O and e2.type = O
        else if (s1s2 == OO)
        {
            calculateSampleCost_O_O<num_states>(*seq1_region, *seq2_region, blength, total_factor, model);
        }
        // 3.2. e1.type = O and e2.type = R or A/C/G/T
        else if (seq1_region->type == TYPE_O)
        {
            calculateSampleCost_O_RACGT<num_states>(*seq1_region, *seq2_region, blength, end_pos, total_factor, aln, model);
        }
        // 4. e1.type = A/C/G/T
        // 4.1. e1.type =  e2.type
        else if (seq1_region->type == seq2_region->type)
        {
            calculateSampleCost_identicalACGT(*seq1_region, blength, lh_cost, model);
        }
        // e1.type = A/C/G/T and e2.type = O/A/C/G/T
        // 4.2. e1.type = A/C/G/T and e2.type = O
        else if (seq2_region->type == TYPE_O)
        {
            calculateSampleCost_ACGT_O<num_states>(*seq1_region, *seq2_region, blength, lh_cost, total_factor, model);
        }
        // 4.3. e1.type = A/C/G/T and e2.type = R or A/C/G/T
        else
        {
            calculateSampleCost_ACGT_RACGT(*seq1_region, *seq2_region, blength, end_pos, total_factor, aln, model);
        }
         
        // avoid underflow on total_factor
        // approximately update lh_cost and total_factor
        if (total_factor <= MIN_CARRY_OVER)
        {
            if (total_factor < MIN_POSITIVE)
                return MIN_NEGATIVE;
            
            //lh_cost += log(total_factor);
            //total_factor = 1.0;
            total_factor *= MAX_POSITIVE;
            lh_cost -= LOG_MAX_POSITIVE;
        }
        
        // update pos
        pos = end_pos + 1;
    }
    
    return lh_cost + log(total_factor);
}

template <const StateType num_states>
void cmaple::TreeBase::updateZeroBlength(const Index index, PhyloNode& node, std::stack<Index> &node_stack)
{
    // get the top node in the phylo-node
    /*Node* top_node = node->getTopNode();
    ASSERT(top_node);
    SeqRegions* upper_left_right_regions = top_node->neighbor->getPartialLhAtNode(aln, model, threshold_prob);
    SeqRegions* lower_regions = top_node->getPartialLhAtNode(aln, model, threshold_prob);*/
    const NumSeqsType node_vec_index = index.getVectorIndex();
    const Index parent_index = node.getNeighborIndex(TOP);
    PhyloNode& parent_node = nodes[parent_index.getVectorIndex()];
    const std::unique_ptr<SeqRegions>& upper_left_right_regions = parent_node.getPartialLh(parent_index.getMiniIndex());
    const std::unique_ptr<SeqRegions>& lower_regions = node.getPartialLh(TOP);
    
    RealNumType best_lh = calculateSamplePlacementCost<num_states>(upper_left_right_regions, lower_regions, default_blength);
    RealNumType best_length = default_blength;
    
    // try shorter lengths
    bool found_new_best_length = tryShorterNewBranch<&cmaple::TreeBase::calculateSamplePlacementCost<num_states>>(upper_left_right_regions, lower_regions, best_length, best_lh, min_blength);
    
    // try longer lengths
    if (!found_new_best_length)
        tryLongerNewBranch<&cmaple::TreeBase::calculateSamplePlacementCost<num_states>>(upper_left_right_regions, lower_regions, best_length, best_lh, max_blength);
    
    // update best_length
    /*top_node->length = best_length;
    top_node->neighbor->length = best_length;*/
    node.setUpperLength(best_length);
    
    // add current node and its parent to node_stack to for updating partials further from these nodes
    /*top_node->outdated = true;
    top_node->neighbor->getTopNode()->outdated = true;
    node_stack.push(top_node);
    node_stack.push(top_node->neighbor);*/
    node.setOutdated(true);
    parent_node.setOutdated(true);
    node_stack.push(Index(node_vec_index, TOP));
    node_stack.push(parent_index);
}

std::unique_ptr<SeqRegions>& cmaple::TreeBase::getPartialLhAtNode(const Index index)
{
    // may need assert(index.getVectorIndex() < nodes.size());
    return nodes[index.getVectorIndex()].getPartialLh(index.getMiniIndex());
}

template <const StateType num_states>
void cmaple::TreeBase::updateLowerLh(RealNumType& total_lh, std::unique_ptr<SeqRegions>& new_lower_lh, PhyloNode& node, const std::unique_ptr<SeqRegions>& lower_lh_1, const std::unique_ptr<SeqRegions>& lower_lh_2, const Index neighbor_1_index, PhyloNode& neighbor_1, const Index neighbor_2_index, PhyloNode& neighbor_2, const PositionType& seq_length)
{
    lower_lh_1->mergeTwoLowers<num_states>(new_lower_lh, neighbor_1.getUpperLength(), *lower_lh_2, neighbor_2.getUpperLength(), aln, model, params->threshold_prob);
    
    // NHANLT: LOGS FOR DEBUGGING
    /*if (params->debug)
        std::cout << "new_lower_lh " << new_lower_lh->size() << std::endl;*/
     
    // if new_lower_lh is NULL -> we need to update the branch lengths connecting the current node to its children
    if (!new_lower_lh)
    {
        if (neighbor_1.getUpperLength() <= 0) // next_node_1->length <= 0)
        {
            stack<Index> node_stack;
            // NHANLT: note different from original maple
            // updateBLen(nodeList,node,mutMatrix) -> the below codes update from next_node_1 instead of top_node
            updateZeroBlength<num_states>(neighbor_1_index, neighbor_1, node_stack); // next_node_1->neighbor, node_stack, params->threshold_prob);
            updatePartialLh<num_states>(node_stack);
        }
        else if (neighbor_2.getUpperLength() <= 0) // next_node_2->length <= 0)
        {
            stack<Index> node_stack;
            updateZeroBlength<num_states>(neighbor_2_index, neighbor_2, node_stack); // updateZeroBlength<num_states>(next_node_2->neighbor, node_stack, params->threshold_prob);
            updatePartialLh<num_states>(node_stack);
        }
        else
            outError("Strange, branch lengths > 0 but inconsistent lower lh creation in refreshAllLowerLhs()");
    }
    // otherwise, everything is good -> update the lower lh of the current node
    else
        node.setPartialLh(TOP, std::move(new_lower_lh));
}

template <const StateType num_states>
void cmaple::TreeBase::updateLowerLhAvoidUsingUpperLRLh(RealNumType& total_lh, std::unique_ptr<SeqRegions>& new_lower_lh, PhyloNode& node, const std::unique_ptr<SeqRegions>& lower_lh_1, const std::unique_ptr<SeqRegions>& lower_lh_2, const Index neighbor_1_index, PhyloNode& neighbor_1, const Index neighbor_2_index, PhyloNode& neighbor_2, const PositionType& seq_length)
{
    lower_lh_1->mergeTwoLowers<num_states>(new_lower_lh, neighbor_1.getUpperLength(), *lower_lh_2, neighbor_2.getUpperLength(), aln, model, params->threshold_prob);
     
    // if new_lower_lh is NULL -> we need to update the branch lengths connecting the current node to its children
    if (!new_lower_lh)
    {
        if (cmaple::verbose_mode >= cmaple::VB_DEBUG)
         std::cout << "Set a zero branch length to the minmimum branch length " + convertDoubleToString(min_blength) + " to avoid computation error." << std::endl;
        
        // set zero-length to the min_blength
        if (neighbor_1.getUpperLength() <= 0) // next_node_1->length <= 0)
            neighbor_1.setUpperLength(min_blength);
        else if (neighbor_2.getUpperLength() <= 0) // next_node_2->length <= 0)
            neighbor_2.setUpperLength(min_blength);
        else
            outError("Strange, branch lengths > 0 but inconsistent lower lh creation in refreshAllLowerLhs()");
        
        // recompute lower_lh
        lower_lh_1->mergeTwoLowers<num_states>(node.getPartialLh(TOP), neighbor_1.getUpperLength(), *lower_lh_2, neighbor_2.getUpperLength(), aln, model, params->threshold_prob);

    }
    // otherwise, everything is good -> update the lower lh of the current node
    else
        node.setPartialLh(TOP, std::move(new_lower_lh));
}

template <const StateType num_states>
void cmaple::TreeBase::computeLhContribution(RealNumType& total_lh, std::unique_ptr<SeqRegions>& new_lower_lh, PhyloNode& node, const std::unique_ptr<SeqRegions>& lower_lh_1, const std::unique_ptr<SeqRegions>& lower_lh_2, const Index neighbor_1_index, PhyloNode& neighbor_1, const Index neighbor_2_index, PhyloNode& neighbor_2, const PositionType& seq_length)
{
    RealNumType lh_contribution = lower_lh_1->mergeTwoLowers<num_states>(new_lower_lh, neighbor_1.getUpperLength(), *lower_lh_2, neighbor_2.getUpperLength(), aln, model, params->threshold_prob, true);
    total_lh += lh_contribution;
    
    // record the likelihood contribution at this node
    // if likelihood contribution of this node has not yet existed -> add a new one
    if (node.getNodelhIndex() == 0)
    {
        node_lhs.emplace_back(lh_contribution);
        node.setNodeLhIndex(node_lhs.size() - 1);
    }
    // otherwise, update it
    else
        node_lhs[node.getNodelhIndex()].setLhContribution(lh_contribution);
    
    // if new_lower_lh is NULL
    // ASSERT(params.has_value());
    ASSERT(params);
    if (!new_lower_lh)
        outError("Strange, inconsistent lower genome list creation in calculateTreeLh(); old list, and children lists");
    // otherwise, everything is good -> update the lower lh of the current node
    else if (new_lower_lh->areDiffFrom(*node.getPartialLh(TOP), seq_length, num_states, *params))
        outError("Strange, while calculating tree likelihood encountered non-updated lower likelihood!");
}

template <void(cmaple::TreeBase::*task)(RealNumType&, std::unique_ptr<SeqRegions>&, PhyloNode&, const std::unique_ptr<SeqRegions>&, const std::unique_ptr<SeqRegions>&, const Index, PhyloNode&, const Index, PhyloNode&, const PositionType&)>
RealNumType cmaple::TreeBase::performDFS()
{
    // dummy variables
    RealNumType total_lh = 0;
    const PositionType seq_length = aln->ref_seq.size();
    
    // start from root
    Index node_index = Index(root_vector_index, TOP);
    Index last_node_index;
    
    // traverse to the deepest tip, calculate the likelihoods upward from the tips
    while (node_index.getMiniIndex() != UNDEFINED) //node)
    {
        PhyloNode& node = nodes[node_index.getVectorIndex()];
        // we reach a top node by a downward traversing
        if (node_index.getMiniIndex() == TOP) //node->is_top)
        {
            // if the current node is a leaf -> we reach the deepest tip -> traversing upward to calculate the lh of its parent
            if (!node.isInternal()) // node->isLeave())
            {
                /*last_node = node;
                 node = node->neighbor;*/
                last_node_index = node_index;
                node_index = node.getNeighborIndex(TOP);
            }
            // otherwise, keep traversing downward to find the deepest tip
            else
                // node = node->next->neighbor;
                node_index = node.getNeighborIndex(RIGHT);
        }
        // we reach the current node by an upward traversing from its children
        else
        {
            // if we reach the current node by an upward traversing from its first children -> traversing downward to its second children
            if (node.getNeighborIndex(RIGHT) == last_node_index) // node->getTopNode()->next->neighbor == last_node)
            {
                // node = node->getTopNode()->next->next->neighbor;
                node_index = node.getNeighborIndex(LEFT);
            }
            // otherwise, all children of the current node are updated -> update the lower lh of the current node
            else
            {
                // calculate the new lower lh of the current node from its children
                /*Node* top_node = node->getTopNode();
                 Node* next_node_1 = top_node->next;
                 Node* next_node_2 = next_node_1->next;*/
                Index neighbor_1_index = node.getNeighborIndex(RIGHT);
                Index neighbor_2_index = node.getNeighborIndex(LEFT);
                PhyloNode& neighbor_1 = nodes[neighbor_1_index.getVectorIndex()];
                PhyloNode& neighbor_2 = nodes[neighbor_2_index.getVectorIndex()];
                
                std::unique_ptr<SeqRegions> new_lower_lh = nullptr;
                const std::unique_ptr<SeqRegions>& lower_lh_1 = neighbor_1.getPartialLh(TOP); // next_node_1->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob);
                const std::unique_ptr<SeqRegions>& lower_lh_2 = neighbor_2.getPartialLh(TOP); // next_node_2->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob);
                // lower_lh_1->mergeTwoLowers<num_states>(new_lower_lh, next_node_1->length, *lower_lh_2, next_node_2->length, aln, model, params->threshold_prob);
                
                (this->*task)(total_lh, new_lower_lh, node, lower_lh_1, lower_lh_2, neighbor_1_index, neighbor_1, neighbor_2_index, neighbor_2, seq_length);
                
                last_node_index = Index(node_index.getVectorIndex(), TOP);
                node_index = node.getNeighborIndex(TOP);
            }
        }
    }
    
    return total_lh;
}

template <const StateType num_states>
void cmaple::TreeBase::calculate_aRLT(const bool allow_replacing_ML_tree)
{
    // set all nodes outdated
    resetSPRFlags(false, true, true);
    
    // reseve space for node_lhs
    node_lhs.reserve(nodes.size() * 0.5);
    
    // compute the likelihood at root
    RealNumType lh_at_root = nodes[root_vector_index].getPartialLh(TOP)->computeAbsoluteLhAtRoot<num_states>(model);
    
    // LT1 = tree_total_lh = likelihood at root + total likelihood contribution at all internal nodes
    RealNumType tree_total_lh = lh_at_root + performDFS<&cmaple::TreeBase::computeLhContribution<num_states>>();;
    
    // traverse tree to calculate aLRT-SH for each internal branch
    PhyloNode& root = nodes[root_vector_index];
    std::stack<Index> node_stack;
    if (root.isInternal())
    {
        node_stack.push(root.getNeighborIndex(RIGHT));
        node_stack.push(root.getNeighborIndex(LEFT));
    }
    
    while (!node_stack.empty())
    {
        const Index node_index = node_stack.top();
        node_stack.pop();
        const NumSeqsType node_vec = node_index.getVectorIndex();
        PhyloNode& node = nodes[node_vec];
        
        // only consider internal branches
        if (node.isInternal() && node.isOutdated())
        {
            node.setOutdated(false);
            
            const Index child_1_index = node.getNeighborIndex(RIGHT);
            const Index child_2_index = node.getNeighborIndex(LEFT);
            
            // NHANLT: Debug aLRT
            /* if (child_1_index.getVectorIndex() == 3186 || child_1_index.getVectorIndex() == 2841)
                std::cout << "Update" << std::endl;*/
            
            // add its children into node_stack for further traversal
            node_stack.push(child_1_index);
            node_stack.push(child_2_index);
            
            // only compute the aLRT for internal non-zero branches
            if (node.getUpperLength() > 0)
            {
                // show the likelihood contribution at this branch/node
                // std::cout << "Lh contribution of node (node_vec) " << node_vec << " : " << node_lhs[node.getNodelhIndex()].lh_contribution_ << std::endl;
                
                // caculate aLRT
                PhyloNode& child_1 = nodes[child_1_index.getVectorIndex()];
                PhyloNode& child_2 = nodes[child_2_index.getVectorIndex()];
                const Index parent_index = node.getNeighborIndex(TOP);
                PhyloNode& parent = nodes[parent_index.getVectorIndex()];
                const Index sibling_index = parent.getNeighborIndex(parent_index.getFlipMiniIndex());
                PhyloNode& sibling = nodes[sibling_index.getVectorIndex()];
                
                // compute the likelihood differences between each nni neighbor and the current tree
                RealNumType neighbor_2_lh_diff = 0, neighbor_3_lh_diff = 0;
                // if calculateNNILh return false => the current tree was replaced by a newly found ML tree
                if (!calculateNNILh<num_states>(node_stack, neighbor_2_lh_diff, node, child_1, child_2, sibling, parent, parent_index, lh_at_root, allow_replacing_ML_tree))
                {
                    tree_total_lh += neighbor_2_lh_diff;
                    continue;
                }
                // if calculateNNILh return false => the current tree was replaced by a newly found ML tree
                if (!calculateNNILh<num_states>(node_stack, neighbor_3_lh_diff, node, child_2, child_1, sibling, parent, parent_index, lh_at_root, allow_replacing_ML_tree))
                {
                    tree_total_lh += neighbor_3_lh_diff;
                    continue;
                }
                
                // record neighbor_2_lh_diff, neighbor_2_lh_diff
                NodeLh& node_lh = node_lhs[node.getNodelhIndex()];
                node_lh.setLhDiff2(neighbor_2_lh_diff);
                node_lh.setLhDiff3(neighbor_3_lh_diff);
            }
            // return zero for zero-length internal branches
            else
            {
                NodeLh& node_lh = node_lhs[node.getNodelhIndex()];
                node_lh.setLhDiff2(0);
                node_lh.setLhDiff3(0);
            }
        }
    }
    
    // std::cout << std::setprecision(20) << "Expected tree_total_lh: " << tree_total_lh << std::endl;
}

template <const StateType num_states>
void cmaple::TreeBase::calSiteLhDiffRoot(std::vector<RealNumType>& site_lh_diff, std::vector<RealNumType>& site_lh_root_diff, const std::vector<RealNumType>& site_lh_root, std::unique_ptr<SeqRegions>& parent_new_lower_lh, const RealNumType& child_2_new_blength, PhyloNode& current_node, PhyloNode& child_1, PhyloNode& child_2, PhyloNode& sibling, PhyloNode& parent, const Index parent_index)
{
    const RealNumType threshold_prob = params->threshold_prob;
    const RealNumType child_1_blength = child_1.getUpperLength(); // ~new_branch_length
    const std::unique_ptr<SeqRegions>& child_1_lower_regions = child_1.getPartialLh(TOP); // ~subtree_regions
    const std::unique_ptr<SeqRegions>& sibling_lower_lh = sibling.getPartialLh(TOP);
    const std::unique_ptr<SeqRegions>& child_2_lower_lh = child_2.getPartialLh(TOP);
    std::unique_ptr<SeqRegions> new_parent_new_lower_lh = nullptr;
    std::unique_ptr<SeqRegions> tmp_lower_lh = nullptr;
    const PositionType seq_length = aln->ref_seq.size();
    std::vector<RealNumType> site_lh_diff_old(seq_length, 0);
    
    // 2. estimate x ~ the length of the new branch connecting the parent and the new_parent nodes
    RealNumType parent_new_blength = default_blength;
    // merge 2 lower vector into one
    RealNumType best_parent_lh = parent_new_lower_lh->mergeTwoLowers<num_states>(new_parent_new_lower_lh, parent_new_blength, *child_1_lower_regions, child_1_blength, aln, model, threshold_prob, true);
    best_parent_lh += new_parent_new_lower_lh->computeAbsoluteLhAtRoot<num_states>(model);
    // Try shorter branch lengths at root
    tryShorterBranchAtRoot<num_states>(child_1_lower_regions, parent_new_lower_lh, new_parent_new_lower_lh, parent_new_blength, best_parent_lh, child_1_blength);
    
    // 3. estimate the length of the new branch connecting child_1 to the new_parent node
    RealNumType child_1_new_blength = child_1_blength;
    estimateLengthNewBranchAtRoot<num_states>(child_1_lower_regions, parent_new_lower_lh, new_parent_new_lower_lh, child_1_new_blength, best_parent_lh, parent_new_blength, double_min_blength, child_1_blength <= 0);
    
    // 4. compute the new lower_lh of the new_parent node
    parent_new_lower_lh->mergeTwoLowers<num_states>(new_parent_new_lower_lh, parent_new_blength, *child_1_lower_regions, child_1_new_blength, aln, model, threshold_prob, false);
    
    // 5. compute the new upper_left/right lh for the parent and the new parent nodes
    // 5.1. for the new parent node
    std::unique_ptr<SeqRegions> new_parent_new_upper_lr_1 = nullptr;
    child_1_lower_regions->computeTotalLhAtRoot<num_states>(new_parent_new_upper_lr_1, model, child_1_new_blength);
    std::unique_ptr<SeqRegions> new_parent_new_upper_lr_2 = nullptr;
    parent_new_lower_lh->computeTotalLhAtRoot<num_states>(new_parent_new_upper_lr_2, model, parent_new_blength);
    // 5.2. for the parent node
    std::unique_ptr<SeqRegions> parent_new_upper_lr_1 = nullptr;
    new_parent_new_upper_lr_1->mergeUpperLower<num_states>(parent_new_upper_lr_1, parent_new_blength, *sibling_lower_lh, sibling.getUpperLength(), aln, model, threshold_prob);
    std::unique_ptr<SeqRegions> parent_new_upper_lr_2 = nullptr;
    new_parent_new_upper_lr_1->mergeUpperLower<num_states>(parent_new_upper_lr_2, parent_new_blength, *child_2_lower_lh, child_2_new_blength, aln, model, threshold_prob);
    
    // 6. re-optimize the lengths of the upper branches of
    // 6.1. child_2
    const RealNumType child_2_best_blength = estimateBranchLengthWithCheck<num_states>(parent_new_upper_lr_1, child_2_lower_lh, child_2_new_blength);
    // 6.2. sibling
    const RealNumType sibling_best_blength = estimateBranchLengthWithCheck<num_states>(parent_new_upper_lr_2, sibling_lower_lh, sibling.getUpperLength());
    // 6.3. parent
    const RealNumType parent_best_blength = estimateBranchLengthWithCheck<num_states>(new_parent_new_upper_lr_1, parent_new_lower_lh, parent_new_blength);
    // 6.4. child_1
    const RealNumType child_1_best_blength = estimateBranchLengthWithCheck<num_states>(new_parent_new_upper_lr_2, child_1_lower_regions, child_1_new_blength);
    
    // 7. caculate the likelihood contribution changed at
    // 7.1. the parent node
    child_2_lower_lh->calculateSiteLhContributions<num_states>(site_lh_diff, parent_new_lower_lh, child_2_best_blength, *sibling_lower_lh, sibling_best_blength, aln, model, threshold_prob);
    child_2_lower_lh->calculateSiteLhContributions<num_states>(site_lh_diff_old, tmp_lower_lh, child_2.getUpperLength(), *child_1_lower_regions, child_1.getUpperLength(), aln, model, threshold_prob);
    // 7.2. the new_parent node
    parent_new_lower_lh->calculateSiteLhContributions<num_states>(site_lh_diff, new_parent_new_lower_lh, parent_best_blength, *child_1_lower_regions, child_1_best_blength, aln, model, threshold_prob);
    current_node.getPartialLh(TOP)->calculateSiteLhContributions<num_states>(site_lh_diff_old, tmp_lower_lh, current_node.getUpperLength(), *sibling_lower_lh, sibling.getUpperLength(), aln, model, threshold_prob);
    // 7.3 the absolute likelihood at root
    new_parent_new_lower_lh->computeSiteLhAtRoot<num_states>(site_lh_root_diff, model);
    // update site_lh_root_diff
    for (PositionType j = 0; j < seq_length; ++j)
        site_lh_root_diff[j] -= site_lh_root[j];
    
    // calculate the site-lh differences between the new and the old ones
    for (PositionType j = 0; j < seq_length; ++j)
        site_lh_diff[j] -= site_lh_diff_old[j];
}

template <const StateType num_states>
void cmaple::TreeBase::calSiteLhDiffNonRoot(std::vector<RealNumType>& site_lh_diff, std::vector<RealNumType>& site_lh_root_diff, const std::vector<RealNumType>& site_lh_root, std::unique_ptr<SeqRegions>& parent_new_lower_lh, const RealNumType& child_2_new_blength, PhyloNode& current_node, PhyloNode& child_1, PhyloNode& child_2, PhyloNode& sibling, PhyloNode& parent, const Index parent_index)
{
    // dummy variables
    const RealNumType threshold_prob = params->threshold_prob;
    const RealNumType child_1_blength = child_1.getUpperLength(); // ~new_branch_length
    const std::unique_ptr<SeqRegions>& child_1_lower_regions = child_1.getPartialLh(TOP); // ~subtree_regions
    const std::unique_ptr<SeqRegions>& sibling_lower_lh = sibling.getPartialLh(TOP);
    const std::unique_ptr<SeqRegions>& child_2_lower_lh = child_2.getPartialLh(TOP);
    std::unique_ptr<SeqRegions> new_parent_new_lower_lh = nullptr;
    std::unique_ptr<SeqRegions> tmp_lower_lh = nullptr;
    const PositionType seq_length = aln->ref_seq.size();
    std::vector<RealNumType> site_lh_diff_old(seq_length, 0);
    
    const std::unique_ptr<SeqRegions>& grand_parent_upper_lr = getPartialLhAtNode(parent.getNeighborIndex(TOP));
    std::unique_ptr<SeqRegions> best_parent_regions = nullptr;
    RealNumType parent_blength = parent.getUpperLength();
    // update parent_blength if it's <= 0
    if (parent_blength <= 0)
    {
        // dummy variables
        RealNumType best_lh = MIN_NEGATIVE;
        bool blength_changed = false;
        optimizeBlengthBeforeSeekingSPR<num_states>(parent, parent_blength, best_lh, blength_changed, grand_parent_upper_lr, parent_new_lower_lh);
    }
    // because mid_branch_lh is outdated! => compute a new one
    const RealNumType parent_mid_blength = 0.5 * parent_blength;
    std::unique_ptr<SeqRegions> parent_new_mid_branch_lh = nullptr;
    grand_parent_upper_lr->mergeUpperLower<num_states>(parent_new_mid_branch_lh, parent_mid_blength, *parent_new_lower_lh, parent_mid_blength, aln, model, threshold_prob);
    RealNumType best_parent_lh = calculateSubTreePlacementCost<num_states>(parent_new_mid_branch_lh, child_1_lower_regions, child_1_blength);
    RealNumType best_parent_blength_split = parent_mid_blength;
    
    // 2. estimate x (~the position) to regraft child_1 in the branch connecting the parent and grand parent nodes
    // try with a shorter split
    bool found_new_split = tryShorterBranch<num_states, &cmaple::TreeBase::calculateSubTreePlacementCost<num_states>>(parent_blength, best_parent_regions, child_1_lower_regions, grand_parent_upper_lr, parent_new_lower_lh, best_parent_lh, best_parent_blength_split, child_1_blength, false);
    // try with a longer split
    if (!found_new_split)
    {
        // try on the second half of the branch
        found_new_split = tryShorterBranch<num_states, &cmaple::TreeBase::calculateSubTreePlacementCost<num_states>>(parent_blength, best_parent_regions, child_1_lower_regions, grand_parent_upper_lr, parent_new_lower_lh, best_parent_lh, best_parent_blength_split, child_1_blength, true);
        
        if (found_new_split)
            best_parent_blength_split = parent_blength - best_parent_blength_split;
    }
    
    // Delay cloning SeqRegions
    if (!best_parent_regions)
        best_parent_regions = std::make_unique<SeqRegions>(SeqRegions(parent_new_mid_branch_lh));
    
    // 3. estimate new l5 ~ the length for the new branch re-connecting child_1 to the tree
    RealNumType child_1_new_blength = child_1_blength;
    estimateLengthNewBranch<&cmaple::TreeBase::calculateSubTreePlacementCost<num_states>>(best_parent_lh, best_parent_regions, child_1_lower_regions, child_1_new_blength, child_1_blength * 10, double_min_blength, (child_1_blength <= 0));
    
    // finalize x and (l_0 -x)
    RealNumType parent_new_blength = best_parent_blength_split;
    RealNumType new_parent_new_blength = parent_blength - parent_new_blength;
    if (best_parent_blength_split <= 0)
    {
        parent_new_blength = -1;
        new_parent_new_blength = parent_blength;
    }
    
    // 4. recompute the lower_lh of the new_parent
    parent_new_lower_lh->mergeTwoLowers<num_states>(new_parent_new_lower_lh, parent_new_blength, *child_1_lower_regions, child_1_new_blength, aln, model, params->threshold_prob, false);
    
    // 5. compute the new upper_left/right lh for the parent and the new parent nodes
    // 5.1. for the new parent node
    std::unique_ptr<SeqRegions> new_parent_new_upper_lr_1 = nullptr;
    grand_parent_upper_lr->mergeUpperLower<num_states>(new_parent_new_upper_lr_1, new_parent_new_blength, *child_1_lower_regions, child_1_new_blength, aln, model, threshold_prob);
    std::unique_ptr<SeqRegions> new_parent_new_upper_lr_2 = nullptr;
    grand_parent_upper_lr->mergeUpperLower<num_states>(new_parent_new_upper_lr_2, new_parent_new_blength, *parent_new_lower_lh, parent_new_blength, aln, model, threshold_prob);
    // 5.2. for the parent node
    std::unique_ptr<SeqRegions> parent_new_upper_lr_1 = nullptr;
    new_parent_new_upper_lr_1->mergeUpperLower<num_states>(parent_new_upper_lr_1, parent_new_blength, *sibling_lower_lh, sibling.getUpperLength(), aln, model, threshold_prob);
    std::unique_ptr<SeqRegions> parent_new_upper_lr_2 = nullptr;
    new_parent_new_upper_lr_1->mergeUpperLower<num_states>(parent_new_upper_lr_2, parent_new_blength, *child_2_lower_lh, child_2_new_blength, aln, model, threshold_prob);
    
    // 6. re-optimize the lengths of the upper branches of
    // 6.1. child_2
    const RealNumType child_2_best_blength = estimateBranchLengthWithCheck<num_states>(parent_new_upper_lr_1, child_2_lower_lh, child_2_new_blength);
    // 6.2. sibling
    const RealNumType sibling_best_blength = estimateBranchLengthWithCheck<num_states>(parent_new_upper_lr_2, sibling_lower_lh, sibling.getUpperLength());
    // 6.3. parent
    const RealNumType parent_best_blength = estimateBranchLengthWithCheck<num_states>(new_parent_new_upper_lr_1, parent_new_lower_lh, parent_new_blength);
    // 6.4. child_1
    const RealNumType child_1_best_blength = estimateBranchLengthWithCheck<num_states>(new_parent_new_upper_lr_2, child_1_lower_regions, child_1_new_blength);
    // 6.5. new_parent
    const RealNumType new_parent_best_blength = estimateBranchLengthWithCheck<num_states>(grand_parent_upper_lr, new_parent_new_lower_lh, new_parent_new_blength);
    
    // 7. caculate the likelihood contribution changed at
    // 7.1. the parent node
    child_2_lower_lh->calculateSiteLhContributions<num_states>(site_lh_diff, parent_new_lower_lh, child_2_best_blength, *sibling_lower_lh, sibling_best_blength, aln, model, threshold_prob);
    child_2_lower_lh->calculateSiteLhContributions<num_states>(site_lh_diff_old, tmp_lower_lh, child_2.getUpperLength(), *child_1_lower_regions, child_1.getUpperLength(), aln, model, threshold_prob);
    // 7.2. the new_parent node
    RealNumType prev_lh_diff = parent_new_lower_lh->calculateSiteLhContributions<num_states>(site_lh_diff, new_parent_new_lower_lh, parent_best_blength, *child_1_lower_regions, child_1_best_blength, aln, model, threshold_prob) - node_lhs[parent.getNodelhIndex()].getLhContribution() ;
    // 7.3. other ancestors on the path from the new_parent to root (stop when the change is insignificant)
    Index node_index = parent_index;
    RealNumType bk_tmp_blength = parent_best_blength;
    RealNumType bk_tmp_sibling_blength = child_1_best_blength;
    const Index current_node_index = child_1.getNeighborIndex(TOP);
    NumSeqsType bk_sibling_vec = current_node.getNeighborIndex(current_node_index.getMiniIndex()).getVectorIndex();
    std::unique_ptr<SeqRegions> bk_new_lower_lh = std::move(parent_new_lower_lh);
    std::unique_ptr<SeqRegions> new_lower_lh = std::move(new_parent_new_lower_lh);
    std::unique_ptr<SeqRegions> tmp_new_lower_lh = nullptr;
    RealNumType tmp_blength = new_parent_best_blength;
    while (true)
    {
        NumSeqsType node_vec = node_index.getVectorIndex();
        PhyloNode& node = nodes[node_vec];
        
        // if the new lower lh is different from the old one -> traverse upward further
        if (node.getPartialLh(TOP)->areDiffFrom(*new_lower_lh, seq_length, num_states, *params) || fabs(prev_lh_diff) > threshold_prob)
        {
            // update lh_diff
            PhyloNode& tmp_child_1 = nodes[node.getNeighborIndex(RIGHT).getVectorIndex()];
            PhyloNode& tmp_child_2 = nodes[node.getNeighborIndex(LEFT).getVectorIndex()];
            tmp_child_1.getPartialLh(TOP)->calculateSiteLhContributions<num_states>(site_lh_diff_old, tmp_lower_lh, tmp_child_1.getUpperLength(), *tmp_child_2.getPartialLh(TOP), tmp_child_2.getUpperLength(), aln, model, threshold_prob);
            
            // cases when node is non-root
            if (root_vector_index != node_vec)
            {
                // re-caculate the lower likelihood (likelihood contribution) at the grand-parent node
                const Index tmp_parent_index = node.getNeighborIndex(TOP);
                const NumSeqsType tmp_parent_vec = tmp_parent_index.getVectorIndex();
                PhyloNode& tmp_parent = nodes[tmp_parent_vec];
                const NumSeqsType tmp_sibling_vec = tmp_parent.getNeighborIndex(tmp_parent_index.getFlipMiniIndex()).getVectorIndex();
                PhyloNode& tmp_sibling = nodes[tmp_sibling_vec];
                
                prev_lh_diff = new_lower_lh->calculateSiteLhContributions<num_states>(site_lh_diff, tmp_new_lower_lh, tmp_blength, *(tmp_sibling.getPartialLh(TOP)), tmp_sibling.getUpperLength(), aln, model, params->threshold_prob) - node_lhs[tmp_parent.getNodelhIndex()].getLhContribution();
                
                // move a step upwards
                node_index = tmp_parent_index;
                bk_sibling_vec = tmp_sibling_vec;
                bk_tmp_sibling_blength = tmp_sibling.getUpperLength();
                bk_tmp_blength = tmp_blength;
                tmp_blength = tmp_parent.getUpperLength();
                bk_new_lower_lh = std::move(new_lower_lh);
                new_lower_lh = std::move(tmp_new_lower_lh);
            }
            // case when node is root
            else
            {
                // re-calculate likelihood at root
                new_lower_lh->computeSiteLhAtRoot<num_states>(site_lh_root_diff, model);
                
                // update site_lh_root_diff
                for (PositionType j = 0; j < seq_length; ++j)
                    site_lh_root_diff[j] -= site_lh_root[j];
                
                // stop traversing further
                break;
            }
        }
        // otherwise, stop traversing further
        else
        {
            // cancel the contribution of the last merging in site_lh_diff by adding it into site_lh_diff_old
            PhyloNode& tmp_sibling = nodes[bk_sibling_vec];
            bk_new_lower_lh->calculateSiteLhContributions<num_states>(site_lh_diff_old, tmp_new_lower_lh, bk_tmp_blength, *(tmp_sibling.getPartialLh(TOP)), bk_tmp_sibling_blength, aln, model, params->threshold_prob);
            break;
        }
    }
    
    // calculate the site-lh differences between the new and the old ones
    for (PositionType j = 0; j < seq_length; ++j)
        site_lh_diff[j] -= site_lh_diff_old[j];
}

template <const StateType num_states>
void cmaple::TreeBase::calSiteLhDiff(std::vector<RealNumType>& site_lh_diff, std::vector<RealNumType>& site_lh_root_diff, const std::vector<RealNumType>& site_lh_root, PhyloNode& current_node, PhyloNode& child_1, PhyloNode& child_2, PhyloNode& sibling, PhyloNode& parent, const Index parent_index)
{
    // 1. recompute the lowerlh at the parent node after swaping child_1 and sibling
    std::unique_ptr<SeqRegions> parent_new_lower_lh = nullptr;
    RealNumType child_2_new_blength = child_2.getUpperLength();
    
    if (child_2_new_blength > 0)
    {
        if (current_node.getUpperLength() > 0)
            child_2_new_blength += current_node.getUpperLength();
    }
    else
        child_2_new_blength = current_node.getUpperLength();
    child_2.getPartialLh(TOP)->mergeTwoLowers<num_states>(parent_new_lower_lh, child_2_new_blength, *(sibling.getPartialLh(TOP)), sibling.getUpperLength(), aln, model, params->threshold_prob, false);
    
    // if the (old) parent is root
    // for more information, pls see https://tinyurl.com/5n8m5c8y
    if (root_vector_index == parent_index.getVectorIndex())
        calSiteLhDiffRoot<num_states>(site_lh_diff, site_lh_root_diff, site_lh_root, parent_new_lower_lh, child_2_new_blength, current_node, child_1, child_2, sibling, parent, parent_index);
    // otherwise, the (old) parent node is non-root
    // for more information, pls see https://tinyurl.com/ymr49jy8
    else
        calSiteLhDiffNonRoot<num_states>(site_lh_diff, site_lh_root_diff, site_lh_root, parent_new_lower_lh, child_2_new_blength, current_node, child_1, child_2, sibling, parent, parent_index);
}

void findTwoLargest(const RealNumType a, const RealNumType b, const RealNumType c, RealNumType& largest, RealNumType& second_largest) {
    // only consider a and b
    largest = a;
    second_largest = b;
    
    if (b > largest) {
        largest = b;
        second_largest = a;
    }
    
    // take into account c
    if (c > largest) {
        second_largest = largest;
        largest = c;
    } else if (c > second_largest) {
        second_largest = c;
    }
}

template <const StateType num_states>
PositionType cmaple::TreeBase::count_aRLT_SH_branch(std::vector<RealNumType>& site_lh_contributions, std::vector<RealNumType>& site_lh_root, PhyloNode& node, const RealNumType& LT1)
{
    // caculate aLRT
    const Index child_1_index = node.getNeighborIndex(RIGHT);
    const Index child_2_index = node.getNeighborIndex(LEFT);
    PhyloNode& child_1 = nodes[child_1_index.getVectorIndex()];
    PhyloNode& child_2 = nodes[child_2_index.getVectorIndex()];
    const Index parent_index = node.getNeighborIndex(TOP);
    PhyloNode& parent = nodes[parent_index.getVectorIndex()];
    const Index sibling_index = parent.getNeighborIndex(parent_index.getFlipMiniIndex());
    PhyloNode& sibling = nodes[sibling_index.getVectorIndex()];
    PositionType sh_count{0};
    const PositionType seq_length = aln->ref_seq.size();
    const NodeLh& nodelh = node_lhs[node.getNodelhIndex()];
    const RealNumType LT2 = LT1 + nodelh.getLhDiff2();
    const RealNumType LT3 = LT1 + nodelh.getLhDiff3();
    
    // debug
    /*if (!child_1.isInternal() && (aln->data[child_1.getSeqNameIndex()].seq_name == "52"))
        std::cout << "sfdfds " <<std::endl;
    
    if (!child_2.isInternal() && (aln->data[child_2.getSeqNameIndex()].seq_name == "52"))
        std::cout << "sfdfds " <<std::endl;*/
    
    // calculate site_lh differences
    // neighbor 2
    std::vector<RealNumType> site_lh_diff_2(seq_length, 0);
    std::vector<RealNumType> site_lh_root_diff_2(seq_length, 0);
    calSiteLhDiff<num_states>(site_lh_diff_2, site_lh_root_diff_2, site_lh_root, node, child_1, child_2, sibling, parent, parent_index);
    // neighbor 3
    std::vector<RealNumType> site_lh_diff_3(seq_length, 0);
    std::vector<RealNumType> site_lh_root_diff_3(seq_length, 0);
    calSiteLhDiff<num_states>(site_lh_diff_3, site_lh_root_diff_3, site_lh_root, node, child_2, child_1, sibling, parent, parent_index);
    
    // validate the results
    RealNumType lh_diff_2{0}, lh_diff_3{0};
    for (PositionType j = 0; j < seq_length; ++j)
    {
        lh_diff_2 += site_lh_diff_2[j];
        lh_diff_2 += site_lh_root_diff_2[j];
        
        lh_diff_3 += site_lh_diff_3[j];
        lh_diff_3 += site_lh_root_diff_3[j];
    }
    ASSERT(isinf(lh_diff_2) || fabs(lh_diff_2 - nodelh.getLhDiff2()) < 1e-3);
    ASSERT(isinf(lh_diff_3) || fabs(lh_diff_3 - nodelh.getLhDiff3()) < 1e-3);
    
    // iterate a number of replicates
    #pragma omp parallel reduction(+:sh_count)
    {
      std::vector<PositionType> selected_sites(seq_length);
      int thread_id = 0;
      #ifdef _OPENMP
      thread_id = omp_get_thread_num();
      #endif

      // init random generators
      std::default_random_engine gen(params->ran_seed + thread_id);
      std::uniform_int_distribution<> rng_distrib(0, seq_length - 1);
      #pragma omp for
      for (PositionType i = 0; i < params->aLRT_SH_replicates; ++i)
      {
        // reset all LTX*
        RealNumType LT1_star{0}, LT2_star{0}, LT3_star{0};
        
        // generate a vector of sampled sites
        for (PositionType j = 0; j < seq_length; ++j)
            selected_sites[j] = rng_distrib(gen);
        
        // compute LT1*, LT2*, LT3*
        // NOTES: LT2*, LT3* are now the lh differences between the actual LT2*, LT3* and LT1*
        for (PositionType j = 0; j < seq_length; ++j)
        {
            const PositionType site_index = selected_sites[j];
            
            LT1_star += site_lh_contributions[site_index];
            LT1_star += site_lh_root[site_index];
            
            LT2_star += site_lh_diff_2[site_index];
            LT2_star += site_lh_root_diff_2[site_index];
            
            LT3_star += site_lh_diff_3[site_index];
            LT3_star += site_lh_root_diff_3[site_index];
        }
        
        // compute the actual LT2* and LT3*
        LT2_star += LT1_star;
        LT3_star += LT1_star;
        
        // compute the centered sums CS1*, CS2*, CS3*
        // where CSX* = LTX* - LTX
        const RealNumType CS1 = LT1_star - LT1;
        const RealNumType CS2 = LT2_star - LT2;
        const RealNumType CS3 = LT3_star - LT3;
        
        // find CS_first and CS_second which are the highest and the second highest among CSX* values
        RealNumType CS_first, CS_second;
        findTwoLargest(CS1, CS2, CS3, CS_first, CS_second);
        
        // increase sh_count if the condition (aLRT > 2(CS_first - CS_second) + epsilon) is satisfied
        // <=> half_aLRT > CS_first - CS_second + half_epsilon
        if (nodelh.getHalf_aLRT() > (CS_first - CS_second + params->aLRT_SH_half_epsilon))
            ++sh_count;
      } // for aLRT_SH_replicates
    } // omp parallel
    
    return sh_count;
}

template <const StateType num_states>
void cmaple::TreeBase::calculate_aRLT_SH(std::vector<RealNumType>& site_lh_contributions, std::vector<RealNumType>& site_lh_root, const RealNumType& LT1)
{
    const RealNumType replicate_inverse = 100.0 / params->aLRT_SH_replicates;
    
    // traverse tree to calculate aLRT-SH for each internal branch
    PhyloNode& root = nodes[root_vector_index];
    
    // aLRT-SH at root branch is zero
    node_lhs[root.getNodelhIndex()].set_aLRT_SH(0);
    
    std::stack<Index> node_stack;
    if (root.isInternal())
    {
        node_stack.push(root.getNeighborIndex(RIGHT));
        node_stack.push(root.getNeighborIndex(LEFT));
    }
    
    while (!node_stack.empty())
    {
        const Index node_index = node_stack.top();
        node_stack.pop();
        const NumSeqsType node_vec = node_index.getVectorIndex();
        PhyloNode& node = nodes[node_vec];
        
        // only consider internal branches
        if (node.isInternal())
        {
            // add its children into node_stack for further traversal
            node_stack.push(node.getNeighborIndex(RIGHT));
            node_stack.push(node.getNeighborIndex(LEFT));
            
            // only compute the aLRT for internal non-zero branches
            if (node.getUpperLength() > 0)
            {
                node_lhs[node.getNodelhIndex()].set_aLRT_SH(replicate_inverse * count_aRLT_SH_branch<num_states>(site_lh_contributions, site_lh_root, node, LT1));
                
                // print out the aLRT (for debugging only)
                // std::cout << std::setprecision(3) << "aLRT-SH (node_vec: " << node_vec << "): " << node_lhs[node.getNodelhIndex()].get_aLRT_SH() << std::endl;
            }
            // return zero for zero-length internal branches
            else
            {
                node_lhs[node.getNodelhIndex()].set_aLRT_SH(0);
                // print out the aLRT (for debugging only)
                // std::cout << std::setprecision(3) << "aLRT-SH (node_vec: " << node_vec << "): 0 (*zero-length branch)" << std::endl;
            }
        }
    }
}

template <const StateType num_states>
bool cmaple::TreeBase::calculateNNILh(std::stack<Index>& node_stack_aLRT, RealNumType& lh_diff, PhyloNode& current_node, PhyloNode& child_1, PhyloNode& child_2, PhyloNode& sibling, PhyloNode& parent, const Index parent_index, RealNumType& lh_at_root, const bool allow_replacing_ML_tree)
{
    // 1. recompute the lowerlh at the parent node after swaping child_1 and sibling
    std::unique_ptr<SeqRegions> parent_new_lower_lh = nullptr;
    RealNumType child_2_new_blength = child_2.getUpperLength();
    
    if (child_2_new_blength > 0)
    {
        if (current_node.getUpperLength() > 0)
            child_2_new_blength += current_node.getUpperLength();
    }
    else
        child_2_new_blength = current_node.getUpperLength();
    child_2.getPartialLh(TOP)->mergeTwoLowers<num_states>(parent_new_lower_lh, child_2_new_blength, *(sibling.getPartialLh(TOP)), sibling.getUpperLength(), aln, model, params->threshold_prob, false);
    
    // if the (old) parent is root
    // for more information, pls see https://tinyurl.com/5n8m5c8y
    if (root_vector_index == parent_index.getVectorIndex())
    {
        return calculateNNILhRoot<num_states>(node_stack_aLRT, lh_diff, parent_new_lower_lh, child_2_new_blength, current_node, child_1, child_2, sibling, parent, parent_index, lh_at_root, allow_replacing_ML_tree);
    }
    // otherwise, the (old) parent node is non-root
    // for more information, pls see https://tinyurl.com/ymr49jy8
    else
    {
        return calculateNNILhNonRoot<num_states>(node_stack_aLRT, lh_diff, parent_new_lower_lh, child_2_new_blength, current_node, child_1, child_2, sibling, parent, parent_index, lh_at_root, allow_replacing_ML_tree);
    }
    
    return true;
}

template <const StateType num_states>
bool cmaple::TreeBase::calculateNNILhRoot(std::stack<Index>& node_stack_aLRT, RealNumType& lh_diff, std::unique_ptr<SeqRegions>& parent_new_lower_lh, const RealNumType& child_2_new_blength, PhyloNode& current_node, PhyloNode& child_1, PhyloNode& child_2, PhyloNode& sibling, PhyloNode& parent, const Index parent_index, RealNumType& lh_at_root, const bool allow_replacing_ML_tree)
{
    const RealNumType threshold_prob = params->threshold_prob;
    const RealNumType child_1_blength = child_1.getUpperLength(); // ~new_branch_length
    const std::unique_ptr<SeqRegions>& child_1_lower_regions = child_1.getPartialLh(TOP); // ~subtree_regions
    const std::unique_ptr<SeqRegions>& sibling_lower_lh = sibling.getPartialLh(TOP);
    const std::unique_ptr<SeqRegions>& child_2_lower_lh = child_2.getPartialLh(TOP);
    std::unique_ptr<SeqRegions> new_parent_new_lower_lh = nullptr;
    
    // 2. estimate x ~ the length of the new branch connecting the parent and the new_parent nodes
    RealNumType parent_new_blength = default_blength;
    // merge 2 lower vector into one
    RealNumType best_parent_lh = parent_new_lower_lh->mergeTwoLowers<num_states>(new_parent_new_lower_lh, parent_new_blength, *child_1_lower_regions, child_1_blength, aln, model, threshold_prob, true);
    best_parent_lh += new_parent_new_lower_lh->computeAbsoluteLhAtRoot<num_states>(model);
    // Try shorter branch lengths at root
    tryShorterBranchAtRoot<num_states>(child_1_lower_regions, parent_new_lower_lh, new_parent_new_lower_lh, parent_new_blength, best_parent_lh, child_1_blength);
    
    // 3. estimate the length of the new branch connecting child_1 to the new_parent node
    RealNumType child_1_new_blength = child_1_blength;
    estimateLengthNewBranchAtRoot<num_states>(child_1_lower_regions, parent_new_lower_lh, new_parent_new_lower_lh, child_1_new_blength, best_parent_lh, parent_new_blength, double_min_blength, child_1_blength <= 0);
    
    // 4. compute the new lower_lh of the new_parent node
    parent_new_lower_lh->mergeTwoLowers<num_states>(new_parent_new_lower_lh, parent_new_blength, *child_1_lower_regions, child_1_new_blength, aln, model, threshold_prob, false);
    
    // 5. compute the new upper_left/right lh for the parent and the new parent nodes
    // 5.1. for the new parent node
    std::unique_ptr<SeqRegions> new_parent_new_upper_lr_1 = nullptr;
    child_1_lower_regions->computeTotalLhAtRoot<num_states>(new_parent_new_upper_lr_1, model, child_1_new_blength);
    std::unique_ptr<SeqRegions> new_parent_new_upper_lr_2 = nullptr;
    parent_new_lower_lh->computeTotalLhAtRoot<num_states>(new_parent_new_upper_lr_2, model, parent_new_blength);
    // 5.2. for the parent node
    std::unique_ptr<SeqRegions> parent_new_upper_lr_1 = nullptr;
    new_parent_new_upper_lr_1->mergeUpperLower<num_states>(parent_new_upper_lr_1, parent_new_blength, *sibling_lower_lh, sibling.getUpperLength(), aln, model, threshold_prob);
    std::unique_ptr<SeqRegions> parent_new_upper_lr_2 = nullptr;
    new_parent_new_upper_lr_1->mergeUpperLower<num_states>(parent_new_upper_lr_2, parent_new_blength, *child_2_lower_lh, child_2_new_blength, aln, model, threshold_prob);
    
    // 6. re-optimize the lengths of the upper branches of
    // 6.1. child_2
    const RealNumType child_2_best_blength = estimateBranchLengthWithCheck<num_states>(parent_new_upper_lr_1, child_2_lower_lh, child_2_new_blength);
    // 6.2. sibling
    const RealNumType sibling_best_blength = estimateBranchLengthWithCheck<num_states>(parent_new_upper_lr_2, sibling_lower_lh, sibling.getUpperLength());
    // 6.3. parent
    const RealNumType parent_best_blength = estimateBranchLengthWithCheck<num_states>(new_parent_new_upper_lr_1, parent_new_lower_lh, parent_new_blength);
    // 6.4. child_1
    const RealNumType child_1_best_blength = estimateBranchLengthWithCheck<num_states>(new_parent_new_upper_lr_2, child_1_lower_regions, child_1_new_blength);
    
    // 7. caculate the likelihood contribution changed at
    // 7.1. the parent node
    lh_diff += child_2_lower_lh->mergeTwoLowers<num_states>(parent_new_lower_lh, child_2_best_blength, *sibling_lower_lh, sibling_best_blength, aln, model, threshold_prob, true) - node_lhs[current_node.getNodelhIndex()].getLhContribution();
    // 7.2. the new_parent node
    lh_diff += parent_new_lower_lh->mergeTwoLowers<num_states>(new_parent_new_lower_lh, parent_best_blength, *child_1_lower_regions, child_1_best_blength, aln, model, threshold_prob, true) - node_lhs[parent.getNodelhIndex()].getLhContribution();
    // 7.3 the absolute likelihood at root
    lh_diff += new_parent_new_lower_lh->computeAbsoluteLhAtRoot<num_states>(model) - lh_at_root;
    
    // if we found an NNI neighbor with higher lh => replace the ML tree
    if (lh_diff > 0)
    {
        // Check if we could replace the ML tree
        if (allow_replacing_ML_tree)
        {
            if (cmaple::verbose_mode >= cmaple::VB_DEBUG)
                std::cout << std::setprecision(10) << "Replace the ML tree by a newly found NNI neighbor tree (root), improving tree loglh by: " << lh_diff << std::endl;
            // std::cout << "Tree lh (before replacing): " << calculateTreeLh() << std::endl;
            replaceMLTreebyNNIRoot<num_states>(node_stack_aLRT, lh_diff, current_node, child_1, child_2, sibling, parent, lh_at_root, child_1_best_blength, child_2_best_blength, sibling_best_blength, parent_best_blength);
            // std::cout << "Tree lh (after replacing): " << calculateTreeLh() << std::endl;
            
            // return false to let us know that we found a new ML tree
            return false;
        }
        else if (cmaple::verbose_mode >= cmaple::VB_MED)
            outWarning("Found an NNI neighbor tree with a higher likelihood (by " + convertDoubleToString(lh_diff) + ") than the current ML tree");
    }
    
    return true;
}

// NHANLT: Debug aLRT
/* void cmaple::TreeBase::log_current(std::stack<Index>& node_stack_aLRT)
{
    // NHANLT: Debug
    {
        // check if node 1249 exists in the stack
        stack<Index> tmp_stack = node_stack_aLRT;
        while (!tmp_stack.empty())
        {
            if (tmp_stack.top().getVectorIndex() == 3185)
            {
                std::cout << "Found 3185 " << (nodes[3185].isOutdated() ? "Outdated" : "Updated") << std::endl;
                break;
            }
            tmp_stack.pop();
        }
        std::stack<cmaple::Index> node_stack;
        RealNumType lh_at_root = 0;
        RealNumType tree_total_lh = 0;
        PhyloNode& tmp_node = nodes[3185];
            // show the likelihood contribution at this branch/node
            // std::cout << "Lh contribution of node (node_vec) " << node_vec << " : " << node_lhs[node.getNodelhIndex()].lh_contribution_ << std::endl;
            
            // caculate aLRT
            PhyloNode& child_1 = nodes[2841];
            PhyloNode& child_2 = nodes[3186];
            const Index parent_index = tmp_node.getNeighborIndex(TOP);
            PhyloNode& parent = nodes[parent_index.getVectorIndex()];
            const Index sibling_index = parent.getNeighborIndex(parent_index.getFlipMiniIndex());
            PhyloNode& sibling = nodes[sibling_index.getVectorIndex()];
            
            // compute the likelihood differences between each nni neighbor and the current tree
            RealNumType neighbor_2_lh_diff = 0, neighbor_3_lh_diff = 0;
            // if calculateNNILh return false => the current tree was replaced by a newly found ML tree
            if (!calculateNNILh<4>(node_stack, neighbor_2_lh_diff, tmp_node, child_1, child_2, sibling, parent, parent_index, lh_at_root))
            {
                tree_total_lh += neighbor_2_lh_diff;
                //continue;
            }
            // if calculateNNILh return false => the current tree was replaced by a newly found ML tree
            if (!calculateNNILh<4>(node_stack, neighbor_3_lh_diff, tmp_node, child_2, child_1, sibling, parent, parent_index, lh_at_root))
            {
                tree_total_lh += neighbor_3_lh_diff;
                //continue;
            }
            
            // record neighbor_2_lh_diff, neighbor_2_lh_diff
            NodeLh& node_lh = node_lhs[tmp_node.getNodelhIndex()];
        if ( fabs(neighbor_2_lh_diff -  node_lh.getLhDiff2()) > 1e-3)
        {
            // node_lh.setLhDiff2(neighbor_2_lh_diff);
            // node_lh.setLhDiff3(neighbor_3_lh_diff);
            if (!nodes[3185].isOutdated())
                std::cout << neighbor_2_lh_diff << " " << node_lh.getLhDiff2() << std::endl;
        }
    }
}*/

template <const StateType num_states>
bool cmaple::TreeBase::calculateNNILhNonRoot(std::stack<Index>& node_stack_aLRT, RealNumType& lh_diff, std::unique_ptr<SeqRegions>& parent_new_lower_lh, const RealNumType& child_2_new_blength, PhyloNode& current_node, PhyloNode& child_1, PhyloNode& child_2, PhyloNode& sibling, PhyloNode& parent, const Index parent_index, RealNumType& lh_at_root, const bool allow_replacing_ML_tree)
{
    // dummy variables
    const RealNumType threshold_prob = params->threshold_prob;
    const RealNumType child_1_blength = child_1.getUpperLength(); // ~new_branch_length
    const std::unique_ptr<SeqRegions>& child_1_lower_regions = child_1.getPartialLh(TOP); // ~subtree_regions
    const std::unique_ptr<SeqRegions>& sibling_lower_lh = sibling.getPartialLh(TOP);
    const std::unique_ptr<SeqRegions>& child_2_lower_lh = child_2.getPartialLh(TOP);
    std::unique_ptr<SeqRegions> new_parent_new_lower_lh = nullptr;
    
    const std::unique_ptr<SeqRegions>& grand_parent_upper_lr = getPartialLhAtNode(parent.getNeighborIndex(TOP));
    std::unique_ptr<SeqRegions> best_parent_regions = nullptr;
    RealNumType parent_blength = parent.getUpperLength();
    // update parent_blength if it's <= 0
    if (parent_blength <= 0)
    {
        // dummy variables
        RealNumType best_lh = MIN_NEGATIVE;
        bool blength_changed = false;
        optimizeBlengthBeforeSeekingSPR<num_states>(parent, parent_blength, best_lh, blength_changed, grand_parent_upper_lr, parent_new_lower_lh);
    }
    // because mid_branch_lh is outdated! => compute a new one
    const RealNumType parent_mid_blength = 0.5 * parent_blength;
    std::unique_ptr<SeqRegions> parent_new_mid_branch_lh = nullptr;
    grand_parent_upper_lr->mergeUpperLower<num_states>(parent_new_mid_branch_lh, parent_mid_blength, *parent_new_lower_lh, parent_mid_blength, aln, model, threshold_prob);
    RealNumType best_parent_lh = calculateSubTreePlacementCost<num_states>(parent_new_mid_branch_lh, child_1_lower_regions, child_1_blength);
    RealNumType best_parent_blength_split = parent_mid_blength;
    
    // 2. estimate x (~the position) to regraft child_1 in the branch connecting the parent and grand parent nodes
    // try with a shorter split
    bool found_new_split = tryShorterBranch<num_states, &cmaple::TreeBase::calculateSubTreePlacementCost<num_states>>(parent_blength, best_parent_regions, child_1_lower_regions, grand_parent_upper_lr, parent_new_lower_lh, best_parent_lh, best_parent_blength_split, child_1_blength, false);
    // try with a longer split
    if (!found_new_split)
    {
        // try on the second half of the branch
        found_new_split = tryShorterBranch<num_states, &cmaple::TreeBase::calculateSubTreePlacementCost<num_states>>(parent_blength, best_parent_regions, child_1_lower_regions, grand_parent_upper_lr, parent_new_lower_lh, best_parent_lh, best_parent_blength_split, child_1_blength, true);
        
        if (found_new_split)
            best_parent_blength_split = parent_blength - best_parent_blength_split;
    }
    
    // Delay cloning SeqRegions
    if (!best_parent_regions)
        best_parent_regions = std::make_unique<SeqRegions>(SeqRegions(parent_new_mid_branch_lh));
    
    // 3. estimate new l5 ~ the length for the new branch re-connecting child_1 to the tree
    RealNumType child_1_new_blength = child_1_blength;
    estimateLengthNewBranch<&cmaple::TreeBase::calculateSubTreePlacementCost<num_states>>(best_parent_lh, best_parent_regions, child_1_lower_regions, child_1_new_blength, child_1_blength * 10, double_min_blength, (child_1_blength <= 0));
    
    // finalize x and (l_0 -x)
    RealNumType parent_new_blength = best_parent_blength_split;
    RealNumType new_parent_new_blength = parent_blength - parent_new_blength;
    if (best_parent_blength_split <= 0)
    {
        parent_new_blength = -1;
        new_parent_new_blength = parent_blength;
    }
    
    // 4. recompute the lower_lh of the new_parent
    parent_new_lower_lh->mergeTwoLowers<num_states>(new_parent_new_lower_lh, parent_new_blength, *child_1_lower_regions, child_1_new_blength, aln, model, params->threshold_prob, false);
    
    // 5. compute the new upper_left/right lh for the parent and the new parent nodes
    // 5.1. for the new parent node
    std::unique_ptr<SeqRegions> new_parent_new_upper_lr_1 = nullptr;
    grand_parent_upper_lr->mergeUpperLower<num_states>(new_parent_new_upper_lr_1, new_parent_new_blength, *child_1_lower_regions, child_1_new_blength, aln, model, threshold_prob);
    std::unique_ptr<SeqRegions> new_parent_new_upper_lr_2 = nullptr;
    grand_parent_upper_lr->mergeUpperLower<num_states>(new_parent_new_upper_lr_2, new_parent_new_blength, *parent_new_lower_lh, parent_new_blength, aln, model, threshold_prob);
    // 5.2. for the parent node
    std::unique_ptr<SeqRegions> parent_new_upper_lr_1 = nullptr;
    new_parent_new_upper_lr_1->mergeUpperLower<num_states>(parent_new_upper_lr_1, parent_new_blength, *sibling_lower_lh, sibling.getUpperLength(), aln, model, threshold_prob);
    std::unique_ptr<SeqRegions> parent_new_upper_lr_2 = nullptr;
    new_parent_new_upper_lr_1->mergeUpperLower<num_states>(parent_new_upper_lr_2, parent_new_blength, *child_2_lower_lh, child_2_new_blength, aln, model, threshold_prob);
    
    // 6. re-optimize the lengths of the upper branches of
    // 6.1. child_2
    const RealNumType child_2_best_blength = estimateBranchLengthWithCheck<num_states>(parent_new_upper_lr_1, child_2_lower_lh, child_2_new_blength);
    // 6.2. sibling
    const RealNumType sibling_best_blength = estimateBranchLengthWithCheck<num_states>(parent_new_upper_lr_2, sibling_lower_lh, sibling.getUpperLength());
    // 6.3. parent
    const RealNumType parent_best_blength = estimateBranchLengthWithCheck<num_states>(new_parent_new_upper_lr_1, parent_new_lower_lh, parent_new_blength);
    // 6.4. child_1
    const RealNumType child_1_best_blength = estimateBranchLengthWithCheck<num_states>(new_parent_new_upper_lr_2, child_1_lower_regions, child_1_new_blength);
    // 6.5. new_parent
    const RealNumType new_parent_best_blength = estimateBranchLengthWithCheck<num_states>(grand_parent_upper_lr, new_parent_new_lower_lh, new_parent_new_blength);
    
    // 7. caculate the likelihood contribution changed at
    // 7.1. the parent node
    // std::cout << "lh_contribution (before): " << node_lhs[current_node.getNodelhIndex()].getLhContribution() << std::endl;
    lh_diff += child_2_lower_lh->mergeTwoLowers<num_states>(parent_new_lower_lh, child_2_best_blength, *sibling_lower_lh, sibling_best_blength, aln, model, threshold_prob, true) - node_lhs[current_node.getNodelhIndex()].getLhContribution();
    // std::cout << "lh_contribution (after): " << child_2_lower_lh->mergeTwoLowers<num_states>(parent_new_lower_lh, child_2_best_blength, *sibling_lower_lh, sibling_best_blength, aln, model, threshold_prob, true) << std::endl;
    // 7.2. the new_parent node
    RealNumType prev_lh_diff = parent_new_lower_lh->mergeTwoLowers<num_states>(new_parent_new_lower_lh, parent_best_blength, *child_1_lower_regions, child_1_best_blength, aln, model, threshold_prob, true) - node_lhs[parent.getNodelhIndex()].getLhContribution();
    // 7.3. other ancestors on the path from the new_parent to root (stop when the change is insignificant)
    const PositionType seq_length = aln->ref_seq.size();
    
    NumSeqsType node_vec = parent_index.getVectorIndex();
    std::unique_ptr<SeqRegions> new_lower_lh = std::move(new_parent_new_lower_lh);
    std::unique_ptr<SeqRegions> tmp_new_lower_lh = nullptr;
    RealNumType tmp_blength = new_parent_best_blength;
    while (true)
    {
        PhyloNode& node = nodes[node_vec];
    
        // if new_lower_lh == null -> bad NNI -> stop
        if (!new_lower_lh)
        {
            lh_diff = MIN_NEGATIVE;
            break;
        }
        
        // if the new lower lh is different from the old one -> traverse upward further
        if (node.getPartialLh(TOP)->areDiffFrom(*new_lower_lh, seq_length, num_states, *params) || fabs(prev_lh_diff) > threshold_prob)
        {
            // update lh_diff
            // std::cout << "lh_contribution (before): " << node_lhs[node.getNodelhIndex()].getLhContribution() << std::endl;
            lh_diff += prev_lh_diff;
            // std::cout << "lh_contribution (after): " << node_lhs[node.getNodelhIndex()].getLhContribution() + prev_lh_diff << std::endl;
            
            // cases when node is non-root
            if (root_vector_index != node_vec)
            {
                // re-caculate the lower likelihood (likelihood contribution) at the grand-parent node
                const Index tmp_parent_index = node.getNeighborIndex(TOP);
                const NumSeqsType tmp_parent_vec = tmp_parent_index.getVectorIndex();
                PhyloNode& tmp_parent = nodes[tmp_parent_vec];
                PhyloNode& tmp_sibling = nodes[tmp_parent.getNeighborIndex(tmp_parent_index.getFlipMiniIndex()).getVectorIndex()];
                
                prev_lh_diff = new_lower_lh->mergeTwoLowers<num_states>(tmp_new_lower_lh, tmp_blength, *(tmp_sibling.getPartialLh(TOP)), tmp_sibling.getUpperLength(), aln, model, params->threshold_prob, true) - node_lhs[tmp_parent.getNodelhIndex()].getLhContribution();
                
                // move a step upwards
                node_vec = tmp_parent_vec;
                tmp_blength = tmp_parent.getUpperLength();
                new_lower_lh = std::move(tmp_new_lower_lh);
            }
            // case when node is root
            else
            {
                // re-calculate likelihood at root
                // std::cout << "lh at root (before): " << lh_at_root << std::endl;
                lh_diff += new_lower_lh->computeAbsoluteLhAtRoot<num_states>(model) - lh_at_root;
                // std::cout << "lh at root (after): " << new_lower_lh->computeAbsoluteLhAtRoot(num_states, model) << std::endl;
                
                // stop traversing further
                break;
            }
        }
        // otherwise, stop traversing further
        else
            break;
    }
    
    // if we found an NNI neighbor with higher lh => replace the ML tree
    if (lh_diff > 0)
    {
        // Check if we can replace the ML tree
        if (allow_replacing_ML_tree)
        {
            // NHANLT: Debug aLRT
            // log_current(node_stack_aLRT);
            if (cmaple::verbose_mode >= cmaple::VB_DEBUG)
                std::cout << "Replace the ML tree by a newly found NNI neighbor tree (non-root), improving tree loglh by: " << lh_diff << std::endl;
            // std::cout << "Tree lh (before replacing): " <<  std::setprecision(20)<< calculateTreeLh() << std::endl;
            replaceMLTreebyNNINonRoot<num_states>(node_stack_aLRT, lh_diff, current_node, child_1, child_2, sibling, parent, lh_at_root, child_1_best_blength, child_2_best_blength, sibling_best_blength, parent_best_blength, new_parent_best_blength);
            // std::cout << "Tree lh (after replacing): " <<  std::setprecision(20)<< calculateTreeLh() << std::endl;
            
            // NHANLT: Debug aLRT
            // log_current(node_stack_aLRT);
            
            // return false to let us know that we found a new ML tree
            return false;
        }
        else if (cmaple::verbose_mode >= cmaple::VB_MED)
            outWarning("Found an NNI neighbor tree with a higher likelihood (by " + convertDoubleToString(lh_diff) + ") than the current ML tree");
    }
    
    return true;
}

template <const StateType num_states>
void cmaple::TreeBase::replaceMLTreebyNNIRoot(std::stack<Index>& node_stack_aLRT, RealNumType& lh_diff, PhyloNode& current_node, PhyloNode& child_1, PhyloNode& child_2, PhyloNode& sibling, PhyloNode& parent, RealNumType& lh_at_root, const RealNumType child_1_best_blength, const RealNumType child_2_best_blength, const RealNumType sibling_best_blength, const RealNumType parent_best_blength)
{
    // get index of related nodes
    const RealNumType threshold_prob = params->threshold_prob;
    const Index current_node_index = child_1.getNeighborIndex(TOP);
    const MiniIndex current_node_mini = current_node_index.getMiniIndex();
    const Index child_1_index = current_node.getNeighborIndex(current_node_mini);
    const Index child_2_index = current_node.getNeighborIndex(current_node_index.getFlipMiniIndex());
    const Index parent_index = sibling.getNeighborIndex(TOP);
    const MiniIndex parent_mini = parent_index.getMiniIndex();
    const Index sibling_index = parent.getNeighborIndex(parent_mini);
    const std::unique_ptr<SeqRegions>& child_1_lower_lh = child_1.getPartialLh(TOP);
    const std::unique_ptr<SeqRegions>& child_2_lower_lh = child_2.getPartialLh(TOP);
    const std::unique_ptr<SeqRegions>& sibling_lower_lh = sibling.getPartialLh(TOP);
    
    // move child_1 to parent
    child_1.setNeighborIndex(TOP, parent_index);
    parent.setNeighborIndex(parent_mini, child_1_index);
    
    // move sibling to current_node
    sibling.setNeighborIndex(TOP, current_node_index);
    current_node.setNeighborIndex(current_node_mini, sibling_index);
    
    // update 5 blengths
    // recompute aLRT of child_2 if we change its blength from zero to non-zero => don't need to do so because child_2 should be added before, only need to set it outdated
    child_2.setOutdated(true);
    updateBlengthReplaceMLTree<num_states>(node_stack_aLRT, lh_diff, child_2, child_2_index, child_2_best_blength);
    updateBlengthReplaceMLTree<num_states>(node_stack_aLRT, lh_diff, sibling, sibling_index, sibling_best_blength);
    // we don't need to add current_node to node_stack_aLRT since they will be added later
    current_node.setUpperLength(parent_best_blength);
    // we don't need to add child_1 to node_stack_aLRT because child_1 should be added before, only need to set it outdated
    child_1.setOutdated(true);
    updateBlengthReplaceMLTree<num_states>(node_stack_aLRT, lh_diff, child_1, child_1_index, child_1_best_blength);
    
    // stack to update upper_lr_lh later
    std::stack<Index> node_stack_update_upper_lr;
    node_stack_update_upper_lr.push(child_2_index);
    node_stack_update_upper_lr.push(sibling_index);
    
    // update lower_lh of the current node
    std::unique_ptr<SeqRegions>& current_node_lower_lh = current_node.getPartialLh(TOP);
    node_lhs[current_node.getNodelhIndex()].setLhContribution(child_2_lower_lh->mergeTwoLowers<num_states>(current_node_lower_lh, child_2_best_blength, *sibling_lower_lh, sibling_best_blength, aln, model, threshold_prob, true));
    // because the lower_lh of the current node was updated
    // => we need to re-compute aLRT of that node
    current_node.setOutdated(true);
    node_stack_aLRT.push(Index(current_node_index.getVectorIndex(), TOP));
    // => we need to re-compute aLRT of its sibling
    sibling.setOutdated(true);
    node_stack_aLRT.push(Index(sibling_index.getVectorIndex(), TOP));
    // => we need to re-compute aLRT of its parent node
    parent.setOutdated(true);
    node_stack_aLRT.push(Index(parent_index.getVectorIndex(), TOP));
    // => we need to re-compute aLRT of its sibling node => but in this case, no need to add its sibling because child_1 has been added before
    /*child_1.setOutdated(true);
    node_stack_aLRT.push(child_1_index);*/
    // => we need to update the upper_lr of its sibling
    node_stack_update_upper_lr.push(child_1_index);
    
    // update the upper_lr of the parent node
    const MiniIndex parent_current_node_mini = current_node.getNeighborIndex(TOP).getMiniIndex();
    std::unique_ptr<SeqRegions>& parent_current_node_upper_lr = parent.getPartialLh(parent_current_node_mini);
    child_1_lower_lh->computeTotalLhAtRoot<num_states>(parent_current_node_upper_lr, model, child_1_best_blength);
    current_node_lower_lh->computeTotalLhAtRoot<num_states>(parent.getPartialLh(child_1.getNeighborIndex(TOP).getMiniIndex()), model, current_node.getUpperLength());
    // the upper_lr of the parent node is changed -> the aLRT of all its grand_children will be computed later '
    // the parent node has two children: child_1 and the current node
    // recompute aLRT of the children of child_1
    recompute_aLRT_GrandChildren(child_1, node_stack_aLRT);
    // recompute aLRT of the children of the current node -> will be done later (don't need to explicitly add here)
    // update the upper_lr of the current node
    const MiniIndex current_node_child_2_mini = child_2.getNeighborIndex(TOP).getMiniIndex();
    parent_current_node_upper_lr->mergeUpperLower<num_states>(current_node.getPartialLh(current_node_child_2_mini), current_node.getUpperLength(), *sibling_lower_lh, sibling_best_blength, aln, model, threshold_prob);
    const MiniIndex current_node_sibling_mini = sibling.getNeighborIndex(TOP).getMiniIndex();
    parent_current_node_upper_lr->mergeUpperLower<num_states>(current_node.getPartialLh(current_node_sibling_mini), current_node.getUpperLength(), *child_2_lower_lh, child_2_best_blength, aln, model, threshold_prob);
    // because upper_lr of the current node is changed -> we need to recompute aLRT of its grand_children
    // current node has two children: child_2 and sibling (after applying NNI)
    // recompute aLRT of the children of child_2
    recompute_aLRT_GrandChildren(child_2, node_stack_aLRT);
    // recompute aLRT of the children of sibling
    recompute_aLRT_GrandChildren(sibling, node_stack_aLRT);
    
    // update the lower_lh of the parent node
    std::unique_ptr<SeqRegions>& parent_new_lower = parent.getPartialLh(TOP);
    node_lhs[parent.getNodelhIndex()].setLhContribution(current_node.getPartialLh(TOP)->mergeTwoLowers<num_states>(parent_new_lower, current_node.getUpperLength(), *child_1_lower_lh, child_1_best_blength, aln, model, threshold_prob, true));
    // update the absolute likelihood at root
    lh_at_root = parent_new_lower->computeAbsoluteLhAtRoot<num_states>(model);
    
    // traverse downward to update the upper_left/right_region until the changes is insignificant
    updateUpperLR<num_states>(node_stack_update_upper_lr, node_stack_aLRT);
}

template <const StateType num_states>
void cmaple::TreeBase::replaceMLTreebyNNINonRoot(std::stack<Index>& node_stack_aLRT, RealNumType& lh_diff, PhyloNode& current_node, PhyloNode& child_1, PhyloNode& child_2, PhyloNode& sibling, PhyloNode& parent, RealNumType& lh_at_root, const RealNumType child_1_best_blength, const RealNumType child_2_best_blength, const RealNumType sibling_best_blength, const RealNumType parent_best_blength, const RealNumType new_parent_best_blength)
{
    // get index of related nodes
    const RealNumType threshold_prob = params->threshold_prob;
    const Index current_node_index = child_1.getNeighborIndex(TOP);
    const MiniIndex current_node_mini = current_node_index.getMiniIndex();
    const Index child_1_index = current_node.getNeighborIndex(current_node_mini);
    const Index child_2_index = current_node.getNeighborIndex(current_node_index.getFlipMiniIndex());
    const Index parent_index = sibling.getNeighborIndex(TOP);
    const MiniIndex parent_mini = parent_index.getMiniIndex();
    const Index sibling_index = parent.getNeighborIndex(parent_mini);
    const std::unique_ptr<SeqRegions>& child_1_lower_lh = child_1.getPartialLh(TOP);
    const std::unique_ptr<SeqRegions>& child_2_lower_lh = child_2.getPartialLh(TOP);
    const std::unique_ptr<SeqRegions>& sibling_lower_lh = sibling.getPartialLh(TOP);
    const PositionType seq_length = aln->ref_seq.size();
    
    // move child_1 to parent
    child_1.setNeighborIndex(TOP, parent_index);
    parent.setNeighborIndex(parent_mini, child_1_index);
    
    // move sibling to current_node
    sibling.setNeighborIndex(TOP, current_node_index);
    current_node.setNeighborIndex(current_node_mini, sibling_index);
    
    // update 5 blengths
    // recompute aLRT of child_2 if we change its blength from zero to non-zero => don't need to do so because child_2 should be added before, only need to set it outdated
    child_2.setOutdated(true);
    updateBlengthReplaceMLTree<num_states>(node_stack_aLRT, lh_diff, child_2, child_2_index, child_2_best_blength);
    updateBlengthReplaceMLTree<num_states>(node_stack_aLRT, lh_diff, sibling, sibling_index, sibling_best_blength);
    // we don't need to add current_node and parent to node_stack_aLRT since they will be added later
    current_node.setUpperLength(parent_best_blength);
    parent.setUpperLength(new_parent_best_blength);
    // we don't need to add child_1 to node_stack_aLRT because child_1 should be added before, only need to set it outdated
    child_1.setOutdated(true);
    updateBlengthReplaceMLTree<num_states>(node_stack_aLRT, lh_diff, child_1, child_1_index, child_1_best_blength);
    
    // stack to update upper_lr_lh later
    std::stack<Index> node_stack_update_upper_lr;
    node_stack_update_upper_lr.push(child_2_index);
    node_stack_update_upper_lr.push(sibling_index);
    
    // update lower_lh of the current node
    std::unique_ptr<SeqRegions>& current_node_lower_lh = current_node.getPartialLh(TOP);
    // std::cout << "lh_contribution (before): " << node_lhs[current_node.getNodelhIndex()].getLhContribution() << std::endl;
    node_lhs[current_node.getNodelhIndex()].setLhContribution(child_2_lower_lh->mergeTwoLowers<num_states>(current_node_lower_lh, child_2_best_blength, *sibling_lower_lh, sibling_best_blength, aln, model, threshold_prob, true));
    // std::cout << "lh_contribution (after): " << node_lhs[current_node.getNodelhIndex()].getLhContribution() << std::endl;
    // because the lower_lh of the current node was updated
    // => we need to re-compute aLRT of that node
    current_node.setOutdated(true);
    node_stack_aLRT.push(Index(current_node_index.getVectorIndex(), TOP));
    // => we need to re-compute aLRT of the sibling
    sibling.setOutdated(true);
    node_stack_aLRT.push(Index(sibling_index.getVectorIndex(), TOP));
    // => we need to re-compute aLRT of its parent node
    parent.setOutdated(true);
    node_stack_aLRT.push(Index(parent_index.getVectorIndex(), TOP));
    // => we need to re-compute aLRT of its sibling node => but in this case, no need to add its sibling because child_1 has been added before
    /*child_1.setOutdated(true);
    node_stack_aLRT.push(child_1_index);*/
    // => we need to update the upper_lr of its sibling
    node_stack_update_upper_lr.push(child_1_index);
    
    // update the upper_lr of the parent node
    const std::unique_ptr<SeqRegions>& grand_parent_upper_lr = getPartialLhAtNode(parent.getNeighborIndex(TOP));
    const MiniIndex parent_current_node_mini = current_node.getNeighborIndex(TOP).getMiniIndex();
    std::unique_ptr<SeqRegions>& parent_current_node_upper_lr = parent.getPartialLh(parent_current_node_mini);
    grand_parent_upper_lr->mergeUpperLower<num_states>(parent_current_node_upper_lr, parent.getUpperLength(), *child_1_lower_lh, child_1_best_blength, aln, model, threshold_prob);
    grand_parent_upper_lr->mergeUpperLower<num_states>(parent.getPartialLh(child_1.getNeighborIndex(TOP).getMiniIndex()), parent.getUpperLength(), *current_node_lower_lh, current_node.getUpperLength(), aln, model, threshold_prob);
    // the upper_lr of the parent node is changed -> the aLRT of all its grand_children will be computed later '
    // the parent node has two children: child_1 and the current node
    // recompute aLRT of the children of child_1
    recompute_aLRT_GrandChildren(child_1, node_stack_aLRT);
    // recompute aLRT of the children of the current node -> will be done later (don't need to explicitly add here)
    // update the upper_lr of the current node
    const MiniIndex current_node_child_2_mini = child_2.getNeighborIndex(TOP).getMiniIndex();
    parent_current_node_upper_lr->mergeUpperLower<num_states>(current_node.getPartialLh(current_node_child_2_mini), current_node.getUpperLength(), *sibling_lower_lh, sibling_best_blength, aln, model, threshold_prob);
    const MiniIndex current_node_sibling_mini = sibling.getNeighborIndex(TOP).getMiniIndex();
    parent_current_node_upper_lr->mergeUpperLower<num_states>(current_node.getPartialLh(current_node_sibling_mini), current_node.getUpperLength(), *child_2_lower_lh, child_2_best_blength, aln, model, threshold_prob);
    // because upper_lr of the current node is changed -> we need to recompute aLRT of its grand_children
    // current node has two children: child_2 and sibling (after applying NNI)
    // recompute aLRT of the children of child_2
    recompute_aLRT_GrandChildren(child_2, node_stack_aLRT);
    // recompute aLRT of the children of sibling
    recompute_aLRT_GrandChildren(sibling, node_stack_aLRT);
    
    // update the lower_lh of the parent node
    std::unique_ptr<SeqRegions> new_lower_lh = nullptr;
    RealNumType new_lh_contribution = current_node_lower_lh->mergeTwoLowers<num_states>(new_lower_lh, current_node.getUpperLength(), *child_1_lower_lh, child_1_best_blength, aln, model, threshold_prob, true);
    
    // parent and other nodes on the path to root
    NumSeqsType node_vec = parent_index.getVectorIndex();
    while (true)
    {
        PhyloNode& node = nodes[node_vec];
        
        // if the new lower lh is different from the old one -> traverse upward further
        if (node.getPartialLh(TOP)->areDiffFrom(*new_lower_lh, seq_length, num_states, *params) || fabs(new_lh_contribution - node_lhs[node.getNodelhIndex()].getLhContribution()) > threshold_prob)
        {
            // update new_lower_lh
            node.setPartialLh(TOP, std::move(new_lower_lh));
            // std::cout << "lh_contribution (before): " << node_lhs[node.getNodelhIndex()].getLhContribution() << std::endl;
            node_lhs[node.getNodelhIndex()].setLhContribution(new_lh_contribution);
            // std::cout << "lh_contribution (after): " << node_lhs[node.getNodelhIndex()].getLhContribution() << std::endl;
            
            // cases when node is non-root
            if (root_vector_index != node_vec)
            {
                // re-caculate the lower likelihood (likelihood contribution) at the parent node
                const Index tmp_parent_index = node.getNeighborIndex(TOP);
                const NumSeqsType tmp_parent_vec = tmp_parent_index.getVectorIndex();
                PhyloNode& tmp_parent = nodes[tmp_parent_vec];
                const MiniIndex parent_sibling_mini = tmp_parent_index.getFlipMiniIndex();
                const Index sibling_index = tmp_parent.getNeighborIndex(parent_sibling_mini);
                PhyloNode& tmp_sibling = nodes[sibling_index.getVectorIndex()];
                const std::unique_ptr<SeqRegions>& node_lower_lh = node.getPartialLh(TOP);
                
                new_lh_contribution = node_lower_lh->mergeTwoLowers<num_states>(new_lower_lh, node.getUpperLength(), *(tmp_sibling.getPartialLh(TOP)), tmp_sibling.getUpperLength(), aln, model, params->threshold_prob, true);
                
                // move a step upwards
                node_vec = tmp_parent_vec;
                
                // because the lower_lh of the current node was updated
                // => we need to re-compute aLRT of that node => in this case, no need to add the current node because it should be added by its child before
                // => we need to re-compute aLRT of its parent node
                tmp_parent.setOutdated(true);
                node_stack_aLRT.push(Index(tmp_parent_vec, TOP));
                // => we need to re-compute aLRT of its sibling node
                tmp_sibling.setOutdated(true);
                node_stack_aLRT.push(sibling_index);
                // => we need to update the upper_lr of its sibling
                node_stack_update_upper_lr.push(sibling_index);
                // => we need to update the upper_lr of its parent (we can do it immediately)
                std::unique_ptr<SeqRegions> new_upper_lr = nullptr;
                std::unique_ptr<SeqRegions>& old_upper_lr = tmp_parent.getPartialLh(parent_sibling_mini);
                // if parent is root
                if (root_vector_index == tmp_parent_vec)
                {
                    node_lower_lh->computeTotalLhAtRoot<num_states>(new_upper_lr, model, node.getUpperLength());
                }
                // if parent is non-root
                else
                {
                    const std::unique_ptr<SeqRegions>& grand_parent_upper_lr = getPartialLhAtNode(tmp_parent.getNeighborIndex(TOP));
                    grand_parent_upper_lr->mergeUpperLower<num_states>(new_upper_lr, tmp_parent.getUpperLength(), *node_lower_lh, node.getUpperLength(), aln, model, threshold_prob);
                }
                // update the upper_lr of its parent if it's been changed
                if (!old_upper_lr || old_upper_lr->areDiffFrom(*new_upper_lr, seq_length, num_states, *params))
                {
                    old_upper_lr = std::move(new_upper_lr);
                    
                    // recompute aLRT of the grand_children of the parent node
                    recompute_aLRT_GrandChildren(tmp_sibling, node_stack_aLRT);
                }
            }
            // case when node is root
            else
            {
                // re-calculate likelihood at root
                // std::cout << "lh at root (before): " << lh_at_root << std::endl;
                lh_at_root = node.getPartialLh(TOP)->computeAbsoluteLhAtRoot<num_states>(model);
                // std::cout << "lh at root (after): " << lh_at_root << std::endl;
                
                // stop traversing further
                break;
            }
        }
        // otherwise, stop traversing further
        else
            break;
    }
    
    // traverse downward to update the upper_left/right_region until the changes is insignificant
    updateUpperLR<num_states>(node_stack_update_upper_lr, node_stack_aLRT);
}

template <const StateType num_states>
RealNumType cmaple::TreeBase::estimateBranchLengthWithCheck(const std::unique_ptr<SeqRegions>& upper_lr_regions, const std::unique_ptr<SeqRegions>& lower_regions, const RealNumType current_blength)
{
    // try to estimate a better blength
    RealNumType new_blength = estimateBranchLength<num_states>(upper_lr_regions, lower_regions);
    if (new_blength > 0 || current_blength > 0)
    {
        RealNumType diff_thresh = 0.01 * new_blength;
        if (new_blength <= 0 || current_blength <= 0 || (current_blength > (new_blength + diff_thresh)) || (current_blength < (new_blength - diff_thresh)))
            return new_blength;
    }
    
    // if not found, return the current blength
    return current_blength;
}

template <const StateType num_states>
void cmaple::TreeBase::updateUpperLR(std::stack<Index>& node_stack, std::stack<Index>& node_stack_aLRT)
{
    const RealNumType threshold_prob = params->threshold_prob;
    const PositionType seq_length = aln->ref_seq.size();
    
    while (!node_stack.empty())
    {
        const Index node_index = node_stack.top();
        node_stack.pop();
        PhyloNode& node = nodes[node_index.getVectorIndex()];
        if (node.isInternal())
        {
            const Index& left_child_index = node.getNeighborIndex(LEFT);
            const Index& right_child_index = node.getNeighborIndex(RIGHT);
            PhyloNode& left_child = nodes[left_child_index.getVectorIndex()];
            PhyloNode& right_child = nodes[right_child_index.getVectorIndex()];
            
            const std::unique_ptr<SeqRegions>& parent_upper_lr = getPartialLhAtNode(node.getNeighborIndex(TOP));
            std::unique_ptr<SeqRegions> tmp_upper_lr = nullptr;
            parent_upper_lr->mergeUpperLower<num_states>(tmp_upper_lr, node.getUpperLength(), *(left_child.getPartialLh(TOP)), left_child.getUpperLength(), aln, model, threshold_prob);
            
            // if upper_lr_lh of a node is changed -> we need to re-compute aLRT of that child and the grand-children
            if (updateNewPartialIfDifferent(node, RIGHT, tmp_upper_lr, node_stack, seq_length) && right_child.isInternal())
            {
                // re-compute aLRT of that child
                right_child.setOutdated(true);
                node_stack_aLRT.push(right_child_index);
                
                // re-compute aLRT of the grand-children
                recompute_aLRT_GrandChildren(right_child, node_stack_aLRT);
            }
            
            parent_upper_lr->mergeUpperLower<num_states>(tmp_upper_lr, node.getUpperLength(), *(right_child.getPartialLh(TOP)), right_child.getUpperLength(), aln, model, threshold_prob);
            
            // if upper_lr_lh of a node is changed -> we need to re-compute aLRT of that child and the grand-children
            if (updateNewPartialIfDifferent(node, LEFT, tmp_upper_lr, node_stack, seq_length) && left_child.isInternal())
            {
                // re-compute aLRT of that child
                left_child.setOutdated(true);
                node_stack_aLRT.push(left_child_index);
                
                // re-compute aLRT of the grand-children
                recompute_aLRT_GrandChildren(left_child, node_stack_aLRT);
            }
        }
    }
}

void cmaple::TreeBase::recompute_aLRT_GrandChildren(PhyloNode& parent, std::stack<Index>& node_stack_aLRT)
{
    const Index grand_child_1_index = parent.getNeighborIndex(LEFT);
    const Index grand_child_2_index = parent.getNeighborIndex(RIGHT);
    PhyloNode& grand_child_1 = nodes[grand_child_1_index.getVectorIndex()];
    PhyloNode& grand_child_2 = nodes[grand_child_2_index.getVectorIndex()];
    
    grand_child_1.setOutdated(true);
    node_stack_aLRT.push(grand_child_1_index);
    grand_child_2.setOutdated(true);
    node_stack_aLRT.push(grand_child_2_index);
}

template <const StateType num_states>
RealNumType cmaple::TreeBase::calculateSiteLhs(std::vector<RealNumType>& site_lh_contributions, std::vector<RealNumType>& site_lh_root)
{
    // initialize the total_lh by the likelihood from root
    const PositionType seq_length = aln->ref_seq.size();
    const RealNumType threshold_prob = params->threshold_prob;
    if (site_lh_contributions.size() != seq_length)
    {
        site_lh_contributions.resize(seq_length, 0);
        site_lh_root.resize(seq_length, 0);
    }
    RealNumType total_lh = nodes[root_vector_index].getPartialLh(TOP)->computeSiteLhAtRoot<num_states>(site_lh_root, model);

    // start from root
    Index node_index = Index(root_vector_index, TOP);
    Index last_node_index;
    
    // traverse to the deepest tip, calculate the likelihoods upward from the tips
    while (node_index.getMiniIndex() != UNDEFINED) //node)
    {
        PhyloNode& node = nodes[node_index.getVectorIndex()];
        // we reach a top node by a downward traversing
        if (node_index.getMiniIndex() == TOP) //node->is_top)
        {
            // if the current node is a leaf -> we reach the deepest tip -> traversing upward to calculate the lh of its parent
            if (!node.isInternal()) // node->isLeave())
            {
                /*last_node = node;
                 node = node->neighbor;*/
                last_node_index = node_index;
                node_index = node.getNeighborIndex(TOP);
            }
            // otherwise, keep traversing downward to find the deepest tip
            else
                // node = node->next->neighbor;
                node_index = node.getNeighborIndex(RIGHT);
        }
        // we reach the current node by an upward traversing from its children
        else
        {
            // if we reach the current node by an upward traversing from its first children -> traversing downward to its second children
            if (node.getNeighborIndex(RIGHT) == last_node_index) // node->getTopNode()->next->neighbor == last_node)
            {
                // node = node->getTopNode()->next->next->neighbor;
                node_index = node.getNeighborIndex(LEFT);
            }
            // otherwise, all children of the current node are updated -> update the lower lh of the current node
            else
            {
                // calculate the new lower lh of the current node from its children
                /*Node* top_node = node->getTopNode();
                 Node* next_node_1 = top_node->next;
                 Node* next_node_2 = next_node_1->next;*/
                Index neighbor_1_index = node.getNeighborIndex(RIGHT);
                Index neighbor_2_index = node.getNeighborIndex(LEFT);
                PhyloNode& neighbor_1 = nodes[neighbor_1_index.getVectorIndex()];
                PhyloNode& neighbor_2 = nodes[neighbor_2_index.getVectorIndex()];
                
                std::unique_ptr<SeqRegions> new_lower_lh = nullptr;
                const std::unique_ptr<SeqRegions>& lower_lh_1 = neighbor_1.getPartialLh(TOP); // next_node_1->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob);
                const std::unique_ptr<SeqRegions>& lower_lh_2 = neighbor_2.getPartialLh(TOP); // next_node_2->neighbor->getPartialLhAtNode(aln, model, params->threshold_prob);
                
                total_lh += lower_lh_1->calculateSiteLhContributions<num_states>(site_lh_contributions, new_lower_lh, neighbor_1.getUpperLength(), *lower_lh_2, neighbor_2.getUpperLength(), aln, model, threshold_prob);
                
                // if new_lower_lh is NULL
                if (!new_lower_lh)
                    outError("Strange, inconsistent lower genome list creation in calculateTreeLh(); old list, and children lists");
                // otherwise, everything is good -> update the lower lh of the current node
                else if (new_lower_lh->areDiffFrom(*node.getPartialLh(TOP), seq_length, num_states, *params))
                {
                    // outError("Strange, while calculating tree likelihood encountered non-updated lower likelihood!");
                    // non-updated lower likelihood may be due to a replacement of ML tree by an NNI neighbor
                    node.setPartialLh(TOP, std::move(new_lower_lh));
                }
                
                last_node_index = Index(node_index.getVectorIndex(), TOP);
                node_index = node.getNeighborIndex(TOP);
            }
        }
    }
    
    return total_lh;
}

NumSeqsType cmaple::TreeBase::parseFile(std::istream &infile, char& ch, RealNumType& branch_len, PositionType& in_line, PositionType& in_column, const std::map<std::string, NumSeqsType>& map_seqname_index, bool& missing_blengths)
{
    int maxlen = 1000;
    string seqname;
    int seqlen;
    RealNumType brlen = -1;

    // create new node
    NumSeqsType tmp_root_vec;
    // start with "(" -> an internal node
    if (ch == '(')
    {
        createAnInternalNode();
        tmp_root_vec = nodes.size() - 1;
    }
    // otherwise, it's a leaf
    else
    {
        createALeafNode(0);
        tmp_root_vec = nodes.size() - 1;
    }
    
    // start with "(" -> an internal node
    if (ch == '(')
    {
        MiniIndex child_mini = RIGHT;
        
        ch = readNextChar(infile, in_line, in_column);
        while (ch != ')' && !infile.eof())
        {
            const NumSeqsType tmp_node_vec = parseFile(infile, ch, brlen, in_line, in_column, map_seqname_index, missing_blengths);
            
            if (child_mini == UNDEFINED)
            {
                if (cmaple::verbose_mode > cmaple::VB_QUIET)
                    std::cout << "Converting a mutifurcating to a bifurcating tree" << std::endl;
                
                // create a new parent node
                createAnInternalNode();
                const NumSeqsType new_tmp_root_vec = nodes.size() - 1;
                // connect the current root node to the new parent node
                nodes[new_tmp_root_vec].setNeighborIndex(RIGHT, Index(tmp_root_vec, TOP));
                PhyloNode& current_root = nodes[tmp_root_vec];
                current_root.setNeighborIndex(TOP, Index(new_tmp_root_vec, RIGHT));
                current_root.setUpperLength(0);
                
                // the new parent becomes the current root node -> new child will be added as the left child of the (new) root node
                tmp_root_vec = new_tmp_root_vec;
                child_mini = LEFT;
            }
            
            PhyloNode& node = nodes[tmp_node_vec];
            nodes[tmp_root_vec].setNeighborIndex(child_mini, Index(tmp_node_vec, TOP));
            node.setNeighborIndex(TOP, Index(tmp_root_vec, child_mini));
            // If the branch length is not specify -> set it to default_blength and mark the tree with missing blengths so that we can re-estimate the blengths later
            if (brlen == -1)
            {
                brlen = default_blength;
                missing_blengths = true;
            }
            node.setUpperLength(brlen);
            
            // change to the second child
            child_mini = (child_mini == RIGHT) ? LEFT : UNDEFINED;
                
            if (infile.eof())
                throw "Expecting ')', but end of file instead";
            if (ch == ',')
                ch = readNextChar(infile, in_line, in_column);
            else if (ch != ')') {
                string err = "Expecting ')', but found '";
                err += ch;
                err += "' instead";
                throw err;
            }
        }
        if (!infile.eof()) ch = readNextChar(infile, in_line, in_column);
    }
    
    // now read the node name
    seqlen = 0;
    char end_ch = 0;
    if (ch == '\'' || ch == '"') end_ch = ch;
    seqname = "";

    while (!infile.eof() && seqlen < maxlen)
    {
        if (end_ch == 0) {
            if (is_newick_token(ch) || controlchar(ch)) break;
        }
        seqname += ch;
        seqlen++;
        ch = infile.get();
        in_column++;
        if (end_ch != 0 && ch == end_ch) {
            seqname += ch;
            seqlen++;
            break;
        }
    }
    if ((controlchar(ch) || ch == '[' || ch == end_ch) && !infile.eof())
        ch = readNextChar(infile, in_line, in_column, ch);
    if (seqlen == maxlen)
        throw "Too long name ( > 1000)";
    PhyloNode& tmp_root = nodes[tmp_root_vec];
    if (!tmp_root.isInternal())
        renameString(seqname);
//    seqname[seqlen] = 0;
    if (seqlen == 0 && !tmp_root.isInternal())
        throw "Redundant double-bracket ‘((…))’ with closing bracket ending at";
    if (seqlen > 0 && !tmp_root.isInternal())
    {
        auto it = map_seqname_index.find(seqname);
        if (it == map_seqname_index.end())
            outError("Leaf " + seqname + " is not found in the alignment. Please check and try again!");
        else
        {
            const NumSeqsType sequence_index = it->second;
            tmp_root.setSeqNameIndex(sequence_index);
            tmp_root.setPartialLh(TOP, aln->data[sequence_index].getLowerLhVector(aln->ref_seq.size(), aln->num_states, aln->getSeqType()));
            
            // mark the sequece as added (to the tree)
            sequence_added[sequence_index] = true;
        }        
    }
    /*if (!tmp_root.isInternal()) {
        // is a leaf, assign its ID
        if (leafNum == 0)
            cmaple::TreeBase::root_vector_index = root;
    }*/

    if (ch == ';' || infile.eof())
        return tmp_root_vec;
    // parse branch length
    if (ch == ':')
    {
        // string saved_comment = in_comment;
        ch = readNextChar(infile, in_line, in_column);
        /*if (in_comment.empty())
            in_comment = saved_comment;*/
        seqlen = 0;
        seqname = "";
        while (!is_newick_token(ch) && !controlchar(ch) && !infile.eof() && seqlen < maxlen)
        {
//            seqname[seqlen] = ch;
            seqname += ch;
            seqlen++;
            ch = infile.get();
            in_column++;
        }
        if ((controlchar(ch) || ch == '[') && !infile.eof())
            ch = readNextChar(infile, in_line, in_column, ch);
        if (seqlen == maxlen || infile.eof())
            throw "branch length format error.";
        branch_len = convert_real_number(seqname.c_str());
    }
    // handle the case of multiple-length branches but lack of the average branch length
    // e.g A[0.1/0.2/0.3/0.4] instead of A[0.1/0.2/0.3/0.4]:0.25
    /*else if (in_comment.length())
    {
        // set default average branch length at 0
        string default_avg_length = "0";
        parseBranchLength(default_avg_length, branch_len);
    }*/
    
    return tmp_root_vec;
}

const char cmaple::TreeBase::readNextChar(std::istream& in, PositionType& in_line, PositionType& in_column, const char& current_ch) const
{
    char ch;
    if (current_ch == '[')
        ch = current_ch;
    else {
        in.get(ch);
        in_column++;
        
        if (ch == 10)
        {
            in_line++;
            in_column = 1;
        }
    }
    // check if ch is a control character (ascii <= 32)
    while (controlchar(ch) && !in.eof())
    {
        in.get(ch);
        in_column++;
        
        if (ch == 10)
        {
            in_line++;
            in_column = 1;
        }
    }
    string in_comment = "";
    // ignore comment
    while (ch=='[' && !in.eof()) {
        while (ch!=']' && !in.eof()) {
            in.get(ch);
            if (ch != ']')
                in_comment += ch;
            in_column++;
            if (ch == 10) {
                in_line++;
                in_column = 1;
            }
        }
        if (ch != ']') throw "Comments not ended with ]";
        in_column++;
        in.get(ch);
        if (ch == 10) {
            in_line++;
            in_column = 1;
        }
        while (controlchar(ch) && !in.eof()) {
            in_column++;
            in.get(ch);
            if (ch == 10) {
                in_line++;
                in_column = 1;
            }
        }
        if (in_comment.length() && cmaple::verbose_mode > cmaple::VB_QUIET)
            std::cout << "Ignore [" + in_comment +"]" << std::endl;
    }
    return ch;
}

std::map<std::string, NumSeqsType> cmaple::TreeBase::initMapSeqNameIndex()
{
    ASSERT(aln);
    
    // create the map
    std::map<std::string, NumSeqsType> map_seqname_index;
    const std::vector<Sequence>& sequences = aln->data;
    for (NumSeqsType i = 0; i < sequences.size(); ++i)
        map_seqname_index.emplace(sequences[i].seq_name, i);
    
    return map_seqname_index;
}

NumSeqsType cmaple::TreeBase::markAnExistingSeq(const std::string& seq_name, const std::map<std::string, NumSeqsType>& map_name_index)
{
    NumSeqsType new_seq_index = 0;
    
    // Find the sequence name
    auto iter = map_name_index.find(seq_name);
    // If it's found -> mark it as added
    if (iter != map_name_index.end())
    {
        sequence_added[iter->second] = true;
        new_seq_index = iter->second;
    }
    // otherwise, return an error
    else
        outError("Taxon " + seq_name + " is not found in the new alignment!");
    
    // return new_seq_index
    return new_seq_index;
}

void cmaple::TreeBase::remarkExistingSeqs()
{
    ASSERT(aln);
    
    // init a mapping between sequence names and its index in the alignment
    std::map<std::string, NumSeqsType> map_name_index = initMapSeqNameIndex();
    
    // reset all marked sequences
    resetSeqAdded();
    
    // browse all nodes to mark existing leave as added
    for (auto i = 0; i < nodes.size(); ++i)
    {
        // only consider leave
        if (!nodes[i].isInternal())
        {
            PhyloNode& node = nodes[i];
            
            // mark the leaf itself
            node.setSeqNameIndex(markAnExistingSeq(seq_names[node.getSeqNameIndex()], map_name_index));
            
            // mark its less-info sequences
            std::vector<NumSeqsType>& less_info_seqs = node.getLessInfoSeqs();
            for (auto j = 0; j < less_info_seqs.size(); ++j)
                less_info_seqs[j] = markAnExistingSeq(seq_names[less_info_seqs[j]], map_name_index);
        }
    }
}

bool cmaple::TreeBase::readTree(std::istream& in)
{
    // Flag to check whether the tree contains missing branch length
    bool missing_blengths = false;
    
    // create a map between leave and sequences in the alignment
    std::map<std::string, NumSeqsType> map_seqname_index = initMapSeqNameIndex();
    
    if (cmaple::verbose_mode >= cmaple::VB_MED)
        std::cout << "Reading a tree" << std::endl;
    
    // Read tree from the stream
    PositionType in_line = 1;
    PositionType in_column = 1;
    //std::string in_comment{};
    
    try {
        char ch;
        ch = readNextChar(in, in_line, in_column);
        if (ch != '(') {
            cout << in.rdbuf() << endl;
            throw "Tree file does not start with an opening-bracket '('";
        }

        RealNumType branch_len;
        const NumSeqsType tmp_node_vec = parseFile(in, ch, branch_len, in_line, in_column, map_seqname_index, missing_blengths);
        
        // set root
        if (nodes[tmp_node_vec].isInternal())
            root_vector_index = tmp_node_vec;
        else
        {
            for (NumSeqsType i = 0; i < nodes.size(); ++i)
                if (nodes[i].isInternal())
                {
                    root_vector_index = i;
                    break;
                }
        }
        // 2018-01-05: assuming rooted tree if root node has two children
        /*if (is_rooted || (branch_len != 0.0) || node->degree() == 2) {
            if (branch_len == -1.0) branch_len = 0.0;
            if (branch_len < 0.0)
                throw ERR_NEG_BRANCH;
            root = newNode(leafNum, ROOT_NAME);
            root->addNeighbor(node, branch_len);
            node->addNeighbor(root, branch_len);
            leafNum++;
            rooted = true;
            
            // parse key/value from comment
            string KEYWORD="&";
            bool in_comment_contains_key_value = in_comment.length() > KEYWORD.length()
                                                  && !in_comment.substr(0, KEYWORD.length()).compare(KEYWORD);
            if (in_comment_contains_key_value)
                parseKeyValueFromComment(in_comment, root, node);
            
            
        } else { // assign root to one of the neighbor of node, if any
            FOR_NEIGHBOR_IT(node, NULL, it)
            if ((*it)->node->isLeaf()) {
                root = (*it)->node;
                break;
            }
        }
        // make sure that root is a leaf
        ASSERT(root->isLeaf());*/

        if (in.eof() || ch != ';')
            throw "Tree file must be ended with a semi-colon ';'";
    } catch (bad_alloc) {
        outError(ERR_NO_MEMORY);
    } catch (const char *str) {
        outError(str, " (line " + convertIntToString(in_line) + " column " + convertIntToString(in_column - 1) + ")");
    } catch (string str) {
        outError(str.c_str(), " (line " + convertIntToString(in_line) + " column " + convertIntToString(in_column - 1) + ")");
    } catch (ios::failure) {
        outError(ERR_READ_INPUT, " (line " + convertIntToString(in_line) + " column " + convertIntToString(in_column - 1) + ")");
    } catch (...) {
        // anything else
        outError(ERR_READ_ANY, " (line " + convertIntToString(in_line) + " column " + convertIntToString(in_column - 1) + ")");
    }
    
    if (cmaple::verbose_mode >= cmaple::VB_DEBUG)
        std::cout << "Collapsing zero-branch-length leaves into its sibling's vector of less-info-seqs..." << std::endl;
    collapseAllZeroLeave();
    
    return missing_blengths;
}

void cmaple::TreeBase::collapseAllZeroLeave()
{
    // count the number of collapsed nodes -> to make sure we have reseved enough space to store all nodes when we expand (i.e., adding less-informative sequences back to) the tree
    NumSeqsType num_collapsed_nodes = 0;
    
    // the default min_blength in CMaple is not small enough -> I set it at min_blength * 0.1 for a higher accuracy when calculating aLRT-SH
    const RealNumType new_min_blength = (params->fixed_min_blength == -1) ? min_blength * 0.1 : min_blength;
    
    // start from root
    Index node_index = Index(root_vector_index, TOP);
    Index last_node_index;
    
    // traverse to the deepest tip
    while (node_index.getMiniIndex() != UNDEFINED)
    {
        PhyloNode& node = nodes[node_index.getVectorIndex()];
        // we reach a top node by a downward traversing
        if (node_index.getMiniIndex() == TOP)
        {
            // if the current node is a leaf -> we reach the deepest tip -> traversing upward
            if (!node.isInternal())
            {
                last_node_index = node_index;
                node_index = node.getNeighborIndex(TOP);
            }
            // otherwise, keep traversing downward to find the deepest tip
            else
                node_index = node.getNeighborIndex(RIGHT);
        }
        // we reach the current node by an upward traversing from its children
        else
        {
            // if we reach the current node by an upward traversing from its first children -> traversing downward to its second children
            if (node.getNeighborIndex(RIGHT) == last_node_index)
            {
                node_index = node.getNeighborIndex(LEFT);
            }
            // otherwise, all children of the current node are updated
            else
            {
                // calculate the new lower lh of the current node from its children
                Index neighbor_1_index = node.getNeighborIndex(RIGHT);
                Index neighbor_2_index = node.getNeighborIndex(LEFT);
                PhyloNode& neighbor_1 = nodes[neighbor_1_index.getVectorIndex()];
                PhyloNode& neighbor_2 = nodes[neighbor_2_index.getVectorIndex()];
                
                if (neighbor_2.getUpperLength() <= 0 && neighbor_1.getUpperLength() <= 0)
                {
                    // only consider collapsing zero-branch-length leave into its sibling's less-info-seqs if they're both leave
                    if (!neighbor_1.isInternal() && !neighbor_2.isInternal())
                    {
                        if (root_vector_index == node_index.getVectorIndex() || node.getUpperLength() <= 0)
                        {
                            collapseOneZeroLeaf(node, node_index, neighbor_1, neighbor_1_index, neighbor_2);
                            ++num_collapsed_nodes;
                        }
                    }
                    // NHANLT - temporarily solution to keep writing and re-reading tree consistently - need more testing
                    /*else
                    {
                        if (cmaple::verbose_mode >= cmaple::VB_DEBUG)
                            std::cout << "Detect two zero-length branches connecting to an internal node. Reset one of those branches at the minmimum branch length " + convertDoubleToString(new_min_blength) + " to avoid potential error" << std::endl;
                        neighbor_1.setUpperLength(new_min_blength);
                    }*/
                }
                
                last_node_index = Index(node_index.getVectorIndex(), TOP);
                node_index = node.getNeighborIndex(TOP);
            }
        }
    }
    
    // update the capacity of nodes
    nodes.reserve(nodes.capacity() + num_collapsed_nodes + num_collapsed_nodes);
}

void cmaple::TreeBase::collapseOneZeroLeaf(PhyloNode& node, Index& node_index, PhyloNode& neighbor_1, const Index neighbor_1_index, PhyloNode& neighbor_2)
{
    if (cmaple::verbose_mode >= cmaple::VB_DEBUG)
        std::cout << "Collapse " << seq_names[neighbor_2.getSeqNameIndex()] << " into the vector of less-info_seqs of " << seq_names[neighbor_1.getSeqNameIndex()]  << std::endl;
    
    // add neighbor_2 and its less-info-seqs into that of neigbor_1
    neighbor_1.addLessInfoSeqs(neighbor_2.getSeqNameIndex());
    std::vector<NumSeqsType>& neighbor_1_vector = neighbor_1.getLessInfoSeqs();
    std::vector<NumSeqsType>& neighbor_2_vector = neighbor_2.getLessInfoSeqs();
    neighbor_1_vector.insert(neighbor_1_vector.end(), neighbor_2_vector.begin(), neighbor_2_vector.end());
    
    // if node is root -> neighbor_1 becomes the new root
    if (root_vector_index == node_index.getVectorIndex())
        root_vector_index = neighbor_1_index.getVectorIndex();
    // otherwise connect child_1 to the parent of node (~the grandparent of child_1)
    else
    {
        // the length of the new branch (child_1 - its grandparent) = sum(child_1->upperlength + node->upperlength)
        RealNumType new_blength = neighbor_1.getUpperLength();
        if (node.getUpperLength() > 0)
            new_blength = new_blength > 0 ? new_blength + node.getUpperLength() : node.getUpperLength();
        
        // connect child_1 and its grandparent
        const Index parent_index = node.getNeighborIndex(TOP);
        const MiniIndex parent_mini = parent_index.getMiniIndex();
        nodes[parent_index.getVectorIndex()].setNeighborIndex(parent_mini, neighbor_1_index);
        neighbor_1.setNeighborIndex(TOP, parent_index);
        neighbor_1.setUpperLength(new_blength);
    }
    
    node_index = neighbor_1_index;
}

template <const StateType num_states>
void cmaple::TreeBase::updateModelParams()
{
    // perform a DFS -> at each leaf, update the pesudoCount of the model based on the sequence of that leaf
    performDFSAtLeave<&cmaple::TreeBase::updatePesudoCountModel<num_states>>();
    
    // update model params based on the pseudo count
    model->updateMutationMatEmpirical(aln);
}

template <const StateType num_states>
void cmaple::TreeBase::updatePesudoCountModel(PhyloNode& node, const Index node_index, const Index parent_index)
{
    std::unique_ptr<SeqRegions>& upper_lr_regions = getPartialLhAtNode(parent_index);
    std::unique_ptr<SeqRegions>& lower_regions = node.getPartialLh(TOP);
    if (upper_lr_regions && lower_regions)
        model->updatePesudoCount(aln, *upper_lr_regions, *lower_regions);
}

template <const StateType num_states>
void cmaple::TreeBase::expandTreeByOneLessInfoSeq(PhyloNode& node, const Index node_index, const Index parent_index)
{
    // don't expand the tree if the internal branch length is zero or this node doesn't have less-info-sequences
    if (node.getLessInfoSeqs().empty() || node.getUpperLength() <= 0 || parent_index.getMiniIndex() == UNDEFINED) return;
    std::vector<NumSeqsType>& less_info_seqs = node.getLessInfoSeqs();
    const NumSeqsType seq_name_index = less_info_seqs.back();
    less_info_seqs.pop_back();
    
    // debug
    if (cmaple::verbose_mode >= cmaple::VB_DEBUG)
        std::cout << "Add less-info-seq " + seq_names[seq_name_index] + " into the tree" << std::endl;
    
    // dummy variables
    std::unique_ptr<SeqRegions> lower_regions = aln->data[seq_name_index].getLowerLhVector(aln->ref_seq.size(), num_states, aln->getSeqType());
    const std::unique_ptr<SeqRegions>& upper_left_right_regions = getPartialLhAtNode(parent_index);
    std::unique_ptr<SeqRegions> best_child_regions = nullptr;
    const RealNumType top_distance = node.getUpperLength();
    upper_left_right_regions->mergeUpperLower<num_states>(best_child_regions, top_distance, *(node.getPartialLh(TOP)), 0, aln, model, params->threshold_prob);
    
    // add a new node representing the less-info-seq
    // the default min_blength in CMaple is not small enough -> I set it at min_blength * 0.1 for a higher accuracy when calculating aLRT-SH
    const RealNumType new_min_blength = (params->fixed_min_blength == -1) ? min_blength * 0.1 : min_blength;
    connectNewSample2Branch<num_states>(lower_regions, seq_name_index, node_index, node, top_distance, 0, new_min_blength, best_child_regions, upper_left_right_regions);
}

template <void(cmaple::TreeBase::*task)(PhyloNode&, const Index, const Index)>
void cmaple::TreeBase::performDFSAtLeave()
{
    // start from root
    Index node_index = Index(root_vector_index, TOP);
    Index last_node_index;
    
    // traverse to the deepest tip
    while (node_index.getMiniIndex() != UNDEFINED)
    {
        PhyloNode& node = nodes[node_index.getVectorIndex()];
        // we reach a top node by a downward traversing
        if (node_index.getMiniIndex() == TOP)
        {
            // if the current node is a leaf -> we reach the deepest tip -> traversing upward
            if (!node.isInternal())
            {
                // do some tasks at a leaf
                const Index parent_index = node.getNeighborIndex(TOP);
                (this->*task)(node, node_index, parent_index);
                
                last_node_index = node_index;
                node_index = parent_index;
            }
            // otherwise, keep traversing downward to find the deepest tip
            else
                node_index = node.getNeighborIndex(RIGHT);
        }
        // we reach the current node by an upward traversing from its children
        else
        {
            // if we reach the current node by an upward traversing from its first children -> traversing downward to its second children
            if (node.getNeighborIndex(RIGHT) == last_node_index)
            {
                node_index = node.getNeighborIndex(LEFT);
            }
            // otherwise, all children of the current node are reached
            else
            {
                // traverse upward
                last_node_index = Index(node_index.getVectorIndex(), TOP);
                node_index = node.getNeighborIndex(TOP);
            }
        }
    }
}

template <const StateType num_states>
void cmaple::TreeBase::updateBlengthReplaceMLTree(std::stack<Index>& node_stack_aLRT, RealNumType& lh_diff, PhyloNode& node, const Index node_index, const RealNumType best_blength)
{
    RealNumType old_blength = node.getUpperLength();
    
    // update blength
    node.setUpperLength(best_blength);
    
    if (!node.isInternal() && old_blength <= 0 && best_blength > 0 && !node.getLessInfoSeqs().empty())
        addLessInfoSeqReplacingMLTree<num_states>(node_stack_aLRT, lh_diff, node, node_index, node.getNeighborIndex(TOP));
}

template <const StateType num_states>
void cmaple::TreeBase::addLessInfoSeqReplacingMLTree(std::stack<Index>& node_stack_aLRT, RealNumType& lh_diff, PhyloNode& node, const Index node_index, const Index parent_index)
{
    // add one less-info-seq of the current leaf into the tree
    expandTreeByOneLessInfoSeq<num_states>(node, node_index, parent_index);
    
    const Index new_internal_from_node_index = node.getNeighborIndex(TOP);
    const NumSeqsType new_internal_vec = new_internal_from_node_index.getVectorIndex();
    PhyloNode& new_internal = nodes[new_internal_vec];
    const Index sibling_index = new_internal.getNeighborIndex(new_internal_from_node_index.getFlipMiniIndex());
    PhyloNode& sibling = nodes[sibling_index.getVectorIndex()];
    
    // compute the likelihood contribution at the new_internal and update the total_lh
    std::unique_ptr<SeqRegions> new_lower_lh = nullptr;
    computeLhContribution<num_states>(lh_diff, new_lower_lh, new_internal, node.getPartialLh(TOP), sibling.getPartialLh(TOP), node_index, node, sibling_index, sibling, aln->ref_seq.size());
    
    // add newly created internal node to the queue to compute the aLRT later
    new_internal.setOutdated(true);
    node_stack_aLRT.push(Index(new_internal_vec, TOP));
}

void cmaple::TreeBase::resetSPRFlags(const bool n_SPR_applied, const bool update_outdated, const bool n_outdated)
{
    // browse all nodes
    for (auto i = 0; i < nodes.size(); ++i)
    {
        PhyloNode& node = nodes[i];
        
        // update SPR_applied
        node.setSPRApplied(n_SPR_applied);
        
        // update outdated if necessary
        if (update_outdated)
            node.setOutdated(n_outdated);
    }
}

bool cmaple::TreeBase::isComplete()
{
    // make sure aln is not null
    if (aln)
    {
        // browse sequences in the alignment one by one
        for (auto i = 0; i < aln->data.size(); ++i)
            // if any of sequence has yet added -> this tree is incomplete
            if (!sequence_added[i])
                return false;
    }
    
    // Return complete (by default)
    return true;
}
