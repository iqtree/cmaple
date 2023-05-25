#include "../tree/tree.h"
#include "../utils/timeutil.h"
#include "../alignment/seqregions.h"

#pragma once

/** CMaple class contains all main methods*/
class CMaple {
private:
    
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
    
    /**
        Template of doInference()
     */
    template <const cmaple::StateType  num_states>
    void doInferenceTemplate();
    
    /**
        Template of postInference()
     */
    template <const cmaple::StateType  num_states>
    void postInferenceTemplate();
    
    /**
        Setup function pointers
     */
    void setupFuncPtrs(const cmaple::StateType  num_states);
    
    /**
        Build an Initial Tree
        @return TRUE if there is any node added, thus, we need to optimize the tree later
     */
    template <const cmaple::StateType  num_states>
    bool buildInitialTree();
    
    /**
        Optimize the current tree
     */
    template <const cmaple::StateType  num_states>
    void optimizeTree(const bool new_sequences_added);
    
    /**
        Optimize the tree topology
     */
    template <const cmaple::StateType  num_states>
    void optimizeTreeTopology(bool short_range_search = false);
    
    /**
        Optimize the branch lengths of the current tree
     */
    template <const cmaple::StateType  num_states>
    void optimizeBranchLengthsOfTree();
    
    /**
        export output files
     */
    void exportOutput(const std::string &filename, const bool show_branch_support = false);
    
    /**
        calculate branch supports
     */
    template <const cmaple::StateType  num_states>
    void calculateBranchSupports();
    
    /**
        Load input tree
     */
    template <const cmaple::StateType  num_states>
    void loadInputTree();
    
public:
    
    /** The phylogenetic tree */
    cmaple::Tree tree;
    
    /**
    *  CMaple constructor
    */
    CMaple() = default;
    
    /**
    *  CMaple constructor
     @param params user-specified parameters
    */
    // CMaple(cmaple::Params params):tree(std::move(params)){};
    CMaple(cmaple::Params params):tree(params){};
    
    /**
    *  CMaple destructor
    */
    ~CMaple() = default;
    
    /**
        Load input data
     */
    void loadInput();
    
    /**
        Prepare for the inference
     */
    void preInference();
    
    /**
        Do the inference
     */
    void doInference();
    
    /**
        Complete the inference
     */
    void postInference();
};

/** Method to run CMaple */
void runCMaple(cmaple::Params &params);
