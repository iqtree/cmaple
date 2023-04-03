#include "tree/tree.h"
#include "utils/timeutil.h"
#include "alignment/seqregions.h"

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
    template <const StateType num_states>
    void doInferenceTemplate();
    
    /**
        Template of postInference()
     */
    template <const StateType num_states>
    void postInferenceTemplate();
    
    /**
        Setup function pointers
     */
    void setupFuncPtrs(const StateType num_states);
    
    /**
        Build an Initial Tree
     */
    template <const StateType num_states>
    void buildInitialTree();
    
    /**
        Optimize the current tree
     */
    template <const StateType num_states>
    void optimizeTree();
    
    /**
        Optimize the tree topology
     */
    template <const StateType num_states>
    void optimizeTreeTopology(bool short_range_search = false);
    
    /**
        Optimize the branch lengths of the current tree
     */
    template <const StateType num_states>
    void optimizeBranchLengthsOfTree();
    
    /**
        export output files
     */
    void exportOutput(const std::string &filename, const bool show_branch_support = false);
    
    /**
        calculate branch supports
     */
    template <const StateType num_states>
    void calculateBranchSupports();
    
public:
    
    /** The phylogenetic tree */
    Tree tree;
    
    /**
    *  CMaple constructor
    */
    CMaple() = default;
    
    /**
    *  CMaple constructor
     @param params user-specified parameters
    */
    CMaple(Params params):tree(std::move(params)){};
    
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
void runCMaple(Params &params);
