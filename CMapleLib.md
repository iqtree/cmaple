# What's CMaple library?

Introduction to CMaple library...


<br>
# How to use?
Include CMaple as a submodule...

### An example of usage

    	#include “cmaple.h”
    	
    	// Create an alignment from a file
    	cmaple::Alignment aln("alignment.maple");
    	
    	// Create an auto-selected model according to the data type from the alignment (i.e., GTR for DNA, and LG for protein data)
    	Model model(cmaple::ModelBase::MODEL_AUTO, aln.getSeqType());
    	
    	// Create a tree, attach the alignment and model to the tree
    	cmaple::Tree tree(&aln, &model);
    	
    	// Infer a phylogenetic tree from the alignment and the model using MAPLE algorithm
    	cout << tree.autoProceedMAPLE() << endl;
    	
    	// Compute the branch supports for the inferred tree
    	cout << tree.computeBranchSupport() << endl;
    	
    	// Compute the likelihood of the tree
    	cout << "- Tree log likelihood: " << tree.computeLh() << endl;
    	
    	// Export the tree (with branch supports) in NEWICK format
    	cout << "- Tree: " << tree.exportNewick("cmaple::Tree::BIN_TREE, true) << endl;

### Tips for Debugging
CMaple outputs debugging messages to the standard output `std::cout`. One could control the amount of those messages via setting `cmaple::verbose_mode` to one of the following values.
<br> - `VB_QUIET`: no messages except errors.
<br> - `VB_MED` (default): common messages (showing the processing progress).
<br> - `VB_DEBUG`: as many messages as possible (useful for debugging).


<br>
# More APIs?
CMaple APIs are exposed in the following classes:
<br>[**Alignment**](classcmaple_1_1_alignment.html)
<br>[**Model**](classcmaple_1_1_model.html)
<br>[**Tree**](classcmaple_1_1_tree.html)


<br>
# How to cite CMaple?
CMaple manuscript...


