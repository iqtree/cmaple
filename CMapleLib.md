# What's CMaple library?

Introduction to CMaple library...


<br>
# How to use?
In the following, we give an instruction of how to include CMaple library into a project using CMake build system (instructions for other build systems will be updated soon). Here, we'll take [IQ-TREE](https://github.com/iqtree/iqtree2) as example.

### Include CMaple as a submodule of your project

In your project directory, run
    	
    	git submodule add https://github.com/trongnhanuit/cmaple.git
    	
   
### Update CMakeList to include CMaple

    	project(iqtree)
    	
    	# Step 1: Add cmaple directory to the build
    	add_subdirectory(cmaple)
    	# add the main and other directories of your project
    	add_subdirectory(main)
    	
    	# Step 2: Build two executables for DNA and protein data
    	# build an executable for DNA data
    	add_executable(iqtree2
			main.cpp main.h
		)
		# build another executable for protein data
    	add_executable(iqtree2-aa
			main.cpp main.h
		)
		
    	# Step 3: Add the binary tree to the search path for include files so that we will find cmaple/cmaple_config.h
    	include_directories("${PROJECT_BINARY_DIR}/cmaple")
    	
    	# Step 4: Add linking libraries
    	target_link_libraries(iqtree2 main maple cmaple_alignment)
    	target_link_libraries(iqtree2-aa main maple-aa cmaple_alignment-aa)

### An example of APIs usage

    	#include “cmaple.h”
    	
    	// Create an alignment from a file
    	cmaple::Alignment aln("alignment.maple");
    	
    	// Check if the input alignment is suitable for using [C]Maple method
    	if (checkMapleSuitability(aln))
    	{
    		// Create a default model according to the data type from the alignment (i.e., GTR for DNA, and LG for protein data)
    		Model model(cmaple::ModelBase::DEFAULT, aln.getSeqType());
    	
    		// Create a tree, attach the alignment and model to the tree
    		cmaple::Tree tree(&aln, &model);
    	
    		// Infer a phylogenetic tree from the alignment and the model using [C]Maple algorithm
    		cout << tree.autoProceedMAPLE() << endl;
    	
    		// Compute the branch supports for the inferred tree
    		cout << tree.computeBranchSupport() << endl;
    	
    		// Compute the likelihood of the tree
    		cout << "- Tree log likelihood: " << tree.computeLh() << endl;
    	
    		// Export the tree (with branch supports) in NEWICK format
    		cout << "- Tree: " << tree.exportNewick("cmaple::Tree::BIN_TREE, true) << endl;
    	}
    	else
    		std::cout << "The input sequences are too divergent, which are inappropriate for [C]Maple method. We highly recommend users to use other methods (e.g., IQ-TREE, RAXML) to analyse this alignment!" << std::endl;

### Tips for Debugging
CMaple outputs debugging messages to the standard output `std::cout`. One could control the amount of those messages via setting `cmaple::verbose_mode` to one of the following values.
<br> - `VB_QUIET`: no messages except errors.
<br> - `VB_MED` (default): common messages (showing the processing progress).
<br> - `VB_DEBUG`: as many messages as possible (useful for debugging).


<br>
# More APIs?
CMaple APIs are exposed in the following.
<br>[**CMaple**](group__cmaple.html)
<br>[**Alignment**](classcmaple_1_1_alignment.html)
<br>[**Model**](classcmaple_1_1_model.html)
<br>[**Tree**](classcmaple_1_1_tree.html)


<br>
# How to cite CMaple?
CMaple manuscript...


