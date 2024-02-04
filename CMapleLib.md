# What's CMAPLE library?

CMAPLE is a C++ reimplementation of [MAPLE](https://www.nature.com/articles/s41588-023-01368-0) - a novel likelihood-based phylogenetic inference method for pandemic-scale epidemiological genomic data. CMAPLE is highly optimized for performance and scalability with many new features.

CMAPLE library provides a set of APIs, which allow users to integrate CMAPLE into existing phylogenetic inference methods.

# How to use?
In the following, we give an instruction of how to include CMAPLE library into a project using CMake build system (instructions for other build systems will be updated soon). Here, we'll take [IQ-TREE](https://github.com/iqtree/iqtree2) as an example.

## Include CMAPLE as a submodule of your project

In your project directory, run
    	
    	git submodule add https://github.com/iqtree/cmaple.git
    	
   
## Update CMakeList to include CMAPLE

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
		
    	# Step 3: Add the binary tree to the search path for include files so that we can find cmaple/cmaple_config.h
    	include_directories("${PROJECT_BINARY_DIR}/cmaple")
    	
    	# Step 4: Add linking libraries
    	target_link_libraries(iqtree2 main maple)
    	target_link_libraries(iqtree2-aa main maple-aa)

## An example of APIs usage

    	#include “cmaple.h”
    	
    	// Read an alignment from a FASTA file
    	cmaple::Alignment aln("alignment.fasta");
    	
    	// Check if the [C]Maple algorithm is effective to analyse the input alignment
    	if (cmaple::isEffective(aln))
    	{
    		// Initialise a substitution model as GTR
    		cmaple::Model model(cmaple::ModelBase::GTR);
    	
    		// Declare a tree, attaching it to the alignment and model 
			cmaple::Tree tree(&aln, &model);
    	
    		// Infer a maximum likelihood tree
    		tree.infer();
    	
    		// Compute aLRT-SH branch supports of the inferred tree
			tree.computeBranchSupport();
    	
    		// Compute the likelihood of the tree
    		cout << "Log-likelihood: " << tree.computeLh() << endl;
    	
    		// Export the tree in NEWICK format
    		cout << "Tree: " << tree.exportNewick() << endl;
    	}
    	else
    	{
    		// Execute existing methods (e.g., IQ-TREE, RAXML, PHYML) to analyse this alignment
    	}

## Tips for Debugging
CMAPLE outputs debugging messages to the standard output `std::cout`. One could control the amount of those messages via setting `cmaple::verbose_mode` to one of the following values.
<br> - `VB_QUIET`: no messages except errors.
<br> - `VB_MED` (default): common messages (showing the processing progress).
<br> - `VB_DEBUG`: as many messages as possible (useful for debugging).


<br>
# More APIs?
CMAPLE APIs are exposed in the following.
<br>[**CMaple**](group__cmaple.html)
<br>[**Alignment**](classcmaple_1_1_alignment.html)
<br>[**Model**](classcmaple_1_1_model.html)
<br>[**Tree**](classcmaple_1_1_tree.html)


<br>
# How to cite CMAPLE?

To be updated...


