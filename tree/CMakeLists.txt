# DNA data
add_library(cmaple_tree-dna
tree.h tree.cpp
updatingnode.h updatingnode.cpp
traversingnode.h traversingnode.cpp
phylonode.h phylonode.cpp
leaf.h
internal.h
)
target_link_libraries(cmaple_tree-dna cmaple_model-dna cmaple_alignment-dna)

# Protein data
add_library(cmaple_tree-aa
tree.h tree.cpp
updatingnode.h updatingnode.cpp
traversingnode.h traversingnode.cpp
phylonode.h phylonode.cpp
leaf.h
internal.h
)
target_link_libraries(cmaple_tree-aa cmaple_model-aa cmaple_alignment-aa)
