# DNA data
add_library(tree-dna
tree.h tree.cpp
updatingnode.h updatingnode.cpp
traversingnode.h traversingnode.cpp
phylonode.h phylonode.cpp
leaf.h
internal.h
)
target_link_libraries(tree-dna model-dna alignment-dna)

# Protein data
add_library(tree-aa
tree.h tree.cpp
updatingnode.h updatingnode.cpp
traversingnode.h traversingnode.cpp
phylonode.h phylonode.cpp
leaf.h
internal.h
)
target_link_libraries(tree-aa model-aa alignment-aa)
