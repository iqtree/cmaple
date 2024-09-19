# DNA data
add_library(cmaple_tree
tree.h tree.cpp
updatingnode.h updatingnode.cpp
traversingnode.h traversingnode.cpp
phylonode.h phylonode.cpp
leaf.h
internal.h
altbranch.h
)
target_link_libraries(cmaple_tree cmaple_model cmaple_alignment cmaple_utils)

# Protein data
add_library(cmaple_tree-aa
tree.h tree.cpp
updatingnode.h updatingnode.cpp
traversingnode.h traversingnode.cpp
phylonode.h phylonode.cpp
leaf.h
internal.h
altbranch.h
)
target_link_libraries(cmaple_tree-aa cmaple_model-aa cmaple_alignment-aa cmaple_utils)
