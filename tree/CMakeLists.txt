add_library(tree
tree.h tree.cpp
updatingnode.h updatingnode.cpp
traversingnode.h traversingnode.cpp
phylonode.h phylonode.cpp
leaf.h
internal.h
)

target_link_libraries(tree model alignment)
