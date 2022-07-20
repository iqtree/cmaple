add_library(tree
tree.h tree.cpp
node.h node.cpp
)

target_link_libraries(tree model alignment)
