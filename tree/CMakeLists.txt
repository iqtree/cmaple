add_library(tree
node.h node.cpp
)

target_link_libraries(tree model alignment)
