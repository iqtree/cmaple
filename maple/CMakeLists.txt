# DNA data
add_library(maple
cmaple.h cmaple.cpp
)
target_link_libraries(maple cmaple_tree cmaple_alignment cmaple_model cmaple_utils)

# Protein data
add_library(maple-aa
cmaple.h cmaple.cpp
)
target_link_libraries(maple-aa cmaple_tree-aa cmaple_alignment-aa cmaple_model-aa cmaple_utils)

