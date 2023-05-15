# DNA data
add_library(maple-dna
cmaple.h cmaple.cpp
cmaple_api.h cmaple_api.cpp
)
target_link_libraries(maple-dna tree-dna alignment-dna)

# Protein data
add_library(maple-aa
cmaple.h cmaple.cpp
)
target_link_libraries(maple-aa tree-aa alignment-aa)

