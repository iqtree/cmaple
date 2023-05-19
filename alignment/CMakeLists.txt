# DNA data
add_library(cmaple_alignment-dna
mutation.h mutation.cpp
seqregion.h seqregion.cpp
seqregions.h seqregions.cpp
sequence.h sequence.cpp
alignment.h alignment.cpp
)
target_compile_definitions(cmaple_alignment-dna PUBLIC NUM_STATES=4)
target_link_libraries(cmaple_alignment-dna cmaple_utils)

# Protein data
add_library(cmaple_alignment-aa
mutation.h mutation.cpp
seqregion.h seqregion.cpp
seqregions.h seqregions.cpp
sequence.h sequence.cpp
alignment.h alignment.cpp
)
target_compile_definitions(cmaple_alignment-aa PUBLIC NUM_STATES=20)
target_link_libraries(cmaple_alignment-aa cmaple_utils)
