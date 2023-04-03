# DNA data
add_library(alignment-dna
mutation.h mutation.cpp
seqregion.h seqregion.cpp
seqregions.h seqregions.cpp
sequence.h sequence.cpp
alignment.h alignment.cpp
)
target_compile_definitions(alignment-dna PUBLIC NUM_STATES=4)
target_link_libraries(alignment-dna)

# Protein data
add_library(alignment-aa
mutation.h mutation.cpp
seqregion.h seqregion.cpp
seqregions.h seqregions.cpp
sequence.h sequence.cpp
alignment.h alignment.cpp
)
target_compile_definitions(alignment-aa PUBLIC NUM_STATES=20)
target_link_libraries(alignment-aa)
