add_library(alignment
mutation.h mutation.cpp
mutationindel.h mutationindel.cpp
seqregion.h seqregion.cpp
seqregions.h seqregions.cpp
sequence.h sequence.cpp
alignment.h alignment.cpp
)

target_link_libraries(alignment model)
