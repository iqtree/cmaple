add_library(alignment
mutation.h mutation.cpp
mutationindel.h mutationindel.cpp
region.h region.cpp
regions.h regions.cpp
sequence.h sequence.cpp
alignment.h alignment.cpp
)

target_link_libraries(alignment model)
