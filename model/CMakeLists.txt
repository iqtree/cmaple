# DNA data
add_library(model-dna
model.h model.cpp
model_dna.h model_dna.cpp
model_aa.h model_aa.cpp
)
target_link_libraries(model-dna alignment-dna utils)

# Protein data
add_library(model-aa
model.h model.cpp
model_dna.h model_dna.cpp
model_aa.h model_aa.cpp
)
target_link_libraries(model-aa alignment-aa utils)
