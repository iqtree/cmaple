# DNA data
add_library(cmaple_model-dna
model.h model.cpp
model_dna.h model_dna.cpp
model_aa.h model_aa.cpp
)
target_link_libraries(cmaple_model-dna cmaple_alignment-dna cmaple_utils)

# Protein data
add_library(cmaple_model-aa
model.h model.cpp
model_dna.h model_dna.cpp
model_aa.h model_aa.cpp
)
target_link_libraries(cmaple_model-aa cmaple_alignment-aa cmaple_utils)
