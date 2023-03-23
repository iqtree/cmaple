add_library(model
model.h model.cpp
model_dna.h model_dna.cpp
model_aa.h model_aa.cpp
)

target_link_libraries(model utils)
