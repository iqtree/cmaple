# DNA data
add_library(cmaple_model
model.h model.cpp
modelbase.h modelbase.cpp
model_dna.h model_dna.cpp
model_aa.h model_aa.cpp
model_dna_rate_variation.h model_dna_rate_variation.cpp
)
target_link_libraries(cmaple_model cmaple_alignment cmaple_utils ncl nclextra)

if (USE_CMAPLE_AA)
    # Protein data
    add_library(cmaple_model-aa
    model.h model.cpp
    modelbase.h modelbase.cpp
    model_dna.h model_dna.cpp
    model_aa.h model_aa.cpp
    model_dna_rate_variation.h model_dna_rate_variation.cpp
    )
    target_link_libraries(cmaple_model-aa cmaple_alignment-aa cmaple_utils ncl nclextra)
endif()
