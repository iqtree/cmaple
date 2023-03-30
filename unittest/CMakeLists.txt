enable_testing()

file(COPY "${PROJECT_SOURCE_DIR}/example" DESTINATION "${PROJECT_BINARY_DIR}")

add_executable(
  maintest          ## todo: should be called 'SeqRegions_test'
  #tree_test.h tree_test.cpp
  #phylonode_test.cpp
  #seqregions_test.cpp
  #model_test.cpp
  #alignment_test.cpp
  #sequence_test.cpp
  #seqregion_test.cpp
  mutation_test.cpp
)
target_link_libraries(
  maintest
  GTest::gtest_main    ## contains the main() function. No need to write one yourself
  utils
  model-dna
  alignment-dna
  tree-dna
  maple-dna
)

include(GoogleTest)
gtest_discover_tests(maintest)
