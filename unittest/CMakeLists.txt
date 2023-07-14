enable_testing()

file(COPY "${PROJECT_SOURCE_DIR}/example" DESTINATION "${PROJECT_BINARY_DIR}")

add_executable(
  maintest          ## todo: should be called 'SeqRegions_test'
  #tree_test.h tree_test.cpp
  alignment_test.cpp
  model_test.cpp
  matrix_test.cpp
  mutation_test.cpp
  phylonode_test.cpp
  seqregion_test.cpp
  seqregions_test.cpp
  sequence_test.cpp
)
target_link_libraries(
  maintest
  GTest::gtest_main    ## contains the main() function. No need to write one yourself
  utils
  model-aa
  alignment-aa
  tree-aa
  maple-aa
)

include(GoogleTest)
gtest_discover_tests(maintest)
