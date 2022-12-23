enable_testing()

add_executable(
  maintest          ## todo: should be called 'SeqRegions_test'
  seqregions_test.h seqregions_test.cpp
)
target_link_libraries(
  maintest
  GTest::gtest_main    ## contains the main() function. No need to write one yourself
  utils
  model
  alignment
  tree
  maple
)

include(GoogleTest)
gtest_discover_tests(maintest)
