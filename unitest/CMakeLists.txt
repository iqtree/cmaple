enable_testing()

add_executable(
  maintest
  maintest.cpp
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