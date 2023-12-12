#include "gtest/gtest.h"
#include "../maple/cmaple.h"

using namespace cmaple;

TEST(CMapleTest, GetVersionTest) {
    // Call the method under test
    std::string version = cmaple::getVersion();

    // Check the result
    EXPECT_EQ(version, "CMAPLE " + convertIntToString(cmaple_VERSION_MAJOR) + "." + convertIntToString(cmaple_VERSION_MINOR) + cmaple_VERSION_PATCH);
}

TEST(CMapleTest, GetCitationsTest) {
    // Call the method under test
    std::string citations = cmaple::getCitations();

    // Check the result
    EXPECT_TRUE(citations.length() > 0);
}

TEST(CMapleTest, checkMapleSuitability) {
    
    // detect the path to the example directory
    std::string example_dir = "../../example/";
    if (!fileExists(example_dir + "example.maple"))
        example_dir = "../example/";
    
    Alignment aln(example_dir + "test_100.maple");
    // Check the result
    EXPECT_TRUE(cmaple::checkMapleSuitability(aln));
    
    aln.read(example_dir + "test_5K.maple");
    // Check the result
    EXPECT_TRUE(cmaple::checkMapleSuitability(aln));
    
    aln.read(example_dir + "input.fa");
    // Check the result
    EXPECT_FALSE(cmaple::checkMapleSuitability(aln));
    
    aln.read(example_dir + "input.phy");
    // Check the result
    EXPECT_FALSE(cmaple::checkMapleSuitability(aln));;
}
