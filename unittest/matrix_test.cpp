#include "gtest/gtest.h"
#include "utils/matrix.h"

#include <algorithm> // for std::for_each
#include <array>
#include <numeric> // for std::iota

/*
 Test template <int length, typename RealType> inline RealType dotProduct(const RealType* p1, const RealType* p2)
 */

template <typename FloatT> void testDotProduct()
{
  std::array<FloatT, 30> data;
  std::iota(std::begin(data), std::end(data), 0); // Fill with 0...29.

  std::array<FloatT, 30> data_neg = data;
  // make numbers negative
  std::for_each(data_neg.begin(), data_neg.end(), [](FloatT& n) { n = -n; });

  // test the generic function
  ASSERT_FLOAT_EQ(dotProduct<10>(&data[0], &data[0]), 385);
  ASSERT_FLOAT_EQ(dotProduct<10>(&data[0], &data_neg[0]), -385);


  // test SIMD/AVX overloads
  ASSERT_FLOAT_EQ(dotProduct<20>(&data[0], &data[0]), 2870);
  ASSERT_FLOAT_EQ(dotProduct<20>(&data[0], &data_neg[0]), -2870);
}

TEST(Model, dotProduct)
{
  testDotProduct<float>();
  testDotProduct<double>();
}
