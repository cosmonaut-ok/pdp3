#include <gtest/gtest.h>
#include "Constant.h"

TEST(Universe, isOk)
{
  EXPECT_TRUE(constant::PI == 3.1415926535897932);
  EXPECT_TRUE(constant::EPSILON0 == 8.85E-12);
  EXPECT_TRUE(constant::EL_MASS == 9.1E-31);
  EXPECT_TRUE(constant::EL_CHARGE == 1.6E-19);
  EXPECT_TRUE(constant::LIGHT_SPEED == 3.0E8);
  EXPECT_TRUE(constant::MAGN_CONST == 1.26E-6);
}

int main(int ac, char* av[])
{
  testing::InitGoogleTest(&ac, av);
  return RUN_ALL_TESTS();
}
