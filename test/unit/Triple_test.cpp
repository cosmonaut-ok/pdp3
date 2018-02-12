#include <gtest/gtest.h>
#include "Triple.h"
#include "Triple.cpp"

double first_component = 0.1;
double second_component = 1.1;
double third_component = 1.0;

Triple test_triple(first_component, second_component, third_component);

TEST(Triple, init)
{
  EXPECT_TRUE(test_triple.first == first_component);
  EXPECT_TRUE(test_triple.second == second_component);
  EXPECT_TRUE(test_triple.third == third_component);
}
