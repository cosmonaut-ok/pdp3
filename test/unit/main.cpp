#include <gtest/gtest.h>

#pragma once
#include "tinyvec3d.cpp"

int main(int ac, char* av[])
{
  testing::InitGoogleTest(&ac, av);
  return RUN_ALL_TESTS();
}
