#include <gtest/gtest.h>
#include "Geometry.h"
#include "Geometry.cpp"

namespace geometry {
  double comp_l_1 = 0.0;
  double comp_l_2 = 0.0;
  double comp_l_3 = 1.1;
  double sigma_1_t = 1e-5;
  double sigma_2_t = 7e-2;

  double r_size = 0.1;
  double z_size = 0.4;
  double n_grid_r = 16;
  double n_grid_z = 64;

  Geometry * geometry = new Geometry(r_size, z_size, n_grid_r, n_grid_z)

  TEST(Geometry, init)
  {
    ASSERT_EQ (geometry->first_size, r_size);
    ASSERT_EQ (geometry->second_size, z_size);
    ASSERT_EQ (geometry->n_grid_1, n_grid_r);
    ASSERT_EQ (geometry->n_grid_2, n_grid_z);
  }
}
