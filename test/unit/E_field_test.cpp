#include <gtest/gtest.h>
#include "E_field.h"
#include "E_field.cpp"
#include "Geometry.h"

namespace e_field {

  double r_size = 0.1;
  double z_size = 0.4;
  double n_grid_r = 16;
  double n_grid_z = 64;

  Geometry * geometry = new Geometry(r_size, z_size, n_grid_r, n_grid_z);

  E_field * efield = new E_field(geometry);

  TEST(E_field, DUMMY)
  {
    EXPECT_TRUE(2 == 2);
  }
}
