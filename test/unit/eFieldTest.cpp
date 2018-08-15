#include <gtest/gtest.h>
#pragma once
#include "field.h"
// #include "tinyvec3d.h"
// #include "tinyvec3d.cpp"
#include "eField.h"
#include "eField.cpp"
#include "geometry.h"

namespace e_field {

  double r_size = 0.1;
  double z_size = 0.4;
  double n_grid_r = 16;
  double n_grid_z = 64;

  Geometry * geometry = new Geometry(r_size, z_size, n_grid_r, n_grid_z);

  EField * efield = new EField(geometry);

  TEST(E_field, object_created)
  {
    EXPECT_TRUE(efield);
  }
}
