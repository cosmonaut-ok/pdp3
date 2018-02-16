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

  Geometry * geometry = new Geometry(r_size, z_size, n_grid_r, n_grid_z);
  Geometry * geometry_default = new Geometry(r_size, z_size, n_grid_r, n_grid_z);

  TEST(Geometry, init)
  {
    ASSERT_EQ (geometry->first_size, r_size);
    ASSERT_EQ (geometry->second_size, z_size);
    ASSERT_EQ (geometry->n_grid_1, n_grid_r);
    ASSERT_EQ (geometry->n_grid_2, n_grid_z);
  }

  TEST(Geometry, set_epsilon)
  {
    geometry->set_epsilon();
    for(int i=0;i<(n_grid_r);i++)
      for(int k=0;k<(n_grid_z);k++)
      {
        ASSERT_EQ (geometry->epsilon[i][k], 1);
        ASSERT_NE (geometry_default->epsilon[i][k], 1);
      }
  }

  TEST(Geometry, set_pml)
  {
    geometry->set_pml(comp_l_1, comp_l_2, comp_l_3, sigma_1_t, sigma_2_t);

     for(int i=0; i < n_grid_r; i++)
      for(int k=0; k < n_grid_z; k++)
        ASSERT_NE (geometry->sigma[i][k], geometry_default->sigma[i][k]);
  }

  TEST(Geometry, get_dr)
  {
    double r_ratio = r_size/(n_grid_r-1);
    ASSERT_EQ (geometry->get_dr(), r_ratio);
  }

  TEST(Geometry, get_dz)
  {
    double z_ratio = z_size/(n_grid_z-1);
    ASSERT_EQ (geometry->get_dz(), z_ratio);
  }
}
