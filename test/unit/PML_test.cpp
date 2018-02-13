#include <gtest/gtest.h>
#include "PML.h"
#include "PML.cpp"
#include "Geometry.h"

// positive case
double comp_l_1 = 0.0;
double comp_l_2 = 0.0;
double comp_l_3 = 1.1;
double sigma_1_t = 1e-5;
double sigma_2_t = 7e-2;

double r_size = 0.1;
double z_size = 0.4;
double n_grid_r = 16;
double n_grid_z = 64;

PML * pml_default = new PML(comp_l_1, comp_l_2, comp_l_3, sigma_1_t, sigma_2_t);

Geometry * geometry = new Geometry(r_size, z_size, n_grid_r, n_grid_z, pml_default);
Geometry * geometry_zero = new Geometry(r_size, z_size, n_grid_r, n_grid_z);

TEST(PML, init)
{  
  ASSERT_EQ(pml_default->comparative_l_1, comp_l_1);
  ASSERT_EQ(pml_default->comparative_l_2, comp_l_2);
  ASSERT_EQ(pml_default->comparative_l_3, comp_l_3);
  ASSERT_EQ(pml_default->sigma1, sigma_1_t);
  ASSERT_EQ(pml_default->sigma2, sigma_2_t);
}

TEST(PML, calc_sigma)
{
  // calc_sigma used directly in geometry constructor,
  // so can not be separated from it
  // Comparing already initialized geometry objects
  // TODO: make more advanced testing
  ASSERT_NE(geometry->sigma, geometry_zero->sigma);
}
