#include <gtest/gtest.h>
#include "Geometry.h"
#include "Particles.h"
#include "Particles.cpp"

namespace particles {
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

	char* p_name = (char*)"Mockytron";
	double charge = 1;
	double mass = 1;
	double number = 1e6;
	double left_density = 2e11;
	double right_density = 2.01e11;
	double temperature = 10;

	Particles * prtls = new Particles(p_name, charge, mass, number, geometry);
	Particles * prtls_clean = new Particles(p_name, charge, mass, number, geometry);

	TEST(Particles, init)
	{
		EXPECT_TRUE(prtls);
		ASSERT_EQ(prtls->name, p_name);
		ASSERT_EQ(prtls->charge, charge*constant::EL_CHARGE);
		ASSERT_EQ(prtls->mass, mass*constant::EL_MASS);
		ASSERT_EQ(prtls->number, number);

		for (int i = 0; i < number; i++)
			prtls->is_alive[i] = 1;
	}

	TEST(Particles, load_spatial_distribution_with_variable_mass)
	{
		prtls->load_spatial_distribution_with_variable_mass(left_density, right_density, 0,0);
		for (int i = 0; i<number;i++)
		{
			ASSERT_NE(prtls->x1[i], prtls_clean->x1[i]);
			ASSERT_NE(prtls->x3[i], prtls_clean->x3[i]);
		}
	}

	TEST(Particles, velocity_distribution)
	{
		prtls->velocity_distribution(temperature);
	}

	TEST(Particles, get_gamma)
	{
		double gamma = prtls->get_gamma(101);
		double gamma_preset = 1.0000334481148347; // just result of small experiment :-P
		ASSERT_EQ(gamma, gamma_preset);
	}

	TEST(Particles, get_gamma_inv)
	{
		double gamma = prtls->get_gamma_inv(101);
		double gamma_preset = 1.0000334458774316; // just result of small experiment :-P
		ASSERT_EQ(gamma, gamma_preset);
	}

}

// set_v_0
// set_x_0
// step_v
// half_step_coord
// charge_weighting
// velocity_distribution
// random_reverse
// load_spatial_distribution
// simple_j_weighting
// simple_constrho_j_weighting
// j_weighting
// azimuthal_j_weighting
// strict_motion_weighting
