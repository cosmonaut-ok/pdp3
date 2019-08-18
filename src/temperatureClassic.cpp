#include "temperatureClassic.h"

// use classical calculations, if velocity lower, than minimal
#ifndef REL_LIMIT
#define REL_LIMIT 5e7
#endif
// define REL_LIMIT^2 to decrease number of operations
const double REL_LIMIT_POW_2 = pow (REL_LIMIT, 2);

// C^2 define c^2 to decrease number of operations
const double LIGHT_SPEED_POW_2 = pow (LIGHT_SPEED, 2);

//// constructor
TemperatureClassic::TemperatureClassic(Geometry *geom1):geom(geom1)
{
  t_src = new double *[geom->n_grid_r - 1];
  t = new double *[geom->n_grid_r - 1];
  t_sum = new double *[geom->n_grid_r];
  t_vec_r = new double *[geom->n_grid_r];
  t_vec_phi = new double *[geom->n_grid_r];
  t_vec_z = new double *[geom->n_grid_r];
  count = new double *[geom->n_grid_r];
  count_sum = new double *[geom->n_grid_r];

  // filling second demension
  // for potential
  for (int i = 0; i < (geom->n_grid_r); i++)
  {
    t_src[i] = new double[geom->n_grid_z - 1];
    t[i] = new double[geom->n_grid_z - 1];
    t_sum[i] = new double[geom->n_grid_z];
    t_vec_r[i] = new double[geom->n_grid_z];
    t_vec_phi[i] = new double[geom->n_grid_z];
    t_vec_z[i] = new double[geom->n_grid_z];
    count[i] = new double[geom->n_grid_z];
    count_sum[i] = new double[geom->n_grid_z];
  }
}

void TemperatureClassic::reset(void)
{
// #pragma omp parallel for
  for (unsigned int i = 0; i < geom->n_grid_r; i++)
    for (unsigned int j = 0; j < geom->n_grid_z; j++)
    {
      t_src[i][j] = 0;
      t[i][j] = 0;
      t_sum[i][j] = 0;
      t_vec_r[i][j] = 0;
      t_vec_phi[i][j] = 0;
      t_vec_z[i][j] = 0;
      count[i][j] = 0;
      count_sum[i][j] = 0;
    }
}

void TemperatureClassic::calc_t(Particles *prtls)
{
  reset();

  double dr = geom->dr;
  double dz = geom->dz;
  double r1, r2, r3; // temp variables for calculation
  double dz1, dz2; // temp var.: width of k and k + 1 cell

  double ro_v = 0; // charge density Q/V, V - volume of particle
  double ro_v_r = 0;
  double ro_v_phi = 0;
  double ro_v_z = 0;
  double ro_p = 0; // charge density Q/V, V - volume of particle
  double v_0 = 0; // charge density Q/V, V - volume of particle
  double v_1 = 0; // volume of [i][k] cell
  double v_2 = 0; // volume of [i + 1][k] cell
  double vel = 0;
  double vel_r = 0;
  double vel_phi = 0;
  double vel_z = 0;

  double value = 0;
  // double **temp = ro1->get_rho();

  for(unsigned int i = 0; i < prtls->number; i++)
  {
    // finding number of i and k cell. example: dr = 0.5; r = 0.4; i =0
    unsigned int r_i = CELL_NUMBER(prtls->pos[i][0], geom->dr);
    unsigned int z_k = CELL_NUMBER(prtls->pos[i][2], geom->dz);

    vel_r = prtls->vel[i][0];
    vel_phi = prtls->vel[i][1];
    vel_z = prtls->vel[i][2];

    vel = lib::sq_rt(vel_r * vel_r + vel_phi * vel_phi + vel_z * vel_z);

    inc_tmpr(r_i, z_k, vel);
    inc_vec_r(r_i, z_k, vel_r);
    inc_vec_phi(r_i, z_k, vel_phi);
    inc_vec_z(r_i, z_k, vel_z);
    inc_count(r_i, z_k, 1);

  }
}

void TemperatureClassic::inc_tmpr(unsigned int r, unsigned int z, double value)
{
  t_sum[r][z] += value;
}

void TemperatureClassic::inc_vec_r(unsigned int r, unsigned int z, double value)
{
  t_vec_r[r][z] += value;
}

void TemperatureClassic::inc_vec_phi(unsigned int r, unsigned int z, double value)
{
  t_vec_phi[r][z] += value;
}

void TemperatureClassic::inc_vec_z(unsigned int r, unsigned int z, double value)
{
  t_vec_z[r][z] += value;
}

void TemperatureClassic::inc_count(unsigned int r, unsigned int z, double value)
{
  count_sum[r][z] += value;
}

void TemperatureClassic::normalize(Particles *prtls)
{
  for (unsigned int r = 0; r < geom->n_grid_r - 1; r++)
    for (unsigned int z = 0; z < geom->n_grid_z - 1; z++)
    {
      double v_vec_sum_2 = t_vec_r[r][z] * t_vec_r[r][z]
        + t_vec_phi[r][z] * t_vec_phi[r][z] + t_vec_z[r][z] * t_vec_z[r][z];

      double v_sc_sum_2 = (t_sum[r][z] * t_sum[r][z] - v_vec_sum_2) / (count_sum[r][z] * count_sum[r][z]);

      if (v_sc_sum_2 < REL_LIMIT * REL_LIMIT)
	t_src[r][z] = prtls->mass * v_sc_sum_2 / 2;
      else
	{
	  double gamma = lib::get_gamma(v_sc_sum_2);
	  t_src[r][z] = prtls->mass * LIGHT_SPEED_POW_2 * (gamma - 1);
	}
      

      // convert joules to electronvolts
      t_src[r][z] /= abs(prtls->charge);
    }

#if defined TEMPERATURE_POSTPROC_BILINEAR
  lib::bilinear_interpolation(t_src, t, geom->n_grid_r, geom->n_grid_z);
#elif defined TEMPERATURE_POSTPROC_BICUBIC
  lib::bicubic_interpolation(t_src, t, geom->n_grid_r, geom->n_grid_z);
#else
  t = t_src;
#endif

}
