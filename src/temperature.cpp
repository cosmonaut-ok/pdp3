#include "temperature.h"

//// constructor
Temperature::Temperature(Geometry *geom1):geom(geom1)
{
  t = new double *[geom->n_grid_r];
  t_sum = new double *[geom->n_grid_r];
  count = new double *[geom->n_grid_r];

  // filling second demension
  // for potential
  for (int i = 0; i < (geom->n_grid_r); i++)
  {
    t[i] = new double[geom->n_grid_z];
    t_sum[i] = new double[geom->n_grid_z];
    count[i] = new double[geom->n_grid_z];
  }
}

void Temperature::reset(void)
{
// #pragma omp parallel for
  for (unsigned int i = 0; i < geom->n_grid_r; i++)
    for (unsigned int j = 0; j < geom->n_grid_z; j++)
    {
      t[i][j] = 0;
      t_sum[i][j] = 0;
      count[i][j] = 0;
    }
}

void Temperature::calc_t_r(Particles *prtls)
{
  reset();

  for (unsigned int i = 0; i < prtls->number; i++)
  {
    unsigned int i_r = CELL_NUMBER(prtls->pos[i][0], geom->dr);
    unsigned int k_z = CELL_NUMBER(prtls->pos[i][2], geom->dz);

    t_sum[i_r][k_z] = t_sum[i_r][k_z] + VEL_TO_TEMPR(abs(prtls->vel[i][0]), prtls->mass);
    count[i_r][k_z]++;
  }

// #pragma omp parallel for
  for (unsigned int i = 0; i < geom->n_grid_r; i++)
    for (unsigned int j = 0; j < geom->n_grid_z; j++)
      {
	if (count[i][j] == 0) count[i][j] = 1;
	t[i][j] = t_sum[i][j] / count[i][j];
      }
}

void Temperature::calc_t_phi(Particles *prtls)
{
  reset();

  for (unsigned int i = 0; i < prtls->number; i++)
  {
    unsigned int i_r = CELL_NUMBER(prtls->pos[i][0], geom->dr);
    unsigned int k_z = CELL_NUMBER(prtls->pos[i][2], geom->dz);

    t_sum[i_r][k_z] = t_sum[i_r][k_z] + VEL_TO_TEMPR(abs(prtls->vel[i][1]), prtls->mass);
    count[i_r][k_z]++;
  }

// #pragma omp parallel for
  for (unsigned int i = 0; i < geom->n_grid_r; i++)
    for (unsigned int j = 0; j < geom->n_grid_z; j++)
      {
	if (count[i][j] == 0) count[i][j] = 1;
	t[i][j] = t_sum[i][j] / count[i][j];
      }
}

void Temperature::calc_t_z(Particles *prtls)
{
  reset();

  for (unsigned int i = 0; i < prtls->number; i++)
  {
    unsigned int i_r = CELL_NUMBER(prtls->pos[i][0], geom->dr);
    unsigned int k_z = CELL_NUMBER(prtls->pos[i][2], geom->dz);

    t_sum[i_r][k_z] = t_sum[i_r][k_z] + VEL_TO_TEMPR(abs(prtls->vel[i][2]), prtls->mass);
    count[i_r][k_z]++;
  }

// #pragma omp parallel for
  for (unsigned int i = 0; i < geom->n_grid_r; i++)
    for (unsigned int j = 0; j < geom->n_grid_z; j++)
      {
	if (count[i][j] == 0) count[i][j] = 1;
	t[i][j] = t_sum[i][j] / count[i][j];
      }
}

void Temperature::calc_t(Particles *prtls)
{
  reset();

  for (unsigned int i = 0; i < prtls->number; i++)
  {
    unsigned int i_r = CELL_NUMBER(prtls->pos[i][0], geom->dr);
    unsigned int k_z = CELL_NUMBER(prtls->pos[i][2], geom->dz);

    double t_r_cmp = VEL_TO_TEMPR(prtls->vel[i][0], prtls->mass);
    double t_phi_cmp = VEL_TO_TEMPR(prtls->vel[i][1], prtls->mass);
    double t_z_cmp = VEL_TO_TEMPR(prtls->vel[i][2], prtls->mass);

    t_sum[i_r][k_z] = t_sum[i_r][k_z] + lib::sq_rt(t_r_cmp * t_r_cmp + t_phi_cmp * t_phi_cmp + t_z_cmp * t_z_cmp);
    count[i_r][k_z]++;
  }

// #pragma omp parallel for
  for (unsigned int i = 0; i < geom->n_grid_r; i++)
    for (unsigned int j = 0; j < geom->n_grid_z; j++)
      {
	if (count[i][j] == 0) count[i][j] = 1;
	t[i][j] = t_sum[i][j] / count[i][j];
      }
}
