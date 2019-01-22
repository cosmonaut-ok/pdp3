#include "temperature.h"


//// constructor
Temperature::Temperature(Geometry *geom1):geom(geom1)
{
  t_r = new double *[geom->n_grid_r];
  t_phi = new double *[geom->n_grid_r];
  t_z = new double *[geom->n_grid_r];

  t_r_sum = new double *[geom->n_grid_r];
  t_phi_sum = new double *[geom->n_grid_r];
  t_z_sum = new double *[geom->n_grid_r];

  count = new double *[geom->n_grid_r];

  // filling second demension
  // for potential
  for (int i = 0; i < (geom->n_grid_r); i++)
  {
    t_r[i] = new double[geom->n_grid_z];
    t_phi[i] = new double[geom->n_grid_z];
    t_z[i] = new double[geom->n_grid_z];

    t_r_sum[i] = new double[geom->n_grid_z];
    t_phi_sum[i] = new double[geom->n_grid_z];
    t_z_sum[i] = new double[geom->n_grid_z];

    count[i] = new double[geom->n_grid_z];
  }
}

void Temperature::calc_t_r(Particles *prtls)
{
  for (unsigned int i = 0; i < prtls->number; i++)
  {
    unsigned int i_r = CELL_NUMBER(prtls->pos[i][0], geom->dr);
    unsigned int k_z = CELL_NUMBER(prtls->pos[i][2], geom->dz);


    t_r_sum[i_r][k_z] = t_r_sum[i_r][k_z] + prtls->vel[i][0];
    count[i_r][k_z] = count[i_r][k_z]++;
  }

// #pragma omp parallel for
  for (unsigned int i = 0; i < geom->n_grid_r; i++)
    for (unsigned int j = 0; j < geom->n_grid_z; j++)
      t_r[i][j] = t_r_sum[i][j] / count[i][j];
}

void Temperature::calc_t_phi(Particles *prtls)
{
  for (unsigned int i = 0; i < prtls->number; i++)
  {
    unsigned int i_r = CELL_NUMBER(prtls->pos[i][0], geom->dr);
    unsigned int k_z = CELL_NUMBER(prtls->pos[i][2], geom->dz);

    t_phi_sum[i_r][k_z] = t_phi_sum[i_r][k_z] + prtls->vel[i][1];
    count[i_r][k_z] = count[i_r][k_z]++;
  }

// #pragma omp parallel for
  for (unsigned int i = 0; i < geom->n_grid_r; i++)
    for (unsigned int j = 0; j < geom->n_grid_z; j++)
      t_phi[i][j] = t_phi_sum[i][j] / count[i][j];
}

void Temperature::calc_t_z(Particles *prtls)
{
  for (unsigned int i = 0; i < prtls->number; i++)
  {
    unsigned int i_r = CELL_NUMBER(prtls->pos[i][0], geom->dr);
    unsigned int k_z = CELL_NUMBER(prtls->pos[i][2], geom->dz);

    t_z_sum[i_r][k_z] = t_z_sum[i_r][k_z] + prtls->vel[i][2];
    count[i_r][k_z] = count[i_r][k_z]++;
  }

// #pragma omp parallel for
  for (unsigned int i = 0; i < geom->n_grid_r; i++)
    for (unsigned int j = 0; j < geom->n_grid_z; j++)
      t_z[i][j] = t_z_sum[i][j] / count[i][j];
}



void Temperature::reset(void)
{
// #pragma omp parallel for
  for (unsigned int i = 0; i < geom->n_grid_r; i++)
    for (unsigned int j = 0; j < geom->n_grid_z; j++)
    {
      t_r[i][j] = 0;
      t_phi[i][j] = 0;
      t_z[i][j] = 0;
      t_r_sum[i][j] = 0;
      t_phi_sum[i][j] = 0;
      t_z_sum[i][j] = 0;
      count[i][j] = 0;
    }
}
