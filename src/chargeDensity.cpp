#include "chargeDensity.h"

ChargeDensity::ChargeDensity(void)
{
}

// constructor
ChargeDensity::ChargeDensity(Geometry *geom1_t): geom1(geom1_t)
{
  rho = new double *[geom1->n_grid_r];

#pragma omp parallel
  {
#pragma omp for
    for (int i = 0; i < (geom1->n_grid_r); i++)
      rho[i] = new double[geom1->n_grid_z];

    // initialization
#pragma omp for
    for (int i = 0; i < (geom1->n_grid_r); i++)
      for (int k = 0; k < (geom1->n_grid_z); k++)
        rho[i][k] = 0;
  }
}

ChargeDensity::~ChargeDensity(void)
{
  for (int i = 0; i < (geom1->n_grid_r); i++)
    delete[]rho[i];
  delete[]rho;
}

double **ChargeDensity::get_rho() const
{
  return rho;
}

void ChargeDensity::inc_rho(int i, int k, double value)
{
  rho[i][k] += value;
}

void ChargeDensity::reset_rho()
{
#pragma omp parallel for
  for (int i = 0; i < geom1->n_grid_r; i++)
    for (int k = 0; k < geom1->n_grid_z; k++)
      rho[i][k] = 0.;
}
