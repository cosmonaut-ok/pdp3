#include "charge_density.h"

charge_density::charge_density(void)
{
}

// constructor
charge_density::charge_density(Geometry* geom1_t): geom1(geom1_t)
{
  rho = new double*[geom1->n_grid_1];
#pragma omp parallel
  {
#pragma omp for
    for (int i=0; i<(geom1->n_grid_1);i++)
      rho[i]= new double[geom1->n_grid_2];

    // initialization
#pragma omp for
    for (int i=0; i<(geom1->n_grid_1);i++)
      for (int k=0; k<(geom1->n_grid_2);k++)
        rho[i][k]=0;
  }
}
charge_density::~charge_density(void)
{
  for (int i=0; i<(geom1->n_grid_1);i++)
    delete[]rho[i];
  delete[]rho;
}

double** charge_density::get_ro() const
{
  return rho;
}

void charge_density::set_ro_weighting(int i, int k, double value)
{
  rho[i][k]=rho[i][k]+value;
}

void charge_density::reset_rho()
{
#pragma omp parallel for
  for (int i=0;i<geom1->n_grid_1;i++)
    for (int k=0;k<geom1->n_grid_2;k++)
      rho[i][k] = 0.0;
}
