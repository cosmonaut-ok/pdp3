#include "charge_density.h"

charge_density::charge_density(void)
{
}

// constructor
charge_density::charge_density(Geometry* geom1_t): geom1(geom1_t)
{
  ro = new double*[geom1->n_grid_1];
  for (int i=0; i<(geom1->n_grid_1);i++)
  {
    ro[i]= new double[geom1->n_grid_2];
  }

  // initialization
  for (int i=0; i<(geom1->n_grid_1);i++)
    for (int k=0; k<(geom1->n_grid_2);k++)
    {
      ro[i][k]=0;
    }
}

charge_density::~charge_density(void)
{
  for (int i=0; i<(geom1->n_grid_1);i++)
    delete[]ro[i];
  delete[]ro;
}

double** charge_density::get_ro() const
{
  return ro;
}

void charge_density::set_ro_weighting(int i, int k, double value)
{
  ro[i][k]=ro[i][k]+value;
}

void charge_density::reset_rho()
{
#pragma omp parallel for shared (geom1, ro)
  for (int i=0;i<geom1->n_grid_1;i++)
    for (int k=0;k<geom1->n_grid_2;k++)
      ro[i][k] = 0.0;

}
