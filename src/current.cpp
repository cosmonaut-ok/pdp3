#include "current.h"

Current::Current(void)
{
}

Current::Current(Geometry *geom1_t): geom1(geom1_t)
{
  j_r = new double *[geom1->n_grid_r];
  j_phi = new double *[geom1->n_grid_r];
  j_z = new double *[geom1->n_grid_r];

  j_r_1d = new double[geom1->n_grid_r * geom1->n_grid_z];
  j_phi_1d = new double[geom1->n_grid_r * geom1->n_grid_z];
  j_z_1d = new double[geom1->n_grid_r * geom1->n_grid_z];

#pragma omp parallel
  {
    // jr
#pragma omp for
    for (int i = 0; i < (geom1->n_grid_r); i++)
      j_r[i]= new double[geom1->n_grid_z];

    // jfi
#pragma omp for
    for (int i = 0; i < (geom1->n_grid_r); i++)
      j_phi[i] = new double[geom1->n_grid_z];

    // jz
#pragma omp for
    for (int i = 0; i < (geom1->n_grid_r); i++)
      j_z[i] = new double[geom1->n_grid_z];
  }

  // initialization
#pragma omp parallel
  {
#pragma omp for
    for (int i = 0; i < (geom1->n_grid_r-1); i++)
      for (int k = 0; k < (geom1->n_grid_z-1); k++)
      {
        j_r[i][k]=0.;
        j_phi[i][k]=0.;
        j_z[i][k]=0.;
      }

#pragma omp for
    for (int i = 0; i < (geom1->n_grid_r-1); i++)
      j_r[i][geom1->n_grid_z-1] = 0.;

#pragma omp for
    for(int k = 0; k < (geom1->n_grid_z); k++)
      j_phi[geom1->n_grid_r-1][k] = 0.;

#pragma omp for
    for(int i = 0; i < (geom1->n_grid_r); i++)
      j_phi[i][geom1->n_grid_z-1] = 0.;

#pragma omp for
    for(int k = 0; k < (geom1->n_grid_z-1); k++)
      j_z[geom1->n_grid_r-1][k] = 0.;
  }
}

// Destructor
Current::~Current(void)
{
  for (int i = 0; i < (geom1->n_grid_r-1); i++)
    delete[]j_r[i];
  delete[]j_r;

  for (int i = 0; i < (geom1->n_grid_r); i++)
    delete[]j_phi[i];
  delete[]j_phi;

  for (int i = 0; i < (geom1->n_grid_r-1); i++)
    delete[]j_z[i];
  delete[]j_z;

  delete j_r_1d;
  delete j_phi_1d;
  delete j_z_1d;
}

//
// functions for getting access to j arrays
//

double **Current::get_j_r() const
{
  return this->j_r;
}

double **Current::get_j_phi() const
{
  return this->j_phi;
}

double **Current::get_j_z() const
{
  return this->j_z;
}

//
// functions for changing values of j
//
void Current::inc_j_r(int i, int k, double value)
{
#pragma omp critical
  j_r[i][k] += value;
}

void Current::inc_j_phi(int i, int k, double value)
{
#pragma omp critical
  j_phi[i][k] += value;
}

void Current::inc_j_z(int i, int k, double value)
{
#pragma omp critical
  j_z[i][k] += value;
}

void Current::reset_j()
{
#pragma omp parallel for
  for (int i = 0; i < (geom1->n_grid_r); i++)
    for (int k = 0; k < (geom1->n_grid_z); k++)
    {
      j_r[i][k] = 0;
      j_phi[i][k] = 0;
      j_z[i][k] = 0;
    }
}
