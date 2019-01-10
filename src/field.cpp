#include "field.h"

using namespace std;

Field::Field()
{
}

//// constructor
Field::Field(Geometry *geom1_t): geom1(geom1_t)
{
  // n_grid - number of edges
  field_r = new double *[geom1->n_grid_r];
  field_phi = new double *[geom1->n_grid_r];
  field_z = new double *[geom1->n_grid_r];

#pragma omp parallel for
  // filling second demension
  for (int i=0; i<(geom1->n_grid_r); i++)
  {
    field_r[i] = new double [geom1->n_grid_z];
    field_phi[i] = new double [geom1->n_grid_z];
    field_z[i] = new double [geom1->n_grid_z];
  }
}

// destructor
Field::~Field()
{
  for (int i = 0; i < (geom1->n_grid_r - 1); i++)
  {
    delete[]field_r[i];
    delete[]field_phi[i];
    delete[]field_z[i];
  }

  delete[]field_r;
  delete[]field_phi;
  delete[]field_z;
}
