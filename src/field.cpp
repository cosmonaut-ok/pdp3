#include "field.h"

using namespace std;

Field::Field()
{
}

//// constructor
Field::Field(Geometry *geom1_t): geom1(geom1_t)
{
  // n_grid - number of edges
  //// Field_r
  field_r = new double *[geom1->n_grid_1];
  field_r_1d = new double [(geom1->n_grid_1)*geom1->n_grid_2];
  //// Field_phi
  field_phi = new double *[geom1->n_grid_1];
  field_phi_1d = new double [geom1->n_grid_1*geom1->n_grid_2];
  //// Field_z
  field_z = new double *[geom1->n_grid_1];
  field_z_1d = new double [geom1->n_grid_1*(geom1->n_grid_2)];

#pragma omp parallel for
  // filling second demension
  for (int i=0; i<(geom1->n_grid_1); i++)
  {
    field_r[i]= new double [geom1->n_grid_2];
    field_phi[i]= new double [geom1->n_grid_2];
    field_z[i]= new double [geom1->n_grid_2];
  }
}

// destructor
Field::~Field()
{
  for (int i=0; i<(geom1->n_grid_1-1);i++)
  {
    delete[]field_r[i];
    delete[]field_phi[i];
    delete[]field_z[i];
  }
  delete[]field_r;
  delete[]field_phi;
  delete[]field_z;
}

// Return one dimensional field components
double *Field::get_1d_field_r()
{
  // copy 2d field array into 1d array rowwise
#pragma omp parallel for
  for (int i = 0; i < geom1->n_grid_1 - 1; i++)
    for (int k = 0; k < geom1->n_grid_2; k++)
      field_r_1d[i * geom1->n_grid_2 + k] = field_r[i][k];
  return field_r_1d;
}

double *Field::get_1d_field_phi()
{
  // copy 2d field array into 1d array rowwise
#pragma omp parallel for
  for (int i = 0; i < geom1->n_grid_1; i++)
    for (int k = 0; k < geom1->n_grid_2; k++)
      field_phi_1d[i * geom1->n_grid_2 + k] = field_phi[i][k];
  return field_phi_1d;
}

double *Field::get_1d_field_z()
{
  // copy 2d field array into 1d array rowwise
#pragma omp parallel for
  for (int i = 0; i < geom1->n_grid_1; i++)
    for (int k = 0; k < geom1->n_grid_2 - 1; k++)
      field_z_1d[i * (geom1->n_grid_2 - 1) + k] = field_z[i][k];
  return field_z_1d;
}
