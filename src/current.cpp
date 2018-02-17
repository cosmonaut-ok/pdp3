#include "current.h"

Current::Current(void)
{
}

Current::Current(Geometry* geom1_t): geom1(geom1_t)
{
  //// jr
  // n_grid - number of edges
  j1 = new double*[geom1->n_grid_1-1];
#pragma omp parallel for shared (j1)
  for (int i=0; i<(geom1->n_grid_1-1);i++)
    j1[i]= new double[geom1->n_grid_2];

  //// jfi
  j2 = new double*[geom1->n_grid_1];
#pragma omp parallel for shared (j2)
  for (int i=0; i<(geom1->n_grid_1);i++)
    j2[i]= new double[geom1->n_grid_2];

  //// jz
  j3 = new double*[geom1->n_grid_1];
#pragma omp parallel for shared (j3)
  for (int i=0; i<(geom1->n_grid_1);i++)
    j3[i]= new double[geom1->n_grid_2-1];

  j1_1d = new double[(geom1->n_grid_1-1)*geom1->n_grid_2];

  j2_1d = new double[geom1->n_grid_1*geom1->n_grid_2];

  j3_1d = new double[geom1->n_grid_1*(geom1->n_grid_2-1)];

  // initialization
#pragma omp parallel for shared (j1, j2, j3)
  for (int i=0; i<(geom1->n_grid_1-1);i++)
    for (int k=0; k<(geom1->n_grid_2-1);k++)
    {
      j1[i][k]=0;
      j2[i][k]=0;
      j3[i][k]=0;
    }
#pragma omp parallel for shared (j1)
  for (int i=0; i<(geom1->n_grid_1-1);i++)
    j1[i][geom1->n_grid_2-1]=0;

#pragma omp parallel for shared (j2)
  for(int k=0;k<(geom1->n_grid_2);k++)
    j2[geom1->n_grid_1-1][k]=0;

#pragma omp parallel for shared (j2)
  for(int i=0;i<(geom1->n_grid_1);i++)
    j2[i][geom1->n_grid_2-1]=0;

#pragma omp parallel for shared (j3)
  for(int k=0;k<(geom1->n_grid_2-1);k++)
    j3[geom1->n_grid_1-1][k]=0;
}

// Destructor
Current::~Current(void)
{
  for (int i=0; i<(geom1->n_grid_1-1);i++)
    delete[]j1[i];
  delete[]j1;

  for (int i=0; i<(geom1->n_grid_1);i++)
    delete[]j2[i];
  delete[]j2;

  for (int i=0; i<(geom1->n_grid_1-1);i++)
    delete[]j3[i];
  delete[]j3;

  delete j1_1d;
  delete j2_1d;
  delete j3_1d;
}

//
// functions for getting access to j arrays
//

double** Current::get_j1() const
{
  return this->j1;
}

double** Current::get_j2() const
{
  return this->j2;
}

double** Current::get_j3() const
{
  return this->j3;
}

// Return one dimensional field components
double* Current::get_j1_1d() const
{
  // copy 2d field array into 1d array rowwise
#pragma omp parallel for
  for (int i = 0; i < geom1->n_grid_1 - 1; i++)
    for (int k = 0; k < geom1->n_grid_2; k++)
      j1_1d[i * geom1->n_grid_2 + k] = j1[i][k];
  return j1_1d;
}

double* Current::get_j2_1d() const
{
  // copy 2d field array into 1d array rowwise
#pragma omp parallel for
  for (int i = 0; i < geom1->n_grid_1; i++)
    for (int k = 0; k < geom1->n_grid_2; k++)
      j2_1d[i * geom1->n_grid_2 + k] = j2[i][k];
  return j2_1d;
}

double* Current::get_j3_1d() const
{
  // copy 2d field array into 1d array rowwise
#pragma omp parallel for
  for (int i = 0; i < geom1->n_grid_1; i++)
    for (int k = 0; k < geom1->n_grid_2 - 1; k++)
      j3_1d[i * (geom1->n_grid_2 - 1) + k] = j3[i][k];
  return j3_1d;
}

void Current::j1_add_1d(double *j1_1d)
{
  // copy 2d field array into 1d array rowwise
#pragma omp parallel for
  for (int i = 0; i < geom1->n_grid_1 - 1; i++)
    for (int k = 0; k < geom1->n_grid_2; k++)
      j1[i][k] += j1_1d[i * geom1->n_grid_2 + k];
}

void Current::j2_add_1d(double *j2_1d)
{
  // copy 2d field array into 1d array rowwise
#pragma omp parallel for
  for (int i = 0; i < geom1->n_grid_1; i++)
    for (int k = 0; k < geom1->n_grid_2; k++)
      j2[i][k] += j2_1d[i * geom1->n_grid_2 + k];
}

void Current::j3_add_1d(double *j3_1d)
{
  // copy 2d field array into 1d array rowwise
#pragma omp parallel for
  for (int i = 0; i < geom1->n_grid_1; i++)
    for (int k = 0; k < geom1->n_grid_2 - 1; k++)
      j3[i][k] += j3_1d[i * (geom1->n_grid_2 - 1) + k];
}

//
// functions for changing values of j
//
void Current::set_j1(int i, int k, double value)
{
  j1[i][k] = j1[i][k] + value;
}

void Current::set_j2(int i, int k, double value)
{
  j2[i][k] = j2[i][k] + value;
}

void Current::set_j3(int i, int k, double value)
{
  j3[i][k]=j3[i][k] + value;
}

void Current::reset_j()
{
#pragma omp parallel
  {
#pragma omp for
    for (int i=0; i<(geom1->n_grid_1-1); i++)
      for (int k=0; k<(geom1->n_grid_2-1); k++)
      {
        j1[i][k]=0;
        j3[i][k]=0;
      }
#pragma omp for
    for (int i=0; i<geom1->n_grid_1; i++)
      for (int k=0; k<geom1->n_grid_2; k++)
        j2[i][k] = 0.0;
#pragma omp for
    for (int i=0; i<(geom1->n_grid_1-1); i++)
      j1[i][geom1->n_grid_2-1]=0;

#pragma omp for
    for(int i=0; i<(geom1->n_grid_2-1); i++)
      j3[geom1->n_grid_1-1][i]=0;
  }
}
