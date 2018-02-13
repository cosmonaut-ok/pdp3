#include "Geometry.h"
#include "PML.h"

//constructor with sigma
Geometry::Geometry(double fs, double ss, int ng1, int ng2, PML* pml1_t): pml1(pml1_t)
{
	init_geometry (fs, ss, ng1, ng2);
  // sigma assigning
  pml1->calc_sigma(this);
}

// constructor: sigma = 0
Geometry::Geometry(double fs, double ss,  int ng1, int ng2)
{
	init_geometry (fs, ss, ng1, ng2);
}

Geometry::Geometry()
  {
  }

// generic part of 'Geometry' constructors
void Geometry::init_geometry (double fs, double ss, int ng1, int ng2)
{
  first_size = fs;
  second_size = ss;
  n_grid_1 = ng1;
  n_grid_2 = ng2;
  dr = set_dr();
  dz = set_dz();

  epsilon = new double*[n_grid_1];
#pragma omp parallel for
  for (int i=0; i<(n_grid_1);i++)
    epsilon[i]= new double[n_grid_2];

  sigma = new double*[n_grid_1];
#pragma omp parallel for
  for (int i=0; i<(n_grid_1);i++)
    sigma[i]= new double[n_grid_2];

#pragma omp parallel for
  for (int i=0; i<(n_grid_1);i++)
    for (int k=0; k<(n_grid_2);k++)
      sigma[i][k]=0;
}

void Geometry::set_epsilon()
{
#pragma omp parallel for
  for(int i=0;i<(n_grid_1);i++)
    for(int k=0;k<(n_grid_2);k++)
      epsilon[i][k]=1;
}

Geometry::~Geometry(void)
{
}

double Geometry::set_dr()
{
  return first_size/(n_grid_1-1);
}

double Geometry::set_dz()
{
  return second_size/(n_grid_2-1);
}
