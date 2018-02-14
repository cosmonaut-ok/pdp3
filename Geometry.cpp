#include <math.h>
#include "Geometry.h"

Geometry::Geometry()
  {
  }

// Geometry constructor
Geometry::Geometry (double fs, double ss, int ng1, int ng2)
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

// PML
void Geometry::setPML(double comparative_l_1, double comparative_l_2, double comparative_l_3,
                      double sigma1, double sigma2)
{
  // defining lenght of sigma calculation ion left wall
  double lenght_sigma_left = dz * (floor(n_grid_2*comparative_l_1));

  // defining lenght of sigma calculation on right wall
  double lenght_sigma_right = dz * (floor(n_grid_2*comparative_l_2));

  // defining lenght of sigma calculation on z-wall
  double lenght_sigma_extern = dr * (floor(n_grid_1*comparative_l_3));

  // if pml is only on z walll
  if ((comparative_l_1 == 0) && (comparative_l_2 == 0) && (comparative_l_3 != 0))
  {
    for(int i=0; i<n_grid_1; i++)
      for(int k=0; k<n_grid_2; k++)
        if ((first_size - dr*(i)) <= lenght_sigma_extern)
          sigma[i][k] = sigma1 + (sigma2-sigma1)/(pow(lenght_sigma_extern, 2)) *
            pow((dr*(i+1)-first_size+lenght_sigma_extern), 2);
  }
  else
  {
    // sigma on r=0 and r=r walls
    for(int k=0;k<n_grid_2;k++)
      for(int i=0;i<n_grid_1;i++)
      {
        if (dz*k<lenght_sigma_left)
        {
          sigma[i][k] = sigma1 + (sigma2-sigma1)/(pow(lenght_sigma_left, 2)) *
            pow((lenght_sigma_left-dz*k), 2);
        }
        if ((second_size - dz*(k))<lenght_sigma_right)
        {
          sigma[i][k] = sigma1 + (sigma2-sigma1)/(pow(lenght_sigma_right, 2)) *
            pow((dz*(k+1)-second_size+lenght_sigma_right), 2);
        }
      }
  }

  if (comparative_l_3 == 0)
  {
  }
  else
  {
    // sigma assigning
    // sigma on z=z wall
    for(int i=0;i < n_grid_1; i++)
      for(int k=0;k < n_grid_2; k++)
        if ((first_size - dr*(i)) <= lenght_sigma_extern)
          if (((first_size - dr*(i+1)) < lenght_sigma_extern * dz *
               (k)/lenght_sigma_left) && ( (first_size - dr*(i)) <=
                                           (lenght_sigma_extern/lenght_sigma_right *
                                            (second_size - dz*(k)))))
            sigma[i][k] = sigma1 + (sigma2-sigma1)/(pow(lenght_sigma_extern,2)) *
              pow((dr*(i+1)-first_size+lenght_sigma_extern),2);
  }
  return;
}
