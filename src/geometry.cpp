#include "geometry.h"

Geometry::Geometry()
{
}

// Geometry constructor
Geometry::Geometry (double rs, double zs, int ngr, int ngz)
{
  r_size = rs;
  z_size = zs;
  n_grid_r = ngr;
  n_grid_z = ngz;
  dr = get_dr();
  dz = get_dz();

  epsilon = new double *[n_grid_r];
#pragma omp parallel for
  for (int i = 0; i < n_grid_r; i++)
    epsilon[i]= new double[n_grid_z];

  sigma = new double *[n_grid_r];
#pragma omp parallel for
  for (int i = 0; i < n_grid_r; i++)
    sigma[i]= new double[n_grid_z];

#pragma omp parallel for
  for (int i = 0; i < n_grid_r; i++)
    for (int k = 0; k<(n_grid_z); k++)
      sigma[i][k]=0;
}

void Geometry::set_epsilon()
{
#pragma omp parallel for
  for(int i = 0; i < n_grid_r; i++)
    for(int k = 0; k < n_grid_z; k++)
      epsilon[i][k]=1;
}

Geometry::~Geometry(void)
{
}

double Geometry::get_dr()
{
  return r_size / (n_grid_r - 1);
}

double Geometry::get_dz()
{
  return z_size / (n_grid_z - 1);
}

// PML
void Geometry::set_pml(double comparative_l_1, double comparative_l_2, double comparative_l_3,
                       double sigma1, double sigma2)
{
  // defining lenght of sigma calculation ion left wall
  double lenght_sigma_left = dz * (floor(n_grid_z * comparative_l_1));

  // defining lenght of sigma calculation on right wall
  double lenght_sigma_right = dz * (floor(n_grid_z * comparative_l_2));

  // defining lenght of sigma calculation on z-wall
  double lenght_sigma_extern = dr * (floor(n_grid_r * comparative_l_3));

  // if pml is only on z walll
  if ((comparative_l_1 == 0) && (comparative_l_2 == 0) && (comparative_l_3 != 0))
  {
    for(int i = 0; i < n_grid_r; i++)
      for(int k = 0; k < n_grid_z; k++)
        if (r_size - dr * i <= lenght_sigma_extern)
          sigma[i][k] = sigma1 + (sigma2 - sigma1) / pow(lenght_sigma_extern, 2)
            * pow((dr * (i + 1) - r_size + lenght_sigma_extern), 2);
  }
  else
  {
    // sigma on r=0 and r=r walls
    for(int k = 0; k < n_grid_z; k++)
      for(int i = 0; i < n_grid_r; i++)
      {
        if (dz * k < lenght_sigma_left)
          sigma[i][k] = sigma1 + (sigma2 - sigma1) / pow(lenght_sigma_left, 2)
            * pow(lenght_sigma_left - dz * k, 2);
        if (z_size - dz * k < lenght_sigma_right)
          sigma[i][k] = sigma1 + (sigma2 - sigma1) / pow(lenght_sigma_right, 2)
            * pow(dz * (k + 1) - z_size + lenght_sigma_right, 2);
      }
  }

  if (comparative_l_3 != 0)
  {
    // sigma assigning
    // sigma on z=z wall
    for(int i = 0; i < n_grid_r; i++)
      for(int k = 0; k < n_grid_z; k++)
        if (r_size - dr * i <= lenght_sigma_extern)
          if (
            (r_size - dr * (i + 1) < lenght_sigma_extern * dz * k / lenght_sigma_left)
            && (r_size - dr * i <= lenght_sigma_extern / lenght_sigma_right * (z_size - dz * k))
            )
            sigma[i][k] = sigma1 + (sigma2 - sigma1) / pow(lenght_sigma_extern, 2)
              * pow(dr * (i + 1) - r_size + lenght_sigma_extern, 2);
  }
  return;
}
