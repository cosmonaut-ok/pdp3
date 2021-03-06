#include "eField.h"

using namespace std;
using namespace constant;
using namespace math::fourier;

EField::EField()
{
}

//// constructor
EField::EField(Geometry *geom1_t) : Field (geom1_t)
{
  // fi
  fi = new double *[geom1->n_grid_r];

  // filling second demension
  // for potential
#pragma omp parallel for
  for (int i = 0; i < (geom1->n_grid_r); i++)
    fi[i]= new double[geom1->n_grid_z];

}

// destructor
EField::~EField()
{
  for (int i = 0; i < (geom1->n_grid_r); i++)
    delete[]fi[i];
  delete[]fi;
}

// initial E distribution
void EField::set_homogeneous_efield(double E_r, double E_phi, double E_z)
{
#pragma omp parallel for shared (E_r, E_phi, E_z)
  for(int i = 0; i < (geom1->n_grid_r-1); i++)
    for(int k = 0; k < (geom1->n_grid_z-1); k++)
    {
      field_r[i][k] = E_r;
      field_phi[i][k] = E_phi;
      field_z[i][k] = E_z;
    }
}

void EField::boundary_conditions()
{
  // Border E_r conditions
  // last elements of array [ngrid-1]!!!
#pragma omp parallel
  {
#pragma omp for
    for(int i = 0; i < (geom1->n_grid_r - 1); i++)
    {
      field_r[i][0] = 0;
      field_r[i][geom1->n_grid_z - 1] = 0;
    }

    // Border E_phi conditions
#pragma omp for
    for(int k = 0; k < (geom1->n_grid_z); k++)
    {
      field_phi[0][k] = 0;
      field_phi[geom1->n_grid_r - 1][k] = 0;
    }

#pragma omp for
    for(int i = 0; i < (geom1->n_grid_r); i++)
    {
      field_phi[i][0] = 0;
      field_phi[i][geom1->n_grid_z - 1] = 0;
    }

    // Border E_z conditions
#pragma omp for
    for(int k = 0; k < (geom1->n_grid_z - 1); k++)
    {
      field_z[0][k] = 0;
      field_z[geom1->n_grid_r-1][k] = 0;
    }
  }
}

// Electric field calculation
void EField::calc_field(HField *h_field1,
                         Time *time1,
                         Current *current1)
{
  double **j_r = current1->get_j_r();
  double **j_phi = current1->get_j_phi();
  double **j_z = current1->get_j_z();
  double **h_r = h_field1->field_r;
  double **h_phi = h_field1->field_phi;
  double **h_z = h_field1->field_z;
  double dr = geom1->dr;
  double dz = geom1->dz;

  // Er first[i] value
#pragma omp parallel
  {
#pragma omp for
    for(int k = 1; k < (geom1->n_grid_z - 1); k++)
    {
      int i = 0;
      double epsilonx2 = 2 * geom1->epsilon[i][k] * EPSILON0;
      double sigma_t = geom1->sigma[i][k] * time1->delta_t;

      double koef_e = (epsilonx2 - sigma_t) / (epsilonx2 + sigma_t);

      double koef_h =  2 * time1->delta_t / (epsilonx2 + sigma_t);

      field_r[i][k]=field_r[i][k] * koef_e
        - (j_r[i][k]+(h_phi[i][k] - h_phi[i][k-1]) / dz) * koef_h;
    }

    // Ez=on axis
#pragma omp for
    for(int k = 0; k < (geom1->n_grid_z - 1); k++)
    {
      int i = 0;
      double epsilonx2 = 2 * geom1->epsilon[i][k] * EPSILON0;
      double sigma_t = geom1->sigma[i][k] * time1->delta_t;

      double koef_e = (epsilonx2 - sigma_t) / (epsilonx2 + sigma_t);
      double koef_h = 2 * time1->delta_t / (epsilonx2 + sigma_t);

      field_z[i][k]=field_z[i][k] * koef_e
        - (j_z[i][k] - 4.0 / dr*h_phi[i][k]) * koef_h;
    }

#pragma omp for
    for(int i=1; i < (geom1->n_grid_r - 1); i++)
      for(int k=1; k < (geom1->n_grid_z - 1); k++)
      {
        double epsilonx2 = 2 * geom1->epsilon[i][k] * EPSILON0;
        double sigma_t = geom1->sigma[i][k] * time1->delta_t;

        double koef_e = (epsilonx2 - sigma_t) / (epsilonx2 + sigma_t);
        double koef_h = 2 * time1->delta_t / (epsilonx2 + sigma_t);

        field_r[i][k] = field_r[i][k] * koef_e
          - (j_r[i][k] + (h_phi[i][k] - h_phi[i][k-1]) / dz) * koef_h;

        field_phi[i][k] = field_phi[i][k] * koef_e
          - (j_phi[i][k] - (h_r[i][k] - h_r[i][k-1])
             / dz + (h_z[i][k] - h_z[i-1][k]) / dr) * koef_h;

        field_z[i][k] = field_z[i][k] * koef_e
          - (j_z[i][k] - (h_phi[i][k] - h_phi[i-1][k]) / dr
             - (h_phi[i][k] + h_phi[i-1][k]) / (2.0 * dr * i)
            ) * koef_h;
      }

#pragma omp for
    for(int i=1; i < (geom1->n_grid_r-1); i++)
    {
      int k = 0;
      double epsilonx2 = 2 * geom1->epsilon[i][k] * EPSILON0;
      double sigma_t = geom1->sigma[i][k] * time1->delta_t;

      double koef_e = (epsilonx2 - sigma_t) / (epsilonx2 + sigma_t);
      double koef_h = 2 * time1->delta_t / (epsilonx2 + sigma_t);

      field_z[i][k] = field_z[i][k] * koef_e -
        (j_z[i][k] - (h_phi[i][k] - h_phi[i-1][k]) / dr
         - (h_phi[i][k] + h_phi[i-1][k]) / (2.0 * dr * i)
          ) * koef_h;
    }
  }
}

// set fi on r=z boundary
void EField::set_fi_on_z()
{
#pragma omp parallel for
  for (int k = 0; k < (geom1->n_grid_z); k++)
    fi[geom1->n_grid_r-1][k] = 0;
}

void EField::tridiagonal_solve(const double *a,
                               const double *b,
                               double *c,
                               double *d, double *x,
                               int n)
{
  // Modify the coefficients
  c[0] /= b[0]; // Division by zero risk
  d[0] /= b[0]; // Division by zero would imply a singular matrix

  for(int i = 1; i < n; i++)
  {
    double id = (b[i] - c[i-1] * a[i]); // Division by zero risk
    c[i] /= id; // Last value calculated is redundant
    d[i] = (d[i] - d[i-1] * a[i]) / id;
  }

  // Now, back substitute
  x[n - 1] = d[n - 1];

  for(int i = n - 2; i >= 0; i--)
    x[i] = d[i] - c[i] * x[i+1];
}

double* EField::get_field(double radius, double longitude)
//! function for electric field weighting
{
  int i_r = 0; // number of particle i cell
  int k_z = 0; // number of particle k cell

  double dr = geom1->dr;
  double dz = geom1->dz;
  double r1, r2, r3; // temp variables for calculation
  double dz1, dz2; // temp var.: width of k and k+1 cell
  double er = 0;
  double efi = 0;
  double ez = 0;
  double vol_1 = 0; // volume of i cell; Q/V, V - volume of elementary cell
  double vol_2 = 0; // volume of i+1 cell;

  r1 = radius - 0.5 * dr;
  r3 = radius + 0.5 * dr;

  // weighting of E_r
  // finding number of cell. example dr=0.5, radius = 0.7, i_r =0;!!
  i_r = CELL_NUMBER(radius - 0.5 * dr, dr);
  k_z = CELL_NUMBER(longitude, dz);
  // TODO: workaround: sometimes it gives -1.
  // Just get 0 cell if it happence
  if (i_r < 0) i_r = 0;
  if (k_z < 0) k_z = 0;

  vol_1 = CELL_VOLUME(i_r+1, dr, dz);
  vol_2 = CELL_VOLUME(i_r+3, dr, dz);
  dz1 = (k_z+1) * dz - longitude;
  dz2 = longitude - k_z * dz;
  r2 = (i_r + 1) * dr;

  //weighting Er[i][k]//
  er = er + field_r[i_r][k_z] * CYL_RNG_VOL(dz1, r1, r2) / vol_1;

  //weighting Er[i+1][k]//
  er = er + field_r[i_r+1][k_z] * CYL_RNG_VOL(dz1, r2, r3) / vol_2;

  //weighting Er[i][k+1]//
  er= er + field_r[i_r][k_z+1] * CYL_RNG_VOL(dz2, r1, r2) / vol_1;

  //weighting Er[i+1][k+1]//
  er = er + field_r[i_r+1][k_z+1] * CYL_RNG_VOL(dz2, r2, r3) / vol_2;

  // weighting of E_z
  // finding number of cell. example dz=0.5, longitude = 0.7, z_k =0;!!
  i_r = CELL_NUMBER(radius, dr);
  k_z = CELL_NUMBER(longitude - 0.5 * dz, dz);
  // TODO: workaround: sometimes it gives -1.
  // Just get 0 cell if it happence
  if (i_r < 0) i_r = 0;
  if (k_z < 0) k_z = 0;

  if (radius > dr)
    vol_1 = CELL_VOLUME(i_r, dr, dz);
  else
    vol_1 = CYL_VOL(dz, dr); // volume of first cell

  r2 = (i_r + 0.5) * dr;
  vol_2 = CELL_VOLUME(i_r+2, dr, dz);
  dz1 = (k_z + 1.5) * dz - longitude;
  dz2 = longitude - (k_z + 0.5) * dz;

  // weighting Ez[i][k]
  ez += field_z[i_r][k_z] * CYL_RNG_VOL(dz1, r1, r2) / vol_1;

  // weighting Ez[i+1][k]
  ez += field_z[i_r+1][k_z] * CYL_RNG_VOL(dz1, r2, r3) / vol_2;

  // weighting Ez[i][k+1]
  ez += field_z[i_r][k_z+1] * CYL_RNG_VOL(dz2, r1, r2) / vol_1;

  //weighting Ez[i+1][k+1]//
  ez += field_z[i_r+1][k_z+1] * CYL_RNG_VOL(dz2, r2, r3) / vol_2;

  // weighting of E_fi
  // finding number of cell. example dz=0.5, longitude = 0.7, z_k =1;

  i_r = CELL_NUMBER(radius, dr);
  k_z = CELL_NUMBER(longitude, dz);
  // TODO: workaround: sometimes it gives -1.
  // Just get 0 cell if it happence
  if (i_r < 0) i_r = 0;
  if (k_z < 0) k_z = 0;

  if(radius>dr)
    vol_1 = CELL_VOLUME(i_r, dr, dz);
  else
    vol_1 = CYL_VOL(dz, dr); // volume of first cell

  r2 = (i_r + 0.5) * dr;
  vol_2 = CELL_VOLUME(i_r+2, dr, dz);
  dz1 = (k_z+1)*dz-longitude;
  dz2 = longitude - k_z * dz;

  // weighting Efi[i][k]
  efi += field_phi[i_r][k_z] * CYL_RNG_VOL(dz1, r1, r2) / vol_1;

  // weighting Efi[i+1][k]
  efi += field_phi[i_r+1][k_z] * CYL_RNG_VOL(dz1, r2, r3) / vol_2;

  // weighting Efi[i][k+1]
  efi += field_phi[i_r][k_z+1] * CYL_RNG_VOL(dz2, r1, r2) / vol_1;

  // weighting Efi[i+1][k+1]
  efi += field_phi[i_r+1][k_z+1] * CYL_RNG_VOL(dz2, r2, r3) / vol_2;

  double* components = tinyvec3d::mkvector3d(er, efi, ez);

  return components;
}
