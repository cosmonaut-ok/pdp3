#include "hField.h"

using namespace constant;
using namespace math;

// Constructor
HField::HField(Geometry *geom1_l) : Field(geom1_l)
{
  // Calling parent constructor

  // Hr_half_time
  field_r_half_time = new double *[geom1->n_grid_1];
#pragma omp parallel for
  for (int i=0;i<(geom1->n_grid_1);i++)
    field_r_half_time[i] = new double[geom1->n_grid_2-1];

  // Hf half time
  field_phi_half_time = new double *[geom1->n_grid_1-1];
#pragma omp parallel for
  for (int i=0;i<(geom1->n_grid_1-1);i++)
    field_phi_half_time[i] = new double[geom1->n_grid_2-1];

  // Hz half time
  field_z_half_time = new double *[geom1->n_grid_1-1];
#pragma omp parallel for
  for (int i=0;i<(geom1->n_grid_1-1);i++)
    field_z_half_time[i] = new double[geom1->n_grid_2];

  // Ar
  Ar = new double *[geom1->n_grid_1];
#pragma omp parallel for
  for (int i=0;i<(geom1->n_grid_1);i++)
    Ar[i] = new double[geom1->n_grid_2];

  // Afi
  Afi = new double *[geom1->n_grid_1];
#pragma omp parallel for
  for (int i=0;i<(geom1->n_grid_1);i++)
    Afi[i] = new double[geom1->n_grid_2];

  // Az
  Az = new double *[geom1->n_grid_1];
#pragma omp parallel for
  for (int i=0;i<(geom1->n_grid_1);i++)
  {
    Az[i] = new double[geom1->n_grid_2];
  }

  // initialisation of magnetic potential
#pragma omp parallel
  {
#pragma omp for
    for(int i=0;i<geom1->n_grid_1-1;i++)
      for(int k=0;k<geom1->n_grid_2-1;k++)
      {
        field_r[i][k]=0;
        field_phi[i][k]=0;
        field_z[i][k]=0;
        field_r_half_time[i][k]=0;
        field_phi_half_time[i][k]=0;
        field_z_half_time[i][k]=0;
      }
#pragma omp for
    for(int k=0;k<geom1->n_grid_2-1;k++)
    {
      field_r[geom1->n_grid_1-1][k]=0;
      field_r_half_time[geom1->n_grid_1-1][k]=0;
    }
  }
}
// Destructor
HField::~HField(void)
{
  for (int i=0; i<(geom1->n_grid_1-1);i++)
    delete[]Ar[i];
  delete[]Ar;

  for (int i=0; i<(geom1->n_grid_1-1);i++)
    delete[]Afi[i];
  delete[]Afi;

  for (int i=0; i<(geom1->n_grid_1-1);i++)
    delete[]Az[i];
  delete[]Az;

}

void HField::set_homogeneous_h(double E_r, double E_phi, double E_z)
{
#pragma omp parallel shared (E_r, E_phi, E_z)
  {
#pragma omp for
    for (int i=0;i<geom1->n_grid_1;i++)
      for (int k=0;(k<geom1->n_grid_2-1);k++)
      {
        //if (field_r[i][k]!= NULL)
        field_r[i][k]=E_r;
        field_r_half_time[i][k]=E_r;
      }
#pragma omp for
    for (int i=0;i<(geom1->n_grid_1-1);i++)
      for (int k=0;(k<geom1->n_grid_2-1);k++)
      {
        //if (field_r[i][k]!= NULL)
        field_phi[i][k]=E_phi;
        field_phi_half_time[i][k]=E_phi;
      }
#pragma omp for
    for (int i=0;i<(geom1->n_grid_1-1);i++)
      for (int k=0;(k<geom1->n_grid_2);k++)
      {
        //if (field_r[i][k]!= NULL)
        field_z[i][k]=E_z;
        field_z_half_time[i][k]=E_z;
      }
  }
}

// Field calculation
void HField::calc_field(EField *e_field1, Time *time1)
{
  double **e_r = e_field1->field_r;
  double **e_phi = e_field1->field_phi;
  double **e_z = e_field1->field_z;

  double dr = geom1->dr;
  double dz = geom1->dz;

  // double alpha;

  // Hr - last i value
#pragma omp parallel
  {
#pragma omp for
    for(int k=0; k<(geom1->n_grid_2-1); k++)
    {
      int i=geom1->n_grid_1-1;
      // alpha constant and delta_t production (to optimize calculations)
      double alpha_t = time1->delta_t
        * (e_phi[i][k+1]-e_phi[i][k]) / (dz * MAGN_CONST);

      field_r_half_time[i][k] = field_r[i][k] + alpha_t / 2;
      field_r[i][k] = field_r[i][k] + alpha_t;
    }

#pragma omp for
    for(int i=0; i<(geom1->n_grid_1 - 1); i++)
      for(int k=0; k<(geom1->n_grid_2 - 1); k++)
      {
        double alpha_t = time1->delta_t
          * (e_phi[i][k+1]-e_phi[i][k]) / (dz * MAGN_CONST);

        field_r_half_time[i][k] = field_r[i][k] + alpha_t / 2;
        field_r[i][k] = field_r[i][k] + alpha_t;
      }

#pragma omp for
    for(int i=0; i<(geom1->n_grid_1 - 1); i++)
      for(int k=0; k<(geom1->n_grid_2 - 1); k++)
      {
        double alpha_t = time1->delta_t
          * ((e_z[i+1][k] - e_z[i][k]) / dr
             - (e_r[i][k+1] - e_r[i][k]) / dz)
          / MAGN_CONST;

        field_phi_half_time[i][k] = field_phi[i][k] + alpha_t / 2;
        field_phi[i][k] = field_phi[i][k] + alpha_t;
      }

#pragma omp for
    for(int i=0; i<(geom1->n_grid_1 - 1); i++)
      for(int k=0; k<(geom1->n_grid_2 - 1); k++)
      {
        double alpha_t = time1->delta_t
          * ((e_phi[i+1][k] + e_phi[i][k]) / (2.0 * dr * (i + 0.5))
             + (e_phi[i+1][k] - e_phi[i][k]) / dr)
          / MAGN_CONST;

        field_z_half_time[i][k] = field_z[i][k] - alpha_t / 2;
        field_z[i][k] = field_z[i][k] - alpha_t;
      }
  }
}

double* HField::get_field(double radius, double longitude)
//! function for magnetic field weighting
{
  int i_r=0; // number of particle i cell
  int k_z=0; // number of particle k cell

  double dr = geom1->dr;
  double dz = geom1->dz;
  double r1, r2, r3; // temp variables for calculation
  double dz1, dz2; // temp var.: width of k and k+1 cell
  double hr =0;
  double hfi =0;
  double hz =0;
  double vol_1 =0; // volume of i cell; Q/V, V - volume of elementary cell
  double vol_2 =0; // volume of i+1 cell;

  r1 = radius-0.5*dr;
  r3 = radius+0.5*dr;

  //// weighting of H_z
  //finding number of cell. example dr=0.5, radius = 0.7, i_r =0;!!
  i_r = CELL_NUMBER(radius-0.5*dr, dr);
  k_z = CELL_NUMBER(longitude, dz);
  // TODO: workaround: sometimes it gives -1.
  // Just get 0 cell if it happence
  if (i_r < 0) i_r = 0;
  if (k_z < 0) k_z = 0;

  vol_1 = CELL_VOLUME(i_r+1, dr, dz);
  vol_2 = CELL_VOLUME(i_r+3, dr, dz);
  dz1 = (k_z+1)*dz-longitude;
  dz2 = longitude - k_z*dz;
  r2 = (i_r+1)*dr;

  // weighting Hz[i][k]
  hz = hz + field_z_half_time[i_r][k_z] * CYL_RNG_VOL(dz1, r1, r2) / vol_1;

  // weighting Hz[i+1][k]
  hz = hz + field_z_half_time[i_r+1][k_z] * CYL_RNG_VOL(dz1, r2, r3) / vol_2;

  // weighting Hz[i][k+1]
  hz= hz + field_z_half_time[i_r][k_z+1] * CYL_RNG_VOL(dz2, r2, r3) / vol_1;

  // weighting Hz[i+1][k+1]
  hz = hz + field_z_half_time[i_r+1][k_z+1] * CYL_RNG_VOL(dz2, r2, r3) / vol_2;

  //// weighting of Hr
  // finding number of cell. example dz=0.5, longitude = 0.7, z_k =0;!!
  i_r = CELL_NUMBER(radius, dr);
  k_z = CELL_NUMBER(longitude - 0.5 * dz, dz);
  // TODO: workaround: sometimes it gives -1.
  // Just get 0 cell if it happence
  if (i_r < 0) i_r = 0;
  if (k_z < 0) k_z = 0;

  if(radius>dr)
    vol_1 = CELL_VOLUME(i_r, dr, dz);
  else
    vol_1 = CYL_VOL(dz, dr); // volume of first cell

  r2 = (i_r+0.5)*dr;

  vol_2 = PI*dz*dr*dr*(2*i_r+2);
  dz1 = (k_z+1.5)*dz - longitude;
  dz2 = longitude - (k_z+0.5)*dz;

  // weighting Hr[i][k]
  hr += field_r_half_time[i_r][k_z]*(CYL_RNG_VOL(dz1, r1, r2)) / vol_1;

  // weighting Hr[i+1][k]
  hr += field_r_half_time[i_r+1][k_z] * CYL_RNG_VOL(dz1, r2, r3) / vol_2;

  // weighting Hr[i][k+1]
  hr += field_r_half_time[i_r][k_z+1] * CYL_RNG_VOL(dz2, r1, r2) / vol_1;

  // weighting Hr[i+1][k+1]
  hr += field_r_half_time[i_r+1][k_z+1] * CYL_RNG_VOL(dz2, r2, r3) / vol_2;

  //// weighting of H_fi
  // finding number of cell. example dz=0.5, longitude = 0.7, z_k =0;
  i_r = CELL_NUMBER(radius - 0.5 * dr, dr);
  k_z = CELL_NUMBER(longitude - 0.5 * dz, dz);
  // TODO: workaround: sometimes it gives -1.
  // Just get 0 cell if it happence
  if (i_r < 0) i_r = 0;
  if (k_z < 0) k_z = 0;

  r2 = (i_r+1) * dr;
  vol_1 = CELL_VOLUME(i_r+1, dr, dz);
  vol_2 = CELL_VOLUME(i_r+3, dr, dz);
  dz1 = (k_z+1.5) * dz - longitude;
  dz2 = longitude - (k_z+0.5) * dz;

  // weighting Hfi[i][k]
  hfi += field_phi_half_time[i_r][k_z] * CYL_RNG_VOL(dz1, r1, r2) / vol_1;

  // weighting Hfi[i+1][k]
  hfi += field_phi_half_time[i_r+1][k_z] * CYL_RNG_VOL(dz1, r2, r3) / vol_2;

  // weighting Hfi[i][k+1]
  hfi += field_phi_half_time[i_r][k_z+1] * CYL_RNG_VOL(dz2, r1, r2) / vol_1;

  // weighting Hfi[i+1][k+1]
  hfi += field_phi_half_time[i_r+1][k_z+1] * CYL_RNG_VOL(dz2, r2, r3) / vol_2;

  double* components = tinyvec3d::mkvector3d(hr, hfi, hz);

  return components;
}

//// /Return one dimensional field components
double *HField::get_1d_field_r()
{
  // copy 2d field array into 1d array rowwise
#pragma omp parallel for
  for (int i = 0; i < geom1->n_grid_1; i++)
    for (int k = 0; k < geom1->n_grid_2 - 1; k++)
      field_r_1d[i * (geom1->n_grid_2 - 1) + k] = field_r_half_time[i][k];
  return field_r_1d;
}

double *HField::get_1d_field_phi()
{
  // copy 2d field array into 1d array rowwise
#pragma omp parallel for
  for (int i = 0; i < geom1->n_grid_1 - 1; i++)
    for (int k = 0; k < geom1->n_grid_2 - 1; k++)
      field_phi_1d[i * (geom1->n_grid_2 - 1) + k] = field_phi_half_time[i][k];
  return field_phi_1d;
}

double *HField::get_1d_field_z()
{
  // copy 2d field array into 1d array rowwise
#pragma omp parallel for
  for (int i = 0; i < geom1->n_grid_1 - 1; i++)
    for (int k = 0; k < geom1->n_grid_2; k++)
      field_z_1d[i * geom1->n_grid_2 + k] = field_z_half_time[i][k];
  return field_z_1d;
}
