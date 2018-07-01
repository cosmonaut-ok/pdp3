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
  // Calling parent constructor

  //// fi
  fi = new double*[geom1->n_grid_1];
  //// charge_density
  t_charge_density = new double*[geom1->n_grid_1];

#pragma omp parallel for
  // filling second demension
  // for potential and charge density
  for (int i=0; i<(geom1->n_grid_1); i++)
  {
    fi[i]= new double[geom1->n_grid_2];
    t_charge_density[i] = new double[geom1->n_grid_2];
  }

  // set zero initial temp ro
#pragma omp for
  for (int i=0; i<(geom1->n_grid_1); i++)
    for (int k=0; k<(geom1->n_grid_2); k++)
      t_charge_density[i][k] = 0;
}

// destructor
EField::~EField()
{
  for (int i=0; i<(geom1->n_grid_1);i++)
  {
    delete[]fi[i];
    delete[]t_charge_density[i];
  }
  delete[]fi;
  delete[]t_charge_density;
}

// initial E distribution
void EField::set_homogeneous_efield(double E_r, double E_phi, double E_z)
{
#pragma omp parallel for shared (E_r, E_phi, E_z)
  for(int i=0; i<(geom1->n_grid_1-1); i++)
    for(int k=0; k<(geom1->n_grid_2-1); k++)
    {
      // Er
      field_r[i][k] = E_r;
      // Ef
      field_phi[i][k] = E_phi;
      // Ez
      field_z[i][k] = E_z;
    }
}

void EField::boundary_conditions()
{
  //// Border Er conditions
  // last elements of array [ngrid-1]!!!
#pragma omp parallel
  {
#pragma omp for
    for(int i=0; i<(geom1->n_grid_1-1); i++)
    {
      field_r[i][0]=0;
      field_r[i][geom1->n_grid_2-1]=0;
      //field_r[i][0]=limit_conditions();
    }

    //// Border Ef conditions
#pragma omp for
    for(int k=0; k<(geom1->n_grid_2); k++)
    {
      field_phi[0][k]=0;
      field_phi[geom1->n_grid_1-1][k]=0;
    }

#pragma omp for
    for(int i=0; i<(geom1->n_grid_1); i++)
    {
      field_phi[i][0]=0;
      field_phi[i][geom1->n_grid_2-1]=0;
    }

    //// Border Ez conditions
#pragma omp for
    for(int k=0; k<(geom1->n_grid_2-1); k++)
    {
      field_z[0][k]=0;
      field_z[geom1->n_grid_1-1][k]=0;
      // field_r[0][k]=limit_conditions();
    }
  }
}

// Electric field calculation
void EField::calc_field(HField *h_field1,
                         Time *time1,
                         Current *current1)
{
  double **j1 = current1->get_j1();
  double **j2 = current1->get_j2();
  double **j3 = current1->get_j3();
  double koef_e = 0;
  double koef_h = 0;
  double **h_r = h_field1->field_r;
  double **h_phi = h_field1->field_phi;
  double **h_z = h_field1->field_z;
  double dr = geom1->dr;
  double dz = geom1->dz;

  //// Er first[i] value
#pragma omp parallel for
  for(int k=1; k<(geom1->n_grid_2-1); k++)
  {
    int i=0;
    koef_e = (2.0*geom1->epsilon[i][k]*EPSILON0 - geom1->sigma[i][k]*time1->delta_t) /
      (2.0*geom1->epsilon[i][k]*EPSILON0 + geom1->sigma[i][k]*time1->delta_t);

    koef_h =  2*time1->delta_t/(2.0*geom1->epsilon[i][k]*EPSILON0 + geom1->sigma[i][k]*time1->delta_t);

    field_r[i][k]=field_r[i][k] * koef_e  - (j1[i][k]+(h_phi[i][k] - h_phi[i][k-1])/dz)*koef_h;
  }

  //// Ez=on axis// // ???????????
  for(int k=0; k<(geom1->n_grid_2-1); k++)
  {
    int i=0;
    koef_e = (2.0*geom1->epsilon[i][k]*EPSILON0 - geom1->sigma[i][k]*time1->delta_t) /
      (2.0*geom1->epsilon[i][k]*EPSILON0 + geom1->sigma[i][k]*time1->delta_t);
    koef_h =  2*time1->delta_t/(2.0*geom1->epsilon[i][k]*EPSILON0 + geom1->sigma[i][k]*time1->delta_t);

    field_z[i][k]=field_z[i][k]*koef_e - (j3[i][k]-4.0/dr*h_phi[i][k])*koef_h;
  }

  for(int i=1; i<(geom1->n_grid_1-1); i++)
    for(int k=1; k<(geom1->n_grid_2-1); k++)
    {
      koef_e = (2.0*geom1->epsilon[i][k]*EPSILON0 - geom1->sigma[i][k]*time1->delta_t) /
        (2.0*geom1->epsilon[i][k]*EPSILON0 + geom1->sigma[i][k]*time1->delta_t);
      koef_h =  2*time1->delta_t/(2.0*geom1->epsilon[i][k]*EPSILON0 + geom1->sigma[i][k]*time1->delta_t);
      field_r[i][k]=field_r[i][k]*koef_e - (j1[i][k]+(h_phi[i][k]-h_phi[i][k-1])/dz)*koef_h;
      field_phi[i][k]=field_phi[i][k]*koef_e - (j2[i][k]-(h_r[i][k]-h_r[i][k-1])/dz + (h_z[i][k]-h_z[i-1][k])/dr)*koef_h;
      field_z[i][k]=field_z[i][k]*koef_e -(j3[i][k]-(h_phi[i][k]-h_phi[i-1][k])/dr - (h_phi[i][k]+h_phi[i-1][k])/(2.0*dr*i))*koef_h;
    }

  for(int i=1; i<(geom1->n_grid_1-1); i++)
  {
    int k=0;
    koef_e = (2.0*geom1->epsilon[i][k]*EPSILON0 - geom1->sigma[i][k]*time1->delta_t) /
      (2.0*geom1->epsilon[i][k]*EPSILON0 + geom1->sigma[i][k]*time1->delta_t);
    koef_h =  2*time1->delta_t/(2.0*geom1->epsilon[i][k]*EPSILON0 + geom1->sigma[i][k]*time1->delta_t);
    field_z[i][k]=field_z[i][k]*koef_e - (j3[i][k] - (h_phi[i][k]-h_phi[i-1][k])/dr - (h_phi[i][k]+h_phi[i-1][k])/(2.0*dr*i))*koef_h;
  }
}

// set fi on r=z boundary
void EField::set_fi_on_z()
{
#pragma omp parallel for
  for (int k=0; k<(geom1->n_grid_2);k++)
    fi[geom1->n_grid_1-1][k]=0;
}

// poisson equation solving 2
void EField::poisson_equation2(Geometry *geom1, ChargeDensity *ro1)
{
  // const double epsilon0 = EPSILON0;
  double phi0(0.0);
  double *a = new double [geom1->n_grid_1];
  double *b = new double [geom1->n_grid_1];
  double *c = new double [geom1->n_grid_1];
  double *d = new double [geom1->n_grid_1];
  double *c1 = new double [geom1->n_grid_1];
  double *d1 = new double [geom1->n_grid_1];
  double *phi = new double [geom1->n_grid_1];
  Fourier *four1=0;

  // double dr = geom1->dr;
  double dr2 = geom1->dr*geom1->dr;

  double **ro = ro1->get_rho();

  // copy charge_density array in to temp array
#pragma omp parallel for shared (ro)
  for(int i=0; i<(geom1->n_grid_1); i++)
    for(int k=0;k<(geom1->n_grid_2); k++)
      t_charge_density[i][k]= ro[i][k];

  // call function for cosine transform

  int temp=geom1->n_grid_2;

#pragma omp parallel for shared (temp, four1)
  for (int i=0; i<geom1->n_grid_1; i++)
    four1->fast_cosine_transform((double**)t_charge_density, temp, i, false);
    //four1->fast_fourier_transform(t_charge_density, temp, i, false);

  //set coefficients
  b[0] = 1.0;
  c[0] = 1.0;

#pragma omp parallel for shared (a, b, c, d, c1, d1, phi, dr2)
  for (int k = 0; k < geom1->n_grid_2; k++)
  {
    b[0] = 1.0;
    c[0] = -1.0;
    d[0] = dr2/4.0/EPSILON0*t_charge_density[0][k];
    for (int i = 1; i < geom1->n_grid_1 -1; i++)
    {
      a[i] = (1.0 - 1.0/2.0/(double)i);
      b[i] = -2.0 + 2.0*(cos(PI*k/(geom1->n_grid_2-1)) - 1)*geom1->dr*geom1->dr/(geom1->dz*geom1->dz);
      c[i] = (1.0 + 1.0/2.0/(double)i);
      d[i] = t_charge_density[i][k]*dr2/EPSILON0;
      c1[i] = (1.0 - 1.0/2.0/(double)i);
      d1[i] = t_charge_density[i][k]*dr2/EPSILON0;
    }
    a[0] = 0;
    c[geom1->n_grid_1-2] = 0.0;
    d[geom1->n_grid_1-2] -= phi0;
    c1[geom1->n_grid_1-2] = 0.0;
    d1[geom1->n_grid_1-2] -= phi0;
    tridiagonal_solve(a, b, c, d, phi, geom1->n_grid_1-1);
    for (int i = 0; i < geom1->n_grid_1 -1; i++)
      fi[i][k] = phi[i];
  }

// call function for inverse cosine transform
#pragma omp parallel for shared (four1)
  for (int i=0; i<geom1->n_grid_1; i++)
    four1->fast_cosine_transform((double**)fi, geom1->n_grid_2, i, true);

  // calculate electric field
  for (int i=0; i<(geom1->n_grid_1-1); i++)
    for (int k=1; k<(geom1->n_grid_2); k++)
      field_r[i][k]=(fi[i][k]-fi[i+1][k])/geom1->dr;

  for (int i=0; i<(geom1->n_grid_1); i++)
    for (int k=1; k<(geom1->n_grid_2-1); k++)
      field_z[i][k]=(fi[i][k]-fi[i][k+1])/geom1->dz;

  ofstream er("er"), ez("ez"); // TODO: WTF?
#pragma omp parallel for
  for (int i=1; i<geom1->n_grid_1-1; i++)
    for (int k=1; k<geom1->n_grid_2-1; k++)
    {
      er<<field_r[i][k]<<" ";
      ez<<field_z[i][k]<<" ";
    }
  er.close();
  ez.close();

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

  for(int i = 1; i < n; i++){
    double id = (b[i] - c[i-1] * a[i]); // Division by zero risk
    c[i] /= id; // Last value calculated is redundant
    d[i] = (d[i] - d[i-1] * a[i])/id;
  }

  // Now, back substitute
  x[n - 1] = d[n - 1];

  for(int i = n - 2; i >= 0; i--)
    x[i] = d[i] - c[i] * x[i + 1];
}

double* EField::get_field(double x1, double x3)
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

  r1 = x1-0.5*dr;
  r3 = x1+0.5*dr;

  // weighting of E_r
  // finding number of cell. example dr=0.5, x1 = 0.7, i_r =0;!!
  i_r = (int)ceil((x1-0.5*dr)/geom1->dr)-1;
  k_z = (int)ceil((x3)/geom1->dz)-1;
  // TODO: workaround: sometimes it gives -1.
  // Just get 0 cell if it happence
  if (i_r < 0) { i_r = 0; }
  if (k_z < 0) { k_z = 0; }

  vol_1 = PI*dz*dr*dr*(2*i_r+1);
  vol_2 = PI*dz*dr*dr*(2*i_r+3);
  dz1 = (k_z+1)*dz-x3;
  dz2 = x3 - k_z*dz;
  r2 = (i_r+1)*dr;

  //weighting Er[i][k]//
  er = er + field_r[i_r][k_z]*(PI*dz1*(r2*r2-r1*r1))/vol_1;

  //weighting Er[i+1][k]//
  er = er + field_r[i_r+1][k_z]*(PI*dz1*(r3*r3-r2*r2))/vol_2;

  //weighting Er[i][k+1]//
  er= er + field_r[i_r][k_z+1]*(PI*dz2*(r2*r2-r1*r1))/vol_1;

  //weighting Er[i+1][k+1]//
  er = er + field_r[i_r+1][k_z+1]*(PI*dz2*(r3*r3-r2*r2))/vol_2;

  // weighting of E_z//
  // finding number of cell. example dz=0.5, x3 = 0.7, z_k =0;!!
  i_r = (int)ceil((x1)/geom1->dr)-1;
  k_z = (int)ceil((x3-0.5*dz)/geom1->dz)-1;
  // TODO: workaround: sometimes it gives -1.
  // Just get 0 cell if it happence
  if (i_r < 0) { i_r = 0; }
  if (k_z < 0) { k_z = 0; }

  if(x1>dr)
    vol_1 = PI*dz*dr*dr*2*i_r;
  else
    vol_1 = PI*dz*dr*dr/4.0; // volume of first cell

  r2 = (i_r+0.5)*dr;
  vol_2 = PI*dz*dr*dr*(2*i_r+2);
  dz1 = (k_z+1.5)*dz - x3;
  dz2 = x3 - (k_z+0.5)*dz;

  // weighting Ez[i][k]
  ez = ez + field_z[i_r][k_z]*(PI*dz1*(r2*r2-r1*r1))/vol_1;

  // weighting Ez[i+1][k]
  ez = ez + field_z[i_r+1][k_z]*PI*dz1*(r3*r3-r2*r2)/vol_2;

  // weighting Ez[i][k+1]
  ez = ez + field_z[i_r][k_z+1]*PI*dz2*(r2*r2-r1*r1)/vol_1;

  //weighting Ez[i+1][k+1]//
  ez = ez + field_z[i_r+1][k_z+1]*PI*dz2*(r3*r3-r2*r2)/vol_2;

  // weighting of E_fi
  // finding number of cell. example dz=0.5, x3 = 0.7, z_k =1;
  i_r = (int)ceil((x1)/geom1->dr)-1;
  k_z = (int)ceil((x3)/geom1->dz)-1;
  // TODO: workaround: sometimes it gives -1.
  // Just get 0 cell if it happence
  if (i_r < 0) { i_r = 0; }
  if (k_z < 0) { k_z = 0; }

  if(x1>dr)
    vol_1 = PI*dz*dr*dr*2*i_r;
  else
    vol_1 = PI*dz*dr*dr/4.0; // volume of first cell

  r2 = (i_r+0.5)*dr;
  vol_2 = PI*dz*dr*dr*(2*i_r+2);
  dz1 = (k_z+1)*dz-x3;
  dz2 = x3-k_z*dz;

  // weighting Efi[i][k]
  efi += field_phi[i_r][k_z]*PI*dz1*(r2*r2 - r1*r1)/vol_1;

  // weighting Efi[i+1][k]
  efi += field_phi[i_r+1][k_z]*PI*dz1*(r3*r3-r2*r2)/vol_2;

  // weighting Efi[i][k+1]
  efi += field_phi[i_r][k_z+1]*PI*dz2*(r2*r2-r1*r1)/vol_1;

  // weighting Efi[i+1][k+1]
  efi += field_phi[i_r+1][k_z+1]*PI*dz2*(r3*r3-r2*r2)/vol_2;

  double* components = math::tiny_vector3d::mkvector3d(er, efi, ez);

  return components;
}
