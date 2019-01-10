#include "poissonNeumann.h"

using namespace constant;
using namespace math::fourier;

PoissonNeumann::PoissonNeumann(Geometry *cyl_geom):Poisson(cyl_geom)
{
  // charge_density
  t_charge_density = new double*[cyl_geom->n_grid_r];
  for (int i = 0; i < (cyl_geom->n_grid_r); i++)
    t_charge_density[i] = new double[cyl_geom->n_grid_z];
}

PoissonNeumann::~PoissonNeumann(void)
{
  for (int i = 0; i < (cyl_geom->n_grid_r); i++)
    delete[]t_charge_density[i];
  delete[]t_charge_density;
}
void PoissonNeumann::poisson_solve(EField *input_e, ChargeDensity *ro1)
{
  double a = 0;
  double c = 0;
  double b = 0;
  double d = 0;
  double *alpha = new double [cyl_geom->n_grid_r];
  double *beta = new double [cyl_geom->n_grid_r];
  Fourier *four1 = 0;
  double **ro = ro1->get_rho();
  double **e1 = input_e->field_r;
  double **e3 = input_e->field_z;
  double **fi = input_e->fi;
  double dr = cyl_geom->dr;
  double dz = cyl_geom->dz;

#pragma omp parallel shared (alpha, beta, four1, ro, e1, e3, fi, dr, dz)
  {
    // copy charge_density array in to temp array
#pragma omp for
    for (int i = 0; i < (cyl_geom->n_grid_r); i++)
      for(int k = 0; k < (cyl_geom->n_grid_z); k++)
        t_charge_density[i][k]= ro[i][k];

    // call function for cosine transform
#pragma omp for
    for (int i = 0; i < cyl_geom->n_grid_r; i++)
      four1->fast_cosine_transform(t_charge_density, cyl_geom->n_grid_z, i, false);

    // sweep method
#pragma omp for
    for(int k = 0; k < (cyl_geom->n_grid_z); k++)
    {
      b = 2.0 - 2.0*(cos(PI*k/(cyl_geom->n_grid_z-1)) - 1)*dr*dr/(dz*dz);
      d = cyl_geom->dr*cyl_geom->dr*t_charge_density[0][k]/EPSILON0;
      alpha[1]=4.0/(2.0+b);
      beta[1]=d/(2.0+b);

      for (int i=1; i < (cyl_geom->n_grid_r-1); i++)
      {
        //ay-1 - by + cy+1= = -d//

        a = 1.0-1.0/(2.0*(i));
        c = 1.0+1.0/(2.0*(i));
        d = cyl_geom->dr*cyl_geom->dr*t_charge_density[i][k]/EPSILON0;

        alpha[i+1] = c/(b-alpha[i]*a);
        beta[i+1] = (d+beta[i]*a)/(b-alpha[i]*a);
      }

      for(int i=(cyl_geom->n_grid_r-2);i>=0;i--)
        fi[i][k]=beta[i+1]+alpha[i+1]*fi[i+1][k];
    }

    // call function for inverse cosine transform
#pragma omp for
    for (int i = 0; i < cyl_geom->n_grid_r; i++)
      four1->fast_cosine_transform(fi, cyl_geom->n_grid_z, i, true);

    // calculate electric field
    for (int i = 0; i < (cyl_geom->n_grid_r-1); i++)
      for (int k = 0; k < (cyl_geom->n_grid_z); k++)
        e1[i][k]=(fi[i][k]-fi[i+1][k])/cyl_geom->dr;

    for (int i = 0; i < (cyl_geom->n_grid_r); i++)
      for (int k = 0; k < (cyl_geom->n_grid_z-1); k++)
        e3[i][k]=(fi[i][k]-fi[i][k+1])/cyl_geom->dz;

    delete [] alpha;
    delete [] beta;
  }
}
