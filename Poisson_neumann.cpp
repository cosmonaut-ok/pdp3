#include "poisson_neumann.h"
#include "Constant.h"

using namespace constant;

Poisson_neumann::Poisson_neumann(Geometry* cyl_geom):Poisson(cyl_geom)
{
  // charge_density
  t_charge_density = new double*[cyl_geom->n_grid_1];
  for (int i=0; i<(cyl_geom->n_grid_1);i++)
    t_charge_density[i] = new double[cyl_geom->n_grid_2];
}

Poisson_neumann::~Poisson_neumann(void)
{
  for (int i=0; i<(cyl_geom->n_grid_1);i++)
    delete[]t_charge_density[i];
  delete[]t_charge_density;
}
void Poisson_neumann::poisson_solve(E_field* input_e, charge_density* ro1)
{
  double a = 0;
  double c = 0;
  double b = 0;
  double d = 0;
  double* alpha = new double [cyl_geom->n_grid_1];
  double* beta = new double [cyl_geom->n_grid_1];
  Fourier* four1 = 0;
  double** ro = ro1->get_ro();
  double** e1 = input_e->e1;
  double** e3 = input_e->e3;
  double** fi = input_e->fi;
  double dr = cyl_geom->dr;
  double dz = cyl_geom->dz;

#pragma omp parallel shared (alpha, beta, four1, ro, e1, e3, fi, dr, dz)
  {
    // copy charge_density array in to temp array
#pragma omp for
    for (int i=0;i<(cyl_geom->n_grid_1);i++)
      for(int k=0;k<(cyl_geom->n_grid_2);k++)
        t_charge_density[i][k]= ro[i][k];

    // call function for cosine transform
#pragma omp for
    for (int i=0;i<cyl_geom->n_grid_1;i++)
      four1->fast_cosine_transform(t_charge_density, cyl_geom->n_grid_2, i, false);

    // sweep method
#pragma omp for
    for(int k=0;k<(cyl_geom->n_grid_2);k++)
    {
      b = 2.0 - 2.0*(cos(PI*k/(cyl_geom->n_grid_2-1)) - 1)*dr*dr/(dz*dz);
      d = cyl_geom->dr*cyl_geom->dr*t_charge_density[0][k]/EPSILON0;
      alpha[1]=4.0/(2.0+b);
      beta[1]=d/(2.0+b);

      for (int i=1;i<(cyl_geom->n_grid_1-1);i++)
      {
        //ay-1 - by + cy+1= = -d//

        a = 1.0-1.0/(2.0*(i));
        c = 1.0+1.0/(2.0*(i));
        d = cyl_geom->dr*cyl_geom->dr*t_charge_density[i][k]/EPSILON0;

        alpha[i+1] = c/(b-alpha[i]*a);
        beta[i+1] = (d+beta[i]*a)/(b-alpha[i]*a);
      }

      for(int i=(cyl_geom->n_grid_1-2);i>=0;i--)
        fi[i][k]=beta[i+1]+alpha[i+1]*fi[i+1][k];
    }

    // call function for inverse cosine transform
#pragma omp for
    for (int i=0;i<cyl_geom->n_grid_1;i++)
      four1->fast_cosine_transform(fi, cyl_geom->n_grid_2, i, true);

    // calculate electric field
    for (int i=0;i<(cyl_geom->n_grid_1-1);i++)
      for (int k=0;k<(cyl_geom->n_grid_2);k++)
        e1[i][k]=(fi[i][k]-fi[i+1][k])/cyl_geom->dr;

    for (int i=0;i<(cyl_geom->n_grid_1);i++)
      for (int k=0;k<(cyl_geom->n_grid_2-1);k++)
        e3[i][k]=(fi[i][k]-fi[i][k+1])/cyl_geom->dz;

    delete [] alpha;
    delete [] beta;
  }
}
