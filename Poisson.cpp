#include "Poisson.h"
#include "Constant.h"

using namespace constant;

Poisson::Poisson(Geometry*cyl_geom_t):cyl_geom(cyl_geom_t)
{
}

Poisson::~Poisson(void)
{
}

bool Poisson::test_poisson_equation(E_field* input_e ,charge_density *rho)
{
  double dr = cyl_geom->dr;
  double dz = cyl_geom->dz;
  double accur =1e-3;
  double a=0;
  bool res =true;
  double **rho1 = rho->get_ro();
  double** e1 =input_e->e1;
  double** e3 = input_e->e3;
  ofstream delta("delta");
  for (int i=1;i<cyl_geom->n_grid_1-1;i++)
    for (int k=1;k<cyl_geom->n_grid_2-1;k++)
    {
      a = rho1[i][k]/EPSILON0 - (e1[i][k]+e1[i-1][k])/(i*2.0*dr) - (e1[i][k]-e1[i-1][k])/(dr)-(e3[i][k]-e3[i][k-1])/dz;
      if (abs(a)>accur)
        res=false;
      delta<<a<<" ";
    }
  delta.close();
  return res;
}
