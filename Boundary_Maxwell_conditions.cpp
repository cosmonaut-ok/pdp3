#include "Boundary_Maxwell_conditions.h"
#include "Constant.h"

using namespace constant;

Boundary_Maxwell_conditions::Boundary_Maxwell_conditions(E_field* e_fld_t):e_fld(e_fld_t)
{
}

Boundary_Maxwell_conditions::Boundary_Maxwell_conditions(void)
{
}

Boundary_Maxwell_conditions::~Boundary_Maxwell_conditions(void)
{
}

void Boundary_Maxwell_conditions::specify_initial_field(Geometry* cyl_geom,
                                                        double E_fi_upper,
                                                        double E_fi_left,
                                                        double E_fi_right)
{
  int n_grid1 = cyl_geom->n_grid_1;
  int n_grid2 = cyl_geom->n_grid_2;
// setazimuthal component electric field initial value
#pragma omp parallel shared(E_fi_left, E_fi_right, E_fi_upper, n_grid1, n_grid2)
  {
#pragma omp for
    for (int i=0;i<(n_grid1);i++)
    {
      e_fld->e2[i][0]=E_fi_left;
      e_fld->e2[i][n_grid2-1]=E_fi_right;
    }

#pragma omp for
    for(int k=0;k<n_grid2;k++)
    {
      e_fld->e2[n_grid1-1][k]=E_fi_upper;
    }
  }
}
