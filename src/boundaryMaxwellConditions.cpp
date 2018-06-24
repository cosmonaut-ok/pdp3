#include "boundaryMaxwellConditions.h"

using namespace constant;

BoundaryMaxwellConditions::BoundaryMaxwellConditions(EField *e_fld_t):e_fld(e_fld_t)
{
}

BoundaryMaxwellConditions::BoundaryMaxwellConditions(void)
{
}

BoundaryMaxwellConditions::~BoundaryMaxwellConditions(void)
{
}

void BoundaryMaxwellConditions::specify_initial_field(Geometry *cyl_geom,
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
      e_fld->field_phi[i][0]=E_fi_left;
      e_fld->field_phi[i][n_grid2-1]=E_fi_right;
    }

#pragma omp for
    for(int k=0;k<n_grid2;k++)
      e_fld->field_phi[n_grid1-1][k]=E_fi_upper;
  }
}
