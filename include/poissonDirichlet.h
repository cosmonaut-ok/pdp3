#include "poisson.h"

class PoissonDirichlet :
  public Poisson
{
public:
  PoissonDirichlet(Geometry *cyl_geom);
  ~PoissonDirichlet(void);
public:
  double **t_charge_density;
  void poisson_solve(EField *input_e, ChargeDensity *input_rho);
};
