#include "poisson.h"
#include "constant.h"

class PoissonNeumann :
  public Poisson
{
public:
  PoissonNeumann(Geometry *cyl_geom);
  ~PoissonNeumann(void);
public:
  double **t_charge_density;
  void poisson_solve(EField *input_e, ChargeDensity *input_rho);
};
