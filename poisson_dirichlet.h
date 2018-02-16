#pragma once
#include "Poisson.h"

class Poisson_dirichlet :
  public Poisson
{
public:
  Poisson_dirichlet(Geometry* cyl_geom);
  ~Poisson_dirichlet(void);
public:
  double** t_charge_density;
  void poisson_solve(EField* input_e, ChargeDensity* input_rho);
};
