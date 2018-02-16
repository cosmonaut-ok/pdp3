#pragma once
#include "Poisson.h"

class Poisson_neumann :
  public Poisson
{
public:
  Poisson_neumann(Geometry* cyl_geom);
  ~Poisson_neumann(void);
public:
  double** t_charge_density;
  void poisson_solve(EField* input_e, ChargeDensity* input_rho);
};
