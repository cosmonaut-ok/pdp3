#pragma once
#include "Geometry.h"

class charge_density
{
public:
  Geometry* geom1;
  double** get_ro() const;
  void set_ro_weighting(int i, int k, double value);
  charge_density(void);
  charge_density(Geometry* geom1_t);
  void reset_rho();

  ~charge_density(void);
protected:
  double** rho;
};
