#pragma once

#include "geometry.h"

class ChargeDensity
{
public:
  Geometry *geom1;
  double **get_rho() const;
  void inc_rho(int i, int k, double value);
  ChargeDensity(void);
  ChargeDensity(Geometry *geom1_t);
  void reset_rho();

  ~ChargeDensity(void);
protected:
  double **rho;
};
