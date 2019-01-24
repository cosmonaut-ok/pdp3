#pragma once

#include <math.h>
#include "particles.h"
#include "geometry.h"
#include "constant.h"
#include "lib.h"

using namespace constant;

class Temperature
{
public:
  double **t_r;
  double **t_phi;
  double **t_z;
  Geometry *geom;

private:
  double **t_r_sum;
  double **t_phi_sum;
  double **t_z_sum;
  double **count;

public:
  Temperature(Geometry *geom1);
  Temperature(void);

  ~Temperature(void);

  void calc_t_r(Particles *prtls);
  void calc_t_phi(Particles *prtls);
  void calc_t_z(Particles *prtls);
  void reset(void);

private:
  void inc_count(unsigned int i, unsigned int j);
  void inc_sum(unsigned int i, unsigned int j, double t);
};

#define VEL_TO_TEMPR(vel, mass) (mass) * pow(lib::get_gamma(vel) * vel, 2) / 2 / BOLTZMANN
