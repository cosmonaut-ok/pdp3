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
  double **t;

private:
  Geometry *geom;
  double **t_sum;
  double **t_vec_r;
  double **t_vec_z;
  double **t_vec_phi;
  double **masses;
  double **count;
  double **count_sum;

public:
  Temperature(Geometry *geom1);
  Temperature(void);

  ~Temperature(void);

  void calc_t(Particles *prtls);
  void reset(void);
  void inc_tmpr(unsigned int r, unsigned int z, double value);
  void inc_vec_r(unsigned int r, unsigned int z, double value);
  void inc_vec_z(unsigned int r, unsigned int z, double value);
  void inc_vec_phi(unsigned int r, unsigned int z, double value);
  void inc_count(unsigned int r, unsigned int z, double value);
  void normalize(Particles *prtls);

private:
  void inc_count(unsigned int i, unsigned int j);
  void inc_sum(unsigned int i, unsigned int j, double t);
};
