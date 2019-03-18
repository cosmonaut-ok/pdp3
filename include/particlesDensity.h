#pragma once

#include <math.h>
#include "particles.h"
#include "geometry.h"
#include "constant.h"
#include "lib.h"

using namespace constant;

class ParticlesDensity
{
public:
  double **density;

private:
  Geometry *geom;

public:
  ParticlesDensity(Geometry *geom1);
  ParticlesDensity(void);

  ~ParticlesDensity(void);

  void calc_density(Particles *prtls);
  void reset(void);
  void inc_count(unsigned int r, unsigned int z, double value);
};
