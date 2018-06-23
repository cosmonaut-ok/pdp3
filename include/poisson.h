#include <fstream>

#include "math/fourier.h"

#include "constant.h"

#include "eField.h"
#include "geometry.h"

class Poisson
{
public:
  Poisson(Geometry *cyl_geom);
  ~Poisson(void);
public:
  Geometry *cyl_geom;
public:
  virtual void poisson_solve(EField *input_e, ChargeDensity *ro1) = 0;
  bool test_poisson_equation(EField *input_e, ChargeDensity *input_rho);
};
