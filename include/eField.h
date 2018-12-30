#pragma once

#include <fstream>
#include "tinyvec3d.h"

#include "math/fourier.h"

#include "constant.h"

#include "geometry.h"
#include "pdp3Time.h"
#include "chargeDensity.h"
#include "current.h"
#include "field.h"
#include "hField.h"

class HField;

class EField : public Field
{
public:
  double **fi; //potential
  double **t_charge_density;

  EField(Geometry *geom1_t);
  EField();
  ~EField(void);

  void calc_field(HField *h_field1, Time *time1, Current *current1);
  void poisson_equation2(Geometry *geom1, ChargeDensity *ro1);
  void cosine_ftansfrom(double **fi_ro, int lenght_n, int k);
  void set_homogeneous_efield(double E_r, double E_phi, double E_z);
  void set_fi_on_z();
  void boundary_conditions();
  double* get_field(double radius, double longitude);
  bool test_poisson_equation(ChargeDensity *rho);
  void tridiagonal_solve(const double *a,
                         const double *b,
                         double *c,
                         double *d,
                         double *x,
                         int n);
};
