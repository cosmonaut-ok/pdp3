#pragma once
#include "pdp3Time.h"
#include "chargeDensity.h"
#include "current.h"
#include "triple.h"
#include "field.h"

class HField;

class EField : public Field
{
public:
  double** fi; //potential
  double** t_charge_density;

  EField(Geometry* geom1_t);
  EField();
  ~EField(void);

  void calc_field(HField* h_field1, Time* timfield_r, Current* current1);
  void poisson_equation2(Geometry* geom1, ChargeDensity* ro1);
  void cosine_ftansfrom(double** fi_ro, int lenght_n, int k);
  void set_homogeneous_efield(double FIELD_R, double FIELD_PHI, double FIELD_Z);
  void set_fi_on_z();
  void boundary_conditions();
  Triple get_field(double x1, double x3);
  bool test_poisson_equation(ChargeDensity* rho);
  void tridiagonal_solve(const double *a,
                         const double *b,
                         double *c,
                         double *d,
                         double *x,
                         int n);
};
