#pragma once
#include "geometry.h"
#include "pdp3Time.h"
#include "particles.h"
#include"chargeDensity.h"
#include"current.h"
#include"triple.h"

class HField;
class EField
{

public:
  double** e1; //Er
  double** e2; //Ef
  double** e3; //Ez
  double* e1_1d; //Er
  double* e2_1d; //Ef
  double* e3_1d; //Ez
  double** fi; //potential
  double** t_charge_density;
  const Geometry* geom1;

  EField(Geometry* geom1_t);
  EField();
  ~EField(void);

  void calc_field(HField* h_field1, Time* time1, Current* current1);
  void poisson_equation2(Geometry* geom1, ChargeDensity* ro1);
  void cosine_ftansfrom(double** fi_ro, int lenght_n, int k);
  void set_homogeneous_efield(double E1, double E2, double E3);
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
  double* get_1d_e1();
  double* get_1d_e2();
  double* get_1d_e3();
};
