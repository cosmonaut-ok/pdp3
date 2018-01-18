#pragma once
#include "field.h"
#include "Geometry.h"
#include "pdp3_time.h"
//#include "H_field.h"
#include "Particles.h"
#include"PML.h"
#include"charge_density.h"
#include"current.h"
#include"Triple.h"
#include "particles_struct.h"

class H_field;
class E_field
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
  const double epsilon0;
  E_field(Geometry* geom1_t);
  E_field();
  ~E_field(void);
  void calc_field(H_field* h_field1, Time* time1, current* current1, PML* pml1);
  void calc_field(H_field* h_field1, Time* time1, current* current1);
  void poisson_equation2(Geometry* geom1, charge_density* ro1);
  void cosine_ftansfrom(double** fi_ro, int lenght_n, int k);
  void set_homogeneous_efield(double E1, double E2, double E3);
  void set_fi_on_z();
  void boundary_conditions();
  Triple get_field(double x1, double x3);
  bool test_poisson_equation(charge_density* rho);
  void TridiagonalSolve(const double *a,
                        const double *b,
                        double *c,
                        double *d,
                        double *x,
                        int n);
  double* get_1d_e1();
  double* get_1d_e2();
  double* get_1d_e3();
};
