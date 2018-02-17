#pragma once
#include "geometry.h"
//#include "eField.h"
#include "pdp3_time.h"
#include "Particles.h"
#include "Fourier.h"
#include "current.h"
#include "Triple.h"
#include "particles_struct.h"

class EField;

class HField
{
public:
  double** h1;
  double** h2;
  double** h3;
  double* h1_1d;
  double* h2_1d;
  double* h3_1d;
  double** h1_half_time;
  double** h2_half_time;
  double** h3_half_time;
  double** Ar;
  double** Afi;
  double** Az;
  const Geometry* geom1;
  HField(Geometry* geom1);
  HField(void);
  ~HField(void);
  void calc_field(EField* e_field1, Time* time1);
  void set_homogeneous_h(double H1, double H2, double H3);
  Triple get_field(double x1, double x3);
  double* get_1d_h1();
  double* get_1d_h2();
  double* get_1d_h3();
};
