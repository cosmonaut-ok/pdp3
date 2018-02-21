#pragma once
#include "pdp3Time.h"
#include "particles.h"
#include "fourier.h"
#include "current.h"
#include "triple.h"
#include "field.h"

class EField;

class HField : public Field
{
public:
  double** field_r_half_time;
  double** field_phi_half_time;
  double** field_z_half_time;
  double** Ar;
  double** Afi;
  double** Az;

  HField(Geometry* geom1);
  HField(void);
  ~HField(void);
  void calc_field(EField* e_field1, Time* time1);
  void set_homogeneous_h(double E1, double E2, double E3);
  Triple get_field(double x1, double x3);
};
