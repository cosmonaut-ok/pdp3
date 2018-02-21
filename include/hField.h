#pragma once
// #include "geometry.h"
// #include "eField.h"
#include "pdp3Time.h"
#include "particles.h"
#include "fourier.h"
#include "current.h"
#include "triple.h"
#include "field.h"

class EField;
// class Field;

class HField : public Field
{
public:
  double** field_r;
  double** field_phi;
  double** field_z;
  double* field_r_1d;
  double* field_phi_1d;
  double* field_z_1d;
  double** field_r_half_time;
  double** field_phi_half_time;
  double** field_z_half_time;
  double** Ar;
  double** Afi;
  double** Az;
  const Geometry* geom1;
  HField(Geometry* geom1);
  HField(void);
  ~HField(void);
  void calc_field(EField* e_field1, Time* time1);
  void set_homogeneous_h(double FIELD_R, double FIELD_PHI, double FIELD_Z);
  Triple get_field(double x1, double x3);
  double* get_1d_field_r();
  double* get_1d_field_phi();
  double* get_1d_field_z();
};
