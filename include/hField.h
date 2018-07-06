#pragma once

#include "tinyvec3d.h"

#include "constant.h"

#include "field.h"
#include "pdp3Time.h"
#include "eField.h"

class EField;

class HField : public Field
{
public:
  double **field_r_half_time;
  double **field_phi_half_time;
  double **field_z_half_time;
  double **Ar;
  double **Afi;
  double **Az;

  HField(Geometry *geom1);
  HField(void);
  ~HField(void);
  void calc_field(EField *e_field1, Time *time1);
  void set_homogeneous_h(double E_r, double E_phi, double E_z);
  double* get_field(double x1, double x3);

  double *get_1d_field_r();
  double *get_1d_field_phi();
  double *get_1d_field_z();
};
