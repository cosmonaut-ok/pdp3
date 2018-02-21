#pragma once
#include "geometry.h"

class Field
{
public:
  double **field_r;
  double **field_phi;
  double **field_z;

  double *field_r_1d;
  double *field_phi_1d;
  double *field_z_1d;

	const Geometry *geom1;

  Field(Geometry* geom1);
  Field(void);

  ~Field(void);

  double* get_1d_field_r();
  double* get_1d_field_phi();
  double* get_1d_field_z();
};
