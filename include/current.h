#pragma once

#include <string.h>
#include "geometry.h"

class Current
{
public:
  Geometry *geom1;
  Current(void);
  Current(Geometry *geom1);
  ~Current(void);
  double **get_j_r() const;
  double **get_j_phi() const;
  double **get_j_z() const;

  void inc_j_r(int i, int k, double value);
  void inc_j_phi(int i, int k, double value);
  void inc_j_z(int i, int k, double value);
  void reset_j();

protected:
  double **j_r;
  double **j_phi;
  double **j_z;
  double *j_r_1d;
  double *j_phi_1d;
  double *j_z_1d;
};
