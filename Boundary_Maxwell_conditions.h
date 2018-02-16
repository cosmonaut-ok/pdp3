#pragma once
#include "E_field.h"
#include "current.h"

class Boundary_Maxwell_conditions
{
public:
  E_field* e_fld;

public:
  Boundary_Maxwell_conditions(void);
  Boundary_Maxwell_conditions(E_field* e_fld);
  ~Boundary_Maxwell_conditions(void);
  void specify_initial_field(Geometry* cyl_geom,double E_fi_upper, double E_fi_left, double E_fi_right);
};
