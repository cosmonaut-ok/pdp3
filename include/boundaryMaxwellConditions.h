#pragma once
#include "eField.h"
#include "current.h"

class BoundaryMaxwellConditions
{
public:
  EField *e_fld;

public:
  BoundaryMaxwellConditions(void);
  BoundaryMaxwellConditions(EField *e_fld);
  ~BoundaryMaxwellConditions(void);
  void specify_initial_field(Geometry *cyl_geom,double E_fi_upper, double E_fi_left, double E_fi_right);
};
