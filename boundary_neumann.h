#pragma once
#include "Boundary_conditions.h"
#include "charge_density.h"

class boundary_neumann :
  public Boundary_conditions
{
public:
  boundary_neumann(E_field* e_f ,charge_density* rho, Geometry *cyl_geom);
  boundary_neumann(void);
  ~boundary_neumann(void);
public:
  Geometry *cyl_geom;
  E_field * e_f;
  charge_density *rho;
  void specify_boundary_conditions(double E_fi_upper,
                                   double E_fi_left,
                                   double E_fi_right,
                                   double fi_upper_wall);
};
