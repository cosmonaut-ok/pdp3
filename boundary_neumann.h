#pragma once
#include "boundary_conditions.h"
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
	void specify_boundary_conditions(flcuda E_fi_upper, flcuda E_fi_left, flcuda E_fi_right, flcuda fi_upper_wall);
};
