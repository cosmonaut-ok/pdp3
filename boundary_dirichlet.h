#pragma once
#include "boundary_conditions.h"
#include "charge_density.h"

class boundary_dirichlet :
	public Boundary_conditions
{
public:
	boundary_dirichlet(E_field* e_f ,charge_density* rho, Geometry *cyl_geom);
	boundary_dirichlet(void);
	~boundary_dirichlet(void);
public:
	 Geometry *cyl_geom;
	E_field * e_f;
	charge_density *rho;
	void specify_boundary_conditions(flcuda E_fi_upper, flcuda E_fi_left, flcuda E_fi_right, flcuda fi_upper_wall);
};
