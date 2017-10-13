#pragma once
#include "Poisson.h"

class Poisson_neumann :
	public Poisson
{
public:
	Poisson_neumann(Geometry* cyl_geom);
	~Poisson_neumann(void);
public:
	double** t_charge_density;
	void poisson_solve(E_field* input_e, charge_density* input_rho);
};
