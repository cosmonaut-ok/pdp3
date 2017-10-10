#pragma once
#include "Poisson.h"

class Poisson_dirichlet :
	public Poisson
{
public:
	Poisson_dirichlet(Geometry* cyl_geom);
	~Poisson_dirichlet(void);
public:
	flcuda** t_charge_density;
	void poisson_solve(E_field* input_e, charge_density* input_rho);
};
