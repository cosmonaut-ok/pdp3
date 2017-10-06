#pragma once
#include "H_field.h"
#include "E_field.h"
#include "charge_density.h"
class Poisson
{
public:
	Poisson(Geometry* cyl_geom);
	~Poisson(void);
public:
	const flcuda epsilon0;
	Geometry* cyl_geom;
public:
	virtual void poisson_solve(E_field* input_e, charge_density* ro1)=0;
	bool test_poisson_equation(E_field* input_e, charge_density* input_rho);
};
