#pragma once
#include"Geometry.h"

class charge_density
{
public:
	Geometry* geom1;
	flcuda** get_ro() const;
	void set_ro_weighting(int i, int k, flcuda value);
	charge_density(void);
	charge_density(Geometry* geom1_t);
	void reset_rho();
	
	~charge_density(void);
protected:
		flcuda** ro;
};
