#pragma once
#include "Geometry.h"
#include "particles_struct.h"

class Geometry;

class PML
{
public:
	double comparative_l_1;
	double comparative_l_2;
	double comparative_l_3;
	double sigma1;
	double sigma2;

	void calc_sigma(Geometry* geom1);
	PML(double comp_l_1, double comp_l2, double comp_l3, double sigma1_t, double sigma2_t);
	PML(double* pml_params);

	PML(void);
	~PML(void);
};
