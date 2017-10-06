#pragma once
#include"Geometry.h"
#include "particles_struct.h"
class Geometry;
class PML
{
public:
	flcuda comparative_l_1;
	flcuda comparative_l_2;
	flcuda comparative_l_3;
	flcuda sigma1;
	flcuda sigma2;

	void calc_sigma(Geometry* geom1);
	PML(flcuda comp_l_1, flcuda comp_l2, flcuda comp_l3, flcuda sigma1_t, flcuda sigma2_t);
	PML(double* pml_params);

	PML(void);
	~PML(void);  
};
