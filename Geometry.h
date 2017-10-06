#pragma once
#include "PML.h"
#include "particles_struct.h"
class PML;
class Geometry
{
public:
	flcuda first_size;
	flcuda second_size;
	int n_grid_1;
	int n_grid_2;
	flcuda dr;
	flcuda dz;
	flcuda** epsilon;
	flcuda** sigma;
	PML* pml1;
//	Geometry* geom1;
	flcuda set_dr();
	flcuda set_dz();
//	#num geometry_type;

	void set_epsilon();
//	friend void calc_sigma(Geometry *geom1);
	Geometry(flcuda fs, flcuda ss,  int ng1, int ng2, PML* pml1_t);
	Geometry(double * param, PML* pml1_t);
	Geometry(flcuda fs, flcuda ss,  int ng1, int ng2);
	Geometry( PML* pml1_t);
	Geometry();

	~Geometry(void);
};
