#pragma once
#include "PML.h"
#include "particles_struct.h"

class PML;

class Geometry
{
public:
	double first_size;
	double second_size;
	int n_grid_1;
	int n_grid_2;
	double dr;
	double dz;
	double** epsilon;
	double** sigma;
	PML* pml1;
	// Geometry* geom1;
	double set_dr();
	double set_dz();
	// #num geometry_type;

	void set_epsilon();
//	friend void calc_sigma(Geometry *geom1);
	Geometry(double fs, double ss, int ng1, int ng2, PML* pml1_t);
	Geometry(double * param, PML* pml1_t);
	Geometry(double fs, double ss, int ng1, int ng2);
	Geometry( PML* pml1_t);
	Geometry();

	~Geometry(void);
};
