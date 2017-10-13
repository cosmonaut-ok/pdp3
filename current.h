#pragma once
#include "Geometry.h"

class current
{
public:
	Geometry* geom1;
	current(void);
	current(Geometry* geom1);
	~current(void);
	double** get_j1() const;
	double** get_j2() const;
	double** get_j3() const;

	double* get_j1_1d() const;
	double* get_j2_1d() const;
	double* get_j3_1d() const;

	void j1_add_1d(double *input);
	void j2_add_1d(double *input);
	void j3_add_1d(double *input);

	void set_j1(int i, int k, double value);
	void set_j2(int i, int k, double value);
	void set_j3(int i, int k, double value);
	void reset_j();
protected:
	double** j1;
	double** j2;
	double** j3;
	double* j1_1d;
	double* j2_1d;
	double* j3_1d;
};
