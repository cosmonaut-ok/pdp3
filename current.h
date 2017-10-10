#pragma once
#include"Geometry.h"
class current
{
public:
	Geometry* geom1;
	current(void);
	current(Geometry* geom1);
	~current(void);
	flcuda** get_j1() const;
	flcuda** get_j2() const;
	flcuda** get_j3() const;

	flcuda* get_j1_1d() const;
	flcuda* get_j2_1d() const;
	flcuda* get_j3_1d() const;

	void j1_add_1d(flcuda *input);
	void j2_add_1d(flcuda *input);
	void j3_add_1d(flcuda *input);

	void set_j1(int i, int k, flcuda value);
	void set_j2(int i, int k, flcuda value);
	void set_j3(int i, int k, flcuda value);
	void reset_j();
protected:
	flcuda** j1;
	flcuda** j2;
	flcuda** j3;
	flcuda* j1_1d;
	flcuda* j2_1d;
	flcuda* j3_1d;
};
