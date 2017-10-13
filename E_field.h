#pragma once
#include "field.h"
#include "Geometry.h"
#include "pdp3_time.h"
//#include "H_field.h"
#include "Particles.h"
#include"PML.h"
#include"charge_density.h"
#include"current.h"
#include"Triple.h"
#include "particles_struct.h"

class H_field;
class E_field
{

public:
	flcuda** e1; //Er
	flcuda** e2; //Ef
	flcuda** e3; //Ez
	flcuda* e1_1d; //Er
	flcuda* e2_1d; //Ef
	flcuda* e3_1d; //Ez
	flcuda** fi; //potential
	flcuda** t_charge_density;
	const Geometry* geom1;
	const flcuda epsilon0;
	E_field(Geometry* geom1_t);
	E_field();
	~E_field(void);
	void calc_field(H_field* h_field1, Time* time1, current* current1, PML* pml1);
	void calc_field(H_field* h_field1, Time* time1, current* current1);
	void poisson_equation2(Geometry* geom1, charge_density* ro1);
	void cosine_ftansfrom(flcuda** fi_ro, int lenght_n, int k);
	void set_homogeneous_efield(flcuda E1, flcuda E2, flcuda E3);
	void set_fi_on_z();
	void boundary_conditions();
	Triple get_field(flcuda x1, flcuda x3);
	bool test_poisson_equation(charge_density* rho);
	void TridiagonalSolve(const flcuda *a,
                        const flcuda *b,
                        flcuda *c,
                        flcuda *d,
                        flcuda *x,
                        unsigned int n);
	flcuda* get_1d_e1();
	flcuda* get_1d_e2();
	flcuda* get_1d_e3();
};
