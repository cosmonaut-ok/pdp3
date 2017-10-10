#pragma once
#include"Particles.h"
#include <vector>
#include"charge_density.h"
#include"E_field.h"
#include"H_field.h"
using namespace std;
class Particles;
class particles_list
{
public:
	particles_list(int i);
	~particles_list(void);
public:
	vector<Particles*> part_list;
	flcuda** x1_old;
	flcuda** x3_old;

public:
	void charge_weighting(charge_density* ro1);
	void step_v(E_field *e_fld, H_field *h_fld, Time* t);
	void half_step_coord(Time* t);
	//void set_j_0();
	//void set_v_0();
	//void set_x_0();
	//void velocity_distribution(flcuda therm_vel);
	//void load_spatial_distribution(flcuda n1, flcuda n2);
	//void load_velocity_distribution(flcuda v_thermal);
	void j_weighting(Time* time1, current *j1);
	void azimuthal_j_weighting(Time* time1, current *j1);
	void create_coord_arrays(void);
	void copy_coords(void);


};
