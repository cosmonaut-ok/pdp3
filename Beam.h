#pragma once
#include "Particles.h"

class Beam :
	public Particles
{
public:
	Beam(char* p_name,
			 double p_charge,
			 double p_mass,
			 int p_number,
			 Geometry* geom,
			 particles_list* p_list,
			 double b_radius);

	~Beam(void);

public:
	double radius; //beam radius
	double n_beam;// beam density;
	double vel_beam; //beam velocity

public:
	void calc_init_param(Time* time,int particles_in_step,double n_b,double b_vel);
	void beam_inject(Time* time, int particles_in_step,double fi,double koef);
	virtual void	half_step_coord(Time* t);
};
