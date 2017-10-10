#pragma once
#include "E_field.h"
#include"current.h"
class Boundary_Maxwell_conditions
{
public:
E_field* e_fld;

public:
	Boundary_Maxwell_conditions(void );
	Boundary_Maxwell_conditions(E_field* e_fld );
	~Boundary_Maxwell_conditions(void);
	void specify_initial_field(Geometry* cyl_geom,flcuda E_fi_upper, flcuda E_fi_left, flcuda E_fi_right);
	void radiation_source(Geometry* cyl_geom, flcuda region_size, flcuda frequency, int type, flcuda time);
	void probe_mode_exitation(Geometry* cyl_geom, current* j_input, flcuda probe_lenght, flcuda frequency,flcuda time);
};
