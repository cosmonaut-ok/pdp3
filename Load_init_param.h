#pragma once

#include<iostream>
#include "field.h"
#include "E_field.h"
#include "H_field.h"
#include "pdp3_time.h"
#include "Particles.h"
#include "Fourier.h"
#include "Poisson.h"
#include "poisson_neumann.h"
#include "poisson_dirichlet.h"
#include "particles_list.h"
#include <math.h>
#include "Boundary_Maxwell_conditions.h"
#include "input_output_class.h"
#include "Beam.h"
#include "Bunch.h"
#include "time.h"
#include "particles_struct.h"

// #include "tinystr.h"
#include "tinyxml2.h"
#include <fstream>
#include<iostream>

using namespace std;
using namespace tinyxml2;

class Load_init_param
{
public:
	Load_init_param(void);
	Load_init_param(char* xml_file_name);
	~Load_init_param(void);

public:
	char* read_char(char* p_name);
	double* read_double_params(const char* p_name);
	bool save_system_state(void);
	void run(void);

public:
	PML * c_pml;
	Geometry * c_geom;
	Time * c_time;
	Particles* c_part;
	Bunch * c_bunch;
	particles_list* p_list;
	E_field* efield;
	H_field* hfield;
	input_output_class * c_io_class;
	charge_density * c_rho_new;
	charge_density * c_rho_old;
	charge_density * c_rho_beam;
	current * c_current;
	XMLDocument* xml_data;
	int testnumber;

private:
	void read_xml(const char* xml_file_name);
	void init_pml ();
	void init_geometry ();
	void init_fields ();
	void init_time ();
	void init_particles ();
	void init_boundary_maxwell ();
	Bunch* init_bunch();
};
