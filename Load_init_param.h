#pragma once

#include<iostream>
#include "field.h"
#include "E_field.h"
#include "H_field.h"
#include "pdp3_time.h"
#include "Particles.h"
#include "Fourier.h"
#include "Poisson.h"
#include "Poisson_neumann.h"
#include "Poisson_dirichlet.h"
#include "particles_list.h"
#include <math.h>
#include "Boundary_Maxwell_conditions.h"
#include "input_output_class.h"
#include "Beam.h"
#include "Bunch.h"
#include "time.h"
#include "particles_struct.h"
#include "system_host.cuh"

#include "tinystr.h"
#include "tinyxml.h"
#include <fstream>
#include<iostream>
using namespace std;
class Load_init_param
{
public:
	Load_init_param(void);
	Load_init_param(char* xml_file_l);
	~Load_init_param(void);
public:
	void read_xml( );
	int get_int_value(const char* r_str);
	double get_double_value(const char* r_str);
	char* read_char(char* p_name);
	double* read_double_params(const char* p_name);
	void read_load_particles();
	Bunch*  read_load_bunch();
	void load_system();
	bool SaveSystemState(void);
	void Run(void);

public:
	char* xml_file;

		PML * c_pml;
		Geometry * c_geom;
		Time * c_time;
		Particles* c_part;
    	Bunch *     c_bunch;
		particles_list* p_list;
		E_field* efield;
		H_field* hfield;
		input_output_class * c_io_class;
		charge_density * c_rho_new;
		charge_density * c_rho_old;
		charge_density * c_rho_beam;
		current * c_current;
};

