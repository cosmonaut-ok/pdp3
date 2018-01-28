#pragma once

#include <fstream>
#include <iostream>
#include <math.h>

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
#include "Boundary_Maxwell_conditions.h"
#include "input_output_class.h"
#include "Bunch.h"
#include "time.h"
#include "particles_struct.h"

#include "tinyxml2.h"

using namespace std;
using namespace tinyxml2;

class Load_init_param
{
public:
  Load_init_param(void);
  Load_init_param(char* xml_file_name);
  ~Load_init_param(void);

public:
  double* read_double_params(const char* p_name);
  bool save_system_state(double t);
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
  int data_dump_interval;
  int system_state_dump_interval;
  int frames_per_file;
  // dump data
  bool is_dump_e1;
  bool is_dump_e2;
  bool is_dump_e3;
  bool is_dump_h1;
  bool is_dump_h2;
  bool is_dump_h3;
  bool is_dump_rho_beam;

private:
  void read_xml(const char* xml_file_name);
  void init_pml ();
  void init_geometry ();
  void init_fields ();
  void init_time ();
  void init_particles ();
  void init_boundary ();
  void init_bunch();
  void init_file_saving_parameters ();
  bool to_bool (string str);
};
