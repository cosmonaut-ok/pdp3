#pragma once

#include <fstream>
#include <iostream>
#include <math.h>

#include "field.h"
#include "eField.h"
#include "hField.h"
#include "pdp3Time.h"
#include "particles.h"
#include "fourier.h"
#include "poisson.h"
#include "poissonNeumann.h"
#include "poissonDirichlet.h"
#include "particlesList.h"
#include "boundaryMaxwellConditions.h"
#include "inputOutputClass.h"
#include "bunch.h"
#include "time.h"

#include "tinyxml2.h"

using namespace std;
using namespace tinyxml2;

class LoadInitParam
{
public:
  LoadInitParam(void);
  LoadInitParam(char *xml_file_name);
  ~LoadInitParam(void);

public:
  double *read_double_params(const char *p_name);
  void dump_system_state();
  void dump_data(int step_number);
  void run(void);

public:
  Geometry *c_geom;
  Time *c_time;
  Particles *c_part;
  Bunch *c_bunch;
  ParticlesList *p_list;
  EField *efield;
  HField *hfield;
  InputOutputClass *c_io_class;
  ChargeDensity *c_rho_new;
  ChargeDensity *c_rho_old;
  ChargeDensity *c_rho_bunch;
  Current *c_current;
  XMLDocument *xml_data;
  int testnumber;
  int data_dump_interval;
  int system_state_dump_interval;
  int frames_per_file;
  // dump data
  bool is_dump_e_r;
  bool is_dump_e_phi;
  bool is_dump_e_z;
  bool is_dump_h_r;
  bool is_dump_h_phi;
  bool is_dump_h_z;
  bool is_dump_rho_bunch;

private:
  void read_xml(const char *xml_file_name);
  void init_geometry ();
  void init_fields ();
  void init_time ();
  void init_particles ();
  void init_boundary ();
  void init_bunch();
  void init_file_saving_parameters ();
  bool to_bool (string str);
};
