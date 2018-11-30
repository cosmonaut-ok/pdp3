#include <fstream>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <vector>

// enable openmp optional
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

#include "config.h"

#include "math/fourier.h"

#include "lib.h"
#include "boundaryMaxwellConditions.h"
#include "bunch.h"
#include "current.h"
#include "geometry.h"
#include "hField.h"
#include "ioText.h"

#include "writer.h"
#include "writerPlain.h"

#ifdef USE_HDF5
#include "ioHDF5.h"
#endif

#include "particles.h"
#include "particlesList.h"
#include "pdp3Time.h"
#include "poissonDirichlet.h"
#include "parameters.h"

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
  void dump_data(int step_number);
  void print_data(int writer_type, char* component, int step_number, int dump_interval, int* shape);
  void run(void);

public:
  // Geometry *c_geom;
  // Time *c_time;
  Particles *c_part;
  vector<Bunch*> c_bunches;

#ifndef LEGACY
  vector<Writer*> c_writers;
#endif

  unsigned int current_bunch_number;
  ParticlesList *p_list;
  EField *efield;
  HField *hfield;

#ifndef LEGACY
  Writer *c_writer_e_r;
  Writer *c_writer_e_phi;
  Writer *c_writer_e_z;
  Writer *c_writer_h_r;
  Writer *c_writer_h_phi;
  Writer *c_writer_h_z;
  Writer *c_writer_j_r;
  Writer *c_writer_j_phi;
  Writer *c_writer_j_z;
  Writer *c_writer_rho_beam;
#else
  InputOutputClass *c_io_class;
#endif
  ChargeDensity *c_rho_new;
  ChargeDensity *c_rho_old;
  ChargeDensity *c_rho_bunch;
  Current *c_current;
  XMLDocument *xml_data;
  Parameters *params;
  int testnumber;
  int data_dump_interval;
  int system_state_dump_interval;
  int frames_per_file;

  int printed_step = 0;
  // dump data
  bool is_dump_e_r;
  bool is_dump_e_phi;
  bool is_dump_e_z;
  bool is_dump_h_r;
  bool is_dump_h_phi;
  bool is_dump_h_z;
  bool is_dump_rho_bunch;

  bool is_debug;

private:
  void read_xml(const char *xml_file_name);
  void init_geometry ();
  void init_fields ();
  void init_time ();
  void init_particles ();
  void init_boundary ();
  void init_beam();
  void init_file_saving_parameters ();
};
