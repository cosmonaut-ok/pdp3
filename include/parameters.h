//! Class for reading configfile and calculate parameters,
//! requires for modeling

#include <vector>

#include "tinyxml2.h"

#include "lib.h"
#include "geometry.h"
#include "pdp3Time.h"

using namespace std;
using namespace tinyxml2;

struct particle_specie {
  char *name;
  unsigned int mass;
  int charge;
  double number_macro;
  double left_density;
  double right_density;
  double temperature;
};

class Parameters
{
public:
  // Parameters(void);
  Parameters(const char *xml_file_name);
  ~Parameters(void);

public:
  /* double *read_double_params(const char *p_name); */
  /* void dump_system_state(); */
  /* void dump_data(int step_number); */
  /* void run(void); */

public:

  bool debug = false;
  bool use_hdf5 = false;

  Geometry *geom;

  /* <pml> */
  double pml_comparative_l_1; // 0</comparative_l_1>
  double pml_comparative_l_2; // >0</comparative_l_2>
  double pml_comparative_l_3; // >0</comparative_l_3>
  double pml_sigma_1; // >1e-5</sigma_1>
  double pml_sigma_2; // >7e-2</sigma_2>


  /* <time> */
  Time *time;

  vector<particle_specie> particle_species;

  //! use "-particle-" infix for physical particles
  //! use "-macro-" infix for macro- or superparticles
  //! use "beam" for full set of electron bunches
  //! use "bunch" for sinble electron bunch (part of beam)

  char *beam_name;
  //! particles charge in beam
  //! (mean, all bunches in the beam and charge/mass of all particles
  //! in the bunch are equal.
  int beam_particle_charge;
  //! particles mass in beam
  unsigned int beam_particle_mass;
  //! number of particles in beam
  unsigned int beam_number_bunches;
  //! distance between bunches in beam
  double beam_bunches_distance;
  //! initial velocity of bunches in beam
  //! (mean, it is equal for all bucnhes in beam)
  double beam_initial_velocity;
  //! number of macroparticles in bunch
  unsigned int bunch_number_macro;
  //! length of bunch
  double bunch_lenght;
  //! radius of bunch
  double bunch_radius;
  //! particles density in bunch
  double bunch_density;

  double boundary_maxwell_e_phi_upper;
  double boundary_maxwell_e_phi_left;
  double boundary_maxwell_e_phi_right;
  int boundary_conditions;

  char *dump_result_path;
  char *dump_save_state_path;
  unsigned int dump_data_interval;
  unsigned int dump_frames_per_file;
  unsigned int dump_system_state_interval;
  bool dump_compress = false;
  int dump_compress_level = 0;

  // dump different kinds of data
  bool dump_e_r = false;
  bool dump_e_phi = false;
  bool dump_e_z = false;
  bool dump_h_r = false;
  bool dump_h_phi = false;
  bool dump_h_z = false;
  bool dump_position = false;
  bool dump_velocity = false;
  bool dump_rho_beam = false;

  // electrical current dump
  bool dump_current_r = false;
  bool dump_current_phi = false;
  bool dump_current_z = false;

private:
  XMLElement* xml_data;
  XMLElement* try_first_child(XMLElement* element, const char* name);

  void init_particles();
  void init_beam();
  void init_geometry();
  void init_pml();
  void init_time();
  void init_boundary();
  void init_data_dump_parameters();

};
