// #include <sstream>
// #include <string>
// #include <iomanip>
#include <algorithm>
#include <iomanip>
// #include <cctype>
// enable openmp optional
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

#include "Load_init_param.h"
#include "time.h"
#include "tinyxml2.h"

#include <sys/time.h>
#include <ctime>

using namespace std;
using namespace tinyxml2;

#define INITIAL_PARAMS_NAME "Initial_parameters"
#define PML_PARAMS_NAME "PML"
#define GEOMETRY_PARAMS_NAME "geometry"
#define TIME_PARAMS_NAME "Time"
#define PARTICLES_PARAMS_NAME "Particles"
#define BUNCH_PARAMS_NAME "Particles_bunch"
#define BOUNDARY_MAXWELL_PARAMS_NAME "Boundary_Maxwell_conditions"
#define FILE_SAVE_PARAMS_NAME "file_save_parameters"

Load_init_param::Load_init_param(void)
{
}

Load_init_param::Load_init_param(char* xml_file_name)
{
    // NOTE: all system init is too huge,
    // so we use several "subconstructors"
    // to initialise different parts of system

  // define openmp-related options when openmp enabled
#ifdef _OPENMP
  omp_set_dynamic(0); // Explicitly disable dynamic teams
#endif

  // read XML config file
  cout << "Reading configuration file ``" << xml_file_name << "``\n";
  read_xml(xml_file_name);

  // load PML parameters
  cout << "Initialising PML Data\n";
  init_pml();

  // load Geometry parameters
  cout << "Initialising Geometry Parameters\n";
  init_geometry();

  // creating field objects
  cout << "Initialising E/M Fields Data\n";
  init_fields ();

  // load time parameters
  cout << "Initialising Time Data\n";
  init_time ();

  // load particle parameters
  cout << "Initialising Particles Data\n";
  init_particles();

  // load bunch
  cout << "Initialising Particles Bunch Data\n";
  init_bunch();

  cout << "Initialising Bounrary Conditions Data\n";
  init_boundary();

  // load File Path
  cout << "Initialising File System Paths\n";
  init_file_saving_parameters();

  cout << "Initialisation complete\n";
}

Load_init_param::~Load_init_param(void)
{
}

void Load_init_param::read_xml(const char* xml_file_name)
{
  xml_data = new XMLDocument(true, COLLAPSE_WHITESPACE);

  XMLError e_result = xml_data->LoadFile(xml_file_name);
  if (e_result != XML_SUCCESS)
  {
    cerr << "ERROR: Can not read configuration file ``" << xml_file_name << "``\n";
    exit (78);
  }
}

bool Load_init_param::to_bool(string str) {
  std::transform(str.begin(), str.end(), str.begin(), ::tolower);
  std::istringstream is(str);
  bool b;
  is >> std::boolalpha >> b;
  return b;
}

void Load_init_param::init_particles()
{
  const char* p_king_section_name = "particle_kind";
  XMLElement* root = xml_data->FirstChildElement (INITIAL_PARAMS_NAME);
  XMLElement* particle_kind = root
    ->FirstChildElement (PARTICLES_PARAMS_NAME)
    ->FirstChildElement (p_king_section_name);

  char* p_name= new char [50];
  double charge = 0;
  double mass = 0; // TODO: should it be double?
  double number = 0;
  double left_density = 0;
  double right_density = 0;
  double temperature = 0;

  Particles* prtls = 0;

  // initialize particles list
  p_list = new particles_list();

  // creating rho and current arrays
  // WARNING! should be called after geometry initialised
  c_rho_new = new charge_density(c_geom);
  c_rho_old = new charge_density(c_geom);
  c_rho_beam = new charge_density(c_geom);
  c_current = new current(c_geom);

  while(particle_kind)
  {
    strcpy (p_name, particle_kind->FirstChildElement("name")->GetText());
    cout << "Initialising " << p_name << " Data\n";

    charge = atof(particle_kind->FirstChildElement("charge")->GetText());
    mass = atof(particle_kind->FirstChildElement("mass")->GetText());
    number = atof(particle_kind->FirstChildElement("number")->GetText());
    left_density = atof(particle_kind->FirstChildElement("left_density")->GetText());
    right_density = atof(particle_kind->FirstChildElement("right_density")->GetText());
    temperature = atof(particle_kind->FirstChildElement("temperature")->GetText());

    // init and setup particles properties
    prtls = new Particles(strcpy(new char [50], p_name), charge, mass, number, c_geom, p_list);
    prtls->load_spatial_distribution_with_variable_mass(left_density,right_density,0,0);
    prtls->velocity_distribution(temperature);
    particle_kind = particle_kind->NextSiblingElement(p_king_section_name);
  }
}

void Load_init_param::init_bunch()
{
  // initialise particles bunch data
  XMLElement* root = xml_data->FirstChildElement (INITIAL_PARAMS_NAME);
  XMLElement* sub_root =root->FirstChildElement (BUNCH_PARAMS_NAME);

  char* p_name = (char*)sub_root->FirstChildElement("name")->GetText();
  double charge = atof(sub_root
                       ->FirstChildElement("charge")
                       ->GetText());
  double mass = atof(sub_root
                     ->FirstChildElement("mass")
                     ->GetText());
  double number = atof(sub_root
                       ->FirstChildElement("number")
                       ->GetText());
  double duration = atof(sub_root
                         ->FirstChildElement("duration")
                         ->GetText());
  double radius = atof(sub_root
                       ->FirstChildElement("radius")
                       ->GetText());
  double density = atof(sub_root
                        ->FirstChildElement("density")
                        ->GetText());
  double initial_velocity = atof(sub_root
                                 ->FirstChildElement("initial_velocity")
                                 ->GetText());

  Bunch* prtls = new Bunch(p_name, charge, mass,number, c_geom, p_list,
                           duration, radius, density, initial_velocity);

  c_bunch = prtls;
}

void Load_init_param::init_boundary ()
{
  XMLElement* root = xml_data->FirstChildElement (INITIAL_PARAMS_NAME);
  XMLElement* sub_root =root->FirstChildElement (BOUNDARY_MAXWELL_PARAMS_NAME);
  //
  double e_fi_upper = atof(sub_root->
                           FirstChildElement("e_fi_upper")->
                           GetText());
  double e_fi_left = atof(sub_root->
                          FirstChildElement("e_fi_left")->
                          GetText());
  double e_fi_right = atof(sub_root->
                           FirstChildElement("e_fi_right")->
                           GetText());

  int boundary_conditions = atoi(root->FirstChildElement("Boundary_conditions")->GetText());

  // Maxwell initial conditions
  Boundary_Maxwell_conditions maxwell_rad(efield); // TODO: WTF?
  maxwell_rad.specify_initial_field(c_geom, e_fi_upper, e_fi_left, e_fi_right);

  if (boundary_conditions == 0)
  {
    p_list->charge_weighting(c_rho_new);
    Fourier four1();
    Poisson_dirichlet dirih(c_geom);
    dirih.poisson_solve(efield, c_rho_new);
  }
}

void Load_init_param::init_pml ()
{
  XMLElement* root = xml_data->FirstChildElement (INITIAL_PARAMS_NAME);
  XMLElement* sub_root =root->FirstChildElement (PML_PARAMS_NAME);
  //
  double comp_l_1 = atof(sub_root->
                         FirstChildElement("comparative_l_1")->
                         GetText());
  double comp_l_2 = atof(sub_root->
                         FirstChildElement("comparative_l_2")->
                         GetText());
  double comp_l_3 = atof(sub_root->
                         FirstChildElement("comparative_l_3")->
                         GetText());
  double sigma_1_t = atof(sub_root->
                          FirstChildElement("sigma_1")->
                          GetText());
  double sigma_2_t = atof(sub_root->
                          FirstChildElement("sigma_2")->
                          GetText());
  c_pml = new PML(comp_l_1, comp_l_2, comp_l_3, sigma_1_t, sigma_2_t);
}

void Load_init_param::init_geometry ()
{
  XMLElement* root = xml_data->FirstChildElement (INITIAL_PARAMS_NAME);
  XMLElement* sub_root =root->FirstChildElement (GEOMETRY_PARAMS_NAME);
  //
  double r_size = atof(sub_root->
                       FirstChildElement("r_size")->
                       GetText());
  double z_size = atof(sub_root->
                       FirstChildElement("z_size")->
                       GetText());
  int n_grid_r = atoi(sub_root->
                      FirstChildElement("n_grid_r")->
                      GetText());
  int n_grid_z = atoi(sub_root->
                      FirstChildElement("n_grid_z")->
                      GetText());

  // WARNING! should be called after init_pml
  c_geom = new Geometry(r_size, z_size, n_grid_r, n_grid_z, c_pml);
  c_geom->set_epsilon();
}

void Load_init_param::init_fields ()
{
  efield = new E_field(c_geom);
  hfield = new H_field(c_geom);
  efield->boundary_conditions();
  efield->set_homogeneous_efield(0.0, 0.0, 0.0);
  hfield->set_homogeneous_h(0.0, 0.0, 0.0);
  efield->set_fi_on_z();
}

void Load_init_param::init_time ()
{
  XMLElement* root = xml_data->FirstChildElement (INITIAL_PARAMS_NAME);
  XMLElement* sub_root =root->FirstChildElement (TIME_PARAMS_NAME);
  //
  double start_time = atof(sub_root->
                           FirstChildElement("start_time")->
                           GetText());
  double relaxation_time = atof(sub_root->
                                FirstChildElement("relaxation_time")->
                                GetText());
  double current_time = atof(sub_root->
                         FirstChildElement("current_time")->
                         GetText());
  double end_time = atof(sub_root->
                         FirstChildElement("end_time")->
                         GetText());
  double delta_t = atof(sub_root->
                        FirstChildElement("delta_t")->
                        GetText());

  c_time = new Time(current_time, start_time,
                    relaxation_time, end_time, delta_t);
}

void Load_init_param::init_file_saving_parameters ()
{
  XMLElement* root = xml_data->FirstChildElement (INITIAL_PARAMS_NAME);
  XMLElement* sub_root = root->FirstChildElement (FILE_SAVE_PARAMS_NAME);
  XMLElement* dump_data_root = sub_root->FirstChildElement ("dump_data");


  char* path_res = (char*)sub_root->FirstChildElement("path_to_result")->GetText();
  char* path_dump = (char*)sub_root->FirstChildElement("path_to_save_state")->GetText();

  data_dump_interval = atoi(sub_root
                            ->FirstChildElement("data_dump_interval")
                            ->GetText());
  system_state_dump_interval = atoi(sub_root
                                    ->FirstChildElement("system_state_dump_interval")
                                    ->GetText());
  frames_per_file = atoi(sub_root
                         ->FirstChildElement("frames_per_file")
                         ->GetText());

  // dump data configuration
  is_dump_e1 = to_bool(dump_data_root->FirstChildElement("e1")->GetText());
  is_dump_e2 = to_bool(dump_data_root->FirstChildElement("e2")->GetText());
  is_dump_e3 = to_bool(dump_data_root->FirstChildElement("e3")->GetText());
  is_dump_h1 = to_bool(dump_data_root->FirstChildElement("h1")->GetText());
  is_dump_h2 = to_bool(dump_data_root->FirstChildElement("h2")->GetText());
  is_dump_h3 = to_bool(dump_data_root->FirstChildElement("h3")->GetText());
  is_dump_rho_beam = to_bool(dump_data_root->FirstChildElement("rho_beam")->GetText());

  c_io_class = new input_output_class (path_res, path_dump);
}

bool Load_init_param::save_system_state(double timestamp)
{
  cout << "Saving system state at " << timestamp << endl;
  c_io_class->out_field_dump((char*)"e1",efield->e1,c_geom->n_grid_1-1,c_geom->n_grid_2-1);
  c_io_class->out_field_dump((char*)"e2",efield->e2,c_geom->n_grid_1-1,c_geom->n_grid_2-1);
  c_io_class->out_field_dump((char*)"e3",efield->e3,c_geom->n_grid_1-1,c_geom->n_grid_2-1);
  c_io_class->out_field_dump((char*)"h1",hfield->h1,c_geom->n_grid_1-1,c_geom->n_grid_2-1);
  c_io_class->out_field_dump((char*)"h2",hfield->h2,c_geom->n_grid_1-1,c_geom->n_grid_2-1);
  c_io_class->out_field_dump((char*)"h3",hfield->h3,c_geom->n_grid_1-1,c_geom->n_grid_2-1);
  for(unsigned int i=0;i<p_list->part_list.size();i++)
  {
    c_io_class->out_coord_dump(p_list->part_list[i]->name,p_list->part_list[i]->x1, p_list->part_list[i]->x3, p_list->part_list[i]->number);
    c_io_class->out_velocity_dump(p_list->part_list[i]->name,p_list->part_list[i]->v1, p_list->part_list[i]->v2,p_list->part_list[i]->v3, p_list->part_list[i]->number);

  }
  return(true);
}

void Load_init_param::run(void)
{
  cout << "\nLaunch Simulation\n\n";

  this->c_time->current_time = 0.0 ;
  p_list->charge_weighting(c_rho_new);

  // Seems: https://en.wikipedia.org/wiki/Dirichlet_distribution
  Poisson_dirichlet dirih(c_geom);
  dirih.poisson_solve(efield, c_rho_new);

  //variable for out_class function
  p_list->create_coord_arrays();
  int step_number = 0;
  time_t t1 = time(0);

  
  while (c_time->current_time < c_time->end_time)
  {
    c_bunch->bunch_inject(c_time);

    // 1. Calculate H field
    hfield->calc_field(efield, c_time);

    // 2. Calculate v
    c_current->reset_j();
    c_rho_old->reset_rho();
    c_rho_beam->reset_rho();
    p_list->step_v(efield, hfield, c_time);

    // 3. Calculate x, calculate J
    p_list->copy_coords();
    p_list->charge_weighting(c_rho_old);  //continuity equation
    p_list->half_step_coord(c_time);
    p_list->azimuthal_j_weighting(c_time, c_current);
    p_list->half_step_coord(c_time);
    p_list->j_weighting(c_time,c_current);

    // 4. Calculate E
    // maxwell_rad.probe_mode_exitation(&geom1,&current1, 1,7e8, time1.current_time);
    efield->calc_field(hfield,c_time, c_current, c_pml);

    //continuity equation
    c_rho_new->reset_rho();

    p_list->charge_weighting(c_rho_new);  //continuity equation
    //bool res = continuity_equation(c_time, c_geom, c_current, c_rho_old, c_rho_new);

    // print header on every 20 logging steps
    if  ((((int)(c_time->current_time/c_time->delta_t))%(data_dump_interval*20)==0))
    {
      cout << endl
           << left << setw(8) << "Step"
           << left << setw(13) << "Saved Frame"
           << left << setw(18) << "Model Time (sec)"
           << left << setw(32) << "Approx. Step Execution Time (sec)"
           << endl;
    }

    if  ((((int)(c_time->current_time/c_time->delta_t))%data_dump_interval==0))
    {
      cout << left << setw(8) << step_number * data_dump_interval
           << left << setw(13) << step_number
           << left << setw(18) << c_time->current_time
           << left << setw(32) << (double)(time(0) - t1) / data_dump_interval
           << endl;

      c_bunch->charge_weighting(c_rho_beam);
      c_rho_old->reset_rho();
      p_list[0].charge_weighting(c_rho_old);

      if (is_dump_e1) c_io_class->out_data("e1",efield->e1,step_number,frames_per_file,c_geom->n_grid_1-1,c_geom->n_grid_2-1);
      if (is_dump_e2) c_io_class->out_data("e2",efield->e2,step_number,frames_per_file,c_geom->n_grid_1-1,c_geom->n_grid_2-1);
      if (is_dump_e3) c_io_class->out_data("e3",efield->e3,step_number,frames_per_file,c_geom->n_grid_1-1,c_geom->n_grid_2-1);

      if (is_dump_h1) c_io_class->out_data("h1",hfield->h1,step_number,frames_per_file,c_geom->n_grid_1-1,c_geom->n_grid_2-1);
      if (is_dump_h2) c_io_class->out_data("h2",hfield->h2,step_number,frames_per_file,c_geom->n_grid_1-1,c_geom->n_grid_2-1);
      if (is_dump_h3) c_io_class->out_data("h3",hfield->h3,step_number,frames_per_file,c_geom->n_grid_1-1,c_geom->n_grid_2-1);

      if (is_dump_rho_beam) c_io_class->out_data("rho_beam", c_rho_beam->get_ro(),step_number,frames_per_file,c_geom->n_grid_1-1,c_geom->n_grid_2-1);

      step_number += 1;
      t1 = time(0);
      if  ((((int)(c_time->current_time/c_time->delta_t))%system_state_dump_interval==0)&&(step_number!=1))
        this->save_system_state(c_time->current_time);
    }
    c_time->current_time = c_time->current_time + c_time->delta_t;
    //if (!res)
    //  cout<<"Error:"<<c_time->current_time<<"! ";
    cout << "\nSimulation Completed\n\n";    
  }
}
