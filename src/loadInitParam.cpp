#include <algorithm>
#include <iomanip>
#include <ctime>
#include <math.h>

// enable openmp optional
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

#include "loadInitParam.h"
#include "time.h"
#include "tinyxml2.h"

using namespace std;
using namespace tinyxml2;

#define INITIAL_PARAMS_NAME "initial_parameters"
#define PML_PARAMS_NAME "pml"
#define GEOMETRY_PARAMS_NAME "geometry"
#define TIME_PARAMS_NAME "time"
#define PARTICLES_PARAMS_NAME "particles"
#define BUNCH_PARAMS_NAME "particles_bunch"
#define BOUNDARY_MAXWELL_PARAMS_NAME "boundary_maxwell_conditions"
#define FILE_SAVE_PARAMS_NAME "file_save_parameters"

LoadInitParam::LoadInitParam(void)
{
}

LoadInitParam::LoadInitParam(char *xml_file_name)
{
  // NOTE: all system init is too huge,
  // so we use several "subconstructors"
  // to initialize different parts of system

  // define openmp-related options when openmp enabled
#ifdef _OPENMP
  omp_set_dynamic(0); // Explicitly disable dynamic teams
#endif

  // read XML config file
  cout << "Reading configuration file ``" << xml_file_name << "``" << endl;
  read_xml(xml_file_name);

  // load Geometry parameters
  cout << "Initializing Geometry Parameters" << endl;
  init_geometry();

  // creating field objects
  cout << "Initializing E/M Fields Data" << endl;
  init_fields ();

  // load time parameters
  cout << "Initializing Time Data" << endl;
  init_time ();

  // load particle parameters
  cout << "Initializing Particles Data" << endl;
  init_particles();

  // load bunch
  cout << "Initializing Particles Bunch Data" << endl;
  init_bunch();

  cout << "Initializing Bounrary Conditions Data" << endl;
  init_boundary();

  // load File Path
  cout << "Initializing File System Paths" << endl;
  init_file_saving_parameters();

  cout << "Initialization complete" << endl;
}

LoadInitParam::~LoadInitParam(void)
{
}

void LoadInitParam::read_xml(const char *xml_file_name)
{
  xml_data = new XMLDocument(true, COLLAPSE_WHITESPACE);

  XMLError e_result = xml_data->LoadFile(xml_file_name);
  if (e_result != XML_SUCCESS)
  {
    cerr << "ERROR: Can not read configuration file ``" << xml_file_name << "``" << endl;
    exit (78);
  }
}

bool LoadInitParam::to_bool(string str) {
  std::transform(str.begin(), str.end(), str.begin(), ::tolower);
  std::istringstream is(str);
  bool b;
  is >> std::boolalpha >> b;
  return b;
}

void LoadInitParam::init_particles()
{
  const char *p_king_section_name = "particle_kind";
  XMLElement *root = xml_data->FirstChildElement (INITIAL_PARAMS_NAME);
  XMLElement *particle_kind = root
    ->FirstChildElement (PARTICLES_PARAMS_NAME)
    ->FirstChildElement (p_king_section_name);

  char *p_name = new char [50];
  double charge;
  double mass; // TODO: should it be double?
  double number;
  double left_density;
  double right_density;
  double temperature;

  Particles *prtls;

  // initialize particles list
  p_list = new ParticlesList();

  // creating rho and current arrays
  // WARNING! should be called after geometry initialized
  c_rho_new = new ChargeDensity(c_geom);
  c_rho_old = new ChargeDensity(c_geom);
  c_rho_beam = new ChargeDensity(c_geom);
  c_current = new Current(c_geom);

  while(particle_kind)
  {
    strcpy (p_name, particle_kind->FirstChildElement("name")->GetText());
    cout << "  Initializing " << p_name << " Data" << endl;

    charge = atof(particle_kind->FirstChildElement("charge")->GetText());
    mass = atof(particle_kind->FirstChildElement("mass")->GetText());
    number = atof(particle_kind->FirstChildElement("number")->GetText());
    left_density = atof(particle_kind->FirstChildElement("left_density")->GetText());
    right_density = atof(particle_kind->FirstChildElement("right_density")->GetText());
    temperature = atof(particle_kind->FirstChildElement("temperature")->GetText());

    // init and setup particles properties
    prtls = new Particles(strcpy(new char [50], p_name), charge, mass, number, c_geom);
		p_list->part_list.push_back(prtls); // push particles to particles list vector

    prtls->load_spatial_distribution_with_variable_mass(left_density,right_density,0,0);
    prtls->velocity_distribution(temperature);
    particle_kind = particle_kind->NextSiblingElement(p_king_section_name);
  }
}

void LoadInitParam::init_bunch()
{
  // initialize particles bunch data
  XMLElement *root = xml_data->FirstChildElement (INITIAL_PARAMS_NAME);
  XMLElement *sub_root =root->FirstChildElement (BUNCH_PARAMS_NAME);

  char *p_name = (char*)sub_root->FirstChildElement("name")->GetText();
	Bunch *prtls;
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

  prtls = new Bunch(p_name, charge, mass,number, c_geom, // p_list,
										duration, radius, density, initial_velocity);
	p_list->part_list.push_back(prtls); // push bunch to particles list vector

  c_bunch = prtls;
}

void LoadInitParam::init_boundary ()
{
  XMLElement *root = xml_data->FirstChildElement (INITIAL_PARAMS_NAME);
  XMLElement *sub_root =root->FirstChildElement (BOUNDARY_MAXWELL_PARAMS_NAME);
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

  int boundary_conditions = atoi(root->FirstChildElement("boundary_conditions")->GetText());

  // Maxwell initial conditions
  BoundaryMaxwellConditions maxwell_rad(efield); // TODO: WTF?
  maxwell_rad.specify_initial_field(c_geom, e_fi_upper, e_fi_left, e_fi_right);

  if (boundary_conditions == 0)
  {
    p_list->charge_weighting(c_rho_new);
    Fourier four1();
    PoissonDirichlet dirih(c_geom);
    dirih.poisson_solve(efield, c_rho_new);
  }
}

void LoadInitParam::init_geometry ()
{
  XMLElement *root = xml_data->FirstChildElement (INITIAL_PARAMS_NAME);
  XMLElement *sub_root =root->FirstChildElement (GEOMETRY_PARAMS_NAME);
  XMLElement *pml_sub_root =root->FirstChildElement (PML_PARAMS_NAME);

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
  // PML
  double comp_l_1 = atof(pml_sub_root->
                         FirstChildElement("comparative_l_1")->
                         GetText());
  double comp_l_2 = atof(pml_sub_root->
                         FirstChildElement("comparative_l_2")->
                         GetText());
  double comp_l_3 = atof(pml_sub_root->
                         FirstChildElement("comparative_l_3")->
                         GetText());
  double sigma_1_t = atof(pml_sub_root->
                          FirstChildElement("sigma_1")->
                          GetText());
  double sigma_2_t = atof(pml_sub_root->
                          FirstChildElement("sigma_2")->
                          GetText());

  c_geom = new Geometry(r_size, z_size, n_grid_r, n_grid_z);

  // add PML to geometry
  if ((comp_l_1 != 0) || (comp_l_2 != 0) || (comp_l_3 != 0))
    c_geom->set_pml(comp_l_1, comp_l_2, comp_l_3, sigma_1_t, sigma_2_t);

  c_geom->set_epsilon();
}

void LoadInitParam::init_fields ()
{
  efield = new EField(c_geom);
  hfield = new HField(c_geom);
  efield->boundary_conditions();
  efield->set_homogeneous_efield(0.0, 0.0, 0.0);
  hfield->set_homogeneous_h(0.0, 0.0, 0.0);
  efield->set_fi_on_z();
}

void LoadInitParam::init_time ()
{
  XMLElement *root = xml_data->FirstChildElement (INITIAL_PARAMS_NAME);
  XMLElement *sub_root =root->FirstChildElement (TIME_PARAMS_NAME);
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

void LoadInitParam::init_file_saving_parameters ()
{
  XMLElement *root = xml_data->FirstChildElement (INITIAL_PARAMS_NAME);
  XMLElement *sub_root = root->FirstChildElement (FILE_SAVE_PARAMS_NAME);
  XMLElement *dump_data_root = sub_root->FirstChildElement ("dump_data");


  char *path_res = (char*)sub_root->FirstChildElement("path_to_result")->GetText();
  char *path_dump = (char*)sub_root->FirstChildElement("path_to_save_state")->GetText();

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
  is_dump_e_r = to_bool(dump_data_root->FirstChildElement("E_r")->GetText());
  is_dump_e_phi = to_bool(dump_data_root->FirstChildElement("E_phi")->GetText());
  is_dump_e_z = to_bool(dump_data_root->FirstChildElement("E_z")->GetText());
  is_dump_h_r = to_bool(dump_data_root->FirstChildElement("H_r")->GetText());
  is_dump_h_phi = to_bool(dump_data_root->FirstChildElement("H_phi")->GetText());
  is_dump_h_z = to_bool(dump_data_root->FirstChildElement("H_z")->GetText());
  is_dump_rho_beam = to_bool(dump_data_root->FirstChildElement("rho_beam")->GetText());

  c_io_class = new InputOutputClass (path_res, path_dump);
}

bool LoadInitParam::save_system_state(double timestamp)
{
  cout << "Saving system state at " << timestamp << endl;
  c_io_class->out_field_dump((char*)"E_r",efield->field_r,c_geom->n_grid_1-1,c_geom->n_grid_2-1);
  c_io_class->out_field_dump((char*)"E_phi",efield->field_phi,c_geom->n_grid_1-1,c_geom->n_grid_2-1);
  c_io_class->out_field_dump((char*)"E_z",efield->field_z,c_geom->n_grid_1-1,c_geom->n_grid_2-1);
  c_io_class->out_field_dump((char*)"H_r",hfield->field_r,c_geom->n_grid_1-1,c_geom->n_grid_2-1);
  c_io_class->out_field_dump((char*)"H_phi",hfield->field_phi,c_geom->n_grid_1-1,c_geom->n_grid_2-1);
  c_io_class->out_field_dump((char*)"H_z",hfield->field_z,c_geom->n_grid_1-1,c_geom->n_grid_2-1);
  for(unsigned int i=0;i<p_list->part_list.size();i++)
  {
    c_io_class->out_coord_dump(p_list->part_list[i]->name,p_list->part_list[i]->x1, p_list->part_list[i]->x3, p_list->part_list[i]->number);
    c_io_class->out_velocity_dump(p_list->part_list[i]->name,p_list->part_list[i]->v1, p_list->part_list[i]->v2,p_list->part_list[i]->v3, p_list->part_list[i]->number);

  }
  return(true);
}

void LoadInitParam::run(void)
{
  cout << endl << "Launch Simulation" << endl << endl;

  this->c_time->current_time = 0.0 ;
  p_list->charge_weighting(c_rho_new);

  // Seems: https://en.wikipedia.org/wiki/Dirichlet_distribution
  PoissonDirichlet dirih(c_geom);
  dirih.poisson_solve(efield, c_rho_new);

  //variable for out_class function
  p_list->create_coord_arrays();
  int step_number = 0;
  time_t t1 = time(0);
  char avg_step_exec_time[24]; // rounded and formatted average step execution time

  while (c_time->current_time <= c_time->end_time)
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
    efield->calc_field(hfield,c_time, c_current);

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
           << left << setw(32) << "Avg. Step Calculation Time (sec)"
           << endl;
    }

    if  ((((int)(c_time->current_time/c_time->delta_t))%data_dump_interval==0))
    {
      sprintf(avg_step_exec_time, "%.2f", (double)(time(0) - t1) / data_dump_interval);

      cout << left << setw(8) << step_number  *data_dump_interval
           << left << setw(13) << step_number
           << left << setw(18) << c_time->current_time
           << left << setw(32) << avg_step_exec_time
           << endl;

      c_bunch->charge_weighting(c_rho_beam);
      // c_rho_old->reset_rho();
      // p_list[0].charge_weighting(c_rho_old);

      if (is_dump_e_r) c_io_class->out_data("E_r",efield->field_r,step_number,frames_per_file,c_geom->n_grid_1-1,c_geom->n_grid_2-1);
      if (is_dump_e_phi) c_io_class->out_data("E_phi",efield->field_phi,step_number,frames_per_file,c_geom->n_grid_1-1,c_geom->n_grid_2-1);
      if (is_dump_e_z) c_io_class->out_data("E_z",efield->field_z,step_number,frames_per_file,c_geom->n_grid_1-1,c_geom->n_grid_2-1);

      if (is_dump_h_r) c_io_class->out_data("H_e",hfield->field_r,step_number,frames_per_file,c_geom->n_grid_1-1,c_geom->n_grid_2-1);
      if (is_dump_h_phi) c_io_class->out_data("H_phi",hfield->field_phi,step_number,frames_per_file,c_geom->n_grid_1-1,c_geom->n_grid_2-1);
      if (is_dump_h_z) c_io_class->out_data("H_z",hfield->field_z,step_number,frames_per_file,c_geom->n_grid_1-1,c_geom->n_grid_2-1);

      if (is_dump_rho_beam) c_io_class->out_data("rho_beam", c_rho_beam->get_rho(),step_number,frames_per_file,c_geom->n_grid_1-1,c_geom->n_grid_2-1);

      step_number += 1;
      t1 = time(0);
      if  ((((int)(c_time->current_time/c_time->delta_t))%system_state_dump_interval==0)&&(step_number!=1))
        this->save_system_state(c_time->current_time);
    }
    c_time->current_time = c_time->current_time + c_time->delta_t;
    //if (!res)
    //  cout<<"Error:"<<c_time->current_time<<"! ";
  }
  cout << endl << "Simulation Completed" << endl << endl;
}
