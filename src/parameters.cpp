#include "parameters.h"

// using namespace std;
// using namespace tinyxml2;
// using namespace math::fourier;

#define INITIAL_PARAMS_NAME "initial_parameters"
#define PML_PARAMS_NAME "pml"
#define GEOMETRY_PARAMS_NAME "geometry"
#define TIME_PARAMS_NAME "time"
#define PARTICLES_PARAMS_NAME "particles"
#define BEAM_PARAMS_NAME "particle_beam"
#define BOUNDARY_MAXWELL_PARAMS_NAME "boundary_maxwell_conditions"
#define DATA_DUMP_PARAMS_NAME "file_save_parameters"

// Parameters::Parameters(void) // : Parameters("parameters.xml")
// {
// }

Parameters::Parameters(const char *xml_file_name)
{

  //! read given xml file and fill xml_data class variable
  xml_data = new XMLDocument(true, COLLAPSE_WHITESPACE);

  XMLError e_result = xml_data->LoadFile(xml_file_name);
  if (e_result != XML_SUCCESS)
  {
    cerr << "ERROR: Can not read configuration file ``" << xml_file_name << "``" << endl;
    exit (78);
  }

  XMLElement *root = xml_data->FirstChildElement(INITIAL_PARAMS_NAME);

  //! set debug option
  debug = lib::to_bool(root->FirstChildElement("debug")->GetText());

  //! initialize geometry, set PML
  init_geometry();
  init_pml();

  //! initialize time
  init_time ();

  //! initialize particle parameters
  init_particles();

  //! initialize particles beam parameters
  init_beam();

  init_data_dump_parameters();
}






  //! Constructor to initialize start parameters
  //! from parameters.xml and fill particle arrays.
  //! each step is separate method

  //! NOTE: all system init is too huge,
  //! so we use several "subconstructors"
  //! to initialize different parts of system

//   // define openmp-related options when openmp enabled
// #ifdef _OPENMP
//   omp_set_dynamic(1); // Explicitly disable dynamic teams
// #endif

//   //! Steps to initialize:

//   //! 1. read XML config file
//   cout << "Reading configuration file ``" << xml_file_name << "``" << endl;
//   read_xml(xml_file_name);

//   //! 2. load Geometry parameters
//   cout << "Initializing Geometry Parameters" << endl;
//   init_geometry();

//   //! 3. creating field objects
//   cout << "Initializing E/M Fields Data" << endl;
//   init_fields ();

//   //! 4. load time parameters
//   cout << "Initializing Time Data" << endl;
//   init_time ();

//   //! 5. load particle parameters
//   cout << "Initializing Particles Data" << endl;
//   init_particles();

//   //! 6. load bunch
//   cout << "Initializing Particles Bunch Data" << endl;
//   init_bunch();

//   //! 7. load boundary conditions
//   cout << "Initializing Boundary Conditions Data" << endl;
//   init_boundary();

//   //! 8. load File Path
//   cout << "Initializing File System Paths" << endl;
//   init_file_saving_parameters();

//   cout << "Initialization complete" << endl;
// }

// Parameters::~Parameters(void)
// {
// }

// void Parameters::read_xml(const char *xml_file_name)
// {
//   //! read given xml file and fill xml_data class variable
//   xml_data = new XMLDocument(true, COLLAPSE_WHITESPACE);

//   XMLError e_result = xml_data->LoadFile(xml_file_name);
//   if (e_result != XML_SUCCESS)
//   {
//     cerr << "ERROR: Can not read configuration file ``" << xml_file_name << "``" << endl;
//     exit (78);
//   }
//   //! set DEBUG flag
//   XMLElement *root = xml_data->FirstChildElement(INITIAL_PARAMS_NAME);
//   debug = lib::to_bool(root->FirstChildElement("debug")->GetText());
// }

void Parameters::init_particles()
{
  //! initialize particles charge, mass, number,
  //! density and temperature for all particle species
  const char *p_specie_section_name = "particle_specie";
  XMLElement *root = xml_data->FirstChildElement (INITIAL_PARAMS_NAME);
  XMLElement *specie = root
    ->FirstChildElement (PARTICLES_PARAMS_NAME)
    ->FirstChildElement (p_specie_section_name);

  while(specie)
  {
    particle_specie p_specie;

    // set specie name
    // strcpy (p_specie->name, specie->Attribute("name")); // ->GetText());
    // cout << "  Initializing " << p_name << " Data" << endl;

    p_specie.name = (char*)specie->Attribute("name");

    p_specie.charge = atoi(specie->FirstChildElement("charge")->GetText());
    p_specie.mass = atoi(specie->FirstChildElement("mass")->GetText());
    p_specie.macro_count = debug ?
      atof(specie->FirstChildElement("debug_macro_count")->GetText()) :
      atof(specie->FirstChildElement("macro_count")->GetText());
    p_specie.left_density = atof(specie->FirstChildElement("left_density")->GetText());
    p_specie.right_density = atof(specie->FirstChildElement("right_density")->GetText());
    p_specie.temperature = atof(specie->FirstChildElement("temperature")->GetText());

    particle_species.push_back(p_specie);

    specie = specie->NextSiblingElement(p_specie_section_name);
  }
}

void Parameters::init_beam()
{
  //! initialize particles bunch data
  //! particles bunch should be injected to plasma
  XMLElement *root = xml_data->FirstChildElement (INITIAL_PARAMS_NAME);
  XMLElement *sub_root =root->FirstChildElement (BEAM_PARAMS_NAME);

  // strcpy (beam_name, sub_root->FirstChildElement("name")->GetText());
  beam_name = (char*)sub_root->Attribute("name");

  beam_particle_charge = atoi(sub_root
                              ->FirstChildElement("charge")
                              ->GetText());
  beam_particle_mass = atoi(sub_root
                            ->FirstChildElement("mass")
                            ->GetText());
  beam_bunches_count = atoi(sub_root
                            ->FirstChildElement("bunches_count")
                            ->GetText());
  beam_bunches_distance = atof(sub_root
                               ->FirstChildElement("bunches_distance")
                               ->GetText());
  beam_initial_velocity = atof(sub_root
                               ->FirstChildElement("initial_velocity")
                               ->GetText());
  //
  bunch_macro_count = debug ?
    atof(sub_root
         ->FirstChildElement("debug_bunch_macro_count")
         ->GetText()) :
    atof(sub_root
         ->FirstChildElement("bunch_macro_count")
         ->GetText());

  bunch_lenght = atof(sub_root
                      ->FirstChildElement("bunch_length")
                      ->GetText());
  bunch_radius = atof(sub_root
                      ->FirstChildElement("bunch_radius")
                      ->GetText());
  bunch_density = atof(sub_root
                       ->FirstChildElement("bunch_density")
                       ->GetText());
}

void Parameters::init_boundary ()
{
  //! initialize boundaries and boundary conditions
  XMLElement *root = xml_data->FirstChildElement (INITIAL_PARAMS_NAME);
  XMLElement *sub_root =root->FirstChildElement (BOUNDARY_MAXWELL_PARAMS_NAME);
  //
  boundary_maxwell_e_phi_upper = atof(sub_root->
                                     FirstChildElement("e_fi_upper")->
                                     GetText());
  boundary_maxwell_e_phi_left = atof(sub_root->
                                    FirstChildElement("e_fi_left")->
                                    GetText());
  boundary_maxwell_e_phi_right = atof(sub_root->
                                     FirstChildElement("e_fi_right")->
                                     GetText());
  boundary_conditions = atoi(root->FirstChildElement("boundary_conditions")->GetText());
}

void Parameters::init_geometry ()
{
  //! initialize geometry
  XMLElement *root = xml_data->FirstChildElement (INITIAL_PARAMS_NAME);
  XMLElement *sub_root =root->FirstChildElement (GEOMETRY_PARAMS_NAME);

  //
  double r_size = atof(sub_root->
                       FirstChildElement("r_size")->
                       GetText());
  double z_size = atof(sub_root->
                       FirstChildElement("z_size")->
                       GetText());
  int n_grid_r = debug ?
    atoi(sub_root->
         FirstChildElement("debug_n_grid_r")->
         GetText()) :
    atoi(sub_root->
         FirstChildElement("n_grid_r")->
         GetText());
  int n_grid_z = debug ?
    atoi(sub_root->
         FirstChildElement("debug_n_grid_z")->
         GetText()) :
    atoi(sub_root->
         FirstChildElement("n_grid_z")->
         GetText());
  // PML

  geom = new Geometry(r_size, z_size, n_grid_r, n_grid_z);
  geom->set_epsilon();
}

void Parameters::init_pml()
{
  XMLElement *root = xml_data->FirstChildElement (INITIAL_PARAMS_NAME);
  XMLElement *sub_root =root->FirstChildElement (PML_PARAMS_NAME);

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

  //! add PML to geometry
  //! Perfectly Matched Layer description: https://en.wikipedia.org/wiki/Perfectly_matched_layer
  if ((comp_l_1 != 0) || (comp_l_2 != 0) || (comp_l_3 != 0))
    geom->set_pml(comp_l_1, comp_l_2, comp_l_3, sigma_1_t, sigma_2_t);

  pml_comparative_l_1 = comp_l_1;
  pml_comparative_l_2 = comp_l_2;
  pml_comparative_l_3 = comp_l_3;
  pml_sigma_1 = sigma_1_t;
  pml_sigma_2 = sigma_2_t;
}

// void Parameters::init_fields ()
// {
//   //! initialize electrical and magnetic fields
//   efield = new EField(c_geom);
//   hfield = new HField(c_geom);
//   efield->boundary_conditions();
//   efield->set_homogeneous_efield(0.0, 0.0, 0.0);
//   hfield->set_homogeneous_h(0.0, 0.0, 0.0);
//   efield->set_fi_on_z();
// }

void Parameters::init_time ()
{
  //! initialize time parameters: `start_time`, `relaxation_time` `end_time` and `delta_t`
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
  double step_interval = atof(sub_root->
                        FirstChildElement("step_interval")->
                        GetText());

  time = new Time(current_time, start_time,
                  relaxation_time, end_time, step_interval);
}

void Parameters::init_data_dump_parameters ()
{
  //! initialize file saving parameters, like path to computed data files,
  //! path to system state data files max frames number, placed to one file etc.
  XMLElement *root = xml_data->FirstChildElement (INITIAL_PARAMS_NAME);
  XMLElement *sub_root = root->FirstChildElement (DATA_DUMP_PARAMS_NAME);
  XMLElement *dump_data_root = sub_root->FirstChildElement ("dump_data");

  dump_result_path = (char*)sub_root->FirstChildElement("result_path")->GetText();
  dump_save_state_path = (char*)sub_root->FirstChildElement("save_state_path")->GetText();

  dump_data_interval = debug ?
    atoi(sub_root
         ->FirstChildElement("debug_data_dump_interval")
         ->GetText()) :
    atoi(sub_root
         ->FirstChildElement("data_dump_interval")
         ->GetText());

  dump_system_state_interval = atoi(sub_root
                                    ->FirstChildElement("system_state_dump_interval")
                                    ->GetText());
  dump_frames_per_file = atoi(sub_root
                              ->FirstChildElement("frames_per_file")
                              ->GetText());

  //! choose, which parameters should be dumped to data files
  dump_e_r = lib::to_bool(dump_data_root->FirstChildElement("E_r")->GetText());
  dump_e_phi = lib::to_bool(dump_data_root->FirstChildElement("E_phi")->GetText());
  dump_e_z = lib::to_bool(dump_data_root->FirstChildElement("E_z")->GetText());

  dump_h_r = lib::to_bool(dump_data_root->FirstChildElement("H_r")->GetText());
  dump_h_phi = lib::to_bool(dump_data_root->FirstChildElement("H_phi")->GetText());
  dump_h_z = lib::to_bool(dump_data_root->FirstChildElement("H_z")->GetText());

  dump_current_r = lib::to_bool(dump_data_root->FirstChildElement("J_r")->GetText());
  dump_current_phi = lib::to_bool(dump_data_root->FirstChildElement("J_phi")->GetText());
  dump_current_z = lib::to_bool(dump_data_root->FirstChildElement("J_z")->GetText());

  dump_position = lib::to_bool(dump_data_root->FirstChildElement("position")->GetText());
  dump_velocity = lib::to_bool(dump_data_root->FirstChildElement("velocity")->GetText());

  dump_rho_beam = lib::to_bool(dump_data_root->FirstChildElement("rho_bunch")->GetText());



  // c_io_class = new InputOutputClass (path_res, path_dump);
}

// void Parameters::dump_system_state()
// {
//   //! dump system state to file set
//   c_io_class->out_field_dump("E_r", efield->field_r, c_geom->n_grid_1 - 1, c_geom->n_grid_2 - 1);
//   c_io_class->out_field_dump("E_phi", efield->field_phi, c_geom->n_grid_1 - 1, c_geom->n_grid_2 - 1);
//   c_io_class->out_field_dump("E_z", efield->field_z, c_geom->n_grid_1 - 1, c_geom->n_grid_2 - 1);
//   c_io_class->out_field_dump("H_r", hfield->field_r, c_geom->n_grid_1 - 1, c_geom->n_grid_2 - 1);
//   c_io_class->out_field_dump("H_phi", hfield->field_phi, c_geom->n_grid_1 - 1, c_geom->n_grid_2 - 1);
//   c_io_class->out_field_dump("H_z", hfield->field_z, c_geom->n_grid_1 - 1, c_geom->n_grid_2 - 1);
//   for(unsigned int i=0; i<p_list->part_list.size(); i++)
//   {
//     c_io_class->out_pos_dump(p_list->part_list[i]->name,
//                                p_list->part_list[i]->pos,
//                                p_list->part_list[i]->number);

//     c_io_class->out_velocity_dump(p_list->part_list[i]->name,
//                                   p_list->part_list[i]->vel,
//                                   p_list->part_list[i]->number);
//     c_io_class->out_current_time_dump(c_time->current_time);
//   }
// }

// void Parameters::dump_data(int step_number)
// {
//   //! dump claculated data frames (for fields etc.) to files set
//   if (is_dump_e_r)
//     c_io_class->out_data("E_r", efield->field_r, step_number,
//                          frames_per_file, c_geom->n_grid_1 - 1, c_geom->n_grid_2 - 1);

//   if (is_dump_e_phi)
//     c_io_class->out_data("E_phi", efield->field_phi, step_number,
//                          frames_per_file, c_geom->n_grid_1 - 1, c_geom->n_grid_2 - 1);

//   if (is_dump_e_z)
//     c_io_class->out_data("E_z", efield->field_z, step_number,
//                          frames_per_file, c_geom->n_grid_1 - 1, c_geom->n_grid_2 - 1);

//   if (is_dump_h_r)
//     c_io_class->out_data("H_r", hfield->field_r, step_number,
//                          frames_per_file, c_geom->n_grid_1 - 1, c_geom->n_grid_2 - 1);

//   if (is_dump_h_phi)
//     c_io_class->out_data("H_phi", hfield->field_phi, step_number,
//                          frames_per_file, c_geom->n_grid_1 - 1, c_geom->n_grid_2 - 1);

//   if (is_dump_h_z)
//     c_io_class->out_data("H_z", hfield->field_z, step_number,
//                          frames_per_file, c_geom->n_grid_1 - 1, c_geom->n_grid_2 - 1);

//   if (is_dump_rho_bunch)
//     c_io_class->out_data("rho_bunch",  c_rho_bunch->get_rho(), step_number,
//                          frames_per_file, c_geom->n_grid_1 - 1, c_geom->n_grid_2 - 1);
// }

// void Parameters::run(void)
// {
//   //! Launch calculation.
//   //! THIS IS ENTRY POINT TO MAIN PDP3 CALCULATION CYCLE
//   cout << endl << "Launch Simulation" << endl << endl;

//   clock_t the_time = 0;
//   int step_number = 0;
//   time_t t1 = time(0);
//   char avg_step_exec_time[24]; // rounded and formatted average step execution time

//   //! Main calculation cycle
//   while (c_time->current_time <= c_time->end_time)
//   {
//     //! Steps:
// #ifdef PERF_DEBUG
//     the_time = 0;
//     cerr << the_time << " loop iteration begin" << endl;
//     the_time = clock();
// #endif

//     //! 1. inject bunch
//     c_bunch->bunch_inject(c_time);
// #ifdef PERF_DEBUG
//     cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: c_bunch->bunch_inject" << endl;
//     the_time = clock();
// #endif

//     //! 2. Calculate H field
//     hfield->calc_field(efield, c_time);
// #ifdef PERF_DEBUG
//     cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: hfield->calc_field" << endl;
//     the_time = clock();
// #endif

//     //! 3. Calculate v
//     c_current->reset_j();
// #ifdef PERF_DEBUG
//     cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: c_current->reset_j" << endl;
//     the_time = clock();
// #endif
//     c_rho_old->reset_rho();
// #ifdef PERF_DEBUG
//     cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: c_rho_old->reset_rho" << endl;
//     the_time = clock();
// #endif
//     c_rho_bunch->reset_rho();
// #ifdef PERF_DEBUG
//     cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: c_rho_bunch->reset_rho" << endl;
//     the_time = clock();
// #endif
//     p_list->step_v(efield, hfield, c_time);
// #ifdef PERF_DEBUG
//     cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: p_list->step_v" << endl;
//     the_time = clock();
// #endif

//     //! 4. Calculate x, calculate J
//     p_list->dump_position_to_old();
//     // cout << 4 << endl;
//     //! FIXME: for some reason charge_weighting has no effect on result
//     // p_list->charge_weighting(c_rho_old); //continuity equation
// #ifdef PERF_DEBUG
//     cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: p_list->dump_position_to_old" << endl;
//     the_time = clock();
// #endif

//     p_list->half_step_pos(c_time);
// #ifdef PERF_DEBUG
//     cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: p_list->half_step_pos" << endl;
//     the_time = clock();
// #endif
//     p_list->back_position_to_rz();
// #ifdef PERF_DEBUG
//     cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: p_list->back_position_to_rz" << endl;
//     the_time = clock();
// #endif

//     p_list->azimuthal_j_weighting(c_time, c_current);
// #ifdef PERF_DEBUG
//     cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: p_list->azimuthal_j_weighting" << endl;
//     the_time = clock();
// #endif

//     p_list->half_step_pos(c_time);
// #ifdef PERF_DEBUG
//     cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: p_list->half_step_pos" << endl;
//     the_time = clock();
// #endif
//     p_list->back_position_to_rz();
// #ifdef PERF_DEBUG
//     cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: p_list->back_position_to_rz" << endl;
//     the_time = clock();
// #endif

//     p_list->j_weighting(c_time, c_current);
// #ifdef PERF_DEBUG
//     cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: j_weighting" << endl;
//     the_time = clock();
// #endif

//     p_list->back_velocity_to_rz();
// #ifdef PERF_DEBUG
//     cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: back_velocity_to_rz" << endl;
//     the_time = clock();
// #endif

//     //! 5. Calculate E
//     // maxwell_rad.probe_mode_exitation(&geom1,&current1, 1,7e8, time1.current_time);
//     efield->calc_field(hfield, c_time, c_current);
// #ifdef PERF_DEBUG
//     cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: efield->calc_field" << endl;
//     the_time = clock();
// #endif

//     //! 6. Continuity equation
//     c_rho_new->reset_rho(); // TODO: is it 4?
// #ifdef PERF_DEBUG
//     cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: c_rho_new->reset_rho" << endl;
//     the_time = clock();
// #endif

//     //! FIXME: for some reason charge_weighting has no effect on result
//     // p_list->charge_weighting(c_rho_new); // continuity equation

//     //! print header on every 20 logging steps
//     if  ((int)(c_time->current_time / c_time->delta_t) % (data_dump_interval*20) == 0)
//     {
//       cout << endl
//            << left << setw(8) << "Step"
//            << left << setw(13) << "Saved Frame"
//            << left << setw(18) << "Model Time (sec)"
//            << left << setw(32) << "Avg. Step Calculation Time (sec)"
//            << endl;
//     }

//     //! dump data to corresponding files every `parameters.xml->file_save_parameters->data_dump_interval` steps
//     if  ((int)(c_time->current_time / c_time->delta_t) % data_dump_interval == 0)
//     {
//       sprintf(avg_step_exec_time, "%.2f", (double)(time(0) - t1) / data_dump_interval);

//       cout << left << setw(8) << step_number * data_dump_interval
//            << left << setw(13) << step_number
//            << left << setw(18) << c_time->current_time
//            << left << setw(32) << avg_step_exec_time
//            << endl;

//       c_bunch->charge_weighting(c_rho_bunch);
// #ifdef PERF_DEBUG
//       cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: c_bunch->charge_weighting" << endl;
//       the_time = clock();
// #endif
//       // c_rho_old->reset_rho();
//       // p_list[0].charge_weighting(c_rho_old);

//       dump_data(step_number);
// #ifdef PERF_DEBUG
//       cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: dump_data" << endl;
//       the_time = clock();
// #endif

//       step_number += 1;
//       t1 = time(0);
//       if  ((((int)(c_time->current_time/c_time->delta_t))%system_state_dump_interval==0)&&(step_number!=1))
//       {
//         cout << "Saving system state at " << c_time->current_time << endl;
//         dump_system_state();
// #ifdef PERF_DEBUG
// 	cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: dump_system_state" << endl;
// 	the_time = clock();
// #endif
//       }
//     }
//     c_time->current_time = c_time->current_time + c_time->delta_t;
//   }
//   cout << endl << "Simulation Completed" << endl << endl;
// }
  //
