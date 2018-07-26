#include "loadInitParam.h"

using namespace std;
using namespace math::fourier;

LoadInitParam::LoadInitParam(void)
{
}

LoadInitParam::LoadInitParam(char *xml_file_name)
{
  //! Constructor to initialize start parameters
  //! from parameters.xml and fill particle arrays.
  //! each step is separate method

  //! NOTE: all system init is too huge,
  //! so we use several "subconstructors"
  //! to initialize different parts of system

// define openmp-related options when openmp enabled
#ifdef _OPENMP
#ifdef OPENMP_DYNAMIC_THREADS
  omp_set_dynamic(1); // Explicitly disable dynamic teams
#else
  int cores = omp_get_num_procs();
  cout << "CORES: " << cores << endl;
  omp_set_dynamic(0);         // Explicitly disable dynamic teams
  omp_set_num_threads(cores); // Use 4 threads for all consecutive parallel regions
#endif
#endif

  //! Steps to initialize:

  //! 1. read XML config file
  cout << "Reading configuration file ``" << xml_file_name << "``" << endl;
  params = new Parameters(xml_file_name);
  // read_xml(xml_file_name);

  //! 2. load File Saving Paths
  cout << "Initializing File Saving Paths" << endl;
  //! initialize file saving parameters, like path to computed data files,
  //! path to system state data files max frames number, placed to one file etc.
  c_io_class = new InputOutputClass (params->dump_result_path, params->dump_save_state_path);

  //! 2. creating field objects
  cout << "Initializing E/M Fields Data" << endl;
  init_fields ();

  //! 3. load particle parameters
  cout << "Initializing Particles Data" << endl;
  init_particles();

  //! 4. load particles beam
  cout << "Initializing Particles Beam Data" << endl;
  init_beam();

  //! 5. load boundary conditions
  cout << "Initializing Boundary Conditions Data" << endl;
  init_boundary();

  cout << "Initialization complete" << endl;
}

LoadInitParam::~LoadInitParam(void)
{
}

void LoadInitParam::init_particles()
{
  //! initialize particles charge, mass, number,
  //! density and temperature for all particle kinds

  Particles *prtls;

  //! initialize particles list (array)
  p_list = new ParticlesList();

  //! creating rho and current arrays
  // WARNING! should be called after geometry initialized
  c_rho_new = new ChargeDensity(params->geom);
  c_rho_old = new ChargeDensity(params->geom);
  c_rho_bunch = new ChargeDensity(params->geom);
  c_current = new Current(params->geom);

  for (unsigned int i=0; i < params->particle_species.size(); i++)
  {
    particle_specie p_p = params->particle_species[i];

    cout << "  Initializing " << p_p.name << " Data" << endl;

    //! init and setup particles properties
    prtls = new Particles(p_p.name, p_p.charge, p_p.mass, p_p.number_macro, params->geom);
		p_list->part_list.push_back(prtls); // push particles to particles list vector

    // case 2 is for cylindrical distribution
    prtls->load_spatial_distribution(p_p.left_density, p_p.right_density, 0, 2);
    prtls->velocity_distribution(p_p.temperature);
  }

  p_list->charge_weighting(c_rho_new);
}

void LoadInitParam::init_beam()
{
  //! initialize particles bunch data
  //! particles bunch should be injected to plasma
  Bunch *prtls;

  // for (unsigned int i=0; i < params->beam_number_bunches)
  // {
  double duration = params->bunch_lenght / params->beam_initial_velocity;
  double hole_duration = params->beam_bunches_distance / params->beam_initial_velocity;

  for (unsigned int i=0; i < params->beam_number_bunches; i++)
    {
      prtls = new Bunch(params->beam_name,
                        params->beam_particle_charge,
                        params->beam_particle_mass,
                        params->bunch_number_macro,
                        params->geom, // p_list,
                        duration,
                        params->bunch_radius,
                        params->bunch_density,
                        params->beam_initial_velocity,
                        i,
                        hole_duration);
      p_list->part_list.push_back(prtls); // push bunch to particles list vector
      c_bunches.push_back(prtls); // push bunch to bunches list vector
      // }
    }
}

void LoadInitParam::init_boundary ()
{
  // //! initialize boundaries and boundary conditions

  // Maxwell initial conditions
  BoundaryMaxwellConditions maxwell_rad(efield); // TODO: WTF?
  maxwell_rad.specify_initial_field(params->geom,
                                    params->boundary_maxwell_e_phi_upper,
                                    params->boundary_maxwell_e_phi_left,
                                    params->boundary_maxwell_e_phi_right);

  if (params->boundary_conditions == 0)
  {
    p_list->charge_weighting(c_rho_new);
    Fourier four1();
    // Seems: https://en.wikipedia.org/wiki/Dirichlet_distribution
    PoissonDirichlet dirih(params->geom);
    dirih.poisson_solve(efield, c_rho_new);
  }
}

void LoadInitParam::init_fields ()
{
  //! initialize electrical and magnetic fields
  efield = new EField(params->geom);
  hfield = new HField(params->geom);
  efield->boundary_conditions();
  efield->set_homogeneous_efield(0.0, 0.0, 0.0);
  hfield->set_homogeneous_h(0.0, 0.0, 0.0);
  efield->set_fi_on_z();
}

void LoadInitParam::dump_system_state()
{
  //! dump system state to file set
  c_io_class->out_field_dump("E_r", efield->field_r, params->geom->n_grid_1 - 1, params->geom->n_grid_2 - 1);
  c_io_class->out_field_dump("E_phi", efield->field_phi, params->geom->n_grid_1 - 1, params->geom->n_grid_2 - 1);
  c_io_class->out_field_dump("E_z", efield->field_z, params->geom->n_grid_1 - 1, params->geom->n_grid_2 - 1);
  c_io_class->out_field_dump("H_r", hfield->field_r, params->geom->n_grid_1 - 1, params->geom->n_grid_2 - 1);
  c_io_class->out_field_dump("H_phi", hfield->field_phi, params->geom->n_grid_1 - 1, params->geom->n_grid_2 - 1);
  c_io_class->out_field_dump("H_z", hfield->field_z, params->geom->n_grid_1 - 1, params->geom->n_grid_2 - 1);
  for(unsigned int i=0; i<p_list->part_list.size(); i++)
  {
    c_io_class->out_pos_dump(p_list->part_list[i]->name,
                               p_list->part_list[i]->pos,
                               p_list->part_list[i]->number);

    c_io_class->out_velocity_dump(p_list->part_list[i]->name,
                                  p_list->part_list[i]->vel,
                                  p_list->part_list[i]->number);
    c_io_class->out_current_time_dump(params->time->current_time);
  }
}

void LoadInitParam::dump_data(int step_number)
{
  //! dump claculated data frames (for fields etc.) to files set

  // electrical fields
  if (params->dump_e_r)
    c_io_class->out_data("E_r", efield->field_r, step_number,
                         params->dump_frames_per_file, params->geom->n_grid_1 - 1, params->geom->n_grid_2 - 1);

  if (params->dump_e_phi)
    c_io_class->out_data("E_phi", efield->field_phi, step_number,
                         params->dump_frames_per_file, params->geom->n_grid_1 - 1, params->geom->n_grid_2 - 1);

  if (params->dump_e_z)
    c_io_class->out_data("E_z", efield->field_z, step_number,
                         params->dump_frames_per_file, params->geom->n_grid_1 - 1, params->geom->n_grid_2 - 1);

  // magnetic fields
  if (params->dump_h_r)
    c_io_class->out_data("H_r", hfield->field_r, step_number,
                         params->dump_frames_per_file, params->geom->n_grid_1 - 1, params->geom->n_grid_2 - 1);

  if (params->dump_h_phi)
    c_io_class->out_data("H_phi", hfield->field_phi, step_number,
                         params->dump_frames_per_file, params->geom->n_grid_1 - 1, params->geom->n_grid_2 - 1);

  if (params->dump_h_z)
    c_io_class->out_data("H_z", hfield->field_z, step_number,
                         params->dump_frames_per_file, params->geom->n_grid_1 - 1, params->geom->n_grid_2 - 1);

  // currents
  if (params->dump_current_r)
    c_io_class->out_data("J_r", c_current->get_j1(), step_number,
                         params->dump_frames_per_file, params->geom->n_grid_1 - 1, params->geom->n_grid_2 - 1);

  if (params->dump_current_phi)
    c_io_class->out_data("J_phi", c_current->get_j2(), step_number,
                         params->dump_frames_per_file, params->geom->n_grid_1 - 1, params->geom->n_grid_2 - 1);

  if (params->dump_current_z)
    c_io_class->out_data("J_z", c_current->get_j3(), step_number,
                         params->dump_frames_per_file, params->geom->n_grid_1 - 1, params->geom->n_grid_2 - 1);

  // beam density
  if (params->dump_rho_beam)
    c_io_class->out_data("rho_beam", c_rho_bunch->get_rho(), step_number,
                         params->dump_frames_per_file, params->geom->n_grid_1 - 1, params->geom->n_grid_2 - 1);

  // particle positions and velocities
  if (params->dump_position)
    for(unsigned int i=0; i<p_list->part_list.size(); i++)
    {
      char full_name[256];
      strcpy(full_name, p_list->part_list[i]->name);
      strcat(full_name, "_pos");

      c_io_class->out_triple(full_name,
                             p_list->part_list[i]->pos,
                             step_number,
                             params->dump_frames_per_file,
                             p_list->part_list[i]->number);
    }

  if (params->dump_position)
    for(unsigned int i=0; i<p_list->part_list.size(); i++)
    {
      char full_name[256];
      strcpy(full_name, p_list->part_list[i]->name);
      strcat(full_name, "_vel");

      c_io_class->out_triple(full_name,
                             p_list->part_list[i]->vel,
                             step_number,
                             params->dump_frames_per_file,
                             p_list->part_list[i]->number);
    }
}

void LoadInitParam::run(void)
{
  //! Launch calculation.
  //! THIS IS ENTRY POINT TO MAIN PDP3 CALCULATION LOOP
  cout << endl << "Launch Simulation" << endl << endl;

  clock_t the_time = 0;
  int step_number = 0;
  time_t t1 = time(0);
  char avg_step_exec_time[24]; // rounded and formatted average step execution time

  //! Main calculation loop
  while (params->time->current_time < params->time->end_time)
  {
    //! Steps:
#ifdef PERF_DEBUG
    the_time = 0;
    cerr << the_time << " loop iteration begin" << endl;
    the_time = clock();
#endif

    //! 1. inject beam
    for (unsigned int i=0; i < params->beam_number_bunches; i++)
	c_bunches[i]->bunch_inject(params->time);
#ifdef PERF_DEBUG
    cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: c_bunches[i]->bunch_inject" << endl;
    the_time = clock();
#endif

    //! 2. Calculate H field
    hfield->calc_field(efield, params->time);
#ifdef PERF_DEBUG
    cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: hfield->calc_field" << endl;
    the_time = clock();
#endif

    //! 3. Calculate v
    c_current->reset_j();
#ifdef PERF_DEBUG
    cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: c_current->reset_j" << endl;
    the_time = clock();
#endif
    c_rho_old->reset_rho();
#ifdef PERF_DEBUG
    cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: c_rho_old->reset_rho" << endl;
    the_time = clock();
#endif
    c_rho_bunch->reset_rho();
#ifdef PERF_DEBUG
    cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: c_rho_bunch->reset_rho" << endl;
    the_time = clock();
#endif
    p_list->step_v(efield, hfield, params->time);
#ifdef PERF_DEBUG
    cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: p_list->step_v" << endl;
    the_time = clock();
#endif

    //! 4. Calculate x, calculate J
    p_list->dump_position_to_old();
    // cout << 4 << endl;
    //! FIXME: for some reason charge_weighting has no effect on result
    // p_list->charge_weighting(c_rho_old); //continuity equation
#ifdef PERF_DEBUG
    cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: p_list->dump_position_to_old" << endl;
    the_time = clock();
#endif

    p_list->half_step_pos(params->time);
#ifdef PERF_DEBUG
    cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: p_list->half_step_pos" << endl;
    the_time = clock();
#endif
    p_list->back_position_to_rz();
#ifdef PERF_DEBUG
    cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: p_list->back_position_to_rz" << endl;
    the_time = clock();
#endif

    p_list->azimuthal_j_weighting(params->time, c_current);
#ifdef PERF_DEBUG
    cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: p_list->azimuthal_j_weighting" << endl;
    the_time = clock();
#endif

    p_list->half_step_pos(params->time);
#ifdef PERF_DEBUG
    cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: p_list->half_step_pos" << endl;
    the_time = clock();
#endif
    p_list->back_position_to_rz();
#ifdef PERF_DEBUG
    cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: p_list->back_position_to_rz" << endl;
    the_time = clock();
#endif

    p_list->j_weighting(params->time, c_current);
#ifdef PERF_DEBUG
    cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: j_weighting" << endl;
    the_time = clock();
#endif

    p_list->back_velocity_to_rz();
#ifdef PERF_DEBUG
    cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: back_velocity_to_rz" << endl;
    the_time = clock();
#endif

    //! 5. Calculate E
    // maxwell_rad.probe_mode_exitation(&geom1,&current1, 1,7e8, time1.current_time);
    efield->calc_field(hfield, params->time, c_current);
#ifdef PERF_DEBUG
    cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: efield->calc_field" << endl;
    the_time = clock();
#endif

    //! 6. Continuity equation
    c_rho_new->reset_rho(); // TODO: is it 4?
#ifdef PERF_DEBUG
    cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: c_rho_new->reset_rho" << endl;
    the_time = clock();
#endif

    //! FIXME: for some reason charge_weighting has no effect on result
    // p_list->charge_weighting(c_rho_new); // continuity equation

    //! print header on every 20 logging steps
    if  ((int)(params->time->current_time / params->time->delta_t) % (params->dump_data_interval * 20) == 0)
    {
      cout << endl
           << left << setw(8) << "Step"
           << left << setw(13) << "Saved Frame"
           << left << setw(18) << "Model Time (sec)"
           << left << setw(32) << "Avg. Step Calculation Time (sec)"
           << endl;
    }

    //! dump data to corresponding files every `parameters.xml->file_save_parameters->data_dump_interval` steps
    if  ((int)(params->time->current_time / params->time->delta_t) % params->dump_data_interval == 0)
    {
      sprintf(avg_step_exec_time, "%.2f", (double)(time(0) - t1) / params->dump_data_interval);

      cout << left << setw(8) << step_number * params->dump_data_interval
           << left << setw(13) << step_number
           << left << setw(18) << params->time->current_time
           << left << setw(32) << avg_step_exec_time
           << endl;

      for (unsigned int i=0; i < params->beam_number_bunches; i++)
	c_bunches[i]->charge_weighting(c_rho_bunch);
#ifdef PERF_DEBUG
      cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: c_bunches[i]->charge_weighting" << endl;
      the_time = clock();
#endif
      // c_rho_old->reset_rho();
      // p_list[0].charge_weighting(c_rho_old);

      dump_data(step_number);
#ifdef PERF_DEBUG
      cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: dump_data" << endl;
      the_time = clock();
#endif

      step_number += 1;
      t1 = time(0);
      if  ((params->dump_system_state_interval != 0)
           && (((int)(params->time->current_time/params->time->delta_t)) % params->dump_system_state_interval == 0)
           && (step_number != 1))
      {
        cout << "Saving system state at " << params->time->current_time << endl;
        dump_system_state();
#ifdef PERF_DEBUG
	cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: dump_system_state" << endl;
	the_time = clock();
#endif
      }
    }
    params->time->current_time = params->time->current_time + params->time->delta_t;
  }
  cout << endl << "Simulation Completed" << endl << endl;
}
