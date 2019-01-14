#include "loadInitParam.h"

#include <typeinfo>  //for 'typeid' to work

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
  omp_set_dynamic(1); // Explicitly enable dynamic teams
  cerr << "INFO! Number of calculation processors changing dynamically" << endl;
#else
  int cores = omp_get_num_procs();
  cerr << "INFO! Number of calculation processors: " << cores << endl;
  omp_set_dynamic(0); // Explicitly disable dynamic teams
  omp_set_num_threads(cores); // Use 4 threads for all consecutive parallel regions
#endif
#endif
  cerr << endl;

  //! Steps to initialize:

  //! 1. read XML config file
  cerr << "INFO! Reading configuration file ``" << xml_file_name << "``" << endl;
  params = new Parameters(xml_file_name);
  // read_xml(xml_file_name);

  //! 2. load File Saving Paths
  cerr << "INFO! Initializing File Saving Paths" << endl;
  //! initialize file saving parameters, like path to computed data files,
  //! path to system state data files max frames number, placed to one file etc.
  // c_io_class = new InputOutputClass (params->dump_result_path, params->dump_save_state_path);

  if (params->use_hdf5)
  {
#ifdef USE_HDF5
#ifndef LEGACY
    cerr << "ERROR! HDF5 support still not implemented" << endl;
    exit(1);
#else
    c_io_class = new IOHDF5 (params->dump_result_path, params->dump_result_path, params->dump_compress);
#endif
#else
    cerr << "ERROR! PDP3 build without HDF5 support. Please, update your configfile to use plaintext ``<use_hdf5>false</use_hdf5>''" << endl;
    exit(1);
#endif
  }
  else
  {
#ifndef LEGACY
    for (auto i = params->probes.begin(); i != params->probes.end(); ++i)
    {
      char* dump_data_root = new char[100];
      sprintf(dump_data_root, "%s/%s", params->dump_result_path, params->dump_data_root);

      ProbePlain *wp = new ProbePlain(
        dump_data_root, (char*)i->component,
        i->type, i->r_start, i->z_start, i->r_end, i->z_end,
        params->dump_compress, params->dump_compress_level,
        i->schedule);

      c_probes.push_back(wp);
    }
#else
    c_io_class = new IOText (params->dump_result_path, params->dump_result_path, params->dump_compress);
#endif
  }
#ifdef LEGACY
  c_io_class->compress_level = params->dump_compress_level;
#endif

  //! 2. creating field objects
  cerr << "INFO! Initializing E/M Fields Data" << endl;
  init_fields ();

  //! 3. load particle parameters
  cerr << "INFO! Initializing Particles Data" << endl;
  init_particles ();
  current_bunch_number = 0;

  //! 4. load boundary conditions
  cerr << "INFO! Initializing Boundary Conditions Data" << endl;
  init_boundary ();

  cerr << "INFO! Initialization complete" << endl;
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

  for (auto i = params->particle_species.begin(); i != params->particle_species.end(); ++i)
  {
    particle_specie p_p = *i;

    cerr << "      Initializing " << p_p.name << " Data" << endl;

    //! init and setup particles properties
    prtls = new Particles(p_p.name, p_p.charge, p_p.mass, p_p.number_macro, params->geom);
		p_list->part_list.push_back(prtls); // push particles to particles list vector

    // case 2 is for cylindrical distribution
    prtls->load_cylindrical_spatial_distribution(p_p.left_density, p_p.right_density, 0);
    prtls->velocity_distribution(p_p.temperature);
  }

  p_list->charge_weighting(c_rho_new);
}

void LoadInitParam::init_beam()
{
  //! initialize particles bunch data
  //! particles bunch should be injected to plasma

  double duration = params->bunch_lenght / params->beam_initial_velocity;
  double hole_duration = params->beam_bunches_distance / params->beam_initial_velocity;
  double bunch_inject_time = current_bunch_number * (duration + hole_duration);
  char *current_bunch_name =  new char[100];

  char cbn_str_int[100];
  sprintf(cbn_str_int, "%d", current_bunch_number);

  strcpy(current_bunch_name, params->beam_name);
  strcat(current_bunch_name, "#");
  strcat(current_bunch_name, cbn_str_int);

  double estimated_bunch_erase_number = (
    params->time->current_time
    - params->geom->z_size / params->beam_initial_velocity
    - duration )
    / (duration + hole_duration);

  // erase gone bunches
  if (c_bunches.size() > 1) { // strange gcc behavior: size of empty vector of instance references is not zero
    for (auto i = c_bunches.begin(); i != c_bunches.end(); ++i)
    {
      Bunch *cb = *i;
      if (cb->bunch_number < estimated_bunch_erase_number)
      {
        bool is_dead = true;
        for (unsigned int i; i < cb->number; ++i)
          if (cb->is_alive[i])
          {
            is_dead = false;
            break;
          }
        if (is_dead)
	      {
#ifdef DEBUG
          cerr << "Bunch #" << cb->bunch_number << " is out of simulation. Clearing" << endl;
#endif
          for (auto k = p_list->part_list.begin(); k != p_list->part_list.end(); ++k)
            if (*i == *k)
              p_list->part_list.erase(k);
          c_bunches.erase(i);
          delete cb;
	      }
      }
    }
  }
  if (params->time->current_time >= bunch_inject_time // it's not equality, because of distreate simulation nature
      && current_bunch_number < params->beam_number_bunches)
  {
#ifdef DEBUG
    cerr << "Injecting Bunch #" << current_bunch_number << endl;
#endif
    double duration = params->bunch_lenght / params->beam_initial_velocity;
    double hole_duration = params->beam_bunches_distance / params->beam_initial_velocity;
    Bunch *prtls = new Bunch(current_bunch_name,
                             params->beam_particle_charge,
                             params->beam_particle_mass,
                             params->bunch_number_macro,
                             params->geom, // p_list,
                             duration,
                             params->bunch_radius,
                             params->bunch_density,
                             params->beam_initial_velocity,
                             current_bunch_number,
                             hole_duration);
    p_list->part_list.push_back(prtls); // push bunch to particles list vector
    c_bunches.push_back(prtls); // push bunch to bunches list vector

    // calculate next bunch number and inject time
    ++current_bunch_number; // increase current bunch number after creation
  }

  // inject bunch
  for (auto i = c_bunches.begin(); i != c_bunches.end(); ++i)
    (*i)->bunch_inject(params->time);
}

void LoadInitParam::init_boundary ()
{
  //! initialize boundaries and boundary conditions

  // Maxwell initial conditions
  BoundaryMaxwellConditions maxwell_rad(efield); // TODO: WTF?
  maxwell_rad.specify_initial_field(params->geom,
                                    params->boundary_maxwell_e_phi_upper,
                                    params->boundary_maxwell_e_phi_left,
                                    params->boundary_maxwell_e_phi_right);

  if (params->boundary_conditions == 0)
  {
    p_list->charge_weighting(c_rho_new);

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
  efield->set_homogeneous_efield(0., 0., 0.);
  hfield->set_homogeneous_h(0.0, 0., 0.);
  efield->set_fi_on_z();
}

void LoadInitParam::print_data(int probe_type, char* component, int step_number, int dump_interval, int* shape)
{
  int print_header_step = 20;
  if (printed_step % print_header_step == 0)
  {

    cout << endl
         << left << setw(8) << "Step"
         << left << setw(13) << "Saved Frame"
         << left << setw(20) << "Dumping Probe Name"
         << left << setw(21) << "Shape"
         << left << setw(18) << "Model Time (sec)"
         << left << setw(21) << "Simulation Duration"
         << endl;
  }

  char type_comp[100];
  char probe_shape[100];

  switch (probe_type)
  {
  case 0:
    sprintf(type_comp, "frame/%s", component);
    sprintf(probe_shape, "[%i,%i,%i,%i]", shape[0], shape[1], shape[2], shape[3]);
    break;
  case 1:
    sprintf(type_comp, "column/%s", component);
    sprintf(probe_shape, "[%i]", shape[1]);
    break;
  case 2:
    sprintf(type_comp, "row/%s", component);
    sprintf(probe_shape, "[%i]", shape[0]);
    break;
  case 3:
    sprintf(type_comp, "dot/%s", component);
    sprintf(probe_shape, "[%i,%i]", shape[0], shape[1]);
    break;
  case 4:
    sprintf(type_comp, "mpframe/%s", component);
    sprintf(probe_shape, "[%i,%i,%i,%i]", shape[0], shape[1], shape[2], shape[3]);
    break;
  }

  cout << left << setw(8) << step_number * dump_interval
       << left << setw(13) << step_number
       << left << setw(20) << type_comp
       << left << setw(21) << probe_shape
       << left << setw(18) << params->time->current_time
       << left << setw(18) << lib::get_simulation_duration()
       << endl;
  ++printed_step;
}

void LoadInitParam::dump_data(int step_number)
{
  //! dump claculated data frames (for fields etc.) to files set

  // electrical fields
#ifndef LEGACY
  for (auto i = c_probes.begin(); i != c_probes.end(); ++i)
  {
    Probe *w = *i;

    if ((int)(params->time->current_time / params->time->delta_t) % w->schedule == 0)
    {
      int file_number;
      if (w->type == 4)
        file_number = w->step_number;
      else
        file_number = floor(w->step_number / params->dump_frames_per_file);

      char file_name[100];
      sprintf(file_name, "%d", file_number);
      int shape[4] = {w->start_r, w->start_z, w->end_r, w->end_z};

      print_data(w->type, w->component, w->step_number, w->schedule, shape);

      if (strcmp(w->component, "E_r") == 0)
        w->write(file_name, efield->field_r);
      else if (strcmp(w->component, "E_phi") == 0)
        w->write(file_name, efield->field_phi);
      else if (strcmp(w->component, "E_z") == 0)
        w->write(file_name, efield->field_z);
      //
      else if (strcmp(w->component, "H_r") == 0)
        w->write(file_name, hfield->field_r);
      else if (strcmp(w->component, "H_phi") == 0)
        w->write(file_name, hfield->field_phi);
      else if (strcmp(w->component, "H_z") == 0)
        w->write(file_name, hfield->field_z);
      //
      else if (strcmp(w->component, "J_r") == 0)
        w->write(file_name, c_current->get_j_r());
      else if (strcmp(w->component, "J_phi") == 0)
        w->write(file_name, c_current->get_j_phi());
      else if (strcmp(w->component, "J_z") == 0)
        w->write(file_name, c_current->get_j_z());
      //
      else if (strcmp(w->component, "rho_beam") == 0)
        w->write(file_name, c_rho_bunch->get_rho());
      // separate logic for velocity frames
      else if (w->type == 4)
        for (auto i = p_list->part_list.begin(); i != p_list->part_list.end(); ++i)
        {
          if (strcmp((*i)->name, w->component) == 0)
            w->mpwrite(file_name, *i);
        }
      else
        cerr << "WARNING! probe for component ``" << w->component << "'' does not exist. Skipping" << endl;

      ++w->step_number;
    }
  }
#else
  if (params->dump_e_r)
    c_io_class->out_data("E_r", efield->field_r, step_number,
                         params->dump_frames_per_file, params->geom->n_grid_r - 1, params->geom->n_grid_z - 1);

  if (params->dump_e_phi)
    c_io_class->out_data("E_phi", efield->field_phi, step_number,
                         params->dump_frames_per_file, params->geom->n_grid_r - 1, params->geom->n_grid_z - 1);

  if (params->dump_e_z)
    c_io_class->out_data("E_z", efield->field_z, step_number,
                         params->dump_frames_per_file, params->geom->n_grid_r - 1, params->geom->n_grid_z - 1);

  // magnetic fields
  if (params->dump_h_r)
    c_io_class->out_data("H_r", hfield->field_r, step_number,
                         params->dump_frames_per_file, params->geom->n_grid_r - 1, params->geom->n_grid_z - 1);

  if (params->dump_h_phi)
    c_io_class->out_data("H_phi", hfield->field_phi, step_number,
                         params->dump_frames_per_file, params->geom->n_grid_r - 1, params->geom->n_grid_z - 1);

  if (params->dump_h_z)
    c_io_class->out_data("H_z", hfield->field_z, step_number,
                         params->dump_frames_per_file, params->geom->n_grid_r - 1, params->geom->n_grid_z - 1);

  // currents
  if (params->dump_current_r)
    c_io_class->out_data("J_r", c_current->get_j_r(), step_number,
                         params->dump_frames_per_file, params->geom->n_grid_r - 1, params->geom->n_grid_z - 1);

  if (params->dump_current_phi)
    c_io_class->out_data("J_phi", c_current->get_j_phi(), step_number,
                         params->dump_frames_per_file, params->geom->n_grid_r - 1, params->geom->n_grid_z - 1);

  if (params->dump_current_z)
    c_io_class->out_data("J_z", c_current->get_j_z(), step_number,
                         params->dump_frames_per_file, params->geom->n_grid_r - 1, params->geom->n_grid_z - 1);

  // beam density
  if (params->dump_rho_beam)
    c_io_class->out_data("rho_beam", c_rho_bunch->get_rho(), step_number,
                         params->dump_frames_per_file, params->geom->n_grid_r - 1, params->geom->n_grid_z - 1);

#endif
}

void LoadInitParam::run(void)
{
  //! Launch calculation.
  //! THIS IS ENTRY POINT TO MAIN PDP3 CALCULATION LOOP
  cout << endl << "Launch Simulation" << endl << endl;

#ifdef DEBUG
  clock_t the_time = 0;
#endif
#ifdef LEGACY
  time_t t1 = time(0);
  char avg_step_exec_time[24]; // rounded and formatted average step execution time
#endif
  int step_number = 0;

  //! Main calculation loop
  while (params->time->current_time < params->time->end_time)
  {
    //! Steps:
#ifdef DEBUG
    the_time = 0;
    cerr << the_time << " sec. loop iteration begin" << endl;
    the_time = clock();
#endif

    //! 1. inject beam
    init_beam();
#ifdef DEBUG
    cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: c_bunches[i]->bunch_inject" << endl;
    the_time = clock();
#endif

    //! 2. Calculate H field
    hfield->calc_field(efield, params->time);
#ifdef DEBUG
    cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: hfield->calc_field" << endl;
    the_time = clock();
#endif

    //! 3. Calculate v
    c_current->reset_j();
#ifdef DEBUG
    cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: c_current->reset_j" << endl;
    the_time = clock();
#endif

    c_rho_old->reset_rho();
#ifdef DEBUG
    cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: c_rho_old->reset_rho" << endl;
    the_time = clock();
#endif

    c_rho_bunch->reset_rho();
#ifdef DEBUG
    cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: c_rho_bunch->reset_rho" << endl;
    the_time = clock();
#endif

    p_list->step_v(efield, hfield, params->time);
#ifdef DEBUG
    cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: p_list->boris_pusher" << endl;
    the_time = clock();
#endif

#ifdef DEBUG
    cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: full_current_distribution" << endl;
    the_time = clock();
#endif
    p_list->full_current_distribution(c_current, params->time);

    //! 5. Calculate E
    // maxwell_rad.probe_mode_exitation(&geom1,&current1, 1,7e8, time1.current_time);
    efield->calc_field(hfield, params->time, c_current);
#ifdef DEBUG
    cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: efield->calc_field" << endl;
    the_time = clock();
#endif

    //! 6. Continuity equation
    c_rho_new->reset_rho(); // TODO: is it 4?
#ifdef DEBUG
    cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: c_rho_new->reset_rho" << endl;
    the_time = clock();
#endif

    //! FIXME: for some reason charge_weighting has no effect on result
    // p_list->charge_weighting(c_rho_new); // continuity equation

#ifdef LEGACY
    //! print header on every 20 logging steps
    if  ((int)(params->time->current_time / params->time->delta_t) % (params->dump_data_interval * 20) == 0)
    {
      cout << endl
           << left << setw(8) << "Step"
           << left << setw(13) << "Saved Frame"
           << left << setw(18) << "Model Time (sec)"
           << left << setw(18) << "Simulation Time"
           << left << setw(34) << "Avg. Step Calculation Time"
           << endl;
    }

    //! dump data to corresponding files every `parameters.xml->file_save_parameters->data_dump_interval` steps
    if  ((int)(params->time->current_time / params->time->delta_t) % params->dump_data_interval == 0)
    {
#ifdef _OPENMP
      sprintf(avg_step_exec_time, "%.2fs", (double)(time(0) - t1) / params->dump_data_interval / omp_get_num_threads());
#else
      sprintf(avg_step_exec_time, "%.2fs", (double)(time(0) - t1) / params->dump_data_interval);
#endif

      cout << left << setw(8) << step_number * params->dump_data_interval
           << left << setw(13) << step_number
           << left << setw(18) << params->time->current_time
           << left << setw(18) << lib::get_simulation_duration()
           << left << setw(34) << avg_step_exec_time
           << endl;

      for (auto i = c_bunches.begin(); i != c_bunches.end(); ++i)
        (*i)->charge_weighting(c_rho_bunch);

#ifdef DEBUG
      cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: c_bunches[i]->charge_weighting" << endl;
      the_time = clock();
#endif

      dump_data(step_number);

#ifdef DEBUG
      cerr << ((double)(clock() - the_time) / (double)CLOCKS_PER_SEC) << " sec. for: dump_data" << endl;
      the_time = clock();
#endif

      ++step_number;
      t1 = time(0);
    }
#endif

#ifndef LEGACY
    for (auto i = c_bunches.begin(); i != c_bunches.end(); ++i)
      (*i)->charge_weighting(c_rho_bunch);
    dump_data(step_number);
#endif
    params->time->current_time = params->time->current_time + params->time->delta_t;
  }
  cout << endl << "Simulation Completed" << endl << endl;
}
