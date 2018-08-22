#include "parameters.h"

Parameters::Parameters(const char *xml_file_name)
{

  //! read given xml file and fill xml_data class variable
  XMLDocument* xml_document = new XMLDocument(true, COLLAPSE_WHITESPACE);

  XMLError e_result = xml_document->LoadFile(xml_file_name);
  if (e_result != XML_SUCCESS)
  {
    cerr << "ERROR: Can not read configuration file ``" << xml_file_name << "``" << endl;
    exit (78);
  }

  xml_data = xml_document->FirstChildElement("initial_parameters");

  //! set debug option
  debug = lib::to_bool(try_first_child(xml_data, "debug")->GetText());

  //! check if use hdf5 backend
  use_hdf5 = lib::to_bool(try_first_child(xml_data, "use_hdf5")->GetText());

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

XMLElement* Parameters::try_first_child(XMLElement* element, const char* name)
{

  XMLElement* f_child = element->FirstChildElement(name);
  if (f_child == nullptr)
  {
    cerr << "ERROR: Can not get child, named ``"
         << name
         << "'' from XML element ``"
         << element->Name()
         << "''. Seems, configuration file is corrupted, or inconsistent"
         << endl;
    abort();
  }
  return f_child;
}

void Parameters::init_particles()
{
  //! initialize particles charge, mass, number,
  //! density and temperature for all particle species
  const char *p_specie_section_name = "particle_specie";
  XMLElement *specie = try_first_child(try_first_child(xml_data, "particles"), p_specie_section_name);

  while(specie)
  {
    particle_specie p_specie;

    p_specie.name = (char*)specie->Attribute("name");

    p_specie.charge = atoi(try_first_child(specie, "charge")->GetText());
    p_specie.mass = atoi(try_first_child(specie, "mass")->GetText());
    p_specie.number_macro = debug ?
      atof(try_first_child(specie, "debug_number_macro")->GetText()) :
      atof(try_first_child(specie, "number_macro")->GetText());
    p_specie.left_density = atof(try_first_child(specie, "left_density")->GetText());
    p_specie.right_density = atof(try_first_child(specie, "right_density")->GetText());
    p_specie.temperature = atof(try_first_child(specie, "temperature")->GetText());

    particle_species.push_back(p_specie);

    specie = specie->NextSiblingElement(p_specie_section_name);
  }
}

void Parameters::init_beam()
{
  //! initialize particles bunch data
  //! particles bunch should be injected to plasma
  XMLElement *sub_root = try_first_child(xml_data, "particle_beam");

  beam_name = (char*)sub_root->Attribute("name");

  beam_particle_charge = atoi(try_first_child(sub_root, "charge")->GetText());
  beam_particle_mass = atoi(try_first_child(sub_root, "mass")->GetText());
  beam_number_bunches = atoi(try_first_child(sub_root, "number_bunches")->GetText());
  beam_bunches_distance = atof(try_first_child(sub_root, "bunches_distance")->GetText());
  beam_initial_velocity = atof(try_first_child(sub_root, "initial_velocity")->GetText());
  //
  bunch_number_macro = debug ?
    atof(try_first_child(sub_root, "debug_bunch_number_macro")->GetText()) :
    atof(try_first_child(sub_root, "bunch_number_macro")->GetText());

  bunch_lenght = atof(try_first_child(sub_root, "bunch_length")->GetText());
  bunch_radius = atof(try_first_child(sub_root, "bunch_radius")->GetText());
  bunch_density = atof(try_first_child(sub_root, "bunch_density")->GetText());
}

void Parameters::init_boundary ()
{
  //! initialize boundaries and boundary conditions
  XMLElement *sub_root = try_first_child(xml_data, "boundary_maxwell_conditions");
  //
  boundary_maxwell_e_phi_upper = atof(try_first_child(sub_root, "e_fi_upper")->GetText());
  boundary_maxwell_e_phi_left = atof(try_first_child(sub_root, "e_fi_left")->GetText());
  boundary_maxwell_e_phi_right = atof(try_first_child(sub_root, "e_fi_right")->GetText());
  boundary_conditions = atoi(try_first_child(xml_data, "boundary_conditions")->GetText());
}

void Parameters::init_geometry ()
{
  //! initialize geometry
  XMLElement *sub_root = try_first_child(xml_data, "geometry");

  //
  double r_size = atof(try_first_child(sub_root, "r_size")->GetText());
  double z_size = atof(try_first_child(sub_root, "z_size")->GetText());
  int n_grid_r = debug ?
    atoi(try_first_child(sub_root, "debug_n_grid_r")->GetText()) :
    atoi(try_first_child(sub_root, "n_grid_r")->GetText());
  int n_grid_z = debug ?
    atoi(try_first_child(sub_root, "debug_n_grid_z")->GetText()) :
    atoi(try_first_child(sub_root, "n_grid_z")->GetText());
  // PML

  geom = new Geometry(r_size, z_size, n_grid_r, n_grid_z);
  geom->set_epsilon();
}

void Parameters::init_pml()
{
  XMLElement *sub_root = try_first_child(xml_data, "pml");

  double comp_l_1 = atof(try_first_child(sub_root, "comparative_l_1")->GetText());
  double comp_l_2 = atof(try_first_child(sub_root, "comparative_l_2")->GetText());
  double comp_l_3 = atof(try_first_child(sub_root, "comparative_l_3")->GetText());
  double sigma_1_t = atof(try_first_child(sub_root, "sigma_1")->GetText());
  double sigma_2_t = atof(try_first_child(sub_root, "sigma_2")->GetText());

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

void Parameters::init_time ()
{
  //! initialize time parameters: `start_time`, `relaxation_time` `end_time` and `delta_t`
  XMLElement *sub_root = try_first_child(xml_data, "time");
  //
  double start_time = atof(try_first_child(sub_root, "start_time")->GetText());
  double relaxation_time = atof(try_first_child(sub_root, "relaxation_time")->GetText());
  double current_time = atof(try_first_child(sub_root, "current_time")->GetText());
  double end_time = atof(try_first_child(sub_root, "end_time")->GetText());
  double step_interval = atof(try_first_child(sub_root, "step_interval")->GetText());

  time = new Time(current_time, start_time,
                  relaxation_time, end_time, step_interval);
}

void Parameters::init_data_dump_parameters ()
{
  //! initialize file saving parameters, like path to computed data files,
  //! path to system state data files max frames number, placed to one file etc.
  XMLElement *sub_root = try_first_child(xml_data, "file_save_parameters");
  XMLElement *dump_data_root = try_first_child(sub_root, "dump_data");

  dump_result_path = (char*)try_first_child(sub_root, "result_path")->GetText();
  dump_save_state_path = (char*)try_first_child(sub_root, "save_state_path")->GetText();

  dump_data_interval = debug ?
    atoi(try_first_child(sub_root, "debug_data_dump_interval")->GetText()) :
    atoi(try_first_child(sub_root, "data_dump_interval")->GetText());

  dump_system_state_interval = atoi(try_first_child(sub_root, "system_state_dump_interval")->GetText());
  dump_frames_per_file = atoi(try_first_child(sub_root, "frames_per_file")->GetText());

  //! choose, which parameters should be dumped to data files
  dump_e_r = lib::to_bool(try_first_child(dump_data_root, "E_r")->GetText());
  dump_e_phi = lib::to_bool(try_first_child(dump_data_root, "E_phi")->GetText());
  dump_e_z = lib::to_bool(try_first_child(dump_data_root, "E_z")->GetText());

  dump_h_r = lib::to_bool(try_first_child(dump_data_root, "H_r")->GetText());
  dump_h_phi = lib::to_bool(try_first_child(dump_data_root, "H_phi")->GetText());
  dump_h_z = lib::to_bool(try_first_child(dump_data_root, "H_z")->GetText());

  dump_current_r = lib::to_bool(try_first_child(dump_data_root, "J_r")->GetText());
  dump_current_phi = lib::to_bool(try_first_child(dump_data_root, "J_phi")->GetText());
  dump_current_z = lib::to_bool(try_first_child(dump_data_root, "J_z")->GetText());

  dump_position = lib::to_bool(try_first_child(dump_data_root, "position")->GetText());
  dump_velocity = lib::to_bool(try_first_child(dump_data_root, "velocity")->GetText());

  dump_rho_beam = lib::to_bool(try_first_child(dump_data_root, "rho_bunch")->GetText());
}
