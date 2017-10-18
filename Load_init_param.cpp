#include "Load_init_param.h"
#include "time.h"
#include "tinyxml2.h"

// #include "kern_accessor.h"
//extern KernAccessor *kern_access_global;

using namespace tinyxml2;

#define INITIAL_PARAMS_NAME "Initial_parameters"
#define PML_PARAMS_NAME "PML"
#define GEOMETRY_PARAMS_NAME "geometry"
#define TIME_PARAMS_NAME "Time"
#define PARTICLES_PARAMS_NAME "Particles"
#define BUNCH_PARAMS_NAME "Particles_bunch"
#define BOUNDARY_MAXWELL_PARAMS_NAME "Boundary_Maxwell_conditions"

Load_init_param::Load_init_param(void)
{
}

Load_init_param::Load_init_param(char* xml_file_name)
{
		// NOTE: all system init is too huge, so
		// so we use several "subconstructors"
		// to initialise different parts of system

	// read XML config file
	read_xml (xml_file_name);

	// load PML parameters
	init_pml ();

	// load Geometry parameters
	init_geometry ();

	// creating field objects
	init_fields ();

	// load time parameters
	init_time ();

	// load particle parameters
	p_list = new particles_list(0);
	init_particles();

	// load bunch
	c_bunch = init_bunch();

	//Maxwell initial conditions///
	init_boundary_maxwell();

	// creating rho and current arrays
	c_rho_new = new charge_density(c_geom);
	c_rho_old = new charge_density (c_geom);
	c_rho_beam = new charge_density (c_geom);
	c_current = new current(c_geom);

	// boundary conditions
	char* a = read_char("Boundary_conditions");
	if (atoi(a)==0)
	{
		p_list->charge_weighting(c_rho_new);
		Fourier four1(0);
		Poisson_dirichlet dirih(c_geom);
		dirih.poisson_solve(efield, c_rho_new);
	}

	// load File Path
	char* path_res = read_char("PathtoResult");
	char* path_dump = read_char("PathtoSaveState");
	c_io_class	= new input_output_class (path_res , path_dump);

	//////////////////////////////////////////////////

	//	KernAccessor *kern_access = new KernAccessor(c_geom->n_grid_1, c_geom->n_grid_2);
	//	kern_access_global = kern_access;
	//

	cout << "Initialisation complete\n";
}

Load_init_param::~Load_init_param(void)
{
}

void Load_init_param::read_xml(const char* xml_file_name)
{
	cout << "Reading configuration file ``" << xml_file_name << "``\n";

	xml_data = new XMLDocument(true, COLLAPSE_WHITESPACE);

	XMLError e_result = xml_data->LoadFile(xml_file_name);
	if (e_result != XML_SUCCESS)
	{
		cerr << "Can not read configuration file ``" << xml_file_name << "``\n";
		exit (78);
	}
}

char* Load_init_param::read_char(char* p_name)
{
	int number = 0;
	XMLElement* root = xml_data->FirstChildElement (INITIAL_PARAMS_NAME);
	XMLElement* sub_root =root->FirstChildElement (p_name);
	char* vul_arr = new char [50];
	strcpy(vul_arr, sub_root->GetText());

	return vul_arr;
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

	while(particle_kind)
	{
		strcpy (p_name,particle_kind->FirstChildElement("name")->GetText());
		cout << "Initialising " << p_name << " data\n";

		charge = atof(particle_kind->FirstChildElement("charge")->GetText());
		mass = atof(particle_kind->FirstChildElement("mass")->GetText());
		number = atof(particle_kind->FirstChildElement("number")->GetText());
		left_density = atof(particle_kind->FirstChildElement("left_density")->GetText());
		right_density = atof(particle_kind->FirstChildElement("right_density")->GetText());
		temperature = atof(particle_kind->FirstChildElement("temperature")->GetText());

		prtls = new Particles(strcpy(new char [50],p_name), charge, mass, number,
													c_geom, p_list);

		prtls->load_spatial_distribution_with_variable_mass(
			left_density,right_density,0,0);
		// prtls->load_spatial_distribution(params[3],params[4],0,0);
		prtls->velocity_distribution_v2(temperature);
		particle_kind = particle_kind->NextSiblingElement(p_king_section_name);
	}

	return;
}

Bunch* Load_init_param::init_bunch()
{
	// initialise particles bunch data
	XMLElement* root = xml_data->FirstChildElement (INITIAL_PARAMS_NAME);
	XMLElement* sub_root =root->FirstChildElement (BUNCH_PARAMS_NAME);
	cout << "Initialising particle bunch data\n";
	Bunch* prtls = 0;
	char* p_name= new char [50];
	strcpy (p_name, sub_root->FirstChildElement("name")->GetText());
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

	prtls = new Bunch(p_name, charge, mass,number, c_geom, p_list,
										duration, radius, density, initial_velocity);

	return prtls;
}

void Load_init_param::init_boundary_maxwell ()
{
	cout << "Initialising bounrary maxwell conditions data\n";
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

	Boundary_Maxwell_conditions maxwell_rad(efield);
	maxwell_rad.specify_initial_field(c_geom, e_fi_upper, e_fi_left, e_fi_right);
}

void Load_init_param::init_pml ()
{
	cout << "Initialising PML data\n";
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

	// TODO: load_pml should be first. Avoid it
	c_geom = new Geometry(r_size, z_size, n_grid_r, n_grid_z, c_pml);
	c_geom->set_epsilon();
}

void Load_init_param::init_fields ()
{
	efield = new E_field(c_geom);
	hfield = new H_field(c_geom);
	efield->boundary_conditions();
	efield->set_homogeneous_efield(0.0, 0.0, 0); // TODO: unhardcode
	hfield->set_homogeneous_h(0.0, 0.0, 0.0); // TODO: unhardcode
	efield->set_fi_on_z();
}

void Load_init_param::init_time ()
{
	cout << "Initialising time data\n";
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

bool Load_init_param::save_system_state()
{
	c_io_class->out_field_dump("e1",efield->e1,c_geom->n_grid_1-1,c_geom->n_grid_2-1);
	c_io_class->out_field_dump("e2",efield->e2,c_geom->n_grid_1-1,c_geom->n_grid_2-1);
	c_io_class->out_field_dump("e3",efield->e3,c_geom->n_grid_1-1,c_geom->n_grid_2-1);
	c_io_class->out_field_dump("h1",hfield->h1,c_geom->n_grid_1-1,c_geom->n_grid_2-1);
	c_io_class->out_field_dump("h2",hfield->h2,c_geom->n_grid_1-1,c_geom->n_grid_2-1);
	c_io_class->out_field_dump("h3",hfield->h3,c_geom->n_grid_1-1,c_geom->n_grid_2-1);
	for(int i=0;i<p_list->part_list.size();i++)
	{
		c_io_class->out_coord_dump(p_list->part_list[i]->name,p_list->part_list[i]->x1, p_list->part_list[i]->x3, p_list->part_list[i]->number);
		c_io_class->out_velocity_dump(p_list->part_list[i]->name,p_list->part_list[i]->v1, p_list->part_list[i]->v2,p_list->part_list[i]->v3, p_list->part_list[i]->number);

	}
	return true;
}

void Load_init_param::run(void)
{
	cout << "Launching simulation\n";

	this->c_time->current_time = 0.0 ;
	p_list->charge_weighting(c_rho_new);

	Poisson_dirichlet dirih(c_geom);
	dirih.poisson_solve(efield, c_rho_new);

	//variable for out_class function
	p_list->create_coord_arrays();
	int step_number= 0;
	clock_t t1 = clock();

	while (c_time->current_time < c_time->end_time)
	{
		c_bunch->bunch_inject(c_time);

		//1. Calculate H field
		hfield->calc_field(efield, c_time);

		//2. Calculate v
		c_current->reset_j();
		c_rho_old->reset_rho();
		c_rho_beam->reset_rho();
		p_list->step_v(efield, hfield, c_time);

		//3. Calculate x, calculate J
		p_list->copy_coords();
		p_list->charge_weighting(c_rho_old);	//continuity equation
		p_list->half_step_coord(c_time);
		p_list->azimuthal_j_weighting(c_time, c_current);
		p_list->half_step_coord(c_time);
		p_list->j_weighting(c_time,c_current);

		//4. Calculate E
		// maxwell_rad.probe_mode_exitation(&geom1,&current1, 1,7e8, time1.current_time);
		efield->calc_field(hfield,c_time, c_current, c_pml);

		//continuity equation
		c_rho_new->reset_rho();

		p_list->charge_weighting(c_rho_new);	//continuity equation
		//bool res = continuity_equation(c_time, c_geom, c_current, c_rho_old, c_rho_new);

		if	((((int)(c_time->current_time/c_time->delta_t))%5==0))
			//if	((((int)(c_time->current_time/c_time->delta_t)) < 10))
			//if	( abs(time1.current_time - time1.end_time + time1.delta_t) < 1e-13)
		{
			// Logging after step
			cout << "Model step: " << step_number
					 << ", Execution time: "
					 << 1000 * (clock() - t1) / CLOCKS_PER_SEC << " ms"
					 << ", Current time: "
					 << c_time->current_time << "\n";
			c_bunch->charge_weighting(c_rho_beam);
			c_rho_old->reset_rho();
			p_list[0].charge_weighting(c_rho_old);
			c_io_class->out_data("rho_el", c_rho_old->get_ro(),step_number,100,c_geom->n_grid_1-1,c_geom->n_grid_2-1);
			//out_class.out_data("e1",e_field1.e1,100,128,2048);
			c_io_class->out_data("rho_beam", c_rho_beam->get_ro(),step_number,100,c_geom->n_grid_1-1,c_geom->n_grid_2-1);
			c_io_class->out_data("e3",efield->e3,step_number,100,c_geom->n_grid_1-1,c_geom->n_grid_2-1);
			c_io_class->out_data("e1",efield->e1,step_number,100,c_geom->n_grid_1-1,c_geom->n_grid_2-1);
			//c_io_class->out_data("j1",c_current->get_j1(),step_number,100,c_geom->n_grid_1-1,c_geom->n_grid_2-1);
			//c_io_class->out_data("j3",c_current->get_j3(),step_number,100,c_geom->n_grid_1-1,c_geom->n_grid_2-1);
			//c_io_class->out_coord("vels",c_bunch->v1, c_bunch->v3, step_number, 100, c_bunch->number);
			//out_class.out_data("rho",rho_elect.get_ro(),step_number,100,geom1.n_grid_1-1,geom1.n_grid_2-1);
			//out_class.out_coord("vels",electron_bunch.v1, electron_bunch.v3, step_number, 100, electron_bunch.number);
			//out_class.out_coord("coords",electrons.x1, electrons.x3, step_number, 100, electrons.number);
			//out_class.out_coord("vels",electrons.v1, electrons.v3, step_number, 100, electrons.number);
			c_io_class->out_data("h2",hfield->h2,step_number,100,c_geom->n_grid_1-1,c_geom->n_grid_2-1);
			step_number=step_number+1;
			t1 = clock();
			if	((((int)(c_time->current_time/c_time->delta_t))%1000==0)&&(step_number!=1))
				this->save_system_state();
		}
		c_time->current_time = c_time->current_time + c_time->delta_t;
		//if (!res)
		//	cout<<"Error:"<<c_time->current_time<<"! ";

	}
}
