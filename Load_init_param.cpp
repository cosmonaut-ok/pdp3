#include "Load_init_param.h"
#include "time.h"
#include "kern_accessor.h"
//extern KernAccessor *kern_access_global;
Load_init_param::Load_init_param(void)
{
}
Load_init_param::Load_init_param(char* xml_file_l):xml_file(xml_file_l)
{
}


Load_init_param::~Load_init_param(void)
{

}

void Load_init_param::read_xml()
{
	 
	int i=1.2e7;  
	// TiXmlDocument  xml_f("parameters.xml");
	TiXmlDocument  xml_f(xml_file);
	 xml_f.LoadFile();
	 if(xml_f.Error())
	 {
		 cout<< "can't read file";
			 return;
	 }

     TiXmlElement* root = xml_f.FirstChildElement ("Initial_parameters");
	 TiXmlElement* sub_root = root->FirstChildElement ("geometry");
	 char test [30];
	 TiXmlElement* ss_root = sub_root->FirstChildElement();
	const char* a =ss_root->GetText();
	// strcpy( test, sub_root->GetText());
//	 i = get_int_value(sub_root->FirstChildElement( ));
	 
	 cout << i;
	// cin  >> i;
	// TiXmlElement* sub_root = root->FirstChildElement ("geometry");
	// TiXmlElement* value = sub_root->FirstChildElement ();
	// double r_size  = get_double_value(value);
	 //double z_size  = get_double_value(value->NextSiblingElement());


}
int Load_init_param:: get_int_value(const char* r_str)
{
 int i =0;
 	 i = atoi(r_str);
	return i;
}
double Load_init_param:: get_double_value(const char* r_str)
{
 	double i = atof(r_str);
	return i;
}

 char* Load_init_param:: read_char(char* p_name)
{
	TiXmlDocument  xml_f(xml_file);
	 xml_f.LoadFile();
	 if(xml_f.Error())
	 {
		 cout<< "couldn't read file";
			 return 0;
	 }
	 int number= 0;
	 TiXmlElement* root = xml_f.FirstChildElement ("Initial_parameters"); 
    TiXmlElement* sub_root =root->FirstChildElement (p_name); 
	char* vul_arr = new char [50];
 	strcpy(vul_arr, sub_root->GetText());
	 
	   return vul_arr;
}

double* Load_init_param:: read_double_params(const char* p_name)
{
	TiXmlDocument  xml_f(xml_file);
	 xml_f.LoadFile();
	 if(xml_f.Error())
	 {
		 cout<< "couldn't read file";
			 return 0;
	 }
	 int number= 0;
	 TiXmlElement* root = xml_f.FirstChildElement ("Initial_parameters"); 
    TiXmlElement* sub_root =root->FirstChildElement (p_name); 
	TiXmlElement* read_elem =sub_root->FirstChildElement (); 
	  while(read_elem)
        {
       
			read_elem = read_elem->NextSiblingElement();
            number++;
        }
	  double* vul_arr  = new double [number-1];
	 for (int i=0;i<(number);i++)
	 read_elem =sub_root->FirstChildElement (); 

	   for(int i=0;i<(number);i++)
        {
		    vul_arr[i] = get_double_value(read_elem->GetText());
			read_elem=read_elem->NextSiblingElement();
        }	 
	   return vul_arr;
}
///////////////////////////////////////////
void Load_init_param:: read_load_particles()
{
	TiXmlDocument  xml_f(xml_file);
	 xml_f.LoadFile();
	 if(xml_f.Error())
	 {
		 cout<< "couldn't read file";
			 return ;
	 }
	 int number= 0;
	 TiXmlElement* root = xml_f.FirstChildElement ("Initial_parameters"); 
    TiXmlElement* sub_root =root->FirstChildElement ("Particle_Name"); 
	double* params= 0;
	 char* p_name= new char [50];
	 Particles*  prtls =0;
	  while(sub_root)
        {
       
			 const char* test = sub_root->GetText();
			 strcpy (p_name,sub_root->GetText());
			params = read_double_params(p_name);
			 prtls = new Particles(strcpy(new char [50],p_name),params,c_geom, p_list);
			 prtls->load_spatial_distribution_with_variable_mass(params[3],params[4],0,0);
			// prtls->load_spatial_distribution(params[3],params[4],0,0);
			 prtls->velocity_distribution_v2(params[5]);
			 sub_root = sub_root->NextSiblingElement("Particle_Name");      
        }
	 
	   return;
}
/////////////////////////////////////////

Bunch* Load_init_param:: read_load_bunch()
{
	TiXmlDocument  xml_f(xml_file);
	 xml_f.LoadFile();
	 if(xml_f.Error())
	 {
		 cout<< "couldn't read file";
			 return NULL;
	 }
	 int number= 0;
	 TiXmlElement* root = xml_f.FirstChildElement ("Initial_parameters"); 
    TiXmlElement* sub_root =root->FirstChildElement ("Inject_Particles"); 
	double* params= 0;
	 char* p_name= new char [50];
	 Bunch*  prtls =0;

       
			 const char* test = sub_root->GetText();
			 strcpy (p_name,sub_root->GetText());
			params = read_double_params(p_name);
			 prtls = new Bunch(p_name,params[0],params[1],params[2],c_geom, p_list,params[3],params[4]);
			 prtls->calc_init_param(params[5],params[6]);   
	 
	   return  prtls;
}
////////////////////////////////////////////
void Load_init_param:: load_system()
{


	int i=0;
	//load PML parameters///
	///////////////////////////////////
	double * r_params   =  read_double_params("PML");
	c_pml = new PML(r_params);

	//load Geometry parameters///
	/////////////////////////////////////
    r_params   =  read_double_params("geometry");
    c_geom  = new Geometry (r_params ,c_pml);
	c_geom->set_epsilon() ;
	//delete r_params;
	//////////////////////////////////////////////////
	///////////////////////
	///creating field objects

	efield = new E_field(c_geom);
	hfield = new H_field(c_geom);
	efield->boundary_conditions();
	efield->set_homogeneous_efield(0.0, 0.0, 0);
	hfield->set_homogeneous_h(0.0, 0.0, 0.0);
	efield->set_fi_on_z();

	//load time parameters///
	/////////////////////////////////////
	 r_params   =  read_double_params("Time");
	
	 c_time = new Time(r_params);
	 	
	 //load particle parameters///
	/////////////////////////////////////
	 p_list = new particles_list(0); 
	 read_load_particles();
   //////////////////////////////////////

	  //load bunch///
	/////////////////////////////////////
c_bunch = read_load_bunch();
	 
   //////////////////////////////////////

	 //Maxwell initial conditions///

 r_params   =  read_double_params("Boundary_Maxwell_conditions");
  Boundary_Maxwell_conditions maxwell_rad(efield);
maxwell_rad.specify_initial_field(c_geom,r_params[0],r_params[1],r_params[2]);

/////////////////////////////////////
// creating rho and current arrays
	c_rho_new = new charge_density(c_geom);
	c_rho_old = new charge_density (c_geom);
	c_rho_beam = new charge_density (c_geom);
	c_current = new current(c_geom);

////////////////////////////////////
// boundary conditions///
char* a = read_char("Boundary_conditions");
if (get_int_value(a)==0)
{
	p_list->charge_weighting(c_rho_new);
 Fourier four1(0);
 Poisson_dirichlet dirih(c_geom);
	dirih.poisson_solve(efield, c_rho_new);
}


	//load File Path///
	/////////////////////////////////////
    char* path_res  =  read_char("PathtoResult");
	char* path_dump  =  read_char("PathtoSaveState");
    c_io_class  = new input_output_class (path_res , path_dump);
	
	//////////////////////////////////////////////////

//	KernAccessor *kern_access = new KernAccessor(c_geom->n_grid_1, c_geom->n_grid_2);
//	kern_access_global = kern_access;
//
}
bool Load_init_param:: SaveSystemState() 
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
void Load_init_param:: Run(void)
{
	this->c_time->current_time = 0.0 ;
	p_list->charge_weighting(c_rho_new);

	Poisson_dirichlet dirih(c_geom);
	dirih.poisson_solve(efield, c_rho_new);

	//variable for out_class function 
	p_list->create_coord_arrays();
	int step_number= 0;
     
    while (c_time->current_time < c_time->end_time)
	{
		clock_t t1 = clock();
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
			p_list->charge_weighting(c_rho_old);  //continuity equation
			p_list->half_step_coord(c_time);
			p_list->azimuthal_j_weighting(c_time, c_current);
			p_list->half_step_coord(c_time);
			p_list->j_weighting(c_time,c_current);
		

        //4. Calculate E
	  // maxwell_rad.probe_mode_exitation(&geom1,&current1, 1,7e8, time1.current_time);
       efield->calc_field(hfield,c_time, c_current, c_pml);
		
        //continuity equation
		c_rho_new->reset_rho();
		
			p_list->charge_weighting(c_rho_new);  //continuity equation
		//bool res =  continuity_equation(c_time, c_geom, c_current, c_rho_old, c_rho_new); 
		
		cout << "Execution time: " << clock() - t1 << endl;
		
		if  ((((int)(c_time->current_time/c_time->delta_t))%5==0))
		//if  ((((int)(c_time->current_time/c_time->delta_t)) < 10))
		//if  ( abs(time1.current_time - time1.end_time + time1.delta_t) < 1e-13)
		{
			cout<<c_time->current_time<<" ";
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
				if  ((((int)(c_time->current_time/c_time->delta_t))%1000==0)&&(step_number!=1))
				this->SaveSystemState();
		}
	
		c_time->current_time = c_time->current_time + c_time->delta_t;
		//if (!res)
		//	cout<<"Error:"<<c_time->current_time<<"! ";

	}
}