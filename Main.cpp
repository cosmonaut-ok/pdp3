#include<iostream>
#include <Constant.h>
#include "field.h"
#include "E_field.h"
#include "H_field.h"
#include "pdp3_time.h"
#include "Particles.h"
#include "Fourier.h"
#include "Poisson.h"
#include "poisson_neumann.h"
#include "poisson_dirichlet.h"
#include "particles_list.h"
#include <fstream>
#include <math.h>
// #include <tchar.h>

#include "Boundary_Maxwell_conditions.h"
#include "input_output_class.h"
#include "Beam.h"
#include "Bunch.h"
#include "time.h"
#include "particles_struct.h"
#include "system_host.cuh"
#include "Load_init_param.h"

using namespace std;

Particles_struct specie;
#define BUILD_OPENCL

int main(int argc, char **argv)
{
	clock_t start, finish;
	flcuda time_elapsed;
	Load_init_param init_param((char *) "parameters.xml");
	init_param.read_xml();
	init_param.load_system();

	printf("Loaded parameters. Launching Run...\n");
	init_param.Run();

	PML pml1(0.0,0.0, 0.0, 0.000001, 0.07);
	Geometry geom1(0.25,2.0, 255, 2047, &pml1);
	//Geometry geom1(0.2,1.5, 63, 255, &pml1);
	flcuda left_plasma_boundary = geom1.second_size*0.0;

	Time time1(0,0,0,200000e-12,1e-12); // TODO: WTF?
	E_field e_field1(&geom1);
	H_field h_field1(&geom1);
	Fourier four1(0);
	input_output_class out_class((char *) "pdp3_files/_results_ocl/", (char *) "pdp3_result/dump");

	Boundary_Maxwell_conditions maxwell_rad(&e_field1);
	maxwell_rad.specify_initial_field(&geom1,0,0,0);
	bool res = true;
	int k, i;

	current current1(&geom1);
	charge_density rho_new(&geom1);
	charge_density rho_old(&geom1);
	charge_density rho_beam(&geom1);
	e_field1.boundary_conditions();
	e_field1.set_homogeneous_efield(0.0, 0.0, 0);
	h_field1.set_homogeneous_h(0.0, 0.0, 0.0);
	e_field1.set_fi_on_z();

	E_field e_beam(&geom1);
	e_beam.set_homogeneous_efield(0.0, 0.0, 0);

	/////////////////////////////////////////////

	//	////////////////////////////////////////
	//	poisson	 equetion testing//
	//for(k=0;k<geom1.n_grid_2;k++)
	//	rho_new.set_ro_weighting(20,k,1e-7);

	//e_field1.poisson_equation2(&geom1,&rho_new);
	////e_field1.poisson_equation(&geom1,&rho_new);
	//bool res2 = e_field1.test_poisson_equation(&rho_new);
	//////////////////////////////////////////////////////

	////////////////////////////////////////////////
	geom1.set_epsilon() ;
	//	e_field1.set_sigma();
	particles_list p_list(0);
	///////////////////////////////////////////
	// beam part
	//Beam electron_beam("electron_beam", -1, 1, 10e5, &geom1,&p_list,0.01);
	//electron_beam.calc_init_param(&time1,50,5e12,3e7);
	Bunch electron_bunch((char *) "electron_bunch", -1,1,1e6,&geom1,&p_list,2e-8,0.02);
	electron_bunch.calc_init_param(10e12,3.0e7);
	///////////////////////////////////////////
	Particles electrons((char *) "electrons", -1, 1,0e6, &geom1,&p_list);
	Particles ions("ions", 1, 1836, 0e6, &geom1,&p_list);
	p_list.create_coord_arrays();

	//electrons.load_spatial_distribution(0.8e14, 0.81e14, left_plasma_boundary,0);
	electrons.load_spatial_distribution(2e14, 8e14, left_plasma_boundary,0);
	ions.load_spatial_distribution(2e14, 8e14, left_plasma_boundary,0);

	electrons.velocity_distribution_v2(1.0);
	ions.velocity_distribution_v2(0.05);
	//ofstream out_vel("velocities");
	//ofstream out_coords("coords");
	//for (i = 0; i< electrons.number; i++)
	//{
	//	out_vel<<electrons.v1[i]<<" "<<electrons.v2[i]<<" "<<electrons.v3[i]<<" ";
	//	out_coords<<electrons.x1[i]<<" "<<electrons.x3[i]<<" ";
	//}
	//out_vel.close();
	//out_coords.close();

	int cuda_particles_number = 0;
	for (i = 0; i < p_list.part_list.size(); i++)
		if (cuda_particles_number < p_list.part_list[i]->number)
			cuda_particles_number = p_list.part_list[i]->number;

	/*charge_density rho_elect(&geom1);
		electrons.charge_weighting(&rho_elect);
		out_class.out_data("rho",rho_elect.get_ro(),0,100,geom1.n_grid_1-1,geom1.n_grid_2-1);*/

#ifdef BUILD_CUDA
	InitCUDA();
	SetupCUDA(geom1.n_grid_1, geom1.n_grid_2, cuda_particles_number);
#endif

	/////////////////////////////////
	//0. Half step back

	//p_list.azimuthal_j_weighting(&time1, &current1);
	//p_list.j_weighting(&time1,&current1,);

	p_list.charge_weighting(&rho_new);
	//electrons.charge_weighting(&rho_new);
	//out_class.out_data("rho",rho_new.get_ro(),1,128,2048);

	//weight currents and charges before relaxation period
	//solve Poisson equation


	//Poisson_neumann poisson1(&geom1);

	//poisson1.poisson_solve(&e_field1, &rho_new);
	Poisson_dirichlet dirih(&geom1);
	dirih.poisson_solve(&e_field1, &rho_new);
	bool res1 = dirih.test_poisson_equation(&e_field1, &rho_new);
	//out_class.out_data("e3",e_field1.e3,0,100,geom1.n_grid_1-1,geom1.n_grid_2-1);

	//relaxation period
	while (time1.current_time < time1.relaxation_time)
	{
		//1. Calculate H field
		h_field1.calc_field(&e_field1, &time1);

		//2. Calculate E
		e_field1.calc_field(&h_field1, &time1, &current1);
		time1.current_time = time1.current_time + time1.delta_t;


	}
	time1.current_time = 0.0 ;

	//variable for out_class function
	int step_number= 0;

	while (time1.current_time < time1.end_time)
	{
		//electron_beam.beam_inject(1e14,5e7,&time1);
		//electron_beam.beam_inject(&time1,50,1.6e8,0.5);
		electron_bunch.bunch_inject(&time1);
		//electron_bunch.bunch_inject_calc_E(&geom1, &e_beam, &e_field1, &time1);
		//radiation	 source
		//maxwell_rad.radiation_source(&geom1,0.4,2e9,0,time1.current_time);

		//1. Calculate H field
		h_field1.calc_field(&e_field1, &time1);

		//2. Calculate v
		current1.reset_j();
		rho_old.reset_rho();
		rho_beam.reset_rho();
		p_list.step_v(&e_field1, &h_field1, &time1);

		//3. Calculate x, calculate J
		p_list.copy_coords();
		p_list.charge_weighting(&rho_old);	//continuity equation
		p_list.half_step_coord(&time1);
		p_list.azimuthal_j_weighting(&time1, &current1);
		p_list.half_step_coord(&time1);
		p_list.j_weighting(&time1,&current1);


		//4. Calculate E
		// maxwell_rad.probe_mode_exitation(&geom1,&current1, 1,7e8, time1.current_time);
		e_field1.calc_field(&h_field1, &time1, &current1, &pml1);

		//continuity equation
		rho_new.reset_rho();

		p_list.charge_weighting(&rho_new);	//continuity equation
		res =	 continuity_equation(&time1, &geom1, &current1, &rho_old, &rho_new);



		if	((((int)(time1.current_time/time1.delta_t))%100==0))
			//if	( abs(time1.current_time - time1.end_time + time1.delta_t) < 1e-13)
		{
			cout<<time1.current_time<<" ";
			electron_bunch.charge_weighting(&rho_beam);
			//rho_old.reset_rho();
			//electrons.charge_weighting(&rho_old);
			//out_class.out_data("rho_el", rho_old.get_ro(),step_number,100,geom1.n_grid_1-1,geom1.n_grid_2-1);
			//out_class.out_data("e1",e_field1.e1,100,128,2048);
			out_class.out_data("rho_beam", rho_beam.get_ro(),step_number,100,geom1.n_grid_1-1,geom1.n_grid_2-1);
			out_class.out_data("e3",e_field1.e3,step_number,100,geom1.n_grid_1-1,geom1.n_grid_2-1);
			out_class.out_data("e1",e_field1.e1,step_number,100,geom1.n_grid_1-1,geom1.n_grid_2-1);
			//out_class.out_data("rho",rho_elect.get_ro(),step_number,100,geom1.n_grid_1-1,geom1.n_grid_2-1);
			//out_class.out_coord("vels",electron_bunch.v1, electron_bunch.v3, step_number, 100, electron_bunch.number);
			//out_class.out_coord("coords",electrons.x1, electrons.x3, step_number, 100, electrons.number);
			//out_class.out_coord("vels",electrons.v1, electrons.v3, step_number, 100, electrons.number);
			out_class.out_data("h2",h_field1.h2,step_number,100,geom1.n_grid_1-1,geom1.n_grid_2-1);

			step_number=step_number+1;
		}

		time1.current_time = time1.current_time + time1.delta_t;
		if (!res)
			exit(0); // TODO: is it right?
		//break;
		// cout<<"Error:"<<time1.current_time<<"! ";

	}


};
