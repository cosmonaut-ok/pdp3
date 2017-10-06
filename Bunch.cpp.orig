#include "Bunch.h"
#include "Time.h"
#include "Poisson_dirichlet.h"

Bunch::Bunch(char* p_name, double p_charge, double p_mass, int p_number, Geometry* geom, particles_list* p_list, double b_duration, double b_radius):Particles(p_name,p_charge,p_mass,p_number,geom,p_list)
{
	duration = b_duration;
	radius = b_radius;
};
Bunch::~Bunch(void)
{
}
void Bunch::calc_init_param(double n_b, double b_vel)
{
	n_bunch = n_b;
	vel_bunch = b_vel;
	double b_lenght = duration*vel_bunch;
	double n_in_big = 3.1415*radius*radius*b_lenght*n_bunch/number;
	charge *= n_in_big;
	mass *= n_in_big;
	for(int i=0;i<number;i++)
	{
		is_alive[i]=false;
		v1[i]=0;
		v2[i]=0;
		v3[i]=0;
	}
}



void Bunch::bunch_inject( Time* time)
{
	double dl = vel_bunch*time->delta_t;
	int step_num =  duration/time->delta_t;
	int particles_in_step =  number/step_num;
	int start_number = time->current_time/time->delta_t*particles_in_step;
	double dr = geom1->dr*1.00000001;
	double dz = geom1->dz*1.00000001;
 if (time->current_time<duration)
		for(int i = 0; i <  particles_in_step; i++)
		{
		double	rand_r = random_reverse(i,3);		
		double	rand_z = random_reverse(i,5);
			x1[i+start_number] = (radius)*sqrt(rand_r) + dr/2.0;
		
			x3[i+start_number] = dl*(rand_z)+dz/2.0;
			v3[i+start_number] = vel_bunch;
			v1[i+start_number] = 1e5;
			is_alive[i+start_number] = true;
		}
}
void Bunch::bunch_inject_calc_E(Geometry* geom,E_field*E_beam, E_field*E, Time* time)
{
	double dl = vel_bunch*time->delta_t;
	int step_num =  duration/time->delta_t;
	int particles_in_step =  number/step_num;
	int start_number = time->current_time/time->delta_t*particles_in_step;
	double dr = geom1->dr*1.00000001;
	double dz = geom1->dz*1.00000001;
 if (time->current_time<duration)
		for(int i = 0; i <  particles_in_step; i++)
		{
		double	rand_r = random_reverse(i,3);		
		double	rand_z = random_reverse(i,5);
			x1[i+start_number] = (radius)*sqrt(rand_r) + dr/2.0;
		
			x3[i+start_number] = dl*(rand_z)-dl;
			v3[i+start_number] =vel_bunch;
			v1[i+start_number] = 1e5;
			is_alive[i+start_number] = true;
		}
if (time->current_time==0)
	{
		// calculational field of elementary beam portion
		charge_density rho_beam(geom);
		charge_weighting(&rho_beam);
		Poisson_dirichlet dirih(geom);
		dirih.poisson_solve(E_beam, &rho_beam);
	}
		////Er////
	for(int i=0;i<(geom1->n_grid_1-1);i++)
		for(int k=0;k<(geom1->n_grid_2-1);k++)
		{
			E->e1[i][k]= E->e1[i][k]+E_beam->e1[i][k];
		}
	///Ef////
	for(int i=0;i<(geom1->n_grid_1-1);i++)
		for(int k=0;k<(geom1->n_grid_2-1);k++)
		{
			E->e2[i][k]= E->e2[i][k]+E_beam->e2[i][k];
		}
	///Ez////
	for(int i=0;i<(geom1->n_grid_1-1);i++)
		for(int k=0;k<(geom1->n_grid_2-1);k++)
		{
			E->e3[i][k]= E->e3[i][k]+E_beam->e3[i][k];
		}
	
}

