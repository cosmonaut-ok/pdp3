#include "Beam.h"
#include "Time.h"
#include "Poisson_dirichlet.h"

using namespace std;

Beam::Beam(char* p_name, double p_charge, double p_mass, int p_number, Geometry* geom, particles_list* p_list, double b_radius):Particles(p_name,p_charge,p_mass,p_number,geom,p_list)
{
	radius = b_radius;
}
void Beam::calc_init_param(Time* time,int particles_in_step,double n_b, double b_vel)
{
	n_beam = n_b;
	vel_beam = b_vel;
	double dl = vel_beam*time->delta_t;
	double pi = 3.1415926;
	double n_in_big = 3.1415*radius*radius*dl*n_beam/particles_in_step;
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
void Beam::beam_inject(Time* time,int particles_in_step, double fi, double koef)
{
	flcuda dl = vel_beam*time->delta_t;
	flcuda dr = geom1->dr*1.00000001;
	flcuda dz = geom1->dz*1.00000001;
	double pi = 3.1415926;
	double mod_n = (1 - koef/2.0 + koef/2.0*sin(2*pi*fi*time->current_time));
	int n=0;
	int i=0;
	particles_in_step = particles_in_step*mod_n;
 		while((n<particles_in_step)&&(i<number))
		{
			if(!is_alive[i])
			{
		flcuda	rand_r = random_reverse(i,3);		
		flcuda	rand_z = random_reverse(i,5);
				x1[i] = (radius)*sqrt(rand_r) + dr/2.0;
				x3[i] = dl*(rand_z)-dl;
				v3[i] = vel_beam;
				is_alive[i] = true;
				n=n+1;
			}
		i=i+1;
		}
		for(int i = 0; i <  number; i++)
		if(x3[i]>(geom1->second_size - dz/2.0))
		{
			is_alive[i]=false;
		}

}

void Beam::half_step_coord(Time* t)
{
	int i;
	double dr = geom1->dr;
	double dz = geom1->dz;
	double x1_wall = geom1->first_size - dr/2.0;
	double x3_wall = geom1->second_size - dz/2.0;
	double half_dr = dr/2.0;
	double half_dz = dz/2.0;
	double x1_wallX2 = x1_wall*2.0;
	double x3_wallX2 = x3_wall*2.0;
	double half_dt = t->delta_t/2.0;

	for( i=0;i<number;i++)
		if (is_alive[i])
		{	
			x1[i] = x1[i] + v1[i]*half_dt; 
            x3[i] = x3[i] + v3[i]*half_dt;

			if (x1[i] > x1_wall)
			{
				x1[i] = x1_wallX2 - x1[i];
				v1[i] = -v1[i];
			}

			if (x3[i] > x3_wall)
			{
				x3[i] = x3_wallX2 - x3[i];
				v3[i] = -v3[i];
			}

			if (x1[i] < half_dr)
			{
				x1[i] = dr - x1[i];
				v1[i] = -v1[i];
			}

			if (x3[i] < half_dz)
			{
				//x3[i] = dz - x3[i];
				//v3[i] = -v3[i];
			}
		}
}

//////////////////////////////////////////////////

Beam::~Beam(void)
{
}