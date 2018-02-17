#include "bunch.h"
#include "pdp3Time.h"
#include "poissonDirichlet.h"
#include "constant.h"

using namespace constant;

Bunch::Bunch(char* p_name,
             double p_charge,
             double p_mass,
             int p_number,
             Geometry* geom,
             // particles_list* p_list,
             double b_duration,
             double b_radius,
             double b_density,
             double b_init_velocity):Particles(p_name, p_charge, p_mass,
                                               p_number, geom)
{
  // fill object fields
  duration = b_duration;
  radius = b_radius;
  n_bunch = b_density;
  vel_bunch = b_init_velocity;
  //
  double b_lenght = duration*vel_bunch;
  double n_in_big = PI * pow(radius, 2) * b_lenght * n_bunch / number;
  charge *= n_in_big;
  mass *= n_in_big;

#pragma omp parallel for
  for(int i=0; i<number; i++)
  {
    is_alive[i] = 0;
    v1[i]=0;
    v2[i]=0;
    v3[i]=0;
    mass_array[i] = mass;
    charge_array[i] = charge;
  }
}

Bunch::~Bunch(void)
{
}

void Bunch::bunch_inject(Time* time)
{
  double dl = vel_bunch*time->delta_t;
  int step_num = duration/time->delta_t;
  int particles_in_step = number/step_num;
  int start_number = time->current_time/time->delta_t*particles_in_step;
  // exit (1);
  double dr = geom1->dr*1.00000001; // TODO: WTF?
  double dz = geom1->dz*1.00000001;
  double rand_i;
  double rand_z;

  if (time->current_time<duration)
#pragma omp parallel shared(start_number, dr, dz, dl, time) private(rand_i, rand_z)
  {
#pragma omp for
    for(int i = 0; i <  particles_in_step; i++)
    {
      rand_i = random_reverse(start_number + i, 9);
      rand_z = random_reverse(start_number + i, 11);

      x1[i+start_number] = sqrt(dr * dr / 4.0 + radius * (radius - dr) * rand_i);

      x3[i+start_number] = dl*(rand_z)+dz/2.0;
      v3[i+start_number] = vel_bunch;
      v1[i+start_number] = 0;
      v2[i+start_number] = 0; // fi velocity;
      is_alive[i+start_number] = true;
    }

#pragma omp for
    for(int i = 0; i <  number; i++)
      if(x3[i]>(geom1->second_size - dz/2.0))
      {
        is_alive[i]=false;
      }
  }
}

void Bunch::half_step_coord(Time* t)
{
  double dr = geom1->dr;
  double dz = geom1->dz;
  double x1_wall = geom1->first_size - dr/2.0;
  double x3_wall = geom1->second_size - dz/2.0;
  double half_dr = dr/2.0;
  double half_dz = dz/2.0;
  double x1_wallX2 = x1_wall*2.0;
  // double x3_wallX2 = x3_wall*2.0;
  double half_dt = t->delta_t/2.0;

#pragma omp parallel for shared(dr, dz, x1_wall, x3_wall, half_dr, half_dz, half_dt, x1_wallX2)
  for(int i=0;i<number;i++)
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
        is_alive[i] = false;
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
