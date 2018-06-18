#include "bunch.h"
#include "pdp3Time.h"
#include "poissonDirichlet.h"
#include "constant.h"

using namespace constant;

Bunch::Bunch(char *p_name,
             double p_charge,
             double p_mass,
             int p_number,
             Geometry *geom,
             // particles_list *p_list,
             double b_duration,
             double b_radius,
             double b_density,
             double b_init_velocity):Particles(p_name, p_charge, p_mass,
                                               p_number, geom)
{
  // fill object fields
  duration = b_duration;
  radius = b_radius;
  density = b_density;
  velocity = b_init_velocity;
  //
  double b_lenght = duration*velocity;
  double n_in_big = PI * pow(radius, 2) * b_lenght * density / number;
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

void Bunch::bunch_inject(Time *time)
{
  double dl = velocity * time->delta_t;
  int step_num = duration / time->delta_t;
  int particles_in_step = number / step_num;
  int start_number = time->current_time / time->delta_t * particles_in_step;

  // very local constants
  double half_r_cell_size_pow_2 = pow((geom1->dr / 2), 2);
  double half_z_cell_size = geom1->dz / 2.0;
  double const1 = radius * (radius - geom1->dr); // TODO: what is it? and why?

  if (time->current_time<duration)
#pragma omp parallel shared(start_number, dl, half_r_cell_size_pow_2, half_z_cell_size, const1)
  {
#pragma omp for
    for(int i = 0; i < particles_in_step; i++)
    {
      double rand_i = random_reverse(start_number + i, 9); // TODO: why 9 and 11?
      double rand_z = random_reverse(start_number + i, 11);

      x1[i+start_number] = sqrt(half_r_cell_size_pow_2 + const1 * rand_i);
      x3[i+start_number] = dl * rand_z + half_z_cell_size;
      v3[i+start_number] = velocity;
      v1[i+start_number] = 0;
      v2[i+start_number] = 0; // fi velocity;
      is_alive[i+start_number] = true;
    }

#pragma omp for
    for(int i = 0; i < number; i++)
      if(x3[i]>(geom1->second_size - half_z_cell_size))
      {
        is_alive[i]=false;
      }
  }
}
