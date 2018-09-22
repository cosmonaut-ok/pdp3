#include "bunch.h"

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
             double b_init_velocity,
             int b_number,
             double h_duration):Particles(p_name, p_charge, p_mass,
                                          p_number, geom)
{
  // fill object fields
  duration = b_duration;
  bunch_number = b_number; // bunch number
  hole_duration = h_duration; // bunch number
  radius = b_radius;
  density = b_density;
  velocity = b_init_velocity;
  //
  double b_lenght = duration*velocity;
  double n_in_macro = PI * pow(radius, 2) * b_lenght * density / number;
  charge *= n_in_macro;
  mass *= n_in_macro;

#pragma omp parallel for
  for(unsigned int i=0; i<number; i++)
  {
    is_alive[i] = false;
    vel[i][0]=0;
    vel[i][1]=0;
    vel[i][2]=0;
    mass_array[i] = mass;
    charge_array[i] = charge;
  }
}

// Destructor
Bunch::~Bunch()
{
}

void Bunch::bunch_inject(Time *time)
{
  double bunch_time_begin = duration * bunch_number + hole_duration * bunch_number;
  double bunch_time_end = bunch_time_begin + duration;

  double dl = velocity * time->delta_t;
  int steps_amount = ceil(duration / time->delta_t);
  unsigned int particles_in_step = number / steps_amount;
  int start_number = (time->current_time - bunch_time_begin) / time->delta_t * particles_in_step;

  // very local constants
  double half_r_cell_size_pow_2 = pow((geom1->dr / 2), 2);
  double half_z_cell_size = geom1->dz / 2.0;
  double const1 = radius * (radius - geom1->dr); // TODO: what is it? and why?

  if (time->current_time >= bunch_time_begin && time->current_time < bunch_time_end) // (time->current_time<duration)
#pragma omp parallel shared(start_number, dl, half_r_cell_size_pow_2, half_z_cell_size, const1)
  {
#pragma omp for
    for(unsigned int i = 0; i < particles_in_step; i++)
    {
      double rand_i = lib::random_reverse(start_number + i, 9); // TODO: why 9 and 11?
      double rand_z = lib::random_reverse(start_number + i, 11);

      if (i+start_number < number)
      {
        pos[i+start_number][0] = sqrt(half_r_cell_size_pow_2 + const1 * rand_i);
	pos[i+start_number][1] = 0;
        pos[i+start_number][2] = dl * rand_z + half_z_cell_size;
        vel[i+start_number][2] = velocity;
        vel[i+start_number][0] = 0;
        vel[i+start_number][1] = 0; // fi velocity;
        is_alive[i+start_number] = true;
      }
#ifndef _OPENMP
      else // optimize for singlethread case
        break;
#endif
    }
  }
}

void Bunch::reflection()
{
#pragma omp parallel for
  for(unsigned int i=0; i<number; i++)
    if (is_alive[i])
    {
      double dr = geom1->dr;
      double dz = geom1->dz;
      double x1_wall = geom1->first_size - dr/2.0;
      double x3_wall = geom1->second_size - dz/2.0;
      double half_dr = dr/2.0;
      double half_dz = dz/2.0;

      //! FIXME: fix wall reflections for r-position
      if (pos[i][0] > x1_wall)
        is_alive[i] = false;

      if (pos[i][2] > x3_wall)
        is_alive[i] = false;

      if (pos[i][0] < half_dr)
      {
        pos[i][0] = dr - pos[i][0];
        vel[i][0] = -vel[i][0];
      }

      if (pos[i][2] < half_dz)
      {
        pos[i][2] = dz - pos[i][2];
        vel[i][2] = -vel[i][2];
      }
    }
}
