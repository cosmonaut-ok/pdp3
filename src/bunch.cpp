#include "bunch.h"

using namespace constant;

Bunch::Bunch (char *p_name,
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
  double b_lenght = duration * velocity;
  double v_total = 0.; // total volume of all macroparticles
  double N_total = PI * pow(radius, 2) * b_lenght * density; // total number of regular particles
  double dr = geom1->dr;
  double dz = geom1->dz;

#pragma omp parallel for
  for(unsigned int i = 0; i < number; i++)
  {
    double rand_r = lib::random_reverse(i, 9); // TODO: why 9

    pos[i][0] = (radius - dr) * rand_r + dr / 2;
    pos[i][1] = 0.;
    pos[i][2] = 0.;

    is_alive[i] = false;
    vel[i][0] = 0.;
    vel[i][1] = 0.;
    vel[i][2] = 0.;

    v_total += 2 * PI * pos[i][0] * dr * dz;
  }

  // average volume of single macroparticle
  double v_avg = v_total / number;
  double n_per_macro_avg = N_total / number;

  for (unsigned int i = 0; i < number; i++)
  {
    // coefitient of normalization
    double norm =  2 * PI * pos[i][0] * dr * dz / v_avg;

    // number of real particles per macroparticle
    double n_per_macro = n_per_macro_avg * norm;

    // set charge and mass of macroparticle
    charge_array[i] = charge * n_per_macro;
    mass_array[i] = mass * n_per_macro;
  }

  double charge_total = density * PI * radius * radius * b_lenght * charge;
  double charge_total_macro = 0.;
  int onax = 0;
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
  double half_z_cell_size = geom1->dz / 2.0;

  if (time->current_time >= bunch_time_begin && time->current_time < bunch_time_end) // (time->current_time<duration)
#pragma omp parallel shared(start_number, dl, half_z_cell_size)
  {
#pragma omp for
    for(unsigned int i = 0; i < particles_in_step; i++)
    {
      double rand_z = lib::random_reverse(start_number + i, 11); // TODO: why 11?

      if (i+start_number < number)
      {
        pos[i+start_number][2] = dl * rand_z + half_z_cell_size;
        vel[i+start_number][2] = velocity;
        is_alive[i+start_number] = true;
      }
#ifndef _OPENMP
      else // optimize for singlethread case
        break;
#endif
    }
  }
}

void Bunch::reflection_single(unsigned int i)
{
  double dr = geom1->dr;
  double dz = geom1->dz;
  double radius_wall = geom1->first_size - dr / 2.;
  double longitude_wall = geom1->second_size - dz / 2.;
  double half_dr = dr / 2.;
  double half_dz = dz / 2.;

  //! FIXME: fix wall reflections for r-position
  if (pos[i][0] > radius_wall)
    is_alive[i] = false;

  if (pos[i][2] > longitude_wall)
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
