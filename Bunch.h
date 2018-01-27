#pragma once
#include "Particles.h"
#include "pdp3_time.h"

class Bunch : public Particles
{
public:
  Bunch(char* p_name,
        double p_charge,
        double p_mass,
        int p_number,
        Geometry* geom,
        particles_list* p_list,
        double b_duration,
        double b_radius,
        double b_density,
        double b_init_velocity);

  ~Bunch(void);
public:
  double duration; //bunch duration
  double n_bunch; // bunch density;
  double vel_bunch; // bunch velocity
  double radius; // bunch radius
public:
  void bunch_inject(Time* time);
  void bunch_inject_calc_E(Geometry* geom,
                           E_field * E_beam,
                           E_field* E,
                           Time* time);
  virtual void half_step_coord(Time* t);
};
