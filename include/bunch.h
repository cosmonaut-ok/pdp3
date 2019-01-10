#include "particles.h"

class Bunch : public Particles
{
public:
  Bunch(char *p_name,
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
        double h_duration);

  virtual ~Bunch(void);
public:
  double duration; // bunch duration
  double bunch_number; // bunch number (for series)
  double hole_duration; // duration between two bunches in series
  double density;  // bunch density (n_bunch);
  double velocity; // bunch velocity
  double radius;   // bunch radius
public:
  void bunch_inject(Time *time);
  void bunch_inject_calc_E(Geometry *geom,
                           EField  *E_bunch,
                           EField *E,
                           Time *time);
  void reflection_single(unsigned int i);
  void move_half_reflect(Time *t);
};
