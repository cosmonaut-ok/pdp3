#include <vector>

#include "particles.h"

using namespace std;

class Particles;

class ParticlesList
{
public:
  ParticlesList(void);
  ~ParticlesList(void);
public:
  vector<Particles*> part_list;
  double **x1_old;
  double **x3_old;

public:
  void charge_weighting(ChargeDensity *ro1);
  void step_v(EField *e_fld, HField *h_fld, Time *t);
  void half_step_coord(Time *t);
  void j_weighting(Time *time1, Current *j1);
  void azimuthal_j_weighting(Time *time1, Current *j1);
  void create_coord_arrays(void);
  void copy_coords(void);

  void back_coordinates_to_rz();
  void back_velocity_to_rz();

};
