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

public:
  void back_velocity_to_rz();

  void charge_weighting(ChargeDensity *rho);
  void step_v(EField *e_fld, HField *h_fld, Time *t);
  void full_j_weighting(Current *current, Time *t);
  void move_half_reflect(Time *t);
  void azimuthal_j_weighting(Time *time1, Current *j1);
  void j_weighting(Time *time1, Current *j1);
};
