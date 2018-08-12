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
  void charge_weighting(ChargeDensity *ro1);
  void step_v(EField *e_fld, HField *h_fld, Time *t);
  void half_step_pos(Time *t);
  void reflection();
  void j_weighting(Time *time1, Current *j1);
  void azimuthal_j_weighting(Time *time1, Current *j1);
  void dump_position_to_old(void);

  void back_position_to_rz();
  void back_velocity_to_rz();
};
