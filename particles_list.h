#pragma once
#include "Particles.h"
#include <vector>
#include "chargeDensity.h"
#include "eField.h"
#include "hField.h"

using namespace std;

class Particles;

class particles_list
{
public:
  particles_list(void);
  ~particles_list(void);
public:
  vector<Particles*> part_list;
  double** x1_old;
  double** x3_old;

public:
  void charge_weighting(ChargeDensity* ro1);
  void step_v(EField *e_fld, HField *h_fld, Time* t);
  void half_step_coord(Time* t);
  void j_weighting(Time* time1, current *j1);
  void azimuthal_j_weighting(Time* time1, current *j1);
  void create_coord_arrays(void);
  void copy_coords(void);
};
