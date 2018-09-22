#include "particlesList.h"

using namespace std;

ParticlesList::ParticlesList(void)
{
}

ParticlesList::~ParticlesList(void)
{
}

// sicle for all particles kinds in system
void ParticlesList::charge_weighting(ChargeDensity *rho)
{
  for(size_t i = 0; i < part_list.size(); i++)
    part_list[i]->charge_weighting(rho);
}

void ParticlesList::step_v(EField *e_fld, HField *h_fld, Time *t)
{
  for(std::size_t i = 0; i < part_list.size(); i++)
    part_list[i]->step_v(e_fld, h_fld, t);
}

void ParticlesList::half_step_pos(Time *t)
{
  for(size_t i = 0; i < part_list.size(); i++)
    part_list[i]->half_step_pos(t);
}

void ParticlesList::reflection(void)
{
  for(size_t i = 0; i < part_list.size(); i++)
    part_list[i]->reflection();
}

void ParticlesList::j_weighting(Time *time1, Current *j1)
{
  for(std::size_t i = 0; i < part_list.size(); i++)
    part_list[i]->j_weighting(time1,j1);
}

void ParticlesList::azimuthal_j_weighting(Time *time1, Current *j1)
{
  for(std::size_t i=0; i<part_list.size(); i++)
    part_list[i]->azimuthal_j_weighting(j1);
}

void ParticlesList::dump_position_to_old()
//! fuction for copying particles position
{
  int kinds_number = part_list.size();
  for (int k=0; k < kinds_number; k++)
    part_list[k]->dump_position_to_old();
}

void ParticlesList::back_position_to_rz()
{
  for(size_t i = 0; i < part_list.size(); i++)
    part_list[i]->back_position_to_rz();
}

void ParticlesList::back_velocity_to_rz()
{
  for(size_t i = 0; i < part_list.size(); i++)
    part_list[i]->back_velocity_to_rz();
}
