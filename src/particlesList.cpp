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

void ParticlesList::move_half_reflect(Time *t)
{
  for(std::size_t i = 0; i < part_list.size(); i++)
    part_list[i]->move_half_reflect(t);
}

void ParticlesList::current_distribution(Time *time1, Current *j_r)
{
  for(std::size_t i = 0; i < part_list.size(); i++)
    part_list[i]->current_distribution(time1,j_r);
}

void ParticlesList::azimuthal_current_distribution(Time *time1, Current *j_r)
{
  for(std::size_t i=0; i<part_list.size(); i++)
    part_list[i]->azimuthal_current_distribution(j_r);
}

void ParticlesList::back_velocity_to_rz()
{
  for(size_t i = 0; i < part_list.size(); i++)
    part_list[i]->back_velocity_to_rz();
}

void ParticlesList::full_current_distribution(Current *current, Time *t)
{
  for(std::size_t i = 0; i < part_list.size(); i++)
    part_list[i]->full_current_distribution(current, t);
}
