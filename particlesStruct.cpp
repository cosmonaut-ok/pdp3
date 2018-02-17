#include "particlesStruct.h"

ParticlesStruct CreateParticlesStruct(double charge,
                                      double mass,
                                      int number,
                                      int grid_num1,
                                      int grid_num3,
                                      double dr,
                                      double dz)
{
  ParticlesStruct res;
  res.charge = charge;
  res.mass = mass;
  res.number = number;
  res.grid_num1 = grid_num1;
  res.grid_num3 = grid_num3;
  res.dr = dr;
  res.dz = dz;
  return res;
}
