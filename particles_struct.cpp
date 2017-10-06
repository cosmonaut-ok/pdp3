#include "particles_struct.h"

Particles_struct CreateParticles_struct(flcuda charge, flcuda mass, int number, int grid_num1, int grid_num3, flcuda dr ,flcuda dz)
{
  Particles_struct res;
  res.charge = charge;
  res.mass = mass;
  res.number = number;
  //res.x1 = x1;
  //res.x3 = x3;
  //res.v1 = v1;
  //res.v2 = v2;
  //res.v3 = v3;
  res.grid_num1 = grid_num1;
  res.grid_num3 = grid_num3;
  res.dr = dr;
  res.dz = dz;
  return res;
}