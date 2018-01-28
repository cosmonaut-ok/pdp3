#include "particles_struct.h"

Particles_struct CreateParticles_struct(double charge,
                                        double mass,
                                        int number,
                                        int grid_num1,
                                        int grid_num3,
                                        double dr,
                                        double dz)
{
  Particles_struct res;
  res.charge = charge;
  res.mass = mass;
  res.number = number;
  res.grid_num1 = grid_num1;
  res.grid_num3 = grid_num3;
  res.dr = dr;
  res.dz = dz;
  return res;
}
