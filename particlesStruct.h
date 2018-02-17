#ifndef __PARTICLESSTRUCTH__
#define __PARTICLESSTRUCTH__

// typedef double double;

//#define BUILD_OPENCL

struct ParticlesStruct
{
  // The specie charge
  double charge;
  // The specie mass
  double mass;
  // Number of particles
  int number;
  // Particles' coordinates
  //double* x1;    //r
  //double* x3;    //z
  //// Particles' velocity
  //double* v1; //vr
  //double* v2; //vphi
  //double* v3; //vz

  //indicator if particle is still alive
  bool* is_alive;
  int grid_num1;
  int grid_num3;

  double dr;
  double dz;
};

struct Particle
{
  float pos1;
  float pos2;
  float vel1;
  float vel2;
  float vel3;
  bool is_active;
};

ParticlesStruct CreateParticlesStruct(double charge,
                                      double mass,
                                      int number,
                                      int grid_num1,
                                      int grid_num3,
                                      double dr,
                                      double dz);

#endif
