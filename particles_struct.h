#ifndef __PARTICLESSTRUCTH__
#define __PARTICLESSTRUCTH__

typedef double flcuda;

//#define BUILD_OPENCL


struct Particles_struct
{
	// The specie charge
	flcuda charge;
	// The specie mass
	flcuda mass;
	// Number of particles 
	int number;
	// Particles' coordinates
	//flcuda* x1;  //r  
	//flcuda* x3;  //z
	//// Particles' velocity
	//flcuda* v1; //vr
	//flcuda* v2; //vphi
	//flcuda* v3; //vz

	//indicator if particle is still alive
	bool* is_alive;
	int grid_num1;
	int grid_num3;

	flcuda dr;
	flcuda dz;
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

Particles_struct CreateParticles_struct(flcuda charge, flcuda mass, int number, int grid_num1, int grid_num3, flcuda dr, flcuda dz);


#endif