#ifndef _PARTICLES_KERNEL_H_
#define _PARTICLES_KERNEL_H_

#include <stdio.h>
#include <math.h>
//#include "cutil_math.h"
#include "math_constants.h"
#include "system_host.cuh"
//#include "particles_struct.h"

__global__ void StepV(flcuda* x1, flcuda* x3, flcuda* v1, flcuda* v2, flcuda* v3, 
					  flcuda* e1, flcuda* e2, flcuda* e3, flcuda* h1, flcuda* h2, 
					  flcuda* h3, bool* is_alive, int number, flcuda timestep);
__device__ int get_linear_coord(int index_r, int index_z, int ngrid_r, int component);
__device__ flcuda get_field_e(flcuda x1, flcuda x3, flcuda* e, int component, Particles_struct cuda_specie);
__device__ flcuda get_field_h(flcuda x1, flcuda x3, flcuda* h, int component, Particles_struct cuda_specie);
__device__ flcuda get_gamma(int i, flcuda* v1, flcuda* v2, flcuda* v3);
__device__ flcuda get_gamma_inv(int i, flcuda* v1, flcuda* v2, flcuda* v3);


//__global__ void CopyParticlesToGLBuffer(Particle* ParticlesArray, float3* positions);

#endif