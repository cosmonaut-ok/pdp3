#ifndef __SYSTEMHOST__
#define __SYSTEMHOST__

//#include <vector_types.h>	
//#include <driver_types.h>			// for cudaStream_t

#include "particles_struct.h"

typedef unsigned int		uint;		// should be 4-bytes on CUDA


extern "C"
{

bool InitCUDA(void);
void SetupCUDA(int grid_num1, int grid_num3, int number);
void CopySpecie2Cuda (Particles_struct specie);
void TransferXVToCUDA (flcuda* CPU_x1, flcuda* CPU_x3, flcuda* CPU_v1, flcuda* CPU_v2, flcuda* CPU_v3, bool* CPU_is_alive, int number);
void TransferEHToCUDA (flcuda* CPU_e1, flcuda* CPU_e2, flcuda* CPU_e3, flcuda* CPU_h1, flcuda* CPU_h2, flcuda* CPU_h3, int grid_num1, int grid_num3);
void TransferEHFromCUDA (flcuda* CPU_e1, flcuda* CPU_e2, flcuda* CPU_e3, flcuda* CPU_h1, flcuda* CPU_h2, flcuda* CPU_h3, int grid_num1, int grid_num3);
void TransferXVFromCUDA (flcuda* CPU_x1, flcuda* CPU_x3, flcuda* CPU_v1, flcuda* CPU_v2, flcuda* CPU_v3, int number);
void TransferToCUDA (Particle* CPU_ParticlesArray, int numPoints );
void TransferFromCUDA ( Particle* CPU_ParticlesArray, int numPoints );


void CUDA_StepV(int number, flcuda dt);
//void CUDA_CopyParticlesToGLBuffer(float3* positions);


//void computeNumBlocks (int numPnts, int maxThreads, int &numBlocks, int &numThreads);


}

#endif
