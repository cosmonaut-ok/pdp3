//#include "CL/cl.h"
#include "wrapper.h"
#include "current.h"
#define N_MAX_PART 3000000
//class KernAccessor
//{
//public:
//	KernAccessor(int size1, int size3);
//	~KernAccessor();
//	void halfStepCoord(double *x1, double *x3, double *v1, double *v3, double *params, int *is_alive, int n_particles);
//	void stepV(double *x1, double *x3, double *v1, double *v2, double *v3, 
//		       double *e1, double *e2, double *e3, double *h1, double *h2, double *h3,
//		       double *params, int *is_alive, int n_particles, int grid_num1, int grid_num3);
//
//    void j_weighting(current *j, double *x1, double *x3, double *j1, double *j2, double *j3, 
//		             double *x1_o, double *x3_o, int *is_alive, double *params,
//					 int n_particles, int grid_num1, int grid_num3);
//
//private:
    //cl_context m_context;
    //cl_command_queue m_commandQueue;
    //cl_program m_program, m_program_jw;
    //cl_device_id m_device;
    //cl_kernel m_kernel1, m_kernel2, m_kernel_jw;
    //cl_mem memObjects[16]; 
//};