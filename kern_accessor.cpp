#include "kern_accessor.h"

KernAccessor::KernAccessor(int size1, int size3)
{

    m_context = CreateContext();

	m_commandQueue = CreateCommandQueue(m_context, &m_device);

	m_program = CreateProgram(m_context, m_device, "kern.cl");
	m_program_jw = CreateProgram(m_context, m_device, "j_weighting.cl");

	m_kernel1 = clCreateKernel(m_program, "step_coord_kernel", NULL);
	m_kernel2 = clCreateKernel(m_program, "kernel_stepV", NULL);
	m_kernel_jw = clCreateKernel(m_program_jw, "j_weighting_kernel", NULL);

	// x1
	memObjects[0] = clCreateBuffer(m_context, CL_MEM_READ_WRITE,
									sizeof(double) * N_MAX_PART, NULL, NULL);

	// x2
	memObjects[1] = clCreateBuffer(m_context, CL_MEM_READ_WRITE,
									sizeof(double) * N_MAX_PART, NULL, NULL);

	// v1
	memObjects[2] = clCreateBuffer(m_context, CL_MEM_READ_WRITE,
									sizeof(double) * N_MAX_PART, NULL, NULL);

	// v3
	memObjects[3] = clCreateBuffer(m_context, CL_MEM_READ_WRITE,
									sizeof(double) * N_MAX_PART, NULL, NULL);

	// is_alive
	memObjects[5] = clCreateBuffer(m_context, CL_MEM_READ_ONLY,
									sizeof(int) * N_MAX_PART, NULL, NULL);

	// params
	memObjects[4] = clCreateBuffer(m_context, CL_MEM_READ_ONLY,
									sizeof(double) * 10, NULL, NULL);

	// v2
	memObjects[6] = clCreateBuffer(m_context, CL_MEM_READ_WRITE,
									sizeof(double) * N_MAX_PART, NULL, NULL);

	// e1
	memObjects[7] = clCreateBuffer(m_context, CL_MEM_READ_WRITE,
									sizeof(double) * size1 * size3, NULL, NULL);

	// e2
	memObjects[8] = clCreateBuffer(m_context, CL_MEM_READ_WRITE,
									sizeof(double) * size1 * size3, NULL, NULL);

	// e3
	memObjects[9] = clCreateBuffer(m_context, CL_MEM_READ_WRITE,
									sizeof(double) * size1 * size3, NULL, NULL);

	// h1
	memObjects[10] = clCreateBuffer(m_context, CL_MEM_READ_WRITE,
									sizeof(double) * size1 * size3, NULL, NULL);

	// h2
	memObjects[11] = clCreateBuffer(m_context, CL_MEM_READ_WRITE,
									sizeof(double) * size1 * size3, NULL, NULL);

	// h3
	memObjects[12] = clCreateBuffer(m_context, CL_MEM_READ_WRITE,
									sizeof(double) * size1 * size3, NULL, NULL);

	// j1
	memObjects[13] = clCreateBuffer(m_context, CL_MEM_READ_WRITE,
									sizeof(double) * size1 * size3, NULL, NULL);

	// j2
	memObjects[14] = clCreateBuffer(m_context, CL_MEM_READ_WRITE,
									sizeof(double) * size1 * size3, NULL, NULL);

	// j3
	memObjects[15] = clCreateBuffer(m_context, CL_MEM_READ_WRITE,
									sizeof(double) * size1 * size3, NULL, NULL);

}

KernAccessor::~KernAccessor()
{
}

void KernAccessor::halfStepCoord(double *x1, double *x3, double *v1, double *v3, double *params, int *is_alive, int n_particles)
{
	cl_int errNum;
	errNum = clEnqueueWriteBuffer(m_commandQueue, memObjects[0],
									CL_TRUE, 0, n_particles * sizeof(double),
									x1, 0, NULL, NULL);

	errNum = clEnqueueWriteBuffer(m_commandQueue, memObjects[1],
									CL_TRUE, 0, n_particles * sizeof(double),
									x3, 0, NULL, NULL);

	errNum = clEnqueueWriteBuffer(m_commandQueue, memObjects[2],
									CL_TRUE, 0, n_particles * sizeof(double),
									v1, 0, NULL, NULL);

	errNum = clEnqueueWriteBuffer(m_commandQueue, memObjects[3],
									CL_TRUE, 0, n_particles * sizeof(double),
									v3, 0, NULL, NULL);

	errNum = clEnqueueWriteBuffer(m_commandQueue, memObjects[4],
									CL_TRUE, 0, 10 * sizeof(double),
									params, 0, NULL, NULL);

	errNum = clEnqueueWriteBuffer(m_commandQueue, memObjects[5],
									CL_TRUE, 0, n_particles * sizeof(int),
									is_alive, 0, NULL, NULL);


	errNum =  clSetKernelArg(m_kernel1, 0, sizeof(cl_mem), &memObjects[0]);
	errNum |= clSetKernelArg(m_kernel1, 1, sizeof(cl_mem), &memObjects[1]);
	errNum |= clSetKernelArg(m_kernel1, 2, sizeof(cl_mem), &memObjects[2]);
	errNum |= clSetKernelArg(m_kernel1, 3, sizeof(cl_mem), &memObjects[3]);
	errNum |= clSetKernelArg(m_kernel1, 4, sizeof(cl_mem), &memObjects[4]);
	errNum |= clSetKernelArg(m_kernel1, 5, sizeof(cl_mem), &memObjects[5]);
	errNum |= clSetKernelArg(m_kernel1, 6, sizeof(cl_mem), (void *)(&n_particles));

	size_t globalWorkSize1[1] = { (n_particles / 64 + 1) * 64 };
	size_t localWorkSize1[1] = { 64 };
	// Queue the kernel up for execution across the array
	errNum = clEnqueueNDRangeKernel(m_commandQueue, m_kernel1, 1, NULL,
									globalWorkSize1, localWorkSize1,
									0, NULL, NULL);

	errNum = clEnqueueReadBuffer(m_commandQueue, memObjects[0],
									CL_TRUE, 0, n_particles * sizeof(double),
									x1, 0, NULL, NULL);

	errNum = clEnqueueReadBuffer(m_commandQueue, memObjects[1],
									CL_TRUE, 0, n_particles * sizeof(double),
									x3, 0, NULL, NULL);

	errNum = clEnqueueReadBuffer(m_commandQueue, memObjects[2],
									CL_TRUE, 0, n_particles * sizeof(double),
									v1, 0, NULL, NULL);

	errNum = clEnqueueReadBuffer(m_commandQueue, memObjects[3],
									CL_TRUE, 0, n_particles * sizeof(double),
									v3, 0, NULL, NULL);


}

void KernAccessor::stepV(double *x1, double *x3, double *v1, double *v2, double *v3, 
		                 double *e1, double *e2, double *e3, double *h1, double *h2, double *h3,
		                 double *params, int *is_alive, int n_particles, int grid_num1, int grid_num3)
{

	cl_int errNum;
	errNum = clEnqueueWriteBuffer(m_commandQueue, memObjects[0],
									CL_TRUE, 0, n_particles * sizeof(double),
									x1, 0, NULL, NULL);

	errNum = clEnqueueWriteBuffer(m_commandQueue, memObjects[1],
									CL_TRUE, 0, n_particles * sizeof(double),
									x3, 0, NULL, NULL);

	errNum = clEnqueueWriteBuffer(m_commandQueue, memObjects[2],
									CL_TRUE, 0, n_particles * sizeof(double),
									v1, 0, NULL, NULL);

	errNum = clEnqueueWriteBuffer(m_commandQueue, memObjects[3],
									CL_TRUE, 0, n_particles * sizeof(double),
									v3, 0, NULL, NULL);

	errNum = clEnqueueWriteBuffer(m_commandQueue, memObjects[4],
									CL_TRUE, 0, 10 * sizeof(double),
									params, 0, NULL, NULL);

	errNum = clEnqueueWriteBuffer(m_commandQueue, memObjects[5],
									CL_TRUE, 0, n_particles * sizeof(int),
									is_alive, 0, NULL, NULL);
	//v2
	errNum = clEnqueueWriteBuffer(m_commandQueue, memObjects[6],
									CL_TRUE, 0, n_particles * sizeof(double),
									v2, 0, NULL, NULL);

	//e1
	errNum = clEnqueueWriteBuffer(m_commandQueue, memObjects[7],
									CL_TRUE, 0, (grid_num1 - 1) * grid_num3 * sizeof(double),
									e1, 0, NULL, NULL);

	//e2
	errNum = clEnqueueWriteBuffer(m_commandQueue, memObjects[8],
									CL_TRUE, 0, grid_num1 * grid_num3 * sizeof(double),
									e2, 0, NULL, NULL);

	//e3
	errNum = clEnqueueWriteBuffer(m_commandQueue, memObjects[9],
									CL_TRUE, 0, grid_num1 * (grid_num3 - 1) * sizeof(double),
									e3, 0, NULL, NULL);

	//h1
	errNum = clEnqueueWriteBuffer(m_commandQueue, memObjects[10],
									CL_TRUE, 0, grid_num1 * (grid_num3 - 1) * sizeof(double),
									h1, 0, NULL, NULL);

	//h2
	errNum = clEnqueueWriteBuffer(m_commandQueue, memObjects[11],
									CL_TRUE, 0, (grid_num1 - 1) * (grid_num3 - 1) * sizeof(double),
									h2, 0, NULL, NULL);

	//h3
	errNum = clEnqueueWriteBuffer(m_commandQueue, memObjects[12],
									CL_TRUE, 0, (grid_num1 - 1) * grid_num3 * sizeof(double),
									h3, 0, NULL, NULL);



	errNum =  clSetKernelArg(m_kernel2, 0, sizeof(cl_mem), &memObjects[0]);
	errNum |= clSetKernelArg(m_kernel2, 1, sizeof(cl_mem), &memObjects[1]);
	errNum |= clSetKernelArg(m_kernel2, 2, sizeof(cl_mem), &memObjects[2]);
	errNum |= clSetKernelArg(m_kernel2, 3, sizeof(cl_mem), &memObjects[6]);
	errNum |= clSetKernelArg(m_kernel2, 4, sizeof(cl_mem), &memObjects[3]);
	errNum |= clSetKernelArg(m_kernel2, 5, sizeof(cl_mem), &memObjects[7]);
	errNum |= clSetKernelArg(m_kernel2, 6, sizeof(cl_mem), &memObjects[8]);
	errNum |= clSetKernelArg(m_kernel2, 7, sizeof(cl_mem), &memObjects[9]);
	errNum |= clSetKernelArg(m_kernel2, 8, sizeof(cl_mem), &memObjects[10]);
	errNum |= clSetKernelArg(m_kernel2, 9, sizeof(cl_mem), &memObjects[11]);
	errNum |= clSetKernelArg(m_kernel2, 10, sizeof(cl_mem), &memObjects[12]);

	errNum |= clSetKernelArg(m_kernel2, 11, sizeof(cl_mem), &memObjects[5]);
	errNum |= clSetKernelArg(m_kernel2, 12, sizeof(cl_mem), &memObjects[4]);
	errNum |= clSetKernelArg(m_kernel2, 13, sizeof(cl_mem), (void *)(&n_particles));

	size_t globalWorkSize1[1] = { (n_particles / 64 + 1) * 64 };
	size_t localWorkSize1[1] = { 64 };
	// Queue the kernel up for execution across the array
	errNum = clEnqueueNDRangeKernel(m_commandQueue, m_kernel2, 1, NULL,
									globalWorkSize1, localWorkSize1,
									0, NULL, NULL);


	errNum = clEnqueueReadBuffer(m_commandQueue, memObjects[2],
									CL_TRUE, 0, n_particles * sizeof(double),
									v1, 0, NULL, NULL);

	errNum = clEnqueueReadBuffer(m_commandQueue, memObjects[6],
									CL_TRUE, 0, n_particles * sizeof(double),
									v2, 0, NULL, NULL);

	errNum = clEnqueueReadBuffer(m_commandQueue, memObjects[3],
									CL_TRUE, 0, n_particles * sizeof(double),
									v3, 0, NULL, NULL);
}
/*
__kernel void j_weighting_kernel(double dt, 
                                 __global double *j1, 
				 		         __global double *j2, 
								 __global double *j3, 
								 __global double *x1_o, 
								 __global double *x3_o, 
				                 __global double *x1, 
								 __global double *x3, 
								 __global bool   *is_alive, 
								 __global const double *params,
								 int number)
*/
void KernAccessor::j_weighting(current *j, double *x1, double *x3, double *j1, double *j2, double *j3, 
		                 double *x1_o, double *x3_o, int *is_alive, double *params,
						 int n_particles, int grid_num1, int grid_num3)
{
	
	cl_int errNum;
	double *temp = new double[grid_num1 * grid_num3];

    for (int i = 0; i < grid_num1; i++)
        for (int j = 0; j < grid_num3; j++)
        {
            temp[i * grid_num3 + j] = 0.0;
        }

	errNum = clEnqueueWriteBuffer(m_commandQueue, memObjects[0],
									CL_TRUE, 0, n_particles * sizeof(double),
									x1, 0, NULL, NULL);

	errNum = clEnqueueWriteBuffer(m_commandQueue, memObjects[1],
									CL_TRUE, 0, n_particles * sizeof(double),
									x3, 0, NULL, NULL);

	errNum = clEnqueueWriteBuffer(m_commandQueue, memObjects[2],
									CL_TRUE, 0, n_particles * sizeof(double),
									x1_o, 0, NULL, NULL);

	errNum = clEnqueueWriteBuffer(m_commandQueue, memObjects[3],
									CL_TRUE, 0, n_particles * sizeof(double),
									x3_o, 0, NULL, NULL);


	errNum = clEnqueueWriteBuffer(m_commandQueue, memObjects[4],
									CL_TRUE, 0, 10 * sizeof(double),
									params, 0, NULL, NULL);

	errNum = clEnqueueWriteBuffer(m_commandQueue, memObjects[5],
									CL_TRUE, 0, n_particles * sizeof(int),
									is_alive, 0, NULL, NULL);

	errNum = clEnqueueWriteBuffer(m_commandQueue, memObjects[13],
									CL_TRUE, 0, grid_num1 * grid_num3 * sizeof(double),
									temp, 0, NULL, NULL);

	errNum = clEnqueueWriteBuffer(m_commandQueue, memObjects[14],
									CL_TRUE, 0, grid_num1 * grid_num3 * sizeof(double),
									temp, 0, NULL, NULL);

	errNum = clEnqueueWriteBuffer(m_commandQueue, memObjects[15],
									CL_TRUE, 0, grid_num1 * grid_num3 * sizeof(double),
									temp, 0, NULL, NULL);



	errNum =  clSetKernelArg(m_kernel_jw, 0, sizeof(cl_mem), &memObjects[13]);
	errNum |= clSetKernelArg(m_kernel_jw, 1, sizeof(cl_mem), &memObjects[14]);
	errNum |= clSetKernelArg(m_kernel_jw, 2, sizeof(cl_mem), &memObjects[15]);
	errNum |= clSetKernelArg(m_kernel_jw, 3, sizeof(cl_mem), &memObjects[0]);
	errNum |= clSetKernelArg(m_kernel_jw, 4, sizeof(cl_mem), &memObjects[1]);
	errNum |= clSetKernelArg(m_kernel_jw, 5, sizeof(cl_mem), &memObjects[2]);
	errNum |= clSetKernelArg(m_kernel_jw, 6, sizeof(cl_mem), &memObjects[3]);


	errNum |= clSetKernelArg(m_kernel_jw, 7, sizeof(cl_mem), &memObjects[5]);
	errNum |= clSetKernelArg(m_kernel_jw, 8, sizeof(cl_mem), &memObjects[4]);
	errNum |= clSetKernelArg(m_kernel_jw, 9, sizeof(cl_mem), (void *)(&n_particles));

	size_t globalWorkSize1[1] = { (n_particles / 64 + 1) * 64 };
	size_t localWorkSize1[1] = { 64 };
	// Queue the kernel up for execution across the array
	errNum = clEnqueueNDRangeKernel(m_commandQueue, m_kernel_jw, 1, NULL,
									globalWorkSize1, localWorkSize1,
									0, NULL, NULL);


	errNum = clEnqueueReadBuffer(m_commandQueue, memObjects[13],
									CL_TRUE, 0, grid_num1 * grid_num3 * sizeof(double),
									temp, 0, NULL, NULL);
	j->j1_add_1d(temp);

	errNum = clEnqueueReadBuffer(m_commandQueue, memObjects[14],
									CL_TRUE, 0, grid_num1 * grid_num3 * sizeof(double),
									temp, 0, NULL, NULL);
    j->j2_add_1d(temp);

	errNum = clEnqueueReadBuffer(m_commandQueue, memObjects[15],
									CL_TRUE, 0, grid_num1 * grid_num3 * sizeof(double),
									temp, 0, NULL, NULL);
	j->j3_add_1d(temp);
	delete temp;
}