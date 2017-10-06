#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
//#include <CL/cl.h>
#endif

//cl_context CreateContext();
//cl_command_queue CreateCommandQueue(cl_context context, cl_device_id *device);
//cl_program CreateProgram(cl_context context, cl_device_id device, const char* fileName);
