#include "wrapper.h"

void DisplayPlatformInfo(cl_platform_id id, cl_platform_info name, std::string str)
{
	cl_int errNum;
	std::size_t paramValueSize;
	errNum = clGetPlatformInfo(id, name, 0,	NULL, &paramValueSize);
	if (errNum != CL_SUCCESS)
	{
		std::cerr << "Failed to find OpenCL platform "
		<< str << "." << std::endl;
		return;
	}
	char * info = (char *)alloca(sizeof(char) * paramValueSize);
	errNum = clGetPlatformInfo(id, name, paramValueSize, info,	NULL);
	if (errNum != CL_SUCCESS)
	{
		std::cerr << "Failed to find OpenCL platform "
		<< str << "." << std::endl;
		return;
	}
	std::cout << "\t" << str << ":\t" << info << std::endl;
}

template<typename T> void appendBitfield(T info, T value, std::string name, std::string & str)
{
	if (info & value)
	{
		if (str.length() > 0)
		{
			str.append(" | ");
		}
		str.append(name);
	}
}
template <typename T> class InfoDevice
{
public:
static void display(
cl_device_id id, cl_device_info name, std::string str)
{
	cl_int errNum;
	std::size_t paramValueSize;
	errNum = clGetDeviceInfo(id, name, 0, NULL, &paramValueSize);
	if (errNum != CL_SUCCESS)
	{
	std::cerr << "Failed to find OpenCL device info "
	<< str << "." << std::endl;
	return;
	}

	T * info = (T *)alloca(sizeof(T) * paramValueSize);
	errNum = clGetDeviceInfo(id,name,paramValueSize,info,NULL);
	if (errNum != CL_SUCCESS)
	{
	std::cerr << "Failed to find OpenCL device info "
	<< str << "." << std::endl;
	return;
	}

	switch (name)
	{
		case CL_DEVICE_TYPE:
		{
			std::string deviceType;
			appendBitfield<cl_device_type>(
			*(reinterpret_cast<cl_device_type*>(info)),
			CL_DEVICE_TYPE_CPU, "CL_DEVICE_TYPE_CPU", deviceType);
			appendBitfield<cl_device_type>(
			*(reinterpret_cast<cl_device_type*>(info)),
			CL_DEVICE_TYPE_GPU, "CL_DEVICE_TYPE_GPU", deviceType);
			appendBitfield<cl_device_type>(
			*(reinterpret_cast<cl_device_type*>(info)),
			CL_DEVICE_TYPE_ACCELERATOR,
			"CL_DEVICE_TYPE_ACCELERATOR",
			deviceType);
			appendBitfield<cl_device_type>(
			*(reinterpret_cast<cl_device_type*>(info)),
			CL_DEVICE_TYPE_DEFAULT,
			"CL_DEVICE_TYPE_DEFAULT",
			deviceType);
			std::cout << "\t\t" << str << ":\t"
			<< deviceType << std::endl;
		}
		break;
		case CL_DEVICE_SINGLE_FP_CONFIG:
		{
			std::string fpType;
			appendBitfield<cl_device_fp_config>(
			*(reinterpret_cast<cl_device_fp_config*>(info)),
			CL_FP_DENORM, "CL_FP_DENORM", fpType);
			appendBitfield<cl_device_fp_config>(
			*(reinterpret_cast<cl_device_fp_config*>(info)),
			CL_FP_INF_NAN, "CL_FP_INF_NAN", fpType);
			appendBitfield<cl_device_fp_config>(
			*(reinterpret_cast<cl_device_fp_config*>(info)),
			CL_FP_ROUND_TO_NEAREST,
			"CL_FP_ROUND_TO_NEAREST",
			fpType);
			appendBitfield<cl_device_fp_config>(
			*(reinterpret_cast<cl_device_fp_config*>(info)),
			CL_FP_ROUND_TO_ZERO, "CL_FP_ROUND_TO_ZERO", fpType);
			appendBitfield<cl_device_fp_config>(
			*(reinterpret_cast<cl_device_fp_config*>(info)),
			CL_FP_ROUND_TO_INF, "CL_FP_ROUND_TO_INF", fpType);
			appendBitfield<cl_device_fp_config>(
			*(reinterpret_cast<cl_device_fp_config*>(info)),
			CL_FP_FMA, "CL_FP_FMA", fpType);
			appendBitfield<cl_device_fp_config>(
			*(reinterpret_cast<cl_device_fp_config*>(info)),
			CL_FP_SOFT_FLOAT, "CL_FP_SOFT_FLOAT", fpType);
			std::cout << "\t\t" << str << ":\t" << fpType << std::endl;
		}
		break;
		case CL_DEVICE_GLOBAL_MEM_CACHE_TYPE:
		{
			std::string memType;
			appendBitfield<cl_device_mem_cache_type>(
			*(reinterpret_cast<cl_device_mem_cache_type*>(info)),
			CL_NONE, "CL_NONE", memType);
			appendBitfield<cl_device_mem_cache_type>(
			*(reinterpret_cast<cl_device_mem_cache_type*>(info)),
			CL_READ_ONLY_CACHE, "CL_READ_ONLY_CACHE", memType);
			appendBitfield<cl_device_mem_cache_type>(
			*(reinterpret_cast<cl_device_mem_cache_type*>(info)),
			CL_READ_WRITE_CACHE, "CL_READ_WRITE_CACHE", memType);
			std::cout << "\t\t" << str << ":\t" << memType << std::endl;
		}
		break;
		case CL_DEVICE_LOCAL_MEM_TYPE:
		{
			std::string memType;
			appendBitfield<cl_device_local_mem_type>(
			*(reinterpret_cast<cl_device_local_mem_type*>(info)),
			CL_GLOBAL, "CL_LOCAL", memType);
			appendBitfield<cl_device_local_mem_type>(
			*(reinterpret_cast<cl_device_local_mem_type*>(info)),
			CL_GLOBAL, "CL_GLOBAL", memType);
			std::cout << "\t\t" << str << ":\t" << memType << std::endl;
		}
		break;
		case CL_DEVICE_EXECUTION_CAPABILITIES:
		{
			std::string memType;
			appendBitfield<cl_device_exec_capabilities>(
			*(reinterpret_cast<cl_device_exec_capabilities*>(info)),
			CL_EXEC_KERNEL, "CL_EXEC_KERNEL", memType);
			appendBitfield<cl_device_exec_capabilities>(
			*(reinterpret_cast<cl_device_exec_capabilities*>(info)),
			CL_EXEC_NATIVE_KERNEL, "CL_EXEC_NATIVE_KERNEL", memType);
			std::cout << "\t\t" << str << ":\t" << memType << std::endl;
		}
		break;
		case CL_DEVICE_QUEUE_PROPERTIES:
		{
			std::string memType;
			appendBitfield<cl_device_exec_capabilities>(
			*(reinterpret_cast<cl_device_exec_capabilities*>(info)),
			CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE,
			"CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE", memType);
			appendBitfield<cl_device_exec_capabilities>(
			*(reinterpret_cast<cl_device_exec_capabilities*>(info)),
			CL_QUEUE_PROFILING_ENABLE, "CL_QUEUE_PROFILING_ENABLE",
			memType);
			std::cout << "\t\t" << str << ":\t" << memType << std::endl;
		}
		break;
		default:
			std::cout << "\t\t" << str << ":\t" << *info << std::endl;
			break;
	}
}
};

///
//  Create an OpenCL context on the first available platform using
//  either a GPU or CPU depending on what is available.
//
cl_context CreateContext()
{
    cl_int errNum;
    cl_uint numPlatforms;
    cl_platform_id firstPlatformId;
    cl_context context = NULL;

    // First, select an OpenCL platform to run on.  For this example, we
    // simply choose the first available platform.  Normally, you would
    // query for all available platforms and select the most appropriate one.
    errNum = clGetPlatformIDs(1, &firstPlatformId, &numPlatforms);
    if (errNum != CL_SUCCESS || numPlatforms <= 0)
    {
        std::cerr << "Failed to find any OpenCL platforms." << std::endl;
        return NULL;
    }

    // Next, create an OpenCL context on the platform.  Attempt to
    // create a GPU-based context, and if that fails, try to create
    // a CPU-based context.
    cl_context_properties contextProperties[] =
    {
        CL_CONTEXT_PLATFORM,
        (cl_context_properties)firstPlatformId,
        0
    };
    context = clCreateContextFromType(contextProperties, CL_DEVICE_TYPE_GPU,
                                      NULL, NULL, &errNum);
    if (errNum != CL_SUCCESS)
    {
        std::cout << "Could not create GPU context, trying CPU..." << std::endl;
        context = clCreateContextFromType(contextProperties, CL_DEVICE_TYPE_CPU,
                                          NULL, NULL, &errNum);
        if (errNum != CL_SUCCESS)
        {
            std::cerr << "Failed to create an OpenCL GPU or CPU context." << std::endl;
            return NULL;
        }
    }

	cl_device_id deviceIds;
	errNum = clGetDeviceIDs(firstPlatformId,
							CL_DEVICE_TYPE_GPU,
							1,
							&deviceIds,
							NULL);

	DisplayPlatformInfo(firstPlatformId, CL_PLATFORM_PROFILE, "CL_PLATFORM_PROFILE");
	DisplayPlatformInfo(firstPlatformId, CL_PLATFORM_NAME, "CL_PLATFORM_NAME");
	DisplayPlatformInfo(firstPlatformId, CL_PLATFORM_VENDOR, "CL_PLATFORM_VENDOR");
	DisplayPlatformInfo(firstPlatformId, CL_PLATFORM_VERSION, "CL_PLATFORM_VERSION");
	DisplayPlatformInfo(firstPlatformId, CL_PLATFORM_EXTENSIONS, "CL_PLATFORM_EXTENSIONS");

	InfoDevice<cl_uint>::display(deviceIds,	CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, "DEVICE has max work item dimensions");
	InfoDevice<cl_uint>::display(deviceIds,	CL_DEVICE_MAX_WORK_ITEM_SIZES, "DEVICE has max work item sizes");
	InfoDevice<cl_uint>::display(deviceIds,	CL_DEVICE_MAX_COMPUTE_UNITS, "DEVICE has max compute units");
	InfoDevice<cl_uint>::display(deviceIds,	CL_DEVICE_MAX_WORK_GROUP_SIZE, "DEVICE has max work group size");
	InfoDevice<cl_uint>::display(deviceIds,	CL_DEVICE_GLOBAL_MEM_SIZE, "DEVICE has global mem size");
	InfoDevice<cl_uint>::display(deviceIds,	CL_DEVICE_ADDRESS_BITS, "DEVICE has address space of size");
	InfoDevice<cl_uint>::display(deviceIds,	CL_DEVICE_LOCAL_MEM_SIZE, "DEVICE has local mem of size");
	InfoDevice<cl_uint>::display(deviceIds,	CL_DEVICE_SINGLE_FP_CONFIG, "DEVICE has FP config");

    return context;
}

///
//  Create a command queue on the first device available on the
//  context
//
cl_command_queue CreateCommandQueue(cl_context context, cl_device_id *device)
{
    cl_int errNum;
    cl_device_id *devices;
    cl_command_queue commandQueue = NULL;
    size_t deviceBufferSize = -1;

    // First get the size of the devices buffer
    errNum = clGetContextInfo(context, CL_CONTEXT_DEVICES, 0, NULL, &deviceBufferSize);
    if (errNum != CL_SUCCESS)
    {
        std::cerr << "Failed call to clGetContextInfo(...,GL_CONTEXT_DEVICES,...)";
        return NULL;
    }

    if (deviceBufferSize <= 0)
    {
        std::cerr << "No devices available.";
        return NULL;
    }

    // Allocate memory for the devices buffer
    devices = new cl_device_id[deviceBufferSize / sizeof(cl_device_id)];
    errNum = clGetContextInfo(context, CL_CONTEXT_DEVICES, deviceBufferSize, devices, NULL);
    if (errNum != CL_SUCCESS)
    {
        delete [] devices;
        std::cerr << "Failed to get device IDs";
        return NULL;
    }

    // In this example, we just choose the first available device.  In a
    // real program, you would likely use all available devices or choose
    // the highest performance device based on OpenCL device queries
    commandQueue = clCreateCommandQueue(context, devices[0], 0, NULL);
    if (commandQueue == NULL)
    {
        delete [] devices;
        std::cerr << "Failed to create commandQueue for device 0";
        return NULL;
    }

    *device = devices[0];
    delete [] devices;
    return commandQueue;
}

///
//  Create an OpenCL program from the kernel source file
//
cl_program CreateProgram(cl_context context, cl_device_id device, const char* fileName)
{
    cl_int errNum;
    cl_program program;

    std::ifstream kernelFile(fileName, std::ios::in);
    if (!kernelFile.is_open())
    {
        std::cerr << "Failed to open file for reading: " << fileName << std::endl;
        return NULL;
    }

    std::ostringstream oss;
    oss << kernelFile.rdbuf();

    std::string srcStdStr = oss.str();
    const char *srcStr = srcStdStr.c_str();
    program = clCreateProgramWithSource(context, 1,
                                        (const char**)&srcStr,
                                        NULL, NULL);
    if (program == NULL)
    {
        std::cerr << "Failed to create CL program from source." << std::endl;
        return NULL;
    }

    errNum = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
    if (errNum != CL_SUCCESS)
    {
        // Determine the reason for the error
        char buildLog[16384*5];
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG,
                              sizeof(buildLog), buildLog, NULL);

        std::cerr << "Error in kernel: " << std::endl;
        std::cerr << buildLog;
        clReleaseProgram(program);
        return NULL;
    }

    return program;
}
