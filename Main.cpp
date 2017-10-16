#include<iostream>
#include "Constant.h"
#include "field.h"
#include "E_field.h"
#include "H_field.h"
#include "pdp3_time.h"
#include "Particles.h"
#include "Fourier.h"
#include "Poisson.h"
#include "poisson_neumann.h"
#include "poisson_dirichlet.h"
#include "particles_list.h"
#include <fstream>
#include <math.h>

#include "Boundary_Maxwell_conditions.h"
#include "input_output_class.h"
// #include "Beam.h"
#include "Bunch.h"
#include "time.h"
#include "particles_struct.h"
#include "Load_init_param.h"

using namespace std;

Particles_struct specie;
// #define BUILD_OPENCL

int main(int argc, char **argv)
{
	clock_t start, finish;
	double time_elapsed;
	Load_init_param init_param((char *) "parameters.xml");

	init_param.run();

}
