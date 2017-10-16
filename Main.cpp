// #pragma once

// #include<iostream>
// #include <fstream>
// #include <math.h>

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
