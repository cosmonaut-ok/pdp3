#include "Test_class.h"


Test_class::Test_class(void)
{
}


Test_class::~Test_class(void)
{
}

flcuda Test_class::get_energy(flcuda * vel,flcuda mass, int number)
{
	flcuda energy =0;
	for(int i=0;i<(number-1);i++)
		 energy+= vel[i]*vel[i]*mass/2.0;
	return energy;
}