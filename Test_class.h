#pragma once
#include "particles_struct.h"
class Test_class
{
public:
	Test_class(void);
	~Test_class(void);
public:
	flcuda get_energy(flcuda *,flcuda, int);
};

