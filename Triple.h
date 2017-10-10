#pragma once
#include "particles_struct.h"

class Triple
{
public:
	Triple(flcuda f, flcuda s, flcuda t);
public:
	~Triple(void);

public:
	flcuda first;
	flcuda second;
	flcuda third;
};
