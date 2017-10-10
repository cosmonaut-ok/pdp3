#pragma once
#include"E_field.h"
#include"H_field.h"
#include"Geometry.h"

class Boundary_conditions
{
public:
	Boundary_conditions(void);
	~Boundary_conditions(void);
	
	virtual void specify_boundary_conditions()=0;
};
