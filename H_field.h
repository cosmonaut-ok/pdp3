#pragma once
#include "field.h"
#include "Geometry.h"
//#include "E_field.h"
#include "Time.h"
#include "Particles.h"
#include "Fourier.h"
#include "current.h"
#include "Triple.h"
#include "particles_struct.h"
class E_field;
class H_field
{
public:
	flcuda** h1;
	flcuda** h2;
	flcuda** h3;
	flcuda* h1_1d; 
	flcuda* h2_1d; 
	flcuda* h3_1d; 
	flcuda** h1_half_time;
	flcuda** h2_half_time;
	flcuda** h3_half_time;
	flcuda** Ar;
	flcuda** Afi;
	flcuda** Az;
	const Geometry* geom1;
	H_field(Geometry* geom1);
	H_field(void);
	~H_field(void);
void calc_field(E_field* e_field1, Time* time1);
void set_homogeneous_h(flcuda H1, flcuda H2, flcuda H3);
void magnetostatic_equation(Geometry* geom1);
Triple get_field(flcuda x1, flcuda x3);
	 flcuda* get_1d_h1();
	 flcuda* get_1d_h2();
	 flcuda* get_1d_h3();

};
