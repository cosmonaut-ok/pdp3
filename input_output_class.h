#pragma once
#include <fstream>
#include<iostream>
#include "particles_struct.h"

using namespace std;
class input_output_class
{
private:
	char path_result[100];
	char path_dump[100];

public:

	input_output_class(void);
	input_output_class(char* cpathres,char* cpathdump);
	~input_output_class(void);
	void out_data(char const *comp_name,
                double** out_value,
                int step_number,
                int number,
                int r_step,
                int z_step);
	void out_field_dump(char* comp_name,double** out_value,int r_step,int z_step);

	void out_coord(char* comp_name,
                 double* coord_r,
                 double* coord_z,
                 int step_number,
                 int number,
                 int particles_number);
	void out_coord_dump(char* comp_name,
                      double* coord_r,
                      double* coord_z,
                      int particles_number);
	void out_velocity_dump(char* comp_name,
                         double* v1,
                         double* v2,
                         double* v3,
                         int particles_number);
};
