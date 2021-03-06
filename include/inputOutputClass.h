#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h> /* floor */

using namespace std;
class InputOutputClass
{
private:
  char path_result[100];
  char path_dump[100];

public:
  bool compress = false;
  // char *compress_algo = 'gzip';
  int compress_level = 6;

  virtual ~InputOutputClass(void) {};

  virtual void out_data(char const *comp_name,
                        double **out_value,
                        int step_number,
                        int number,
                        int r_step,
                        int z_step) = 0;
  virtual void out_field_dump(char const *comp_name, double **out_value, int r_step, int z_step) = 0;
  virtual void out_triple(char *comp_name,
                          double **pos,
                          int step_number,
                          int number,
                          int particles_number) = 0;
  virtual void out_pos_dump(char *comp_name,
                            double **pos,
                            int particles_number) = 0;
  virtual void out_velocity_dump(char *comp_name,
                                 double **vel,
                                 int particles_number) = 0;
  virtual void out_current_time_dump(double current_time) = 0;
};
