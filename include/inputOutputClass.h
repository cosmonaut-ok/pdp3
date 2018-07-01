#include <fstream>
#include <iostream>
#include <string>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;
class InputOutputClass
{
private:
  char path_result[100];
  char path_dump[100];

public:

  InputOutputClass(void);
  InputOutputClass(char *cpathres,char *cpathdump);
  ~InputOutputClass(void);
  void out_data(char const *comp_name,
                double **out_value,
                int step_number,
                int number,
                int r_step,
                int z_step);
  void out_field_dump(char *comp_name,double **out_value,int r_step,int z_step);

  void out_pos(char *comp_name,
                 double *pos_r,
                 double *pos_z,
                 int step_number,
                 int number,
                 int particles_number);
  void out_pos_dump(char *comp_name,
                      double **pos,
                      int particles_number);
  void out_velocity_dump(char *comp_name,
                         double **vel,
                         int particles_number);
};
