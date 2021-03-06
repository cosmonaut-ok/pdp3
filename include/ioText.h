#include <fstream>
#include <iostream>
#include <string>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h> /* floor */
#include "inputOutputClass.h"

using namespace std;
class IOText : public InputOutputClass
{
private:
  char path_result[100];

  char path_dump[100];

  void dump_data(char const *path,
                 double **out_value,
                 int r_step,
                 int z_step,
                 bool is_rewrite);

  void dump_components(char const *path,
                       double **components,
                       int size,
                       int dimensions,
                       bool is_rewrite);

  char* get_data_path(char const *comp_name, int number);

  char* get_dump_path(char const *comp_name);

public:
  IOText(void);
  IOText(char *cpathres,char *cpathdump, bool c_compress);
  ~IOText(void);

  void out_data(char const *comp_name,
                double **out_value,
                int step_number,
                int number,
                int r_step,
                int z_step);
  void out_field_dump(char const *comp_name, double **out_value, int r_step, int z_step);
  void out_triple(char *comp_name,
                  double **pos,
                  int step_number,
                  int number,
                  int particles_number);
  void out_pos_dump(char *comp_name,
                    double **pos,
                    int particles_number);
  void out_velocity_dump(char *comp_name,
                         double **vel,
                         int particles_number);
  void out_current_time_dump(double current_time);
};
