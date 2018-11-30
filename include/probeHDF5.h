#include <fstream>
#include <iostream>
#include <string>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h> /* floor */
#include "probe.h"
#include <hdf5.h> // TODO: is it true only for debian-based distros?

using namespace std;
class ProbeHDF5 : public Probe
{
private:
  char path_data[100];
  char path_dump[100];
  char path_result[100];

  char hdf5_file[100];
  char root_group_name[100];
  char result_group_name[100];
  char state_group_name[100];

  herr_t status;

  void dump_h5_dataset(char const *group_name,
                       char const *name,
                       double **out_value,
                       int r_step,
                       int z_step);

  void dump_h5_nd_components(char const *group_name,
                             char const *name,
                             double **components,
                             int size,
                             int dimensions);

  hid_t create_or_open_h5_group(hid_t file_id, char const *group_name);

public:
  ProbeHDF5(void);
  ProbeHDF5(char *cpathres,char *cpathdump, bool c_compress);
  ~ProbeHDF5(void);

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
