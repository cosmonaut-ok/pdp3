#include <fstream>
#include <iostream>
#include <string>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h> /* floor */
#include <hdf5.h>
#include "probe.h"

using namespace std;
class ProbeHDF5 : public virtual Probe
{
 private:
  char hdf5_file[100];
  char probe_path[200];
  char probe_data_set_name[200];

  void write_data(char *out_value, bool is_rewrite);
  char* get_data_path(char *name);

  hid_t create_or_open_h5_group(hid_t file_id, char const *group_name);
  hid_t create_or_open_h5_dataset(hid_t file_id, hid_t group_id, char const *ds_name);

 protected:
  void write_frame(char *name, double **data, bool is_rewrite);
  void write_col(char *name, double **data, bool is_rewrite);
  void write_row(char *name, double **data, bool is_rewrite);
  void write_dot(char *name, double **data, bool is_rewrite);

 public:
  ProbeHDF5(void);

  ProbeHDF5(char *c_path, char *c_data_root, char *c_component, int c_type,
              int c_start_r, int c_start_z, int c_end_r, int c_end_z,
              bool c_compress, int c_compress_level, int schedule);
  ~ProbeHDF5(void);

  void write(char *name, double **out_value);
  void mpwrite(char *name, Particles *p_specie);
};
