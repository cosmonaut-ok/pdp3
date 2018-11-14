#include <fstream>
#include <iostream>
#include <string>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h> /* floor */
#include "writer.h"
/* #include <iterator> */
/* #include <vector> */

using namespace std;
class WriterPlain : public virtual Writer
{
 private:
  void write_data(char *out_value, bool is_rewrite);
  char* get_data_path(char *name);

 protected:
  void write_frame(char *name, double **data, bool is_rewrite);
  void write_col(char *name, double **data, bool is_rewrite);
  void write_row(char *name, double **data, bool is_rewrite);
  void write_dot(char *name, double **data, bool is_rewrite);

 public:
  int upf = 1; // how units (frames/cols/rows/dots) per file
  WriterPlain(void);
  // WriterPlain(char *cpathres,char *cpathdump, bool c_compress);


  WriterPlain(char *c_path, char *c_component, int c_type,
              int c_start_r, int c_start_z, int c_end_r, int c_end_z,
              bool c_compress, int c_compress_level, int schedule);
  ~WriterPlain(void);

  void write(char *name, double **out_value);
};
