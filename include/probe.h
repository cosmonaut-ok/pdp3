#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h> /* floor */

#include "lib.h"
#include "particles.h"

using namespace std;

class Probe
{
 protected:
  char path_result[100];
  char data_root[100];

 protected:
  void write_frame(double **data);
  void write_col(double **data);
  void write_row(double **data);
  void write_dot(double **data);

 public:
  char name[100];
  char path[256];
  char component[10];
  int type; // 0 - frame, 1 - col, 2 - row, 3 - dot
  int schedule = 1;

  int step_number = 0;

  bool compress = false;
  // char *compress_algo = 'gzip';
  int compress_level = 6;

  int start_r = 0;
  int start_z = 0;
  int end_r = 0;
  int end_z = 0;

  virtual ~Probe(void) {};
  virtual void write(char *name, double **out_value) = 0;

  //! mpwrite is a function for MacroParticles frame write output
  //! writes two sets of files:
  //! 3-velocity components and 2-position components
  //! writes to files <name>_vel_[r,phi,z].dat and <name>_pos_[r,z].dat
  virtual void mpwrite(char *name, Particles *p_specie) = 0;
};
