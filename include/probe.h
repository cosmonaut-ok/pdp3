#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h> /* floor */

using namespace std;

class Probe
{
 protected:
  char path_result[100];
  char path_dump[100];

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
  // void write(char *name, double **data);
};
