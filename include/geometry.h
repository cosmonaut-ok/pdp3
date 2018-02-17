#pragma once
#include "particlesStruct.h"

class Geometry
{
public:
  double first_size;
  double second_size;
  int n_grid_1;
  int n_grid_2;
  double dr;
  double dz;
  double** epsilon;
  double** sigma;

  double get_dr();
  double get_dz();

  void set_epsilon();
  void set_pml(double comparative_l_1, double comparative_l_2, double comparative_l_3,
              double sigma1, double sigma2);

	Geometry(double fs, double ss, int ng1, int ng2);
  Geometry();

  ~Geometry(void);

};
