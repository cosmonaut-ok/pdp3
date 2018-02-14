#pragma once
#include "particles_struct.h"

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

  double set_dr();
  double set_dz();

  void set_epsilon();
  void addPML(double comp_l_1, double comp_l_2, double comp_l_3,
              double sigma1_t, double sigma2_t);
  void calcSigma();

	Geometry(double fs, double ss, int ng1, int ng2);
  Geometry();

  ~Geometry(void);

private:
  double comparative_l_1; // length to left wall TODO: orly?
  double comparative_l_2; // length to right wall TODO: orly?
  double comparative_l_3; // length to z-wall
  double sigma1;
  double sigma2;

};
