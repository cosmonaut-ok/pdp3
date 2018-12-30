#pragma once

#include <cmath>

class Geometry
{
public:
  double first_size;
  double second_size;
  int n_grid_1;
  int n_grid_2;
  double dr;
  double dz;
  double **epsilon;
  double **sigma;

  double get_dr();
  double get_dz();

  void set_epsilon();
  void set_pml(double comparative_l_1, double comparative_l_2, double comparative_l_3,
              double sigma1, double sigma2);

	Geometry(double fs, double ss, int ng1, int ng2);
  Geometry();

  ~Geometry(void);

};

// some geometry-related macros
//! \f$ ( \pi \times (dr * (i+0.5))^2 - \pi \times (dr * (i-0.5))^2 ) * dz \f$
#define CELL_VOLUME(i, dr, dz) PI * (dz) * (dr) * (dr) * 2.0 * (i)

//! volume of the cylindrical ring
#define CYL_RNG_VOL(z, r1, r2) PI * (z) * ((r2) * (r2) - (r1) * (r1))

// #define PARTICLE_VOLUME(x,y) (PI * dz * dr * dr * 2.0 * i)

//! get cell number by 'radius'
#define CELL_NUMBER(position, dx) (int)ceil((position) / (dx)) - 1
