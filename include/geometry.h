#pragma once

#include <cmath>

class Geometry
{
public:
  double r_size;
  double z_size;
  int n_grid_r;
  int n_grid_z;
  double dr;
  double dz;
  double **epsilon;
  double **sigma;

  double get_dr();
  double get_dz();

  void set_epsilon();
  void set_pml(double comparative_l_1, double comparative_l_2, double comparative_l_3,
               double sigma1, double sigma2);

	Geometry(double rs, double zs, int ngr, int ngz);
  Geometry();

  ~Geometry(void);

};

// some geometry-related macros
//! \f$ ( \pi \times (dr * (i+0.5))^2 - \pi \times (dr * (i-0.5))^2 ) * dz \f$
#define CELL_VOLUME(i, dr, dz) PI * (dz) * (dr) * (dr) * 2.0 * (i)

//! volume of the cylindrical ring (internal cylinder on r1 is cut out)
#define CYL_RNG_VOL(z, r1, r2) PI * (z) * ((r2) * (r2) - (r1) * (r1))

//! volume of the cylinder
#define CYL_VOL(z, r) PI * (z) * (r) * (r) / 4.

// #define PARTICLE_VOLUME(x,y) (PI * dz * dr * dr * 2.0 * i)

//! get cell number by 'radius'
#define CELL_NUMBER(position, dx) (int)ceil((position) / (dx)) - 1
