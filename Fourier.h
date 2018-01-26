#pragma once
#include "Geometry.h"
#include <complex>

using namespace std;

class Fourier
{
public:
  // Geometry* goem1;
  Fourier(void);
  ~Fourier(void);
  void fast_sine_transform(double** a, int lenght_n, int ir, bool inv);
  void fast_cosine_transform(double** a, int lenght_n, int ir, bool inv);
  void fast_fourier_alg(complex<double>* a, int lenght_n);
  void fast_fourier_transform(double**a,int lenght_n,int ir, bool inv);
  void fastcosinetransform_old(double** a, int tnn, bool inversefct, int k);
  void fastsinetransform_old(double** a, int tnn, bool inversefst, int ir);
  void dct_2(double** a, int lenght_n, int ir, bool inv);

};
