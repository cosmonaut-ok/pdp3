#pragma once
#include "geometry.h"
#include <complex>

using namespace std;

class Fourier
{
public:
  Geometry *goem1;
  Fourier(void);
  ~Fourier(void);
  void fast_sine_transform(double **a, int lenght_n, int ir, bool inv);
  void fast_cosine_transform(double **a, int lenght_n, int ir, bool inv);
  void fast_fourier_alg(complex<double> *a, int lenght_n);
  void fast_fourier_transform(double **a,int lenght_n,int ir, bool inv);
};
