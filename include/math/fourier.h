#pragma once

#include <complex>

#include "constant.h"
#include "geometry.h"

using namespace std;

namespace math
{
  namespace fourier
  {
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
  }
}
