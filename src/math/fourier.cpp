#include "math/fourier.h"

using namespace std;
using namespace constant;

namespace math
{
  namespace fourier
  {

// TODO: https://people.sc.fsu.edu/~jburkardt/cpp_src/fft_openmp/fft_openmp.cpp

    Fourier::~Fourier(void)
    {
    }

    void Fourier::fast_sine_transform(double **a, int lenght_n, int ir, bool inv)
    {
      int d_lenght=2*lenght_n;
      complex<double> im_i (0,1);
      complex<double> *a_sinc = new complex<double> [d_lenght+2];

      a_sinc[0]=(0.0);
      a_sinc[lenght_n]=(0.0);

#pragma omp parallel for shared (lenght_n, d_lenght, a_sinc)
      for (int k=0; k<lenght_n; k++)
      {
        a_sinc[d_lenght-k+1].real(-a[ir][k]);
        a_sinc[d_lenght-k+1].imag(0);
        a_sinc[k+1].real(a[ir][k]);
        a_sinc[k+1].imag(0);
      }

      d_lenght=d_lenght+2;
      Fourier::fast_fourier_alg(a_sinc, d_lenght);

#pragma omp parallel for shared (im_i, a_sinc)
      for (int k=0; k<lenght_n; k++)
      {
        a_sinc[k+1]=-(double)1.0/((double)2.0*im_i)*a_sinc[k+1];
        a[ir][k]=real(a_sinc[k+1]);
      }

      //inverse transform//
      if(inv==true)
      {
#pragma omp parallel for
        for (int k=0; k<lenght_n; k++)
          a[ir][k]=a[ir][k]*2.0/(lenght_n+1);
      }
    }

    void Fourier::fast_cosine_transform(double **a, int lenght_n, int ir, bool inv)
    {
//  void fast_fourier_alg(complex<double> *a, int lenght_n);
      int d_lenght=2*lenght_n;
      complex<double> im_i (0,1);
      complex<double> *a_sinc = new complex<double> [d_lenght-2];

#pragma omp parallel for shared(lenght_n, d_lenght)
      for (int k=0;k<(lenght_n-2);k++)
      {
        a_sinc[lenght_n+k].real(a[ir][lenght_n-k-2]);
        a_sinc[d_lenght-k-2].imag(0);
        a_sinc[k].real(a[ir][k]);
        a_sinc[k].imag(0);
      }
      a_sinc[lenght_n-1].real(a[ir][lenght_n-1]);
      a_sinc[lenght_n-1].imag(0);
      a_sinc[lenght_n-2].real(a[ir][lenght_n-2]);
      a_sinc[lenght_n-2].imag(0);

      d_lenght=d_lenght-2;
      Fourier::fast_fourier_alg(a_sinc, d_lenght);

#pragma omp parallel for shared(a_sinc)
      for (int k=0; k<lenght_n;k++)
      {
        // a_sinc[k]=0.5*a_sinc[k];
        a[ir][k]=0.5*real(a_sinc[k]);
      }

      //inverse transform//
      if(inv==true)
      {
#pragma omp parallel for shared (lenght_n)
        for (int k=0; k<lenght_n; k++)
          a[ir][k]=a[ir][k]*2.0/(lenght_n-1);
      }
    }

    void Fourier::fast_fourier_alg(complex<double> *a, int lenght_n)
    {
#pragma omp parallel for shared (lenght_n)
      for (int k=0; k<lenght_n; k++)
      {
        int bt=lenght_n >> 1;
        int j = 0;
        int inver = k;
        complex<double> temp;

        while (bt>0)
        {
          j+=(inver % 2)*bt;
          inver >>= 1;
          bt >>= 1;
        }
        if (k<j)
        {
          temp = a[k];
          a[k]=a[j];
          a[j]=temp;
        }

      }

      complex<double> im_j (0,1);

      int n_bit = lenght_n-1;
      int n_step = 1; //step of two point calculation
      while (n_bit>0)
      {
#pragma omp parallel for shared (im_j, n_step, lenght_n)
        for (int k=0;k<n_step;k++)
        {
          complex<double> W;

          // int kt=k;
          W=exp((-PI*(double)2.0*k/(n_step*2))*im_j);

          for(int m=k; m<(lenght_n); m+=2*n_step)
          {
            complex<double> t_comp_a;

            t_comp_a=a[m];
            a[m]=a[m]+W*a[m+n_step];
            a[m+n_step]=t_comp_a-W*a[m+n_step];
          }
        }
        n_step=n_step*2;
        n_bit=n_bit>>1;
      }
    }

// use fast_fourier_alg to real data
    void Fourier::fast_fourier_transform(double **a, int lenght_n, int ir, bool inv)
    {
      complex<double> *im_a = new complex<double> [lenght_n];

#pragma omp parallel for shared (im_a, ir)
      for(int k=0;k<lenght_n;k++)
      {
        im_a[k].real(a[ir][k]);
        im_a[k].imag(0);
      }

      Fourier::fast_fourier_alg(im_a, lenght_n);

#pragma omp parallel for shared (im_a, ir)
      for(int k=0;k<lenght_n;k++)
      {
        a[ir][k]=real(im_a[k]);
      }

      if (inv == true)
      {
#pragma omp parallel for shared (ir, lenght_n)
        for(int k=0;k<lenght_n;k++)
          a[ir][k]=1.0/lenght_n*a[ir][k];
      }
    }
  }
}
