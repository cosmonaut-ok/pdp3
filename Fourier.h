#pragma once
#include "Geometry.h"
#include <complex>
using namespace std;
class Fourier
{
public:
	 Geometry* goem1;
	  Fourier(int t);
	~Fourier(void);
	 void fast_sine_transform(flcuda** a, int lenght_n, int ir, bool inv);
	 void fast_cosine_transform(flcuda** a, int lenght_n, int ir, bool inv);
	 void fast_fourier_alg(complex<flcuda>* a, int lenght_n);
	 void fast_fourier_transform(flcuda**a,int lenght_n,int ir, bool inv);
	 void fastcosinetransform_old(flcuda** a, int tnn, bool inversefct, int k);
	 void fastsinetransform_old(flcuda** a, int tnn, bool inversefst, int ir);
	 void dct_2(flcuda** a, int lenght_n, int ir, bool inv);
	
};
