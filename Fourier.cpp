#include "Fourier.h"
#include <complex>
#include "stdio.h"
using namespace std;
Fourier::Fourier(int t=0)
{
	
}

Fourier::~Fourier(void)
{
}


/////////////////////////////////////////////////////
void Fourier::fast_sine_transform(flcuda** a, int lenght_n, int ir, bool inv)
{
	int d_lenght=2*lenght_n;
	complex<flcuda> im_i (0,1);
	complex<flcuda>* a_sinc = new complex<flcuda> [d_lenght+2];


	int k=0;
	for (k=0;k<lenght_n;k++)
	{
		a_sinc[0]=(0,0);
		a_sinc[lenght_n]=(0,0);
		a_sinc[d_lenght-k+1].real(-a[ir][k]);
		a_sinc[d_lenght-k+1].imag(0);
		a_sinc[k+1].real(a[ir][k]);
		a_sinc[k+1].imag(0);
	}

	 d_lenght=d_lenght+2;
	Fourier::fast_fourier_alg(a_sinc, d_lenght);

	  for (k=0;k<lenght_n;k++)
	  {
		  a_sinc[k+1]=-(flcuda)1.0/((flcuda)2.0*im_i)*a_sinc[k+1];
		  a[ir][k]=real(a_sinc[k+1]);
	  }

	//inverse transform//
	  if(inv==true)
	  {
		  for (k=0;k<lenght_n;k++)
			  a[ir][k]=a[ir][k]*2.0/(lenght_n+1);
	  }
 // delete [] a_sinc;
  return;
}

///////////////////////////////////////////////////
void Fourier::dct_2(flcuda** a, int lenght_n, int ir, bool inv)
{
	int k=0, n=0;
	flcuda pi =  4.0*atan(1.0);
	flcuda t_f =0;
	flcuda* a_t = new flcuda [lenght_n];
	for (k=0;k<lenght_n;k++)
	{
		a_t[k]=a[ir][k];
	}
	flcuda sq =  sqrt((flcuda) (lenght_n));
	for (k=0;k<lenght_n;k++)
	{
		t_f=0;
		for (n=0;n<(lenght_n);n++)
			{
				t_f =t_f+ a_t[n]*cos(pi*k*(n+0.5)/(lenght_n));
			}
	 a[ir][k]=t_f;
	}	

	//inverse transform//
	  if(inv==true)
	  {
		  for (k=0;k<lenght_n;k++)
			  a[ir][k]=a[ir][k]*2.0/(lenght_n-1);
	  }
//delete[] a_sinc;
//a_sinc=0;
  return;
 
}

void Fourier::fast_cosine_transform(flcuda** a, int lenght_n, int ir, bool inv)
{
//	void fast_fourier_alg(complex<flcuda>* a, int lenght_n);
	int d_lenght=2*lenght_n;
	complex<flcuda> im_i (0,1);
	complex<flcuda>* a_sinc = new complex<flcuda> [d_lenght-2];


	int k=0;
	for (k=0;k<(lenght_n-2);k++)
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
	
	  for (k=0;k<lenght_n;k++)
	  {
//		  a_sinc[k]=0.5*a_sinc[k];
		  a[ir][k]=0.5*real(a_sinc[k]);
	  }
//	   delete[] a_sinc;
//		a_sinc=0;
	//inverse transform//
	  if(inv==true)
	  {
		  for (k=0;k<lenght_n;k++)
			  a[ir][k]=a[ir][k]*2.0/(lenght_n-1);
	  }
//delete[] a_sinc;
//a_sinc=0;
  return;
 
}
//////////////////////////////////////////////////
void Fourier::fast_fourier_alg(complex<flcuda>* a, int lenght_n)
{
	//bit reversal part//
////////////////////////////////////////////
complex<flcuda> W;
complex<flcuda> t_comp_a;
complex<flcuda> im_j (0,1);

 int bt,k,j,inver;
 complex<flcuda> temp;
for (k=0; k<lenght_n; k++)
	{
		bt=lenght_n >> 1;
		j=0;
		inver=k;
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
//for(int k=0; k<lenght_n; k++)
//printf("\n%d -> %d",k,a[k]);
//////////////////////////////////////////////
	
//fast fourier part//
/////////////////////////////////////////////

///////////////////////////////////////////

	//fast algorithm//
//////////////////////////////////////////

flcuda pi =  4.0*atan(1.0);
int n_bit=lenght_n-1;
int n_step;
int m;
n_step=1; //step of two point calculation
while (n_bit>0)
	{
		
		for (k=0;k<n_step;k++)
		{
			int kt=k;
			W=exp((-pi*(flcuda)2.0*k/(n_step*2))*im_j);

			for(m=k;m<(lenght_n);m+=2*n_step)
			{
				t_comp_a=a[m];
				a[m]=a[m]+W*a[m+n_step];
				a[m+n_step]=t_comp_a-W*a[m+n_step];


			}
		}
		n_step=n_step*2;
		n_bit=n_bit>>1;
	}
//for(int k=0; k<lenght_n; k++);
//	a[ir][k]=real(comp_a[k]);
return;
}
//////////////////////////////////////////////////////



/////////////////////////////////////////////////////
	//use fast_fourier_alg to real data//
////////////////////////////////////////////////////

void Fourier::fast_fourier_transform(flcuda**a,int lenght_n,int ir, bool inv)
{
	int k=0;
	complex<flcuda>* im_a = new complex<flcuda> [lenght_n];
	for(k=0;k<lenght_n;k++)
	{
		im_a[k].real(a[ir][k]);
		im_a[k].imag(0);
	}

	Fourier::fast_fourier_alg(im_a, lenght_n);

    for(k=0;k<lenght_n;k++)
	{
		a[ir][k]=real(im_a[k]);
	}
	if (inv == true)
	{
		for(k=0;k<lenght_n;k++)
			a[ir][k]=1.0/lenght_n*a[ir][k];
	}
//  delete []im_a;
//  im_a=0;
return;
}
	//fast cosine transform//
/////////////////////////////////////////////////////
	// tnn=2^n//
	// array [0..tnn], number of function vulues = 2^n+1 //
void Fourier::fastcosinetransform_old(flcuda** a, int tnn, bool inversefct, int ir)
{
	flcuda pi = 3.14159265;
    int j;
    int n2;
    flcuda sum;
    flcuda y1;
    flcuda y2;
    flcuda theta;
    flcuda wi;
    flcuda wpi;
    flcuda wr;
    flcuda wpr;
    flcuda wtemp;
    flcuda twr;
    flcuda twi;
    flcuda twpr;
    flcuda twpi;
    flcuda twtemp;
    flcuda ttheta;
    int i;
    int i1;
    int i2;
    int i3;
    int i4;
    flcuda c1;
    flcuda c2;
    flcuda h1r;
    flcuda h1i;
    flcuda h2r;
    flcuda h2i;
    flcuda wrs;
    flcuda wis;
    int nn;
    int ii;
    int jj;
    int n;
    int mmax;
    int m;
    int istep;
    int isign;
    flcuda tempr;
    flcuda tempi;

    if( tnn==1 )
    {
        y1 = a[ir][0];
        y2 = a[ir][1];
        a[ir][0] = 0.5*(y1+y2);
        a[ir][1] = 0.5*(y1-y2);
        if( inversefct )
        {
            a[ir][0] = a[ir][0]*2;
            a[ir][0] = a[ir][1]*2;
        }
        return;
    }
    wi = 0;
    wr = 1;
    theta = pi/tnn;
    wtemp = sin(theta*0.5);
    wpr = -2*wtemp*wtemp;
    wpi = sin(theta);
    sum = 0.5*(a[ir][0]-a[ir][tnn]);
    a[ir][0] = 0.5*(a[ir][0]+a[ir][tnn]);
    n2 = tnn+2;
    for(j = 2; j <= tnn/2; j++)
    {
        wtemp = wr;
        wr = wtemp*wpr-wi*wpi+wtemp;
        wi = wi*wpr+wtemp*wpi+wi;
        y1 = 0.5*(a[ir][j-1]+a[ir][n2-j-1]);
        y2 = a[ir][j-1]-a[ir][n2-j-1];
        a[ir][j-1] = y1-wi*y2;
        a[ir][n2-j-1] = y1+wi*y2;
        sum = sum+wr*y2;
    } 
   ttheta = 2*pi/tnn;
    c1 = 0.5;
    c2 = -0.5;
    isign = 1;
    n = tnn;
    nn = tnn/2;
    j = 1;
    for(ii = 1; ii <= nn; ii++)
    {
        i = 2*ii-1;
        if( j>i )
        {
            tempr = a[ir][j-1];
            tempi = a[ir][j];
            a[ir][j-1] = a[ir][i-1];
            a[ir][j] = a[ir][i];
            a[ir][i-1] = tempr;
            a[ir][i] = tempi;
        }
        m = n/2;
        while(m>=2&&j>m)
        {
            j = j-m;
            m = m/2;
        }
        j = j+m;
    }
    mmax = 2; 
   while(n>mmax)
    {
        istep = 2*mmax;
        theta = 2*pi/(isign*mmax);
        wpr = -2.0*pow((sin(0.5*theta)),2);
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        for(ii = 1; ii <= mmax/2; ii++)
        {
            m = 2*ii-1;
            for(jj = 0; jj <= (n-m)/istep; jj++)
            {
                i = m+jj*istep;
                j = i+mmax;
                tempr = wr*a[ir][j-1]-wi*a[ir][j];
                tempi = wr*a[ir][j]+wi*a[ir][j-1];
                a[ir][j-1] = a[ir][i-1]-tempr;
                a[ir][j] = a[ir][i]-tempi;
                a[ir][i-1] = a[ir][i-1]+tempr;
                a[ir][i] = a[ir][i]+tempi;
            }
            wtemp = wr;
            wr = wr*wpr-wi*wpi+wr;
            wi = wi*wpr+wtemp*wpi+wi;
        }
        mmax = istep;
    }
    twpr = -2.0*pow((sin(0.5*ttheta)),2);
    twpi = sin(ttheta);
    twr = 1.0+twpr;
    twi = twpi;
    for(i = 2; i <= tnn/4+1; i++)
    {
        i1 = i+i-2;
        i2 = i1+1;
        i3 = tnn+1-i2;
        i4 = i3+1;
        wrs = twr;
        wis = twi;
        h1r = c1*(a[ir][i1]+a[ir][i3]);
        h1i = c1*(a[ir][i2]-a[ir][i4]);
        h2r = -c2*(a[ir][i2]+a[ir][i4]);
        h2i = c2*(a[ir][i1]-a[ir][i3]);
        a[ir][i1] = h1r+wrs*h2r-wis*h2i;
        a[ir][i2] = h1i+wrs*h2i+wis*h2r;
        a[ir][i3] = h1r-wrs*h2r+wis*h2i;
        a[ir][i4] = -h1i+wrs*h2i+wis*h2r;
        twtemp = twr;
        twr = twr*twpr-twi*twpi+twr;
        twi = twi*twpr+twtemp*twpi+twi;
    }
    h1r = a[ir][0];
    a[ir][0] = h1r+a[ir][1];
    a[ir][1] = h1r-a[ir][1];
    a[ir][tnn] = a[ir][1];
    a[ir][1] = sum;
    j = 4;
    while(j<=tnn)
    {
        sum = sum+a[ir][j-1];
        a[ir][j-1] = sum;
        j = j+2;
    }
    if( inversefct )
    {
        for(j = 0; j <= tnn; j++)
        {
            a[ir][j] = a[ir][j]*2/tnn;
        }
    }
}
////////////////////////////////////////////////////////////////////

	//fast sine transform//
///////////////////////////////////////////////////////////////////
	//tnn=2^n//
// array [0..tnn-1], number of function vulues = 2^n//
void Fourier::fastsinetransform_old(flcuda** a, int tnn, bool inversefst, int ir)
{
	flcuda pi = 3.14159265;
    int jj;
    int j;
    int tm;
    int n2;
    flcuda sum;
    flcuda y1;
    flcuda y2;
    flcuda theta;
    flcuda wi;
    flcuda wr;
    flcuda wpi;
    flcuda wpr;
    flcuda wtemp;
    flcuda twr;
    flcuda twi;
    flcuda twpr;
    flcuda twpi;
    flcuda twtemp;
    flcuda ttheta;
    int i;
    int i1;
    int i2;
    int i3;
    int i4;
    flcuda c1;
    flcuda c2;
    flcuda h1r;
    flcuda h1i;
    flcuda h2r;
    flcuda h2i;
    flcuda wrs;
    flcuda wis;
    int nn;
    int ii;
    int n;
    int mmax;
    int m;
    int istep;
    int isign;
    flcuda tempr;
    flcuda tempi;

    if( tnn==1 )
    {
        a[ir][0] = 0;
        return;
    }
    theta = pi/tnn;
    wr = 1.0;
    wi = 0.0;
    wpr = -2.0*pow((sin(0.5*theta)),2);
    wpi = sin(theta);
    a[ir][0] = 0.0;
    tm = tnn/2;
    n2 = tnn+2;
    for(j = 2; j <= tm+1; j++)
    {
        wtemp = wr;
        wr = wr*wpr-wi*wpi+wr;
        wi = wi*wpr+wtemp*wpi+wi;
        y1 = wi*(a[ir][j-1]+a[ir][n2-j-1]);
        y2 = 0.5*(a[ir][j-1]-a[ir][n2-j-1]);
        a[ir][j-1] = y1+y2;
        a[ir][n2-j-1] = y1-y2;
    }
    ttheta = pi/tnn;
    c1 = 0.5;
    c2 = -0.5;
    isign = 1;
    n = tnn;
	//bit reversal part//
////////////////////////////////////
    nn = tnn/2;
    j = 1;
    for(ii = 1; ii <= nn; ii++)
    {
        i = 2*ii-1;
        if( j>i )
        {
            tempr = a[ir][j-1];
            tempi = a[ir][j];
            a[ir][j-1] = a[ir][i-1];
            a[ir][j] = a[ir][i];
            a[ir][i-1] = tempr;
            a[ir][i] = tempi;
        }
        m = n/2;
        while(m>=2&&j>m)
        {
            j = j-m;
            m = m/2;
        }
        j = j+m;
    }
/////////////////////////////////////
    mmax = 2;
    while(n>mmax)
    {
        istep = 2*mmax;
        theta = 2*pi/(isign*mmax);
        wpr = -2.0*pow((sin(0.5*theta)),2);
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        for(ii = 1; ii <= mmax/2; ii++)
        {
            m = 2*ii-1;
            for(jj = 0; jj <= (n-m)/istep; jj++)
            {
                i = m+jj*istep;
                j = i+mmax;
                tempr = wr*a[ir][j-1]-wi*a[ir][j];
                tempi = wr*a[ir][j]+wi*a[ir][j-1];
                a[ir][j-1] = a[ir][i-1]-tempr;
                a[ir][j] = a[ir][i]-tempi;
                a[ir][i-1] = a[ir][i-1]+tempr;
                a[ir][i] = a[ir][i]+tempi;
            }
            wtemp = wr;
            wr = wr*wpr-wi*wpi+wr;
            wi = wi*wpr+wtemp*wpi+wi;
        }
        mmax = istep;
    }
    twpr = -2.0*pow((sin(0.5*ttheta)),2);
    twpi = sin(ttheta);
    twr = 1.0+twpr;
    twi = twpi;
    for(i = 2; i <= tnn/4+1; i++)
    {
        i1 = i+i-2;
        i2 = i1+1;
        i3 = tnn+1-i2;
        i4 = i3+1;
        wrs = twr;
        wis = twi;
        h1r = c1*(a[ir][i1]+a[ir][i3]);
        h1i = c1*(a[ir][i2]-a[ir][i4]);
        h2r = -c2*(a[ir][i2]+a[ir][i4]);
        h2i = c2*(a[ir][i1]-a[ir][i3]);
        a[ir][i1] = h1r+wrs*h2r-wis*h2i;
        a[ir][i2] = h1i+wrs*h2i+wis*h2r;
        a[ir][i3] = h1r-wrs*h2r+wis*h2i;
        a[ir][i4] = -h1i+wrs*h2i+wis*h2r;
        twtemp = twr;
        twr = twr*twpr-twi*twpi+twr;
        twi = twi*twpr+twtemp*twpi+twi;
    }
    h1r = a[ir][0];
    a[ir][0] = h1r+a[ir][1];
    a[ir][1] = h1r-a[ir][1];
    sum = 0.0;
    a[ir][0] = 0.5*a[ir][0];
    a[ir][1] = 0.0;
    for(jj = 0; jj <= tm-1; jj++)
    {
        j = 2*jj+1;
        sum = sum+a[ir][j-1];
        a[ir][j-1] = a[ir][j];
        a[ir][j] = sum;
    }
    if( inversefst )
    {
        for(j = 1; j <= tnn; j++)
        {
            a[ir][j-1] = a[ir][j-1]*2/tnn;
        }
    }
}



