#include "E_field.h"
#include "H_field.h"
#include "math.h"
#include "Fourier.h"
#include <fstream>
#include <Constant.h>

using namespace std;
// using namespace constant;
const double PI = constant::PI; // reuse to prevent "magic numbers"
const double EPSILON0 = constant::EPSILON0; // reuse to prevent "magic numbers"

E_field::E_field(): epsilon0(EPSILON0)
{
}

////constructor
E_field::E_field(Geometry* geom1_t): geom1(geom1_t), epsilon0(EPSILON0)
{
	//Er//
	///////////////////////////////////
	//n_grid - number of edges
	e1 = new double*[geom1->n_grid_1-1];
	for (int i=0; i<(geom1->n_grid_1-1);i++)
	{
		e1[i]= new double[geom1->n_grid_2];
	}
	e1_1d = new double[(geom1->n_grid_1-1)*geom1->n_grid_2];
	///////////////////////////////////
	//Ef//
	e2 = new double*[geom1->n_grid_1];

	for (int i=0; i<(geom1->n_grid_1);i++)
	{
		e2[i]= new double[geom1->n_grid_2];
	}
	e2_1d = new double[geom1->n_grid_1*geom1->n_grid_2];
	//////////////////////////////////////
	//Ez//
	//////////////////////////////////////
	e3 = new double*[geom1->n_grid_1];

	for (int i=0; i<(geom1->n_grid_1);i++)
	{
		e3[i]= new double[geom1->n_grid_2-1];
	}
	e3_1d = new double[geom1->n_grid_1*(geom1->n_grid_2-1)];
	//////////////////////////////////////////
	//fi//
	/////////////////////////////////////////
	fi = new double*[geom1->n_grid_1];

	for (int i=0; i<(geom1->n_grid_1);i++)
	{
		fi[i]= new double[geom1->n_grid_2];
	}
	////////////////////////////////////////
	//charge_density//
	///////////////////////////////////////
	t_charge_density = new double*[geom1->n_grid_1];
	for (int i=0; i<(geom1->n_grid_1);i++)
	{
		t_charge_density[i] = new double[geom1->n_grid_2];
	}
	/////////////////////////////////////
	/////////////////////////////////////
	//set zero initial temp ro//
	for (int i=0; i<(geom1->n_grid_1);i++)
		for (int k=0; k<(geom1->n_grid_2);k++)
		{
			t_charge_density[i][k]=0;
		}
}

// Destructor. Freeing memory
E_field::~E_field()
{
	for (int i=0; i<(geom1->n_grid_1-1);i++)
		delete[]e1[i];
	delete[]e1;

	for (int i=0; i<(geom1->n_grid_1);i++)
		delete[]e2[i];
	delete[]e2;

	for (int i=0; i<(geom1->n_grid_1);i++)
		delete[]e3[i];
	delete[]e3;

	for (int i=0; i<(geom1->n_grid_1);i++)
		delete[]fi[i];
	delete[]fi;

	for (int i=0; i<(geom1->n_grid_1);i++)
		delete[]t_charge_density[i];
	delete[]t_charge_density;
}
//initial E distribution//
////////////////////////////////////////////
///////////////////////////////////////////
void E_field::set_homogeneous_efield(double E1, double E2, double E3)
{
	////Er////
	for(int i=0;i<(geom1->n_grid_1-1);i++)
		for(int k=0;k<(geom1->n_grid_2-1);k++)
		{
			e1[i][k]= E1;
		}
	///Ef////
	for(int i=0;i<(geom1->n_grid_1-1);i++)
		for(int k=0;k<(geom1->n_grid_2-1);k++)
		{
			e2[i][k]=E2;
		}
	///Ez////
	for(int i=0;i<(geom1->n_grid_1-1);i++)
		for(int k=0;k<(geom1->n_grid_2-1);k++)
		{
			e3[i][k]=E3;
		}
}
/////////////////////////////////////////////////

////////////////////////////////////////////////

void E_field::boundary_conditions()
{
	//border Er conditions//
	////////////////////////////////
	//last elements of array [ngrid-1]!!!//
	///////////////////////////////
	for(int i=0;i<(geom1->n_grid_1-1);i++)
	{
		e1[i][0]=0;
		e1[i][geom1->n_grid_2-1]=0;
		//e1[i][0]=limit_conditions();
	}
	//////////////////////////////

	//Border Ef conditions//
	/////////////////////////////
	for(int k=0;k<(geom1->n_grid_2);k++)
	{
		e2[0][k]=0;
		e2[geom1->n_grid_1-1][k]=0;
	}
	for(int i=0;i<(geom1->n_grid_1);i++)
	{
		e2[i][0]=0;
		e2[i][geom1->n_grid_2-1]=0;
	}

	////////////////////////////
	//Border Ez conditions//
	/////////////////////////////////
	for(int k=0;k<(geom1->n_grid_2-1);k++)
	{
		e3[0][k]=0;
		e3[geom1->n_grid_1-1][k]=0;
		//e1[0][k]=limit_conditions();
	}
	////////////////////////////////
}

//Electric field calculation with absorbing fields on walls//
///////////////////////////////////////////////////////////////////////////////
void E_field::calc_field(H_field* h_field1,
                         Time* time1,
                         current* current1,
                         PML* pml1)
{
	int i=0;
	int k=0;
	double** j1= current1->get_j1();
	double** j2= current1->get_j2();
	double**j3= current1->get_j3();
	double koef_e = 0;
	double koef_h =0;
	double** h1 = h_field1->h1;
	double** h2 = h_field1->h2;
	double** h3 = h_field1->h3;
	double dr = geom1->dr;
	double dz = geom1->dz;

	///////Er first[i] value////
	for(k=1;k<(geom1->n_grid_2-1);k++)
	{
		i=0;
		koef_e = (2.0*geom1->epsilon[i][k]*epsilon0 - geom1->sigma[i][k]*time1->delta_t) / (2.0*geom1->epsilon[i][k]*epsilon0 + geom1->sigma[i][k]*time1->delta_t);
		koef_h =	2*time1->delta_t/(2.0*geom1->epsilon[i][k]*epsilon0 + geom1->sigma[i][k]*time1->delta_t);
		e1[i][k]=e1[i][k] * koef_e	- (j1[i][k]+(h2[i][k] - h2[i][k-1])/dz)*koef_h;
	}
	/////////////////////////////////////
	//Ez=on axis// // ???????????
	for(k=0;k<(geom1->n_grid_2-1);k++)
	{
		i=0;
		koef_e = (2.0*geom1->epsilon[i][k]*epsilon0 - geom1->sigma[i][k]*time1->delta_t) / (2.0*geom1->epsilon[i][k]*epsilon0 + geom1->sigma[i][k]*time1->delta_t);
		koef_h =	2*time1->delta_t/(2.0*geom1->epsilon[i][k]*epsilon0 + geom1->sigma[i][k]*time1->delta_t);

		e3[i][k]=e3[i][k]*koef_e - (j3[i][k]-4.0/dr*h2[i][k])*koef_h;
	}

	for( i=1;i<(geom1->n_grid_1-1);i++)
		for( k=1;k<(geom1->n_grid_2-1);k++)
		{
			koef_e = (2.0*geom1->epsilon[i][k]*epsilon0 - geom1->sigma[i][k]*time1->delta_t) / (2.0*geom1->epsilon[i][k]*epsilon0 + geom1->sigma[i][k]*time1->delta_t);
			koef_h =	2*time1->delta_t/(2.0*geom1->epsilon[i][k]*epsilon0 + geom1->sigma[i][k]*time1->delta_t);
			e1[i][k]=e1[i][k]*koef_e - (j1[i][k]+(h2[i][k]-h2[i][k-1])/dz)*koef_h;
			e2[i][k]=e2[i][k]*koef_e - (j2[i][k]-(h1[i][k]-h1[i][k-1])/dz + (h3[i][k]-h3[i-1][k])/dr)*koef_h;
			e3[i][k]=e3[i][k]*koef_e -(j3[i][k]-(h2[i][k]-h2[i-1][k])/dr - (h2[i][k]+h2[i-1][k])/(2.0*dr*i))*koef_h;
		}

	for(i=1;i<(geom1->n_grid_1-1);i++)
	{
		k=0;
		koef_e = (2.0*geom1->epsilon[i][k]*epsilon0 - geom1->sigma[i][k]*time1->delta_t) / (2.0*geom1->epsilon[i][k]*epsilon0 + geom1->sigma[i][k]*time1->delta_t);
		koef_h =	2*time1->delta_t/(2.0*geom1->epsilon[i][k]*epsilon0 + geom1->sigma[i][k]*time1->delta_t);
		e3[i][k]=e3[i][k]*koef_e - (j3[i][k] - (h2[i][k]-h2[i-1][k])/dr - (h2[i][k]+h2[i-1][k])/(2.0*dr*i))*koef_h;
	}

	return;
}
////////////////////////////////////////////////////////////////////////////////

//Electric field calculation sigma=0//
///////////////////////////////////////////////////////////////////////////////
void E_field::calc_field(H_field* h_field1, Time* time1, current* current1)
{
	int i=0;
	int k=0;
	double** j1= current1->get_j1();
	double** j2= current1->get_j2();
	double**j3= current1->get_j3();
	double** h1 = h_field1->h1;
	double** h2 = h_field1->h2;
	double** h3 = h_field1->h3;
	double dr = geom1->dr;
	double dz = geom1->dz;

	///////Er first[i] value////
	for(k=1;k<(geom1->n_grid_2-1);k++)
	{
		i=0;
		e1[i][k]=e1[i][k]-(j1[i][k])*time1->delta_t/epsilon0;
	}

	//Ez on axis//
	for(k=0;k<(geom1->n_grid_2-1);k++)
	{
		i=0;
		e3[i][k]=e3[i][k]-(j3[i][k]-2.0/dr*h2[i][k])*time1->delta_t/epsilon0;
	}

	for( i=1;i<(geom1->n_grid_1-1);i++)
		for( k=1;k<(geom1->n_grid_2-1);k++)
		{
			e1[i][k]=e1[i][k] - (j1[i][k]+(h2[i][k]-h2[i][k-1])/dz)*time1->delta_t/epsilon0;
			e2[i][k]=e2[i][k] - (j2[i][k]-(h1[i][k]-h1[i][k-1])/dz + (h3[i][k]-h3[i-1][k])/dr)*time1->delta_t/epsilon0;;
			e3[i][k]=e3[i][k] - (j3[i][k]-(h2[i][k]-h2[i-1][k])/dr - (h2[i][k]+h2[i-1][k])/(2.0*dr*i))*time1->delta_t/epsilon0;
		}

	for(i=1;i<(geom1->n_grid_1-1);i++)
	{
		k=0;
		e3[i][k]=e3[i][k] -(j3[i][k]-(h2[i][k]-h2[i-1][k])/dr - (h2[i][k]+h2[i-1][k])/(2.0*dr*i))*time1->delta_t/epsilon0;
	}

	return;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
	//set fi on r=z boundary//
////////////////////////////////////////////////////////////////////////////////
void E_field::set_fi_on_z()
{
	for (int k=0; k<(geom1->n_grid_2);k++)
	{
		fi[geom1->n_grid_1-1][k]=0;
	}
}

//poisson equation solving 2//
///////////////////////////////////////////////////////////////////////////////
void E_field::poisson_equation2(Geometry* geom1, charge_density* ro1)
{
	const double epsilon0 = EPSILON0;
	int i=0;
	int k=0;
	double phi0(0.0);
	double *a = new double [geom1->n_grid_1];
	double *b = new double [geom1->n_grid_1];
	double *c = new double [geom1->n_grid_1];
	double *d = new double [geom1->n_grid_1];
	double *c1 = new double [geom1->n_grid_1];
	double *d1 = new double [geom1->n_grid_1];
	double *phi = new double [geom1->n_grid_1];
	Fourier* four1=0;

	// double dr = geom1->dr;
	double dr2 = geom1->dr*geom1->dr;

	double** ro=ro1->get_ro();

	///////////////////////////////////////////////
	//copy charge_density array in to temp array//
	for( i=0;i<(geom1->n_grid_1);i++)
		for( k=0;k<(geom1->n_grid_2);k++)
		{
			t_charge_density[i][k]= ro[i][k];
		}

	//call function for cosine transform//
////////////////////////////////////////////////
////////////////////////////////////////////////

	int temp=geom1->n_grid_2;

	for (i=0;i<geom1->n_grid_1;i++)
	{
		four1->fast_cosine_transform((double**)t_charge_density, temp, i, false);
		//four1->fast_fourier_transform(t_charge_density, temp, i, false);
	}

	//set coefficients
	b[0] = 1.0;
	c[0] = 1.0;
	for (k = 0; k < geom1->n_grid_2; k++)
	{
		b[0] = 1.0;
		c[0] = -1.0;
		d[0] = dr2/4.0/epsilon0*t_charge_density[0][k];
		for (i = 1; i < geom1->n_grid_1 -1; i++)
		{
			a[i] = (1.0 - 1.0/2.0/(double)i);
			b[i] = -2.0 + 2.0*(cos(PI*k/(geom1->n_grid_2-1)) - 1)*geom1->dr*geom1->dr/(geom1->dz*geom1->dz);
			//b[i] = -2.0 + 2.0*(cos(PI*k/geom1->n_grid_2) - 1)/(geom1->dz*geom1->dz);
			c[i] = (1.0 + 1.0/2.0/(double)i);
			d[i] = t_charge_density[i][k]*dr2/epsilon0;
			c1[i] = (1.0 - 1.0/2.0/(double)i);
			d1[i] = t_charge_density[i][k]*dr2/epsilon0;
		}
		a[0] = 0;
		c[geom1->n_grid_1-2] = 0.0;
		d[geom1->n_grid_1-2] -= phi0;
		c1[geom1->n_grid_1-2] = 0.0;
		d1[geom1->n_grid_1-2] -= phi0;
		TridiagonalSolve(a, b, c, d, phi, geom1->n_grid_1-1);
		for (i = 0; i < geom1->n_grid_1 -1; i++)
		{
			fi[i][k] = phi[i];
		}
	}

// call function for inverse cosine transform//
////////////////////////////////////////////////////////////
	for (i=0;i<geom1->n_grid_1;i++)
	{
		int temp=geom1->n_grid_2;
		four1->fast_cosine_transform((double**)fi, temp,i, true);
		//four1->fast_fourier_transform(fi, temp,i, true);
	}
////////////////////////////////////////////////////////////

	//calculate electric field//
///////////////////////////////////////////////////////////
	for (i=0;i<(geom1->n_grid_1-1);i++)
		for (k=1;k<(geom1->n_grid_2);k++)
		{
			e1[i][k]=(fi[i][k]-fi[i+1][k])/geom1->dr;
		}

	for (i=0;i<(geom1->n_grid_1);i++)
		for (k=1;k<(geom1->n_grid_2-1);k++)
		{
			e3[i][k]=(fi[i][k]-fi[i][k+1])/geom1->dz;
		}

	ofstream er("er"), ez("ez"); // TODO: WTF?
	for (i=1;i<geom1->n_grid_1-1;i++)
		for (k=1;k<geom1->n_grid_2-1;k++)
		{
			er<<e1[i][k]<<" ";
			ez<<e3[i][k]<<" ";
		}
	er.close();
	ez.close();

}
/////////////////////////////////////////////////////////////

void E_field::TridiagonalSolve(const double *a,
                               const double *b,
                               double *c,
                               double *d, double *x,
                               unsigned int n)
{
	int i;

	/* Modify the coefficients. */
	c[0] /= b[0];				/* Division by zero risk. */
	d[0] /= b[0];				/* Division by zero would imply a singular matrix. */
	for(int i = 1; i < n; i++){
		double id = (b[i] - c[i-1] * a[i]);	/* Division by zero risk. */
		c[i] /= id;				/* Last value calculated is redundant. */
		d[i] = (d[i] - d[i-1] * a[i])/id;
	}

	/* Now back substitute. */
	x[n - 1] = d[n - 1];
	for(i = n - 2; i >= 0; i--)
		x[i] = d[i] - c[i] * x[i + 1];
}

////////////////////////////////////////////////////////////
	//function for electric field weighting//
Triple E_field::get_field(double x1, double x3)
{
	int i_r=0;	// number of particle i cell
	int k_z=0;	// number of particle k cell

	// double PI = 3.1415926535897932;
	double dr = geom1->dr;
	double dz = geom1->dz;
	double r1, r2, r3; // temp variables for calculation
	double dz1, dz2;	 //	 temp var.: width of k and k+1 cell
	double er =0;
	double efi =0;
	double ez =0;
	double vol_1 =0; //	 volume of i cell; Q/V, V - volume of elementary cell
	double vol_2 =0; //	 volume of i+1 cell;

	// double value =0;
////////////////////////
	r1 = x1-0.5*dr;
	r3 = x1+0.5*dr;
///////////////////////

	// weighting of E_r//
///////////////////////////////////////////////////
	//finding number of cell. example dr=0.5,	 x1 = 0.7, i_r =0;!!
	i_r = (int)ceil((x1-0.5*dr)/geom1->dr)-1;
	k_z = (int)ceil((x3)/geom1->dz)-1;

	vol_1 = PI*dz*dr*dr*(2*i_r+1);
	vol_2 = PI*dz*dr*dr*(2*i_r+3);
	dz1 = (k_z+1)*dz-x3;
	dz2 = x3 - k_z*dz;
	r2 = (i_r+1)*dr;
	///////////////////////////////////////

	//weighting Er[i][k]//
	er = er + e1[i_r][k_z]*(PI*dz1*(r2*r2-r1*r1))/vol_1;

	//weighting Er[i+1][k]//
	er = er + e1[i_r+1][k_z]*(PI*dz1*(r3*r3-r2*r2))/vol_2;

	//weighting Er[i][k+1]//
	er= er + e1[i_r][k_z+1]*(PI*dz2*(r2*r2-r1*r1))/vol_1;

	//weighting Er[i+1][k+1]//
	er = er + e1[i_r+1][k_z+1]*(PI*dz2*(r3*r3-r2*r2))/vol_2;

///////////////////////////////////////////////////////

	// weighting of E_z//
///////////////////////////////////////////////////////
//finding number of cell. example dz=0.5,	 x3 = 0.7, z_k =0;!!
	i_r = (int)ceil((x1)/geom1->dr)-1;
	k_z = (int)ceil((x3-0.5*dz)/geom1->dz)-1;

///////////////////////////////////

	if(x1>dr)
	{
		vol_1 = PI*dz*dr*dr*2*i_r;
	}
	else
	{
		vol_1 = PI*dz*dr*dr/4.0; //volume of first cell
	}
	r2 = (i_r+0.5)*dr;
	vol_2 = PI*dz*dr*dr*(2*i_r+2);
	dz1 = (k_z+1.5)*dz - x3;
	dz2 = x3 - (k_z+0.5)*dz;
	//////////////////////////////////////

	//weighting Ez[i][k]//
	ez = ez + e3[i_r][k_z]*(PI*dz1*(r2*r2-r1*r1))/vol_1;

	//weighting Ez[i+1][k]//
	ez = ez + e3[i_r+1][k_z]*PI*dz1*(r3*r3-r2*r2)/vol_2;

	//weighting Ez[i][k+1]//
	ez = ez + e3[i_r][k_z+1]*PI*dz2*(r2*r2-r1*r1)/vol_1;

	//weighting Ez[i+1][k+1]//
	ez = ez + e3[i_r+1][k_z+1]*PI*dz2*(r3*r3-r2*r2)/vol_2;

///////////////////////////////////////////////////////

	// weighting of E_fi//
///////////////////////////////////////////////////////
	//finding number of cell. example dz=0.5,	 x3 = 0.7, z_k =1;
	i_r = (int)ceil((x1)/geom1->dr)-1;
	k_z = (int)ceil((x3)/geom1->dz)-1;

	if(x1>dr)
	{
		vol_1 = PI*dz*dr*dr*2*i_r;
	}
	else
	{
		vol_1 = PI*dz*dr*dr/4.0; //volume of first cell
	}

	r2 = (i_r+0.5)*dr;
	vol_2 = PI*dz*dr*dr*(2*i_r+2);
	dz1 = (k_z+1)*dz-x3;
	dz2 = x3-k_z*dz;
	//////////////////////////////////////
	//weighting Efi[i][k]//
	efi = efi + e2[i_r][k_z]*PI*dz1*(r2*r2 - r1*r1)/vol_1;

	//weighting Efi[i+1][k]//
	efi = efi + e2[i_r+1][k_z]*PI*dz1*(r3*r3-r2*r2)/vol_2;

	//weighting Efi[i][k+1]//
	efi = efi + e2[i_r][k_z+1]*PI*dz2*(r2*r2-r1*r1)/vol_1;

	//weighting Efi[i+1][k+1]//
	efi =efi + e2[i_r+1][k_z+1]*PI*dz2*(r3*r3-r2*r2)/vol_2;

	Triple components(er, efi, ez);

	return components;
}
double accur =1e-12;

///////////////////////////////////////////////

///////Return one dimensional field components///////////

double* E_field::get_1d_e1()
{
	//copy 2d field array into 1d array rowwise
	for (int i = 0; i < geom1->n_grid_1 - 1; i++)
		for (int k = 0; k < geom1->n_grid_2; k++)
			e1_1d[i * geom1->n_grid_2 + k] = e1[i][k];
	return e1_1d;
}

double* E_field::get_1d_e2()
{
	//copy 2d field array into 1d array rowwise
	for (int i = 0; i < geom1->n_grid_1; i++)
		for (int k = 0; k < geom1->n_grid_2; k++)
			e2_1d[i * geom1->n_grid_2 + k] = e2[i][k];
	return e2_1d;
}

double* E_field::get_1d_e3()
{
	//copy 2d field array into 1d array rowwise
	for (int i = 0; i < geom1->n_grid_1; i++)
		for (int k = 0; k < geom1->n_grid_2 - 1; k++)
			e3_1d[i * (geom1->n_grid_2 - 1) + k] = e3[i][k];
	return e3_1d;
}
/////////////////////////////////////////////////////
