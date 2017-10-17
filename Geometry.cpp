#include "Geometry.h"
#include "PML.h"

// class constructor with sigma //
/////////////////////////////////////////////
Geometry::Geometry(double fs,
									 double ss,
									 int ng1,
									 int ng2,
									 PML* pml1_t): pml1(pml1_t)//, geom1(geom1_t)
{
	first_size=fs;
	second_size=ss;
	n_grid_1=ng1;
	n_grid_2=ng2;
	dr=set_dr();
	dz=set_dz();
	//////////////////////////////////////
	epsilon = new double*[n_grid_1];
	for (int i=0; i<(n_grid_1);i++)
	{
		epsilon[i]= new double[n_grid_2];
	}

	sigma = new double*[n_grid_1];
	for (int i=0; i<(n_grid_1);i++)
	{
		sigma[i]= new double[n_grid_2];
	}

	for (int i=0; i<(n_grid_1);i++)
	{
		for (int k=0; k<(n_grid_2);k++)
			sigma[i][k]=0;
	}
//////////////////////////////////

	// sigma assigning//
///////////////////////////////
	pml1->calc_sigma(this);
//////////////////////////////
}
///////////////////////////////////////////

////////////////////////////////////
Geometry::Geometry(double* param, PML* pml1_t): pml1(pml1_t)//, geom1(geom1_t)
{
	first_size=param[0];
	second_size=param[1];
	n_grid_1=param[2];
	n_grid_2=param[3];
	dr=set_dr();
	dz=set_dz();
	//////////////////////////////////////
	epsilon = new double*[n_grid_1];
	for (int i=0; i<(n_grid_1);i++)
	{
		epsilon[i]= new double[n_grid_2];
	}

	sigma = new double*[n_grid_1];
	for (int i=0; i<(n_grid_1);i++)
	{
		sigma[i]= new double[n_grid_2];
	}

	for (int i=0; i<(n_grid_1);i++)
	{
		for (int k=0; k<(n_grid_2);k++)
			sigma[i][k]=0;
	}
//////////////////////////////////

	// sigma assigning//
///////////////////////////////
	pml1->calc_sigma(this);
//////////////////////////////
}

// class constructor: sigma =0 //
/////////////////////////////////////////////
Geometry::Geometry(double fs, double ss,	int ng1, int ng2)
{
	first_size=fs;
	second_size=ss;
	n_grid_1=ng1;
	n_grid_2=ng2;
	dr=set_dr();
	dz=set_dz();
	pml1=0;
	//////////////////////////////////////
	epsilon = new double*[n_grid_1];
	for (int i=0; i<(n_grid_1);i++)
	{
		epsilon[i]= new double[n_grid_2];
	}

	sigma = new double*[n_grid_1];
	for (int i=0; i<(n_grid_1);i++)
	{
		sigma[i]= new double[n_grid_2];
	}

	for (int i=0; i<(n_grid_1);i++)
	{
		for (int k=0; k<(n_grid_2);k++)
			sigma[i][k]=0;
	}
}
///////////////////////////////////////////


///////////////////////////////////////////

Geometry::Geometry()
	{
	}
///////////////////////////////////////////
///////////////////////////////////////////
void Geometry::set_epsilon()
{
	for(int i=0;i<(n_grid_1);i++)
		for(int k=0;k<(n_grid_2);k++)
		{
			epsilon[i][k]=1;
		}
}
/////////////////////////////////////////////////

Geometry::~Geometry(void)
{
}

double Geometry::set_dr()
{
	return first_size/(n_grid_1-1);
}

double Geometry::set_dz()
{
	return second_size/(n_grid_2-1);
}
