#include "current.h"

current::current(void)
{
}

current::current(Geometry* geom1_t): geom1(geom1_t)
{
	//jr//
	///////////////////////////////////
	//n_grid - number of edges
	j1 = new double*[geom1->n_grid_1-1];
	for (int i=0; i<(geom1->n_grid_1-1);i++)
	{
		j1[i]= new double[geom1->n_grid_2];
	}
	///////////////////////////////////
	//jfi//
	j2 = new double*[geom1->n_grid_1];

	for (int i=0; i<(geom1->n_grid_1);i++)
	{
		j2[i]= new double[geom1->n_grid_2];
	}
	//////////////////////////////////////
	//jz//
	//////////////////////////////////////
	j3 = new double*[geom1->n_grid_1];

	for (int i=0; i<(geom1->n_grid_1);i++)
	{
		j3[i]= new double[geom1->n_grid_2-1];
	}
	///////////////////////////////////////

	j1_1d = new double[(geom1->n_grid_1-1)*geom1->n_grid_2];
	///////////////////////////////////

	j2_1d = new double[geom1->n_grid_1*geom1->n_grid_2];
	//////////////////////////////////////

	j3_1d = new double[geom1->n_grid_1*(geom1->n_grid_2-1)];
	//////////////////////////////////////
	//initialization//

	for (int i=0; i<(geom1->n_grid_1-1);i++)
		for (int k=0; k<(geom1->n_grid_2-1);k++)
		{
			j1[i][k]=0;
			j2[i][k]=0;
			j3[i][k]=0;
		}
	for (int i=0; i<(geom1->n_grid_1-1);i++)
	{
		j1[i][geom1->n_grid_2-1]=0;
	}

	for(int k=0;k<(geom1->n_grid_2);k++)
	{
		j2[geom1->n_grid_1-1][k]=0;
	}
	for(int i=0;i<(geom1->n_grid_1);i++)
	{
		j2[i][geom1->n_grid_2-1]=0;
	}
	for(int k=0;k<(geom1->n_grid_2-1);k++)
	{
		j3[geom1->n_grid_1-1][k]=0;
	};


}

// Destructor
current::~current(void)
{
	for (int i=0; i<(geom1->n_grid_1-1);i++)
		delete[]j1[i];
	delete[]j1;

	for (int i=0; i<(geom1->n_grid_1);i++)
		delete[]j2[i];
	delete[]j2;

	for (int i=0; i<(geom1->n_grid_1-1);i++)
		delete[]j3[i];
	delete[]j3;

	delete j1_1d;
	delete j2_1d;
	delete j3_1d;
}

//////////////////////////////////////////

//////////////////////////////////////////

/////////////////////////////////////////
//functions for getting access to j arrays//
/////////////////////////////////////////

double** current::get_j1() const
{
	return this->j1;
}

double** current::get_j2() const
{
	return this->j2;
}

double** current::get_j3() const
{
	return this->j3;
}
////////////////////////////////////////////

///////Return one dimensional field components///////////

double* current::get_j1_1d() const
{
	//copy 2d field array into 1d array rowwise
	for (int i = 0; i < geom1->n_grid_1 - 1; i++)
		for (int k = 0; k < geom1->n_grid_2; k++)
			j1_1d[i * geom1->n_grid_2 + k] = j1[i][k];
	return j1_1d;
}

double* current::get_j2_1d() const
{
	//copy 2d field array into 1d array rowwise
	for (int i = 0; i < geom1->n_grid_1; i++)
		for (int k = 0; k < geom1->n_grid_2; k++)
			j2_1d[i * geom1->n_grid_2 + k] = j2[i][k];
	return j2_1d;
}

double* current::get_j3_1d() const
{
	//copy 2d field array into 1d array rowwise
	for (int i = 0; i < geom1->n_grid_1; i++)
		for (int k = 0; k < geom1->n_grid_2 - 1; k++)
			j3_1d[i * (geom1->n_grid_2 - 1) + k] = j3[i][k];
	return j3_1d;
}


void current::j1_add_1d(double *j1_1d)
{
	//copy 2d field array into 1d array rowwise
	for (int i = 0; i < geom1->n_grid_1 - 1; i++)
		for (int k = 0; k < geom1->n_grid_2; k++)
			j1[i][k] += j1_1d[i * geom1->n_grid_2 + k];
	return;
}

void current::j2_add_1d(double *j2_1d)
{
	//copy 2d field array into 1d array rowwise
	for (int i = 0; i < geom1->n_grid_1; i++)
		for (int k = 0; k < geom1->n_grid_2; k++)
			j2[i][k] += j2_1d[i * geom1->n_grid_2 + k];
	return;
}

void current::j3_add_1d(double *j3_1d)
{
	//copy 2d field array into 1d array rowwise
	for (int i = 0; i < geom1->n_grid_1; i++)
		for (int k = 0; k < geom1->n_grid_2 - 1; k++)
			j3[i][k] += j3_1d[i * (geom1->n_grid_2 - 1) + k];
	return;
}

//////////////////////////////////////////
// functions for changing values of j//
void current::set_j1(int i, int k, double value)
{
	j1[i][k]= j1[i][k]+value;
}

void current::set_j2(int i, int k, double value)
{
	j2[i][k]= j2[i][k]+value;
}

void current::set_j3(int i, int k, double value)
{
	j3[i][k]=j3[i][k]+value;
}

void current::reset_j()
{
	int i=0;
	int k=0;
	for (i=0; i<(geom1->n_grid_1-1);i++)
		for (k=0; k<(geom1->n_grid_2-1);k++)
		{
			j1[i][k]=0;
			j3[i][k]=0;
		}

	for (i=0;i<geom1->n_grid_1;i++)
		for (k=0;k<geom1->n_grid_2;k++)
			j2[i][k] = 0.0;

	for (i=0; i<(geom1->n_grid_1-1);i++)
		j1[i][geom1->n_grid_2-1]=0;


	for(k=0;k<(geom1->n_grid_2-1);k++)
		j3[geom1->n_grid_1-1][k]=0;
}
