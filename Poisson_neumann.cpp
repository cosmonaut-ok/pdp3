#include "Poisson_neumann.h"


Poisson_neumann::Poisson_neumann(Geometry* cyl_geom):Poisson(cyl_geom)
{

	//charge_density//
///////////////////////////////////////
   t_charge_density = new flcuda*[cyl_geom->n_grid_1];
	   for (int i=0; i<(cyl_geom->n_grid_1);i++)
		 {
			t_charge_density[i] = new flcuda[cyl_geom->n_grid_2];
		 }
}

Poisson_neumann::~Poisson_neumann(void)
{
	for (int i=0; i<(cyl_geom->n_grid_1);i++)
	    delete[]t_charge_density[i];
    delete[]t_charge_density;
}
void Poisson_neumann::poisson_solve(E_field* input_e, charge_density* ro1)
{	
	const flcuda pi = 3.14159265358979;
	int i=0;
	int k=0;
	flcuda a=0;
	flcuda c=0;
	flcuda b=0;
	flcuda d=0;
	flcuda* alpha = new flcuda [cyl_geom->n_grid_1];
	flcuda* beta = new flcuda [cyl_geom->n_grid_1];
	Fourier* four1=0;
	flcuda** ro=ro1->get_ro();
	flcuda** e1 = input_e->e1;
	flcuda** e3 = input_e->e3;
	flcuda** fi = input_e->fi;
	flcuda dr = cyl_geom->dr;
	flcuda dz = cyl_geom->dz;

///////////////////////////////////////////////
		//copy charge_density array in to temp array//
for( i=0;i<(cyl_geom->n_grid_1);i++)
	for( k=0;k<(cyl_geom->n_grid_2);k++)
	{
		t_charge_density[i][k]= ro[i][k];
	}
	

	//call function for cosine transform//
////////////////////////////////////////////////
////////////////////////////////////////////////

	int temp=cyl_geom->n_grid_2;
	
		for (i=0;i<cyl_geom->n_grid_1;i++)
		{
			four1->fast_cosine_transform(t_charge_density, temp, i, false);
		}

		//sweep method//
//////////////////////////////////////////////////////////////////

	for( k=0;k<(cyl_geom->n_grid_2);k++)
	{
			//b=2.0+ pow((geom1->dr*pi*k/(geom1->dz*geom1->n_grid_2)),2);
	    	b = 2.0 - 2.0*(cos(pi*k/(cyl_geom->n_grid_2-1)) - 1)*dr*dr/(dz*dz);
			d=cyl_geom->dr*cyl_geom->dr*t_charge_density[0][k]/epsilon0;
			alpha[1]=4.0/(2.0+b);
			beta[1]=d/(2.0+b);


			for ( i=1;i<(cyl_geom->n_grid_1-1);i++)
			  {	
			//ay-1 - by + cy+1= = -d//
				
				a = 1.0-1.0/(2.0*(i));
				c = 1.0+1.0/(2.0*(i));
				d = cyl_geom->dr*cyl_geom->dr*t_charge_density[i][k]/epsilon0;

				alpha[i+1] = c/(b-alpha[i]*a);
				beta[i+1]  = (d+beta[i]*a)/(b-alpha[i]*a);

			  }
		   
				//fi[geom1->n_grid_1-2][k]= -d-c*fi[geom1->n_grid_1-1][k]-a*beta[geom1->n_grid_1-2]/(a*alpha[geom1->n_grid_1-2]-b);

				for(i=(cyl_geom->n_grid_1-2);i>=0;i--)
				 {
					fi[i][k]=beta[i+1]+alpha[i+1]*fi[i+1][k];
				 }
	}

	//call function for inverse cosine transform//
////////////////////////////////////////////////////////////
	for (i=0;i<cyl_geom->n_grid_1;i++)
		{
			int temp=cyl_geom->n_grid_2;
			four1->fast_cosine_transform(fi, temp,i, true);

		}
////////////////////////////////////////////////////////////


	//calculate electric field//
///////////////////////////////////////////////////////////
	for (i=0;i<(cyl_geom->n_grid_1-1);i++)
		for (k=0;k<(cyl_geom->n_grid_2);k++)
			{
				e1[i][k]=(fi[i][k]-fi[i+1][k])/cyl_geom->dr;
			}

	for (i=0;i<(cyl_geom->n_grid_1);i++)
		for (k=0;k<(cyl_geom->n_grid_2-1);k++)
			{
				e3[i][k]=(fi[i][k]-fi[i][k+1])/cyl_geom->dz;
			}
////////////////////////////////////////////////////////////

 delete [] alpha;
 delete [] beta;
}
/////////////////////////////////////////////////////////////
