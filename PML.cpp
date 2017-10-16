#include "PML.h"
#include "math.h"

PML::PML(void)
{
}

PML::PML(double comp_l_1, double comp_l_2, double comp_l_3,
				 double sigma1_t, double sigma2_t)
{
	comparative_l_1 = comp_l_1; // lenhgt sigma on left wall
	comparative_l_2 = comp_l_2; //  lenhgt sigma on right_wall
	comparative_l_3 = comp_l_3; // lenght sigma on external wall
	sigma1 = sigma1_t;
	sigma2 = sigma2_t;

}

PML::~PML(void)
{
}
// sigma assignment
void PML::calc_sigma(Geometry* geom1)
{
  // defining lenght of sigma calculation ion left wall
	double lenght_sigma_left =geom1->dz * (floor(geom1->n_grid_2*comparative_l_1));

	// defining lenght of sigma calculation on right wall
	double lenght_sigma_right =geom1->dz * (floor(geom1->n_grid_2*comparative_l_2));

	// defining lenght of sigma calculation on right wall
	double lenght_sigma_extern =geom1->dr * (floor(geom1->n_grid_1*comparative_l_3));

	// if pml is only on z walll
	if ((comparative_l_1==0.0)&&(comparative_l_2==0.0))
	{
		for(int i=0;i<geom1->n_grid_1;i++)
			for(int k=0;k<geom1->n_grid_2;k++)
				if ((geom1->first_size - geom1->dr*(i)) <= lenght_sigma_extern)
					geom1->sigma[i][k] = sigma1 + (sigma2-sigma1)/(pow(lenght_sigma_extern,2))*pow((geom1->dr*(i+1)-geom1->first_size+lenght_sigma_extern),2);
	}
	else
	{
		// sigma on r=0 and r=r walls
		for(int k=0;k<geom1->n_grid_2;k++)
			for(int i=0;i<geom1->n_grid_1;i++)
			{
				if (geom1->dz*k<lenght_sigma_left)
				{
					geom1->sigma[i][k] = sigma1 + (sigma2-sigma1)/(pow(lenght_sigma_left,2))*pow((lenght_sigma_left-geom1->dz*k),2);
				}
				if ((geom1->second_size - geom1->dz*(k))<lenght_sigma_right)
				{
					geom1->sigma[i][k] = sigma1 + (sigma2-sigma1)/(pow(lenght_sigma_right,2))*pow((geom1->dz*(k+1)-geom1->second_size+lenght_sigma_right),2);
				}
			}
	}
/////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////

	if(comparative_l_3==0)
	{
	}
	else
	{
    // sigma assigning//
		//sigma on z=z wall//
//////////////////////////////////////////////
		for(int i=0;i<geom1->n_grid_1;i++)
			for(int k=0;k<geom1->n_grid_2;k++)
				if ((geom1->first_size - geom1->dr*(i)) <= lenght_sigma_extern)
					if (((geom1->first_size - geom1->dr*(i+1)) < lenght_sigma_extern * geom1->dz * (k)/lenght_sigma_left)  &&  ( (geom1->first_size - geom1->dr*(i)) <= (lenght_sigma_extern/lenght_sigma_right*(geom1->second_size - geom1->dz*(k)))))
						geom1->sigma[i][k] = sigma1 + (sigma2-sigma1)/(pow(lenght_sigma_extern,2))*pow((geom1->dr*(i+1)-geom1->first_size+lenght_sigma_extern),2);
/////////////////////////////////////////////
	}
	return;
}
