#include "Boundary_Maxwell_conditions.h"
#include <Constant.h>

const double PI = constant::PI;
const double LIGHT_SPEED = constant::LIGHT_SPEED;

Boundary_Maxwell_conditions::Boundary_Maxwell_conditions(E_field* e_fld_t):e_fld(e_fld_t)
{
}

Boundary_Maxwell_conditions::Boundary_Maxwell_conditions(void)
{
}

Boundary_Maxwell_conditions::~Boundary_Maxwell_conditions(void)
{
}

void Boundary_Maxwell_conditions::specify_initial_field(Geometry* cyl_geom,
                                                        double E_fi_upper,
                                                        double E_fi_left,
                                                        double E_fi_right)
{
	int i=0, k=0;
	int n_grid1=cyl_geom->n_grid_1;
	int n_grid2 = cyl_geom->n_grid_2;
// setazimuthal component electric field initial value
/////////////////////////////////////////////
	for (i=0;i<(n_grid1);i++)
	{
		e_fld->e2[i][0]=E_fi_left;
		e_fld->e2[i][n_grid2-1]=E_fi_right;
	}
	for(k=0;k<n_grid2;k++)
	{
		e_fld->e2[n_grid1-1][k]=E_fi_upper;
	}
/////////////////////////////////////////////

}
void Boundary_Maxwell_conditions::radiation_source(Geometry* cyl_geom,
                                                   double region_size,
                                                   double frequency,
                                                   int wave_type,
                                                   double time)
{
	int i_max =0, i=0;
	i_max = (int) (region_size/cyl_geom->dr);

	switch (wave_type)
	{
		// E electromagnetic wave
	case 0:
	{
		for (i=0;i<i_max;i++)
			e_fld->e1[i][0]=5e5*sin(2*PI*frequency*time);
	}
	break;
	// H electromagnetic wave
	case 1:
	{
		for (i=0;i<i_max;i++)
			e_fld->e2[i][0]=e_fld->e2[i][0]*cos(2*PI*frequency*time);
	}
	break;
	}
}

void Boundary_Maxwell_conditions::probe_mode_exitation(Geometry* cyl_geom,
                                                        current* j_input,
                                                        double n,
                                                        double frequency,
                                                        double time)
{
	int k_max =0, k=0;
	int k_start =0;
	double probe_lenght = LIGHT_SPEED*n/frequency; // TODO: light velicity?
	k_max = (int) (probe_lenght/cyl_geom->dz)+k_start;
	double dz = cyl_geom->dz;
	double **j_p;
	j_p = j_input->get_j3();
	/*if (cyl_geom->dz*k<probe_lenght)
		{
		k_max=3*10e8*time/dz;*/
	for (k=k_start;k<k_max;k++)
	{
		j_p[0][k]=6e7*sin(2*PI*frequency*time)*sin(2.0*n*PI*dz*(k-k_start)/(probe_lenght));
		//e_fld->e1[0][k]=1e3*sin(2*PI*frequency*time)*sin(PI*dz*k/(0.5*3e8/frequency));
	}
	/*}*/

}
