#include "boundary_neumann.h"
#include "Poisson.h"
#include "Poisson_neumann.h"

boundary_neumann::boundary_neumann(E_field* ef_t,charge_density* rho_t,Geometry *geom_t): e_f(ef_t),rho(rho_t),cyl_geom(geom_t)
{
}
boundary_neumann::boundary_neumann(void)
{
}

boundary_neumann::~boundary_neumann(void)
{
}
void boundary_neumann:: specify_boundary_conditions(flcuda E_fi_upper, flcuda E_fi_left, flcuda E_fi_right, flcuda fi_upper_wall)
{
int i=0, k=0;
int n_grid1=cyl_geom->n_grid_1;
int n_grid2 = cyl_geom->n_grid_2;
// setazimuthal component electric field initial value
/////////////////////////////////////////////
for (i=0;i<(n_grid1);i++)
{
	e_f->e2[i][0]=E_fi_left;
	e_f->e2[i][n_grid2-1]=E_fi_right;
}
	for(k=0;k<n_grid2;k++)
	{
		e_f->e2[n_grid1-1][k]=E_fi_upper;
	}
/////////////////////////////////////////////

////////////////////////////////////////////
	//set potential initial value in z wall
	for(k=0;k<n_grid2;k++)
	{
		e_f->fi[n_grid1-1][k]=fi_upper_wall;
	}
///////////////////////////////////////////
	// calculate poisson equation with Neumann boundary conditions
	
	Poisson* solve = new Poisson_neumann(cyl_geom);
	solve->poisson_solve(e_f,rho);
	solve->test_poisson_equation(e_f,rho);
/////////////////////////////////////////
}
