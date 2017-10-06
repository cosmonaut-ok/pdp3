#include <stdio.h>
#include <math.h>
#include "system_kern.cuh"

#define TOTAL_THREADS		65536
#define BLOCK_THREADS		256
	
__device__ Particles_struct cuda_specie;    // simulation data (on device)
__device__ flcuda* x1;
__device__ flcuda* x3;
__device__ flcuda* v1;
__device__ flcuda* v2;
__device__ flcuda* v3;
//electric and magnetic field
__device__ flcuda* e1;
__device__ flcuda* e2;
__device__ flcuda* e3;

__device__ flcuda* h1;
__device__ flcuda* h2;
__device__ flcuda* h3;


//current and charge density
__device__ flcuda* cur1;
__device__ flcuda* cur2;
__device__ flcuda* cur3;
__device__ flcuda* rho;

//is alive array
__device__ bool* is_alive;



__global__ void StepV(flcuda* x1, flcuda* x3, flcuda* v1, flcuda* v2, flcuda* v3, 
					  flcuda* e1, flcuda* e2, flcuda* e3, flcuda* h1, flcuda* h2, 
					  flcuda* h3, bool* is_alive, int number, flcuda timestep)
{
	uint i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;
	flcuda gamma, gamma_inv;
	flcuda e1_val = 0.0, e2_val, e3_val, b1_val, b2_val, b3_val, vv1, vv2, vv3;
	const flcuda mu0 = (flcuda)1.0e-6;
	flcuda const1 = cuda_specie.charge*timestep/(flcuda)2.0/cuda_specie.mass;
	flcuda const2;
	//if (t->current_time == t->start_time) const1 = const1/2.0;
		if (is_alive[i] && (i < number))
		{
	   	      e1_val = get_field_e(x1[i], x3[i], e1, 1, cuda_specie) * const1;
	          e2_val = get_field_e(x1[i], x3[i], e2, 2, cuda_specie) * const1;
	          e3_val = get_field_e(x1[i], x3[i], e3, 3, cuda_specie) * const1;



	        b1_val = get_field_h(x1[i] ,x3[i], h1, 1, cuda_specie)*mu0*const1;
	        b2_val = get_field_h(x1[i], x3[i], h2, 2, cuda_specie)*mu0*const1;
	        b3_val = get_field_h(x1[i], x3[i], h3, 3, cuda_specie)*mu0*const1;



			////1. Multiplication by relativistic factor
			////u(n-1/2) = gamma(n-1/2)*v(n-1/)
			gamma = get_gamma(i, v1, v2, v3);
			v1[i] = gamma*v1[i];
			v2[i] = gamma*v2[i];
			v3[i] = gamma*v3[i];

			//2. Half acceleration in the electric field
			//u'(n) = u(n-1/2) + q*dt/2/m*E(n)

			v1[i] = v1[i] + e1_val;
			v2[i] = v2[i] + e2_val;
			v3[i] = v3[i] + e3_val;

			//3. Rotation in the magnetic field
			//u" = u' + 2/(1+B'^2)[(u' + [u'xB'(n)])xB'(n)]
			//B'(n) = B(n)*q*dt/2/mass/gamma(n)
			gamma_inv = get_gamma_inv(i, v1, v2, v3);
  
			b1_val = b1_val/gamma_inv;
	        b2_val = b2_val/gamma_inv;
	        b3_val = b3_val/gamma_inv;
			const2 = (flcuda)2.0/((flcuda)1.0 + b1_val*b1_val + b2_val*b2_val + b3_val*b3_val);
			vv1 = v1[i];
			vv2 = v2[i];
			vv3 = v3[i];
			v1[i] = vv1 + const2*((vv2 - vv1*b3_val + vv3*b1_val)*b3_val - (vv3 + vv1*b2_val - vv2*b1_val)*b2_val);
			v2[i] = vv2 + const2*(-(vv1 + vv2*b3_val - vv3*b2_val)*b3_val + (vv3 + vv1*b2_val - vv2*b1_val)*b1_val); 
			v3[i] = vv3 + const2*((vv1 + vv2*b3_val - vv3*b2_val)*b2_val - (vv2 - vv1*b3_val + vv3*b1_val)*b1_val);

			//4. Half acceleration in the electric field
			//u(n+1/2) = u(n) + q*dt/2/m*E(n)
			v1[i] = v1[i] + e1_val;
			v2[i] = v2[i] + e2_val;
			v3[i] = v3[i] + e3_val;

			//5. Division by relativistic factor
			gamma = get_gamma_inv(i, v1, v2, v3);
		    v1[i] = v1[i]/gamma;
			v2[i] = v2[i]/gamma;
			v3[i] = v3[i]/gamma;
		}

}


__device__ flcuda get_field_e(flcuda x1, flcuda x3, flcuda* e_input, int component, Particles_struct cuda_specie)
{	
	int i_r=0;  // number of particle i cell 
	int k_z=0;  // number of particle k cell
	int counter;
	
	flcuda pi = 3.1415926535897932;
	flcuda dr = cuda_specie.dr;
	flcuda dz = cuda_specie.dz;
	flcuda r1, r2, r3; // temp variables for calculation
	flcuda dz1, dz2;   //  temp var.: width of k and k+1 cell 
	flcuda er =0;
	flcuda efi =0;
	flcuda ez =0;
	flcuda vol_1 =0; //  volume of i cell; Q/V, V - volume of elementary cell 
	flcuda vol_2 =0; //  volume of i+1 cell;
	
////////////////////////
	r1 = x1-0.5*dr;
	r3 = x1+0.5*dr;
///////////////////////

		switch (component)
	{
      case 1:
      {

    // weighting of E_r//
///////////////////////////////////////////////////
	//finding number of cell. example dr=0.5,  x1 = 0.7, i_r =0;!!
	 i_r = (int)ceil((x1-0.5*dr)/dr)-1;
	 k_z = (int)ceil((x3)/dz)-1;
	
	 vol_1 = pi*dz*dr*dr*(2*i_r+1);
	 vol_2 = pi*dz*dr*dr*(2*i_r+3);
	 dz1 = (k_z+1)*dz-x3;
	 dz2 = x3 - k_z*dz;
	 r2 = (i_r+1)*dr;
    ///////////////////////////////////////
  
    //weighting Er[i][k]//
   counter = get_linear_coord(i_r, k_z, cuda_specie.grid_num3, 1);
   er = e_input[counter];
   er = er*(pi*dz1*(r2*r2-r1*r1))/vol_1;

	//weighting Er[i+1][k]//
   er = er + e_input[get_linear_coord(i_r+1, k_z, cuda_specie.grid_num3, 1)]*(pi*dz1*(r3*r3-r2*r2))/vol_2;

   //weighting Er[i][k+1]//
   er= er + e_input[get_linear_coord(i_r, k_z+1, cuda_specie.grid_num3, 1)]*(pi*dz2*(r2*r2-r1*r1))/vol_1;

   //weighting Er[i+1][k+1]//
   er = er + e_input[get_linear_coord(i_r+1, k_z+1, cuda_specie.grid_num3, 1)]*(pi*dz2*(r3*r3-r2*r2))/vol_2;
   
///////////////////////////////////////////////////////
   return er;
   }
   case 3:
   {
 
	
	     // weighting of E_z//
///////////////////////////////////////////////////////
//finding number of cell. example dz=0.5,  x3 = 0.7, z_k =0;!!
	i_r = (int)ceil((x1)/dr)-1;
	k_z = (int)ceil((x3-0.5*dz)/dz)-1;

///////////////////////////////////

   if(x1>dr)
	{
		vol_1 = pi*dz*dr*dr*2.0*(flcuda)i_r;
    }
   else
   {
	   vol_1 = pi*dz*dr*dr/4.0; //volume of first cell
   }
		  r2 = (i_r+0.5)*dr;
		  vol_2 = pi*dz*dr*dr*(2*i_r+2);
		  dz1 = (k_z+1.5)*dz - x3;
		  dz2 = x3 - (k_z+0.5)*dz;
		  //////////////////////////////////////

		   //weighting Ez[i][k]//
		   ez = ez + e_input[get_linear_coord(i_r, k_z, cuda_specie.grid_num3, 3)]*(pi*dz1*(r2*r2-r1*r1))/vol_1;

		  //weighting Ez[i+1][k]//
		   ez = ez + e_input[get_linear_coord(i_r+1, k_z, cuda_specie.grid_num3, 3)]*pi*dz1*(r3*r3-r2*r2)/vol_2;   

          //weighting Ez[i][k+1]//
		   ez = ez + e_input[get_linear_coord(i_r, k_z+1, cuda_specie.grid_num3, 3)]*pi*dz2*(r2*r2-r1*r1)/vol_1;
   
         //weighting Ez[i+1][k+1]//
		   ez = ez + e_input[get_linear_coord(i_r+1, k_z+1, cuda_specie.grid_num3, 3)]*pi*dz2*(r3*r3-r2*r2)/vol_2;    

	    return ez;
   }
   case 2:
	{

///////////////////////////////////////////////////////

	 // weighting of E_fi//
///////////////////////////////////////////////////////
 //finding number of cell. example dz=0.5,  x3 = 0.7, z_k =1;
	 i_r = (int)ceil((x1)/dr)-1;
     k_z = (int)ceil((x3)/dz)-1;
	
  if(x1>dr)
	{
		vol_1 = pi*dz*dr*dr*2.0*(flcuda)i_r;
    }
  else
  {
	 vol_1 = pi*dz*dr*dr/4.0; //volume of first cell
  }

		  r2 = (i_r+0.5)*dr;
		  vol_2 = pi*dz*dr*dr*(2*i_r+2);
		  dz1 = (k_z+1)*dz-x3;
		  dz2 = x3-k_z*dz;
		  //////////////////////////////////////
		  //weighting Efi[i][k]//
		  efi = efi + e_input[get_linear_coord(i_r, k_z, cuda_specie.grid_num3, 2)]*pi*dz1*(r2*r2 - r1*r1)/vol_1;

		  //weighting Efi[i+1][k]//
		   efi = efi + e_input[get_linear_coord(i_r+1, k_z, cuda_specie.grid_num3, 2)]*pi*dz1*(r3*r3-r2*r2)/vol_2;

          //weighting Efi[i][k+1]//
		   efi = efi + e_input[get_linear_coord(i_r, k_z+1, cuda_specie.grid_num3, 2)]*pi*dz2*(r2*r2-r1*r1)/vol_1;
   
         //weighting Efi[i+1][k+1]//
		   efi =efi + e_input[get_linear_coord(i_r+1, k_z+1, cuda_specie.grid_num3, 2)]*pi*dz2*(r3*r3-r2*r2)/vol_2;
		  return efi;
	}
	}
  		 
	return 0.0;
}

__device__ flcuda get_field_h(flcuda x1, flcuda x3, flcuda* h_input, int component, Particles_struct cuda_specie)
{
  int i_r = 0, k_z = 0;  
  flcuda r1, r2, r3, dz1, dz2, hr = 0.0, hfi = 0.0, hz = 0.0, 
	     vol_1 = 0.0, vol_2 = 0.0; //  volumes of i and i+1 cell; Q/V, V - volume of elementary cell 
	
  flcuda pi = 3.1415926535897932;
  flcuda dr = cuda_specie.dr;
  flcuda dz = cuda_specie.dz;

  ////////////////////////
  r1 = x1-0.5*dr;
  r3 = x1+0.5*dr;
  ///////////////////////

  switch (component)
  {
    case 3:
    {
      // weighting of H_z//
      i_r = (int)ceil((x1-0.5*dr)/dr)-1;
      k_z = (int)ceil((x3)/dz)-1;
      vol_1 = pi*dz*dr*dr*(2*i_r+1);
      vol_2 = pi*dz*dr*dr*(2*i_r+3);
      dz1 = (k_z+1)*dz-x3;
      dz2 = x3 - k_z*dz;
      r2 = (i_r+1)*dr;
      ///////////////////////////////////////  
      //weighting Hz[i][k]//
      hz = hz + h_input[get_linear_coord(i_r, k_z, cuda_specie.grid_num3, 6)]*(pi*dz1*(r2*r2-r1*r1))/vol_1;

      //weighting Hz[i+1][k]//
      hz = hz + h_input[get_linear_coord(i_r+1, k_z, cuda_specie.grid_num3, 6)]*(pi*dz1*(r3*r3-r2*r2))/vol_2;

      //weighting Hz[i][k+1]//
      hz = hz + h_input[get_linear_coord(i_r, k_z+1, cuda_specie.grid_num3, 6)]*(pi*dz2*(r3*r3-r2*r2))/vol_1;

      //weighting Hz[i+1][k+1]//
      hz = hz + h_input[get_linear_coord(i_r+1, k_z+1, cuda_specie.grid_num3, 6)]*(pi*dz2*(r3*r3-r2*r2))/vol_2;
      return hz;
    }
	case 1:
    {
      // weighting of Hr//
      i_r = (int)ceil((x1)/dr)-1;
      k_z = (int)ceil((x3-0.5*dz)/dz)-1;
      if(x1>dr)
      {
        vol_1 = pi*dz*dr*dr*2*i_r;
      }
      else 
      {
        vol_1 = pi*dz*dr*dr/4.0; //volume of first cell
      }
      r2 = (i_r+0.5)*dr;
      vol_2 = pi*dz*dr*dr*(2*i_r+2);
      dz1 = (k_z+1.5)*dz - x3;
      dz2 = x3 - (k_z+0.5)*dz;
      //////////////////////////////////////

      //weighting Hr[i][k]//
      hr = hr + h_input[get_linear_coord(i_r, k_z, cuda_specie.grid_num3, 4)]*(pi*dz1*(r2*r2-r1*r1))/vol_1;

      //weighting Hr[i+1][k]//
      hr = hr + h_input[get_linear_coord(i_r+1, k_z, cuda_specie.grid_num3, 4)]*pi*dz1*(r3*r3-r2*r2)/vol_2;   

      //weighting Hr[i][k+1]//
      hr = hr + h_input[get_linear_coord(i_r, k_z+1, cuda_specie.grid_num3, 4)]*pi*dz2*(r2*r2-r1*r1)/vol_1;
   
      //weighting Hr[i+1][k+1]//
      hr = hr + h_input[get_linear_coord(i_r+1, k_z+1, cuda_specie.grid_num3, 4)]*pi*dz2*(r3*r3-r2*r2)/vol_2; 
      return hr;
	}
	case 2:
    {
      // weighting of H_fi//
      i_r = (int)ceil((x1-0.5*dr)/dr)-1;
      k_z = (int)ceil((x3-0.5*dz)/dz)-1;
      r2 = (i_r+1)*dr;
      vol_1 = pi*dz*dr*dr*(2*i_r+1);
      vol_2 = pi*dz*dr*dr*(2*i_r+3);
      dz1 = (k_z+1.5)*dz-x3;
      dz2 = x3-(k_z+0.5)*dz;
      //weighting Hfi[i][k]//
      hfi = hfi + h_input[get_linear_coord(i_r, k_z, cuda_specie.grid_num3, 5)]*pi*dz1*(r2*r2-r1*r1)/vol_1;
	  //hfi = h_input[get_linear_coord(i_r, k_z, cuda_specie.grid_num3, 5)];
	  //hfi = h_input[1282];
	  //hfi = i_r;

      //weighting Hfi[i+1][k]//
      hfi = hfi + h_input[get_linear_coord(i_r+1, k_z, cuda_specie.grid_num3, 5)]*pi*dz1*(r3*r3-r2*r2)/vol_2;
		   
      //weighting Hfi[i][k+1]//
      hfi = hfi + h_input[get_linear_coord(i_r, k_z+1, cuda_specie.grid_num3, 5)]*dz2*pi*(r2*r2-r1*r1)/vol_1;
   
      //weighting Hfi[i+1][k+1]//
      hfi = hfi + h_input[get_linear_coord(i_r+1, k_z+1, cuda_specie.grid_num3, 5)]*pi*dz2*(r3*r3-r2*r2)/vol_2;  	
      return hfi;
    }
  }
  return 0.0;
}

__device__ int get_linear_coord(int index_r, int index_z, int ngrid_z, int component)
{
	//index components:
	// 1 - er, 2 - e_phi, 3 - e_z
	// 4 - hr, 5 - h_phi, 6 - h_z

	switch (component)
	{
		case 1:
		case 2:
		case 6:
			return (index_r * ngrid_z + index_z);
		case 3:
		case 4:
		case 5:
			return (index_r * (ngrid_z - 1) + index_z);
	}
					 
  return 0;
}

__device__ flcuda get_gamma(int i, flcuda* v1, flcuda* v2, flcuda* v3)
{
	return pow((flcuda)1.0 - (v1[i]*v1[i] + v2[i]*v2[i] + v3[i]*v3[i])/(flcuda)300000000.0/(flcuda)300000000.0,(flcuda)-0.5);
}

__device__ flcuda get_gamma_inv(int i, flcuda* v1, flcuda* v2, flcuda* v3)
{
	return pow((v1[i]*v1[i] + v2[i]*v2[i] + v3[i]*v3[i])/(flcuda)300000000.0/(flcuda)300000000.0 + (flcuda)1.0, (flcuda)0.5);	
}

