#define pi        3.14159265358979

void strict_motion_weighting(double dt, __global double *j1, __global double *j2, __global double *j3, 
                             double x1_new,double x3_new, double x1_old, double x3_old, 
                             __global const double *params, int gid);

void simple_j_weighting(double dt, __global double *j1, __global double *j2, __global double *j3, double x1_new,double x3_new, 
                        double x1_old, double x3_old, int i_n, int k_n, __global const double *params, int gid);

int get_linear_coord(int index_r, int index_z, int ngrid_z, int component)
{
    //index components:
    // 1 - er, 2 - e_phi, 3 - e_z
    // 4 - hr, 5 - h_phi, 6 - h_z

    switch (component)
    {
        case 1:
        case 2:
        case 6:
        case 7:
        case 8:
		{
		    //printf("%d %d %d %d ", index_r, index_z, ngrid_z, index_r * ngrid_z + index_z);
            return (index_r * ngrid_z + index_z);
		}
        case 3:
        case 4:
        case 5:
        case 9:
		{
		    //printf("%d %d %d %d ", index_r, index_z, ngrid_z, index_r * (ngrid_z-1) + index_z);
            return (index_r * (ngrid_z - 1) + index_z);
		}
    }
                     
  return 0;
}

__kernel void j_weighting_kernel(__global double *j1, 
                                 __global double *j2, 
                                 __global double *j3, 
                                 __global double *x1, 
                                 __global double *x3, 
                                 __global double *x1_o, 
                                 __global double *x3_o, 
                                 __global int   *is_alive, 
                                 __global const double *params,
                                             int  number)
{

    double size1 = params[0];
    double size3 = params[1];
    double dr = params[2];
    double dz = params[3];
	double dt = params[4];

    int i, i_n, k_n, i_o, k_o, res_cell;
    //int number = cuda_specie.number;
    double x1_old, x3_old, a, r_boundary, z_boundary, delta_r, delta_z, r1, r2, z1, z2, delta_z1, delta_r2;



    i = get_global_id(0);

    if (i >= number)
        return;
    
    if ((is_alive[i] > 0) && (i < number))
    {

	    //if (i < number)
		//    printf("%d ", i);
        x1_old= x1_o[i];
        x3_old = x3_o[i];
        //finding number new and old cells
        i_n = (int)ceil((x1[i]) /dr) - 1;
        k_n = (int)ceil((x3[i]) /dz) - 1;
        i_o = (int)ceil((x1_old)/dr) - 1;
        k_o = (int)ceil((x3_old)/dz) - 1;
        //j2[get_linear_coord(1, 0, cuda_specie.grid_num3, 8)] = x1[i];
        //j2[get_linear_coord(1, 1, cuda_specie.grid_num3, 8)] = x3[i];
        //j2[get_linear_coord(1, 2, cuda_specie.grid_num3, 8)] = x1_old;
        //j2[get_linear_coord(1, 3, cuda_specie.grid_num3, 8)] = x3_old;
        if (x1_old==(i_o+1)*dr)
            i_o=i_n;
        if(x3_old==(k_o+1)*dz)
            k_o=k_n;
        if (x1[i]==(i_n+1)*dr)
            i_n=i_o;
        if(x3[i]==(k_n+1)*dz)
            k_n=k_o;
        res_cell = abs(i_n-i_o) + abs(k_n-k_o); 
        if ((x1[i]==x1_old)||(x3[i]==x3_old))
        {
            strict_motion_weighting(dt, j1, j2, j3, x1[i],x3[i],x1_old,x3_old, params, i);
			//printf("%d strict motion ", i);
			//printf("%d %e ", i, j3[0]);
        }
        else
        {
            switch (res_cell)
            {
                // 1) charge in four cells
                case 0: 
				{
				    simple_j_weighting(dt, j1, j2, j3, x1[i],x3[i] ,x1_old,x3_old, i_n, k_n, params, i);
					//printf("%d simple weighting, case 0 ", i);
				}
                break;

                // 2) charge in seven cells 
                case 1:
                {
                    // charge in seven cells (i_new != i_old)
                    if ((i_n!=i_o)&&(k_n==k_o))
                    {
                        if (x1_old >(i_n+1)*dr)
                        {
							a = (x1_old - x1[i])/(x3_old - x3[i]);
							r_boundary = (i_n+1)*dr;
							delta_r = r_boundary - x1[i];
							z_boundary = x3[i] + delta_r/a;

							simple_j_weighting(dt, j1, j2, j3, r_boundary,z_boundary ,x1_old,x3_old, i_n+1,k_n, params, i);
							simple_j_weighting(dt, j1, j2, j3, x1[i],x3[i], r_boundary,z_boundary, i_n, k_n, params, i);
                    }
                    else 
                    {
                        a = (x1[i] - x1_old)/(x3[i] - x3_old);
                        r_boundary = (i_n)*dr;
                        delta_r = r_boundary - x1_old;
                        z_boundary = x3_old + delta_r/a;

                        simple_j_weighting(dt, j1, j2, j3, r_boundary,z_boundary ,x1_old,x3_old, i_n-1, k_n, params, i);
                        simple_j_weighting(dt, j1, j2, j3, x1[i],x3[i], r_boundary,z_boundary, i_n, k_n, params, i);
                    }
					}
					//  charge in seven cells (k_new != k_old)
					else if ((i_n==i_o)&&(k_n!=k_o))
					{        
						if (x3_old<k_n*dz)
						{
							z_boundary = k_n*dz;
							delta_z  = z_boundary - x3_old;
							a = (x1[i] - x1_old)/(x3[i] - x3_old);
							r_boundary = x1_old + a*delta_z;
							simple_j_weighting(dt,j1, j2, j3, r_boundary,z_boundary ,x1_old,x3_old, i_n, k_n-1, params,i);
							simple_j_weighting(dt,j1, j2, j3, x1[i],x3[i], r_boundary,z_boundary, i_n, k_n, params,i);
						}
						else 
						{
							z_boundary = (k_n+1)*dz;
							delta_z  = z_boundary - x3[i];
							a = (x1_old - x1[i])/(x3_old - x3[i]);
							r_boundary = x1[i] + a*delta_z;
							simple_j_weighting(dt,j1, j2, j3, r_boundary,z_boundary ,x1_old,x3_old, i_n, k_n+1, params,i);
							simple_j_weighting(dt,j1, j2, j3, x1[i],x3[i], r_boundary,z_boundary,i_n, k_n, params,i);
						}
					}
                }
                break;

				/////////////////////////////////////////
				///////// 3) charge in 10 cells /////////
				/////////////////////////////////////////
                case 2:
                {
					// case, when particle move from [i-1] cell to [i] cell
					if (i_o<i_n)
					/////////////////////////////////////////////////////////////////
					{
						// case, when particle move from [i-1][k-1] -> [i][k] cell
						if(k_o<k_n)
						{
							a = (x1[i] - x1_old)/(x3[i] - x3_old);
							r1 = i_n*dr;
							delta_z1 = (r1 - x1_old)/a;
							z1 = x3_old + delta_z1;
							z2 = k_n*dz;
							delta_r2 = (z2-x3_old)*a;
							r2 = x1_old+ delta_r2;
							if (z1<k_n*dz)
							{
								simple_j_weighting(dt, j1, j2, j3, r1, z1 ,x1_old, x3_old, i_n-1, k_n-1, params,i);
								simple_j_weighting(dt, j1, j2, j3, r2, z2, r1, z1, i_n, k_n-1, params,i);
								simple_j_weighting(dt, j1, j2, j3, x1[i], x3[i], r2, z2, i_n, k_n, params,i);
							}
							else if (z1>k_n*dz)
							{
								simple_j_weighting(dt, j1, j2, j3, r2, z2 ,x1_old, x3_old, i_n-1, k_n-1, params,i);
								simple_j_weighting(dt, j1, j2, j3, r1, z1, r2, z2,i_n-1, k_n, params,i);
								simple_j_weighting(dt, j1, j2, j3, x1[i], x3[i], r1, z1, i_n, k_n, params,i);
							}                
						}
					// case, when particle move from [i-1][k+1] -> [i][k] cell
					else 
					{
							a = (x1[i] - x1_old)/(x3[i] - x3_old);
							r1 = i_n*dr;
							delta_z1 = (r1 - x1_old)/a;
							z1 = x3_old + delta_z1;

							z2 = (k_n+1)*dz;
							delta_r2 = -(x3_old-z2)*a;
							r2 = x1_old+ delta_r2;
							if (z1>(k_n+1)*dz)
							{
								simple_j_weighting(dt, j1, j2, j3, r1, z1 ,x1_old, x3_old, i_n-1, k_n+1, params,i);
								simple_j_weighting(dt, j1, j2, j3, r2, z2, r1, z1, i_n, k_n+1, params,i);
								simple_j_weighting(dt, j1, j2, j3, x1[i], x3[i], r2, z2, i_n, k_n, params,i);
							}
							else if (z1<(k_n+1)*dz)
							{
								simple_j_weighting(dt, j1, j2, j3, r2, z2 ,x1_old, x3_old, i_n-1, k_n+1, params,i);
								simple_j_weighting(dt, j1, j2, j3, r1, z1, r2, z2,i_n-1, k_n, params,i);
								simple_j_weighting(dt, j1, j2, j3, x1[i], x3[i], r1, z1, i_n, k_n, params,i);
							}
						}
					}
					////////////////////////////////////////////////////////////////////////////
					// case, when particle move from [i+1] cell to [i] cell
					else if (i_o>i_n)
					{
						// case, when particle move from [i+1][k-1] -> [i][k] cell
						if(k_o<k_n)
						{
							a = (x1[i] - x1_old)/(x3[i] - x3_old);
							r1 = (i_n+1)*dr;
							delta_z1 = -(x1_old-r1)/a;
							z1 = x3_old + delta_z1;

							z2 = k_n*dz;
							delta_r2 = -(z2-x3_old)*a;
							r2 = x1_old- delta_r2;
                    
							if (z1<(k_n)*dz)
							{
								simple_j_weighting(dt, j1, j2, j3, r1, z1 ,x1_old, x3_old, i_n+1, k_n-1, params,i);
								simple_j_weighting(dt,j1, j2, j3, r2, z2, r1, z1, i_n, k_n-1, params,i);
								simple_j_weighting(dt,j1, j2, j3, x1[i], x3[i], r2, z2, i_n, k_n, params,i);
							}
							else if (z1>(k_n)*dz)
							{
								simple_j_weighting(dt, j1, j2, j3, r2, z2 ,x1_old, x3_old, i_n+1, k_n-1, params,i);
								simple_j_weighting(dt, j1, j2, j3, r1, z1, r2, z2,i_n+1, k_n, params,i);
								simple_j_weighting(dt, j1, j2, j3, x1[i], x3[i], r1, z1, i_n, k_n, params,i);
							}
						}
						// case, when particle move from [i+1][k+1] -> [i][k] cell
						else if (k_o>k_n)
						{
							a = (x1_old-x1[i])/(x3_old-x3[i]);
							r1 = (i_n+1)*dr;
							delta_z1 = (r1-x1[i])/a;
							z1 = x3[i] + delta_z1;

							z2 = (k_n+1)*dz;
							delta_r2 = (z2-x3[i])*a;
							r2 = x1[i] + delta_r2;
                
							if (z1>(k_n+1)*dz)
							{
								simple_j_weighting(dt, j1, j2, j3, r1, z1 ,x1_old,x3_old, i_n+1, k_n+1, params,i);
								simple_j_weighting(dt, j1, j2, j3, r2, z2, r1, z1, i_n, k_n+1, params,i);
								simple_j_weighting(dt, j1, j2, j3, x1[i], x3[i], r2, z2, i_n, k_n, params,i);
							}
							else if (z1<(k_n+1)*dz)
							{
								simple_j_weighting(dt, j1, j2, j3, r2, z2 ,x1_old, x3_old, i_n+1, k_n+1, params,i);
								simple_j_weighting(dt, j1, j2, j3, r1, z1, r2, z2,i_n+1, k_n, params,i);
								simple_j_weighting(dt, j1, j2, j3, x1[i], x3[i], r1, z1, i_n, k_n, params,i);
							}
						}
					}
                } //close case;
                break;
            } //close switch
        } //close condition
    }//close cycle
	//if (is_alive[i])
	//    printf("%d %e ", i, j3[0]);
}

void strict_motion_weighting(double dt, __global double *j1, __global double *j2, __global double *j3, 
                             double x1_new,double x3_new, double x1_old, double x3_old, 
                             __global const double *params, int gid)
{
    double size1 = params[0];
    double size3 = params[1];
    double dr = params[2];
    double dz = params[3];


    double charge = params[5];
    int grid_num3 = (int)size3;
    //defining number of cell
    int i_n = (int)ceil((x1_new)/dr)-1;
    int k_n =(int)ceil((x3_new)/dz)-1;
    int i_o = (int)ceil((x1_old)/dr)-1;
    int k_o =(int)ceil((x3_old)/dz)-1;

	//printf("%d %d %d %d %d", gid, i_n, k_n, i_o, k_o);
	//printf("%d %e ", gid, j3[0]);

    //stirct axis motion
//////////////////////////////////////////
    if (x1_new == x1_old)
    {
      double r1=0, r2=0,r3=0;
      double delta_z = 0.0;
      double value_part = 2.0*pi*x1_new*dr*dz;
      double wj_lower =0;
      r1 = x1_new-0.5*dr;
      r2 = (i_n+0.5)*dr;
      r3 = x1_new+0.5*dr;
      if (i_n==0)
        wj_lower = charge/(dt*pi*dr*dr/4.0) * pi*(r2*r2-r1*r1)/value_part;
      else
        wj_lower = charge/(dt*2*pi*i_n*dr*dr) * pi*(r2*r2-r1*r1)/value_part;

      double wj_upper =  charge/(dt*2*pi*(i_n+1)*dr*dr) *pi*(r3*r3-r2*r2)/value_part;
      double wj=0;
      //j1[get_linear_coord(i_n, k_n, grid_num3, 7)] += 0.0;
      //j1[get_linear_coord(i_n, k_n + 1, grid_num3, 7)] += 0.0;
      int res_k = k_n-k_o;
      switch(res_k)
      {
        case 0: 
        {
          delta_z = x3_new - x3_old;
          wj = wj_lower*delta_z;
		  //printf("%e ", wj);
          j3[get_linear_coord(i_n, k_n, grid_num3, 9)] += wj;


		  //printf("%e ", j3[0]);
          wj = wj_upper*delta_z;
          j3[get_linear_coord(i_n + 1, k_n, grid_num3, 9)] += wj;
        }
        break;

        case 1:
        {
          delta_z = k_n*dz - x3_old;
          wj = wj_lower*delta_z;
          j3[get_linear_coord(i_n, k_n - 1, grid_num3, 9)] += wj_lower;
          wj = wj_upper*delta_z;
          j3[get_linear_coord(i_n + 1, k_n - 1, grid_num3, 9)] += wj_lower;

          delta_z = x3_new - k_n*dz;
          wj = wj_lower*delta_z;
          j3[get_linear_coord(i_n, k_n, grid_num3, 9)] += wj;
          wj = wj_upper*delta_z;
          j3[get_linear_coord(i_n + 1, k_n, grid_num3, 9)] += wj;
        }
        break;

        case -1:
        {        
          delta_z = (k_n+1)*dz - x3_old;
          wj = wj_lower*delta_z;
          j3[get_linear_coord(i_n,k_n + 1, grid_num3, 9)] += wj_lower;
          wj = wj_upper*delta_z;
          j3[get_linear_coord(i_n + 1, k_n + 1, grid_num3, 9)] += wj_lower;

          delta_z = x3_new - (k_n+1)*dz;
          wj = wj_lower*delta_z;
          j3[get_linear_coord(i_n, k_n, grid_num3, 9)] += wj;
          wj = wj_upper*delta_z;
          j3[get_linear_coord(i_n + 1, k_n, grid_num3, 9)] += wj;
        }
        break;
      }
       }
///////////////////////////////////////////////////////
    
    ////stirct radial motion///
//////////////////////////////////////////////////////
    else if (x3_new==x3_old)
    {
        double r0  =(i_n+0.5)*dr;
        double wj= 0;
        double delta_r=0;
        double left_delta_z = 0, right_delta_z = 0;
        double res_j = 0;
        int res_i = i_n - i_o;
        switch(res_i)
        {
        case  0:
            {
                 delta_r = x1_new - x1_old;
                 left_delta_z = (k_n+1)*dz-x3_new;
                 right_delta_z = x3_new - k_n*dz;
                 wj = charge/(pi*4.0*r0*dz*dz*dr*dt)*(delta_r - r0*r0/(x1_old+delta_r) +r0*r0/x1_old + dr*dr/(4.0*(x1_old+delta_r)) - dr*dr/(4.0*x1_old));
                 res_j = wj*left_delta_z;
                 j1[get_linear_coord(i_n, k_n, grid_num3, 7)] += res_j;
                 res_j = wj*right_delta_z;
                 j1[get_linear_coord(i_n, k_n + 1, grid_num3, 7)] += res_j;

            }
            break;
        case 1:
            {
                 delta_r = (i_n)*dr- x1_old;
                 left_delta_z = (k_n+1)*dz-x3_new;
                 right_delta_z = x3_new - k_n*dz;
                 r0 = (i_n-0.5)*dr;
                 wj = charge/(pi*4.0*r0*dz*dz*dr*dt)*(delta_r - r0*r0/(x1_old+delta_r) +r0*r0/x1_old + dr*dr/(4.0*(x1_old+delta_r)) - dr*dr/(4.0*x1_old));
                 res_j = wj*left_delta_z;
                 j1[get_linear_coord(i_n - 1, k_n, grid_num3, 7)] += res_j;
                 res_j = wj*right_delta_z;
                 j1[get_linear_coord(i_n - 1, k_n + 1, grid_num3, 7)] += res_j;
                
                 delta_r = x1_new - i_n*dr;
                 r0 = (i_n+0.5)*dr;
                 wj = charge/(pi*4.0*r0*dz*dz*dr*dt)*(delta_r - r0*r0/(i_n*dr+delta_r) +r0*r0/i_n*dr + dr*dr/(4.0*(i_n*dr+delta_r)) - dr*dr/(4.0*i_n*dr));
                 res_j = wj*left_delta_z;
                 j1[get_linear_coord(i_n, k_n, grid_num3, 7)] += res_j;
                 res_j = wj*right_delta_z;
                 j1[get_linear_coord(i_n, k_n + 1, grid_num3, 7)] += res_j;
            }
            break;
            case -1:
            {
                 delta_r = (i_n+1)*dr - x1_old ;
                 left_delta_z = (k_n+1)*dz-x3_new;
                 right_delta_z = x3_new - k_n*dz;
                 r0 = (i_n+1.5)*dr;
                 wj = charge/(pi*4.0*r0*dz*dz*dr*dt)*(delta_r - r0*r0/(x1_old+delta_r) +r0*r0/x1_old + dr*dr/(4.0*(x1_old+delta_r)) - dr*dr/(4.0*x1_old));
                 res_j = wj*left_delta_z;
                 j1[get_linear_coord(i_n + 1, k_n, grid_num3, 7)] += res_j;
                 res_j = wj*right_delta_z;
                 j1[get_linear_coord(i_n + 1, k_n + 1, grid_num3, 7)] += res_j;
                
                 delta_r = x1_new - (i_n+1)*dr;
                 r0 = (i_n+0.5)*dr;
                 wj = charge/(pi*4.0*r0*dz*dz*dr*dt)*(delta_r - r0*r0/((i_n+1)*dr+delta_r) +r0*r0/(i_n+1)*dr + dr*dr/(4.0*((i_n+1)*dr+delta_r)) - dr*dr/(4.0*(i_n+1)*dr));
                 res_j = wj*left_delta_z;
                 j1[get_linear_coord(i_n, k_n, grid_num3, 7)] += res_j;
                 res_j = wj*right_delta_z;
                 j1[get_linear_coord(i_n, k_n + 1, grid_num3, 7)] += res_j;
            }
         break;
        }

    }
}


void simple_j_weighting(double dt, __global double *j1, __global double *j2, __global double *j3, double x1_new,double x3_new, 
                        double x1_old, double x3_old, int i_n, int k_n, __global const double *params, int gid)
{


    double size1 = params[0];
    double size3 = params[1];
    double dr = params[2];
    double dz = params[3];

    double charge = params[5];
    int grid_num3 = (int)size3;

    //printf("%d %d ", i_n, k_n);

    double wj = 0;


        // distance of particle moving//
        double delta_r = x1_new - x1_old;
        double delta_z = x3_new - x3_old;
        
        if ((delta_r==0)||(delta_z==0))
            return;
        // if i cell is not equal 0 
        if (i_n>=1)
        {
			///////////////////////////////////
			// equation y = k*x+b;//
			// finding k & b//
			double k = delta_r/delta_z;
			double b = x1_old;
			//calculate current jz in [i,k] cell//
			wj = charge/(2*dr*dz*dt*2*pi*i_n*dr*dr) * (dr*delta_z - k*delta_z*delta_z/2.0 - delta_z*b + dr*dr/k * ((i_n+0.5)*(i_n+0.5)-0.25)*log((k*delta_z+b)/b)); 

			// set new weighting current value 
			j3[get_linear_coord(i_n, k_n, grid_num3, 9)] += wj;


			//calculate current in [i+1,k] cell//
			wj = charge/(2*dr*dz*dt*2*pi*(i_n+1)*dr*dr) * (k*delta_z*delta_z/2.0+delta_z*b + delta_z*dr + dr*dr/k * (0.25-(i_n+0.5)*(i_n+0.5)) * log((k*delta_z+b)/b)); 
			// set new weighting current value 
			j3[get_linear_coord(i_n + 1, k_n, grid_num3, 9)] += wj;

        
			///////////////////////////////////
			//calculate current jr in [i,k] cell//
			// equation y = k*x+b;//
			// finding k & b//
			 k = -delta_z/delta_r;
			 double r0 = (i_n+0.5)*dr;
			 double r1 =  x1_old;
			 b= (k_n+1.0)*dz - x3_old;

			//weighting jr in [i][k] cell
			wj = charge/(2*pi*r0*dz*dz*dr*dt) *(r0*k*delta_r+k/2.0 * delta_r*(x1_old+delta_r/2.0)+0.5*delta_r*(b-k*(2*r0+r1))+ delta_r*(b-k*r1)*(4*r0*r0-dr*dr)/(8*x1_old*(x1_old+delta_r)) + (k*(r0*r0/2.0-dr*dr/8.0))*log((x1_old+delta_r)/x1_old));
			j1[get_linear_coord(i_n, k_n, grid_num3, 7)] += wj;

			  b= x3_old- k_n*dz;;
			 //weighting jr in [i][k+1] cell
			wj = charge/(2*pi*r0*dz*dz*dr*dt) * (-r0*k*delta_r - k/2.0*delta_r*(x1_old+delta_r/2.0)+0.5*delta_r*(b+k*(2*r0+r1))+ delta_r*(b+k*r1)*(4*r0*r0-dr*dr)/(8*x1_old*(x1_old+delta_r)) - (k*(r0*r0/2.0-dr*dr/8.0))*log((x1_old+delta_r)/x1_old));
			j1[get_linear_coord(i_n, k_n + 1, grid_num3, 7)] += wj;

        }
        // if i cell is equal 0 
        else
        {
///////////////////////////////////
        // equation y = k*x+b;//
        // finding k & b//
        double k = delta_r/delta_z;
        double b = x1_old;
        //calculate current jz in [i,k] cell//
        wj = charge/(2.0*dr*dz*dt*pi*dr*dr/4.0) * (dr*delta_z - k*delta_z*delta_z/2.0 - delta_z*b ); 
        // set new weighting current value 
        j3[get_linear_coord(i_n, k_n, grid_num3, 9)] += wj;

        //calculate current in [i+1,k] cell//
        wj = charge/(2.0*dr*dz*dt*2.0*pi*dr*dr) * (k*delta_z*delta_z/2.0 + delta_z*dr +delta_z*b); 
        // set new weighting current value 
        j3[get_linear_coord(i_n + 1, k_n, grid_num3, 9)] += wj;

        //j2[get_linear_coord(0, 0, grid_num3, 8)] = charge;
        //j2[get_linear_coord(0, 1, grid_num3, 8)] = dt;
        //j2[get_linear_coord(0, 2, grid_num3, 8)] = delta_z;
        //j2[get_linear_coord(0, 3, grid_num3, 8)] = dr;
        //j2[get_linear_coord(0, 4, grid_num3, 8)] = dz;
        
        ///////////////////////////////////
        //calculate current jr in [i,k] cell//
        // equation y = k*x+b;//
        // finding k & b//
         k = -delta_z/delta_r;
         double r0 = (i_n+0.5)*dr;
         double r1 =  x1_old;
         b= (k_n+1.0)*dz - x3_old;

        //weighting jr in [i][k] cell
        wj = charge/(2*pi*r0*dz*dz*dr*dt) *(r0*k*delta_r+k/2.0 * delta_r*(x1_old+delta_r/2.0)+0.5*delta_r*(b-k*(2*r0+r1))+ delta_r*(b-k*r1)*(4*r0*r0-dr*dr)/(8*x1_old*(x1_old+delta_r)) + (k*(r0*r0/2.0-dr*dr/8.0))*log((x1_old+delta_r)/x1_old));
        j1[get_linear_coord(i_n, k_n, grid_num3, 7)] += wj;

          b= x3_old- k_n*dz;;
         //weighting jr in [i][k+1] cell
         wj = charge/(2*pi*r0*dz*dz*dr*dt) * (-r0*k*delta_r - k/2.0*delta_r*(x1_old+delta_r/2.0)+0.5*delta_r*(b+k*(2*r0+r1))+ delta_r*(b+k*r1)*(4*r0*r0-dr*dr)/(8*x1_old*(x1_old+delta_r)) - (k*(r0*r0/2.0-dr*dr/8.0))*log((x1_old+delta_r)/x1_old));
        j1[get_linear_coord(i_n, k_n + 1, grid_num3, 7)] += wj;
        }
    
    //}
}

