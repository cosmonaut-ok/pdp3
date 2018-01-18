#include "Particles.h"
#include "E_field.h"
#include "H_field.h"
#include "pdp3_time.h"
#include "Triple.h"
#include <math.h>
#include <string.h>
#include "Constant.h"

using namespace std;
using namespace constant;

// C^2 define c^2 to decrease number of operations
const double LIGHT_SPEED_POW_2 = pow (LIGHT_SPEED, 2);

Particles::Particles(void)
{
}

// constructor
Particles::Particles(char const *p_name,
                     double p_charge,
                     double p_mass,
                     int p_number,
                     Geometry* geom,
                     particles_list* t_p_list):geom1(geom),
                                               p_list(t_p_list)
{
  name = (char *) p_name;
  charge = (double)p_charge*(double)EL_CHARGE;
  init_const_mass = p_mass;
  mass = (double)p_mass*(double)EL_MASS;
  number = p_number;

  //allocate memory for coordinates and velocities of particles
  mass_array = new double[number];
  charge_array = new double[number];
  x1 = new double[number];
  x3 = new double[number];
  v1 = new double[number];
  v2 = new double[number];
  v3 = new double[number];
  is_alive = new int[number];

  for (int i = 0; i < number; i++)
  {
    is_alive[i] = 1;
  }

  // insert to particles_lists
  p_list->part_list.push_back(this);
}

// copy constructor
Particles::Particles(Particles &cp_particles)
{
  name = new char[strlen(cp_particles.name)];
  strcpy(name,cp_particles.name);
  charge = cp_particles.charge;
  mass = cp_particles.mass;

  x1 = new double[number];
  x3 = new double[number];
  v1 = new double[number];
  v2 = new double[number];
  v3 = new double[number];
  is_alive = new int[number];
  for (int i=0; i<cp_particles.number; i++)
  {
    x1[i] = cp_particles.x1[i];
    x3[i] = cp_particles.x3[i];
    v1[i] = cp_particles.v1[i];
    v2[i] = cp_particles.v2[i];
    v3[i] = cp_particles.v3[i];
    mass_array[i]=cp_particles.mass_array[i];
    charge_array[i]=cp_particles.charge_array[i];
  }
}

// Destructor
Particles::~Particles()
{
  //delete [] x1;
  //delete [] x3;
  //delete [] v1;
  //delete [] v2;
  //delete [] v3;
  //delete [] is_alive;
}

void Particles::set_v_0()
{
  for(int i=0; i<number; i++)
  {
    v1[i] = 0.0;
    v2[i] = 1.0e5; // TODO: WHY?
    v3[i] = 0.0;
  }
}

void Particles::set_x_0()
{
  for(int i=0; i<number; i++)
  {
    x1[i] = 0.5; // TODO: WHY?
    x3[i] = 0.5; // TODO: WHY?
  }
}

// calculate Lorentz factor
double Particles::get_gamma(int i)
{
  double gamma;
  gamma = pow(1.0 - ((pow(v1[i], 2) + pow(v2[i], 2) + pow(v3[i], 2)) / LIGHT_SPEED_POW_2), -0.5);
  if (isinf(gamma) == 1) { // avoid infinity values
    cerr << "WARNING(get_gamma): gamma (Lorenz factor) girects to infinity for [v1, v2, v3, i]: ["
         << v1[i] << ", " << v2[i] << ", " << v3[i] << ", " << i << "]\n";
    return 1e100; // just return some very big value
  }
  else
  {
    return gamma;
  }
}

// clculate reciprocal Lorentz factor (1/gamma), aka ``alpha''
double Particles::get_gamma_inv(int i) // TODO: it is not alpha
{
  double gamma;
  gamma = pow(1.0 + ((pow(v1[i], 2) + pow(v2[i], 2) + pow(v3[i], 2)) / LIGHT_SPEED_POW_2), 0.5);
  if (isinf(gamma) == 1) { // avoid infinity values
    cerr << "WARNING(get_gamma_inv): reciprocal gamma (Lorenz factor) directs to infinity for [v1, v2, v3, i]: ["
         << v1[i] << ", " << v2[i] << ", " << v3[i] << ", " << i << "]\n";
    return 1e100; // just return some very big value
  }
  else
  {
    return gamma;
  }
}

void Particles::step_v(E_field *e_fld, H_field *h_fld, Time* t)
{
  double gamma, b1, b2, b3, e1, e2, e3, vv1, vv2, vv3;
  Triple E_compon(0.0, 0.0, 0.0), B_compon(0.0, 0.0, 0.0);
  double const1, const2;
  //if (t->current_time == t->start_time) const1 = const1/2.0;
  for(int i=0; i<number; i++)
    if (is_alive[i])
    {
      // check if x1 and x3 are correct
      if (isnan(x1[i]) || isinf(x1[i]) != 0 || isnan(x3[i]) || isinf(x3[i]) != 0)
      {
        cerr << "ERROR(step_v): x1[" << i << "] or x3[" << i << "] is not valid number. Can not continue." << endl;
        exit(1);
      }
      //
      const1 = charge_array[i]*t->delta_t/2.0/mass_array[i];
      E_compon = e_fld->get_field(x1[i],x3[i]);
      B_compon = h_fld->get_field(x1[i],x3[i]);
      e1 = E_compon.first*const1;
      e2 = E_compon.second*const1;
      e3 = E_compon.third*const1;
      b1 = B_compon.first*MAGN_CONST*const1;
      b2 = B_compon.second*MAGN_CONST*const1;
      b3 = B_compon.third*MAGN_CONST*const1;

      // 1. Multiplication by relativistic factor
      // u(n-1/2) = gamma(n-1/2)*v(n-1/2)
      gamma = get_gamma(i);
      v1[i] = gamma*v1[i];
      v2[i] = gamma*v2[i];
      v3[i] = gamma*v3[i];

      // 2. Half acceleration in the electric field
      // u'(n) = u(n-1/2) + q*dt/2/m*E(n)

      v1[i] = v1[i] + e1;
      v2[i] = v2[i] + e2;
      v3[i] = v3[i] + e3;

      // 3. Rotation in the magnetic field
      // u" = u' + 2/(1+B'^2)[(u' + [u'xB'(n)])xB'(n)]
      // B'(n) = B(n)*q*dt/2/mass/gamma(n)
      gamma = get_gamma_inv(i);
      b1 = b1/gamma;
      b2 = b2/gamma;
      b3 = b3/gamma;
      const2 = 2.0/(1.0 + b1*b1 + b2*b2 + b3*b3);
      vv1 = v1[i];
      vv2 = v2[i];
      vv3 = v3[i];
      v1[i] = vv1 + const2*((vv2 - vv1*b3 + vv3*b1)*b3 - (vv3 + vv1*b2 - vv2*b1)*b2);
      v2[i] = vv2 + const2*(-(vv1 + vv2*b3 - vv3*b2)*b3 + (vv3 + vv1*b2 - vv2*b1)*b1);
      v3[i] = vv3 + const2*((vv1 + vv2*b3 - vv3*b2)*b2 - (vv2 - vv1*b3 + vv3*b1)*b1);

      // 4. Half acceleration in the electric field
      // u(n+1/2) = u(n) + q*dt/2/m*E(n)
      v1[i] = v1[i] + e1;
      v2[i] = v2[i] + e2;
      v3[i] = v3[i] + e3;

      // 5. Division by relativistic factor
      gamma = get_gamma_inv(i);
      v1[i] = v1[i]/gamma;
      v2[i] = v2[i]/gamma;
      v3[i] = v3[i]/gamma;
    }
}

void Particles::half_step_coord(Time* t)
{
  double dr = geom1->dr;
  double dz = geom1->dz;
  double x1_wall = geom1->first_size - dr/2.0;
  double x3_wall = geom1->second_size - dz/2.0;
  double half_dr = dr/2.0;
  double half_dz = dz/2.0;
  double x1_wallX2 = x1_wall*2.0;
  double x3_wallX2 = x3_wall*2.0;
  double half_dt = t->delta_t/2.0;

  for(int i=0; i<number; i++)
    if (is_alive[i])
    {
      // check if x1 and x3 are correct
      if (isnan(x1[i]) || isinf(x1[i]) != 0 || isnan(x3[i]) || isinf(x3[i]) != 0)
      {
        cerr << "ERROR(half_step_coord): x1[" << i << "] or x3[" << i << "] is not valid number. Can not continue." << endl;
        exit(1);
      }
      //
      x1[i] = x1[i] + v1[i] * half_dt;
      x3[i] = x3[i] + v3[i] * half_dt;

      if (x1[i] > x1_wall)
      {
        x1[i] = x1_wallX2 - x1[i];
        v1[i] = -v1[i];
      }

      if (x3[i] > x3_wall)
      {
        x3[i] = x3_wallX2 - x3[i];
        v3[i] = -v3[i];
      }

      if (x1[i] < half_dr)
      {
        x1[i] = dr - x1[i];
        v1[i] = -v1[i];
      }

      if (x3[i] < half_dz)
      {
        x3[i] = dz - x3[i];
        v3[i] = -v3[i];
      }
    }
}


// function for charge density weighting
void Particles::charge_weighting(charge_density* ro1)
{
  int r_i=0;  // number of particle i cell
  int z_k=0;  // number of particle k cell

  double dr = geom1->dr;
  double dz = geom1->dz;
  double r1, r2, r3; // temp variables for calculation
  double dz1, dz2; // temp var.: width of k and k+1 cell

  double ro_v = 0; // charge density Q/V, V - volume of particle
  double v_1 = 0; // volume of [i][k] cell
  double v_2 = 0; // volume of [i+1][k] cell
  //double ro_v_2 = 0; // charge density in i+1 cell

  double value =0;
  // double **temp = ro1->get_ro();

  for(int i=0; i<number; i++)
    if (is_alive[i])
    {
      // check if x1 and x3 are correct
      if (isnan(x1[i]) || isinf(x1[i]) != 0 || isnan(x3[i]) || isinf(x3[i]) != 0)
      {
        cerr << "ERROR(charge_weighting): x1[" << i << "] or x3[" << i << "] is not valid number. Can not continue." << endl;
        exit(1);
      }

      // finding number of i and k cell. example: dr = 0.5; r = 0.4; i =0
      r_i = (int)ceil((x1[i])/dr)-1;
      z_k = (int)ceil((x3[i])/dz)-1;
      // TODO: workaround: sometimes it gives -1.
      // Just get 0 cell if it happence
      if (r_i < 0) { r_i = 0; }
      if (z_k < 0) { z_k = 0; }

      // in first cell other alg. of ro_v calc
      if(x1[i]>dr)
      {
        r1 =   x1[i] - 0.5*dr;
        r2 = (r_i+0.5)*dr;
        r3 = x1[i] + 0.5*dr;
        ro_v = charge_array[i]/(2.0*PI*dz*dr*x1[i]);
        v_1 = PI*dz*dr*dr*2.0*(r_i);
        v_2 = PI*dz*dr*dr*2.0*(r_i+1);
        dz1 = (z_k+1)*dz-x3[i];
        dz2 = x3[i] - z_k*dz;

        // weighting in ro[i][k] cell
        value = ro_v*PI*dz1*(r2*r2-r1*r1)/v_1;
        ro1->set_ro_weighting(r_i, z_k, value);

        // weighting in ro[i+1][k] cell
        value = ro_v*PI*dz1*(r3*r3-r2*r2)/v_2;
        ro1->set_ro_weighting(r_i+1,z_k, value);

        // weighting in ro[i][k+1] cell
        value = ro_v*PI*dz2*(r2*r2-r1*r1)/v_1;
        ro1->set_ro_weighting(r_i, z_k+1, value);

        // weighting in ro[i+1][k+1] cell
        value = ro_v*PI*dz2*(r3*r3-r2*r2)/v_2;
        ro1->set_ro_weighting(r_i+1, z_k+1, value);

      }
      else if (x1[i]<=dr/2.0)
      {
        r_i = 0;
        r1 =   0.0;
        r2 = (r_i+0.5)*dr;
        r3 = x1[i]+0.5*dr;
        dz1 = (z_k+1)*dz-x3[i];
        dz2 = x3[i] - z_k*dz;
        ro_v = charge_array[i]/(PI*dz*(2.0*x1[i]*x1[i]+dr*dr/2.0));
        v_1 = PI*dz*dr*dr/4.0;
        v_2 = PI*dz*dr*dr*2.0*(r_i+1);
        ///////////////////////////

        // weighting in ro[i][k] cell
        value = ro_v*PI*dz1*(dr*dr/2.0-x1[i]*dr+x1[i]*x1[i])/v_1;
        ro1->set_ro_weighting(r_i, z_k, value);

        // weighting in ro[i+1][k] cell
        value = ro_v*PI*dz1*(r3*r3-r2*r2)/v_2;
        ro1->set_ro_weighting(r_i+1,z_k, value);

        // weighting in ro[i][k+1] cell
        value = ro_v*PI*dz2*(dr*dr/2.0-x1[i]*dr+x1[i]*x1[i])/v_1;
        ro1->set_ro_weighting(r_i, z_k+1, value);

        // weighting in ro[i+1][k+1] cell
        value = ro_v*PI*dz2*(r3*r3-r2*r2)/v_2;
        ro1->set_ro_weighting(r_i+1, z_k+1, value);

      }
      else
      {
        ///////////////////////////
        r1 =   x1[i] - 0.5*dr;
        r2 = (r_i+0.5)*dr;
        r3 = x1[i]+0.5*dr;
        dz1 = (z_k+1)*dz-x3[i];
        dz2 = x3[i] - z_k*dz;
        ro_v = charge_array[i]/(2.0*PI*dz*dr*x1[i]);
        v_1 = PI*dz*dr*dr/4.0;
        v_2 = PI*dz*dr*dr*2.0*(r_i+1);
        ///////////////////////////

        // weighting in ro[i][k] cell
        value = ro_v*PI*dz1*(r2*r2-r1*r1)/v_1;
        ro1->set_ro_weighting(r_i, z_k, value);

        // weighting in ro[i+1][k] cell
        value = ro_v*PI*dz1*(r3*r3-r2*r2)/v_2;
        ro1->set_ro_weighting(r_i+1,z_k, value);

        // weighting in ro[i][k+1] cell
        value = ro_v*PI*dz2*(r2*r2-r1*r1)/v_1;
        ro1->set_ro_weighting(r_i, z_k+1, value);

        // weighting in ro[i+1][k+1] cell
        value = ro_v*PI*dz2*(r3*r3-r2*r2)/v_2;
        ro1->set_ro_weighting(r_i+1, z_k+1, value);
      }
    }
}


// function for initial (Maxwell) distribution
void Particles::velocity_distribution(double therm_vel)
{
  int j = 0;
  double R =0; // number from [0;1]
  double dv = therm_vel/1e7; // velocity step in calculation integral
  double s =0;
  double ds =0;
  double* v = new double [number]; //velocity vector
  double nr = sqrt(PI/2.0)* pow(therm_vel, 3);

  // double ss = 0;
  double const1 = 2*therm_vel*therm_vel;
  double temp1;

  for(int i=0; i<number; i++)
  {

    R = (double)(i)/(double)number;

    // part of numerical integral calculation
    while (s<R*nr)
    {
      temp1 = dv*j*dv*j;
      ds = temp1 * exp(-temp1 / const1) * dv;
      s = s+ds;
      j = j+1;
    }
    /////////////////////////
    v[i] = dv*j;
    ///////////////////////
    //double R_fi = (number-i)/double(number);
    //v1[i] = v[i]*sin(2.0*PI*R_fi);
    v1[i] = v[i];//*sin(2.0*PI*random_reverse(i,2));
    v3[i] = v[i]*cos(2.0*PI*random_reverse(i,2));
    v2[i] = v[i]*sin(2.0*PI*random_reverse(i,3));
  }

  delete []v;
}

void Particles::velocity_distribution_v2(double tempr_ev)
{
  double therm_vel = sqrt(tempr_ev*2.0*EL_CHARGE/
                          (this->init_const_mass*EL_MASS));
  int j=0;
  // double R =0; // number from [0;1]
  double dv = therm_vel/0.5e7; // velocity step in calculation integral
  double cutoff_vel = 9.0*therm_vel; //cutoff velocity
  int lenght_arr = (int)cutoff_vel/dv;
  double s =0;
  double ds =0;
  double* integ_array = new double [lenght_arr];

  double const1 = 2*therm_vel*therm_vel;
  // double temp1;
  // part of numerical integral calculation//
  int i = 0;
  while (i<lenght_arr)
  {
    /*temp1 = dv*i*dv*i;
      ds = temp1 * exp(-temp1 / const1 ) * dv;*/
    ds = exp(-dv*i*dv*i/const1)*dv;
    s=s+ds;
    integ_array[i] = s;
    i=i+1;
  }

  for(int i_n=0;i_n<number;i_n++)
  {
    double Rr = random_reverse(i_n,3);
    double Rfi = random_reverse(i_n,5);
    double Rz = random_reverse(i_n,7);
    double t_z = sqrt(PI/2.0)*therm_vel;
    // R = rand()/(double)32768;
    double f_vr = Rr*t_z;
    double f_vfi = Rfi*t_z;
    double f_vz = Rz*t_z;
    int sign =1;
    if (i_n%2==1)
      sign =-1;

    // binary search
    int i=0;
    j=lenght_arr;
    int k=0;
    while(i<=j)
    {
      k = i + (j-i)/2;
      if(f_vr>integ_array[k])
        i=k+1;
      else if (f_vr<integ_array[k])
        j=k-1;
      else
        break;
    }
    v1[i_n]=dv*k*sign;

    i=0;
    j=lenght_arr;
    k=0;
    while(i<=j)
    {
      k = i + (j-i)/2;
      if(f_vfi>integ_array[k])
        i=k+1;
      else if (f_vfi<integ_array[k])
        j=k-1;
      else
        break;
    }
    v2[i_n]=dv*k*sign;

    i=0;
    j=lenght_arr;
    k=0;
    while(i<=j)
    {
      k = i + (j-i)/2;
      if(f_vz>integ_array[k])
        i=k+1;
      else if (f_vz<integ_array[k])
        j=k-1;
      else
        break;
    }
    v3[i_n]=dv*k*sign;

    //////////////////////////////////////////////
    //double R_fi = (number-i)/double(number);
    //v1[i] = v[i]*sin(2.0*PI*R_fi);
    //v1[i_n] = v;//*sin(2.0*PI*random_reverse(i,2));
    //v3[i_n] = v*cos(PI*(rand()/32768 - 0.5));//*sin(2.0*PI*random_reverse(i,3));
    ////v3[i_n] = v*cos(2.0*PI*random_reverse(i,2));
    ////v2[i_n] = v*sin(2.0*PI*random_reverse(i,3));
    //v2[i_n] = v*sin(-PI/2 + PI*rand()/32768);

  }
  delete []integ_array;
}

double Particles::random_reverse(double vel, int power)
{
  int int_vel =(int) floor(vel);
  double ost =0;
  double r =0;
  int order =1;
  while(int_vel>=1)
  {
    ost = int_vel % power;

    r = r + ost*pow((double)power,(-order));

    int_vel = (int_vel - ost)/power;
    order = order+1;
  }
  return r;
}

void Particles::load_spatial_distribution(double n1, double n2, double left_plasma_boundary,int type)
{
  // calculate number of electrons in a big particle
  double rand_r;
  double rand_z;
  double dr = geom1->dr*1.00000001;
  double dz = geom1->dz*1.00000001;
  double dn = n2 - n1;
  switch (type) {
  case 0:
  {
    double n_in_big = (PI*geom1->first_size*geom1->first_size*geom1->second_size/number*(n2+n1)/2.0);
    charge *= n_in_big;
    mass *= n_in_big;
    for(int n = 0; n < number; n++)
    {
      rand_r = random_reverse(n,13);
      rand_z = random_reverse(number - 1 - n,11);
      x1[n] = sqrt(rand_r*geom1->first_size*(geom1->first_size-dr)+dr*dr/4.0);
      //x3[n] = (geom1->second_size - dz)*sqrt(rand_z) + dz/2.0;
      x3[n] = (geom1->second_size - left_plasma_boundary - dz)/dn*(sqrt(n1*n1 + rand_z*(2*n1*dn + dn*dn)) - n1) +
        left_plasma_boundary + dz/2.0;
      mass_array[n]=mass;
      charge_array[n]=charge;
    }
  }
  break;
  case 1: //Normal distribution
  {

    double sigma = 0.01;
    double n_in_big = (PI*geom1->second_size*(n2+n1)*sigma*sigma*(1.0-exp(-geom1->first_size*geom1->first_size/(2.0*sigma*sigma)))/number);
    charge *= n_in_big;
    mass *= n_in_big;
    double R_sq= (geom1->first_size - dr/2.0)*(geom1->first_size - dr/2.0 );
    // double tt=0;
    for (int i = 0; i<number;i++)
    {
      rand_r = random_reverse(i,13);
      double int_rd =   exp(-dr*dr/(8.0*sigma*sigma));
      x1[i]=sigma*sqrt(-2.0*log(int_rd - rand_r*(int_rd-exp(-R_sq/(2.0*sigma*sigma)))));

      //x1[i] = (geom1->first_size - dr)*(rand_r)*rand_r + dr/2.0;
      // double tt = exp(-R_sq/(2.0*sigma*sigma));
      rand_z = random_reverse(number - 1 - i,11);
      x3[i] = (geom1->second_size - dz/2.0)*rand_z + dz/2.0;
      //x3[i] = (geom1->second_size - left_plasma_boundary - dz)/dn*(sqrt(n1*n1 + rand_z*(2*n1*dn + dn*dn)) - n1) + left_plasma_boundary + dz/2.0;
    }

//      int i = 0;
//      int j=0;
//      double R =0; // number from [0;1]
//      double dr_int = (geom1->first_size-dr/1.0)/1e7; // velocity step in calculation integral
//      int lenght_arr = 1e7;
//      double s =0;
//      double ds =0;
//      double* integ_array = new double [lenght_arr];
//
//      double const1 = 2*sigma*sigma;
//
//    // part of numerical integral calculation//
//       while (i<lenght_arr)
//      {
//     /*temp1 = dv*i*dv*i;
//     ds = temp1 * exp(-temp1 / const1 ) * dv;*/
//      ds = exp(-dr_int*i*dr_int*i/const1)*dr_int;
//       s=s+ds;
//       integ_array[i] = s;
//       i=i+1;
//       }
// ///////////////////////////////////////
//
// for(int i_n=0;i_n<number;i_n++)
//{
// double Rr = random_reverse(i_n,13);
// double t_z = s;
// //R = rand()/(double)32768;
// double f_vr = Rr*t_z;
//
// //binary search//
// ///////////////////////////////////
// i=0;
// j=lenght_arr;
// int k=0;
//while(i<=j)
//{
//  k = i + (j-i)/2;
//  if(f_vr>integ_array[k])
//    i=k+1;
//  else if (f_vr<integ_array[k])
//    j=k-1;
//  else
//    break;
//}
//x1[i_n]=dr_int*k+dr/2.0;
// }
//////////////////////////////////////////
//
//
//      //x3[n] = (geom1->second_size - dz)*sqrt(rand_z) + dz/2.0;
//        for(n = 0; n < number; n++)
//      {
//        rand_z = random_reverse(number - 1 - n,11);
//        x3[n] = (geom1->second_size - left_plasma_boundary - dz)/dn*(sqrt(n1*n1 + rand_z*(2*n1*dn + dn*dn)) - n1) + left_plasma_boundary + dz/2.0;
//      }

  }
  break;
  }
}

void Particles::load_spatial_distribution_with_variable_mass(double n1,
                                                             double n2,
                                                             double left_plasma_boundary,
                                                             int type)
{
  // calculate number of electrons in a big particle
  double rand_r;
  double rand_z;
  double dr = geom1->dr*1.00000001;
  double dz = geom1->dz*1.00000001;
  double dn = n2 - n1;
  switch (type) {
  case 0:
  {

    double N_big_for_cell=(double) number/( (double) geom1->n_grid_1*geom1->n_grid_2);
    double N_real_i = 8.0*PI*(n2+n1)/2.0*dr*dz;
    double n_in_big =0;
    for(int n = 0; n < number; n++)
    {
      // check if x1 and x3 are correct
      if (isnan(x1[n]) || isinf(x1[n]) != 0 || isnan(x3[n]) || isinf(x3[n]) != 0)
      {
        cerr << "ERROR(load_spatial_distribution_with_variable_mass): x1[" << n << "] or x3[" << n << "] is not valid number. Can not continue." << endl;
        exit(1);
      }

      rand_r = random_reverse(n,13);
      rand_z = random_reverse(number - 1 - n,11);
      x1[n] = (geom1->first_size - dr)*(rand_r) + dr/2.0;
      n_in_big =N_real_i*x1[n]/N_big_for_cell;
      charge_array[n]=charge* n_in_big;
      mass_array[n]= mass*n_in_big;
      //x3[n] = (geom1->second_size - dz)*(rand_z) + dz/2;
      //x3[n] = (geom1->second_size - dz)*sqrt(rand_z) + dz/2.0;

      x3[n] = (geom1->second_size - left_plasma_boundary - dz)/dn*(sqrt(n1*n1 + rand_z*(2*n1*dn + dn*dn)) - n1) +
        left_plasma_boundary + dz/2.0;
    }
  }
  break;
  }
}

void Particles::load_velocity_distribution(double v_thermal)
{
  for (int n=0; n<number; n++)
  {
    v1[n] = 0.0;
    v2[n] = 0.0;
    v3[n] = 0.0;
  }
}

void Particles::simple_j_weighting(Time* time1,
                                   current *j1,
                                   double x1_new,
                                   double x3_new,
                                   double x1_old,
                                   double x3_old,
                                   int i_n,
                                   int k_n,
                                   int p_number)
{
  double dr = geom1->dr;
  double dz = geom1->dz;
  double wj = 0;
  double delta_t = time1->delta_t;

  // distance of particle moving//
  double delta_r = x1_new - x1_old;
  double delta_z = x3_new - x3_old;

  if ((abs(delta_r)<1e-15)||(abs(delta_z)<1e-15))
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
    wj = charge_array[p_number]/(2*dr*dz*delta_t*2*PI*i_n*dr*dr) * (dr*delta_z - k*delta_z*delta_z/2.0 - delta_z*b + dr*dr/k * ((i_n+0.5)*(i_n+0.5)-0.25)*log((k*delta_z+b)/b));
    // set new weighting current value
    j1->set_j3(i_n,k_n, wj);

    //calculate current in [i+1,k] cell//
    wj = charge_array[p_number]/(2*dr*dz*delta_t*2*PI*(i_n+1)*dr*dr) * (k*delta_z*delta_z/2.0+delta_z*b + delta_z*dr + dr*dr/k * (0.25-(i_n+0.5)*(i_n+0.5)) * log((k*delta_z+b)/b));
    // set new weighting current value
    j1->set_j3(i_n+1,k_n, wj);

    ///////////////////////////////////
    //calculate current jr in [i,k] cell//
    // equation y = k*x+b;//
    // finding k & b//
    k = -delta_z/delta_r;
    double r0 = (i_n+0.5)*dr;
    double r1 =   x1_old;
    b= (k_n+1.0)*dz - x3_old;

    //weighting jr in [i][k] cell
    wj = charge_array[p_number]/(2*PI*r0*dz*dz*dr*delta_t) *(r0*k*delta_r+k/2.0 * delta_r*(x1_old+delta_r/2.0)+0.5*delta_r*(b-k*(2*r0+r1))+ delta_r*(b-k*r1)*(4*r0*r0-dr*dr)/(8*x1_old*(x1_old+delta_r)) + (k*(r0*r0/2.0-dr*dr/8.0))*log((x1_old+delta_r)/x1_old));
    j1->set_j1(i_n,k_n, wj);

    b= x3_old- k_n*dz;;
    //weighting jr in [i][k+1] cell
    wj = charge_array[p_number]/(2*PI*r0*dz*dz*dr*delta_t) * (-r0*k*delta_r - k/2.0*delta_r*(x1_old+delta_r/2.0)+0.5*delta_r*(b+k*(2*r0+r1))+ delta_r*(b+k*r1)*(4*r0*r0-dr*dr)/(8*x1_old*(x1_old+delta_r)) - (k*(r0*r0/2.0-dr*dr/8.0))*log((x1_old+delta_r)/x1_old));
    j1->set_j1(i_n, k_n+1, wj);
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
    wj = charge_array[p_number]/(2.0*dr*dz*delta_t*PI*dr*dr/4.0) * (dr*delta_z - k*delta_z*delta_z/2.0 - delta_z*b );
    // set new weighting current value
    j1->set_j3(i_n,k_n, wj);

    //calculate current in [i+1,k] cell//
    wj = charge_array[p_number]/(2.0*dr*dz*delta_t*2.0*PI*dr*dr) * (k*delta_z*delta_z/2.0 + delta_z*dr +delta_z*b);
    // set new weighting current value
    j1->set_j3(i_n+1,k_n, wj);

    ///////////////////////////////////
    //calculate current jr in [i,k] cell//
    // equation y = k*x+b;//
    // finding k & b//
    k = -delta_z/delta_r;
    double r0 = (i_n+0.5)*dr;
    double r1 =   x1_old;
    b= (k_n+1.0)*dz - x3_old;

    //weighting jr in [i][k] cell
    wj = charge_array[p_number]/(2*PI*r0*dz*dz*dr*delta_t) *(r0*k*delta_r+k/2.0 * delta_r*(x1_old+delta_r/2.0)+0.5*delta_r*(b-k*(2*r0+r1))+ delta_r*(b-k*r1)*(4*r0*r0-dr*dr)/(8*x1_old*(x1_old+delta_r)) + (k*(r0*r0/2.0-dr*dr/8.0))*log((x1_old+delta_r)/x1_old));
    j1->set_j1(i_n,k_n, wj);

    b= x3_old- k_n*dz;;
    //weighting jr in [i][k+1] cell
    wj = charge_array[p_number]/(2*PI*r0*dz*dz*dr*delta_t) * (-r0*k*delta_r - k/2.0*delta_r*(x1_old+delta_r/2.0)+0.5*delta_r*(b+k*(2*r0+r1))+ delta_r*(b+k*r1)*(4*r0*r0-dr*dr)/(8*x1_old*(x1_old+delta_r)) - (k*(r0*r0/2.0-dr*dr/8.0))*log((x1_old+delta_r)/x1_old));
    j1->set_j1(i_n, k_n+1, wj);
  }

  //}
}
void Particles::simple_constrho_j_weighting(Time* time1, current *j1, double x1_new,double x3_new, double x1_old, double x3_old, int i_n, int k_n,int p_number)
{
  double dr = geom1->dr;
  double dz = geom1->dz;
  double wj = 0;
  double delta_t = time1->delta_t;

  // distance of particle moving//
  double delta_r = x1_new - x1_old;
  double delta_z = x3_new - x3_old;

  // if i cell is not equal 0
  if (i_n>=1)
  {
    ///////////////////////////////////
    // equation y = k*x+b;//
    // finding k & b//
    // double k = delta_r/delta_z;
    // double b = x1_old-dr/2.0;
//    double r0 = (i_n+0.5)*dr;
    //calculate current jz in [i,k] cell//
    //  wj = particle_rho*PI*delta_z/(delta_t*2*PI*i_n*dr*dr) * (r0*r0 - b*b -a*a*delta_z*delta_z/3.0 - a*b*delta_z);
    // set new weighting current value
    j1->set_j3(i_n,k_n, wj);

    //calculate current in [i+1,k] cell//
    //  wj = particle_rho*PI*delta_z/(delta_t*2*PI*(i_n+1)*dr*dr) * (a*a*delta_z*delta_z/3.0 + a*b*delta_z+b*b-r0*r0);
    // set new weighting current value
    j1->set_j3(i_n+1,k_n, wj);

    ///////////////////////////////////
    // calculate current jr in [i,k] cell
    // equation y = k*x+b;//
    // finding k & b//
    // double k = -delta_z/delta_r;
    // double r0 = (i_n+0.5)*dr;
    // double r1 =   x1_old;
    // double b = (k_n+1.0)*dz - x3_old;

    //weighting jr in [i][k] cell
// wj = charge/(2*PI*r0*dz*dz*dr*delta_t) *();
    j1->set_j1(i_n,k_n, wj);

    // double b = x3_old- k_n*dz;;
    //weighting jr in [i][k+1] cell
    //   wj = charge/(2*PI*r0*dz*dz*dr*delta_t) * ();
    j1->set_j1(i_n, k_n+1, wj);
  }
  // if i cell is equal 0
  else
  {
    // equation y = k*x+b;
    // finding k & b
    double k = delta_r/delta_z;
    double b = x1_old;
    //calculate current jz in [i,k] cell//
    //  wj = charge/(2.0*dr*dz*delta_t*PI*dr*dr/4.0) * (dr*delta_z - k*delta_z*delta_z/2.0 - delta_z*b );
    // set new weighting current value
    j1->set_j3(i_n,k_n, wj);

    //calculate current in [i+1,k] cell//
    //  wj = charge/(2.0*dr*dz*delta_t*2.0*PI*dr*dr) * (k*delta_z*delta_z/2.0 + delta_z*dr +delta_z*b);
    // set new weighting current value
    j1->set_j3(i_n+1,k_n, wj);

    ///////////////////////////////////
    //calculate current jr in [i,k] cell//
    // equation y = k*x+b;//
    // finding k & b//
    k = -delta_z/delta_r;
    double r0 = (i_n+0.5)*dr;
    double r1 =   x1_old;
    b= (k_n+1.0)*dz - x3_old;

    //weighting jr in [i][k] cell
    wj = charge_array[p_number]/(2*PI*r0*dz*dz*dr*delta_t) *(r0*k*delta_r+k/2.0 * delta_r*(x1_old+delta_r/2.0)+0.5*delta_r*(b-k*(2*r0+r1))+ delta_r*(b-k*r1)*(4*r0*r0-dr*dr)/(8*x1_old*(x1_old+delta_r)) + (k*(r0*r0/2.0-dr*dr/8.0))*log((x1_old+delta_r)/x1_old));
    j1->set_j1(i_n,k_n, wj);

    b= x3_old- k_n*dz;;
    //weighting jr in [i][k+1] cell
    wj = charge_array[p_number]/(2*PI*r0*dz*dz*dr*delta_t) * (-r0*k*delta_r - k/2.0*delta_r*(x1_old+delta_r/2.0)+0.5*delta_r*(b+k*(2*r0+r1))+ delta_r*(b+k*r1)*(4*r0*r0-dr*dr)/(8*x1_old*(x1_old+delta_r)) - (k*(r0*r0/2.0-dr*dr/8.0))*log((x1_old+delta_r)/x1_old));
    j1->set_j1(i_n, k_n+1, wj);
  }

  //}
}

void Particles::j_weighting(Time* time1, current *j1, double* x1_o,double* x3_o)
{

  double dr = geom1->dr;
  double dz = geom1->dz;

  //double **J1 = j1->get_j1();
  //double **J2 = j1->get_j2();
  //double **J3 = j1->get_j3();

  for (int i=0;i<number;i++)
    if (is_alive[i])
    {
      double x1_old= x1_o[i];
      double x3_old = x3_o[i];
      //finding number new and old cells
      int i_n = (int)ceil((x1[i])/dr)-1;
      int k_n =(int)ceil((x3[i])/dz)-1;
      int i_o = (int)ceil((x1_old)/dr)-1;
      int k_o =(int)ceil((x3_old)/dz)-1;
      // TODO: workaround: sometimes it gives -1.
      // Just get 0 cell if it happence
      if (i_n < 0) { i_n = 0; }
      if (k_n < 0) { k_n = 0; }
      if (i_o < 0) { i_o = 0; }
      if (k_o < 0) { k_o = 0; }

      if (x1_old==(i_o+1)*dr)
        i_o=i_n;
      if(x3_old==(k_o+1)*dz)
        k_o=k_n;
      if (x1[i]==(i_n+1)*dr)
        i_n=i_o;
      if(x3[i]==(k_n+1)*dz)
        k_n=k_o;
      int res_cell = abs(i_n-i_o) + abs(k_n-k_o);
      if ((abs(x1[i]-x1_old)<1e-15)||(abs(x3[i]-x3_old)<1e-15))
      {
        strict_motion_weighting(time1, j1,x1[i],x3[i],x1_old,x3_old,i);
      }
      else
      {
        switch (res_cell)
        {
          /// 1) charge in four cells
        case 0: simple_j_weighting(time1, j1, x1[i],x3[i] ,x1_old,x3_old, i_n, k_n,i);
          break;

          /// 2) charge in seven cells
        case 1:
        {
          /// charge in seven cells (i_new != i_old)
          if ((i_n!=i_o)&&(k_n==k_o))
          {
            if (x1_old >(i_n+1)*dr)
            {
              double a = (x1_old - x1[i])/(x3_old - x3[i]);
              double r_boundary = (i_n+1)*dr;
              double delta_r = r_boundary - x1[i];
              double z_boundary = x3[i] + delta_r/a;

              simple_j_weighting(time1,j1, r_boundary,z_boundary ,x1_old,x3_old, i_n+1,k_n,i);
              simple_j_weighting(time1,j1, x1[i],x3[i], r_boundary,z_boundary, i_n, k_n,i);
            }
            else
            {
              double a = (x1[i] - x1_old)/(x3[i] - x3_old);
              double r_boundary = (i_n)*dr;
              double delta_r = r_boundary - x1_old;
              double z_boundary = x3_old + delta_r/a;

              simple_j_weighting(time1,j1, r_boundary,z_boundary ,x1_old,x3_old, i_n-1, k_n,i);
              simple_j_weighting(time1,j1, x1[i],x3[i], r_boundary,z_boundary, i_n, k_n,i);
            }

          }
          //  charge in seven cells (k_new != k_old)
          else if ((i_n==i_o)&&(k_n!=k_o))
          {
            if (x3_old<k_n*dz)
            {
              double z_boundary = k_n*dz;
              double delta_z  = z_boundary - x3_old;
              double a = (x1[i] - x1_old)/(x3[i] - x3_old);
              double r_boundary = x1_old + a*delta_z;
              simple_j_weighting(time1,j1, r_boundary,z_boundary ,x1_old,x3_old, i_n, k_n-1,i);
              simple_j_weighting(time1,j1, x1[i],x3[i], r_boundary,z_boundary, i_n, k_n,i);
            }
            else
            {
              double z_boundary = (k_n+1)*dz;
              double delta_z  = z_boundary - x3[i];
              double a = (x1_old - x1[i])/(x3_old - x3[i]);
              double r_boundary = x1[i] + a*delta_z;
              simple_j_weighting(time1,j1, r_boundary,z_boundary ,x1_old,x3_old, i_n, k_n+1,i);
              simple_j_weighting(time1,j1, x1[i],x3[i], r_boundary,z_boundary,i_n, k_n,i);
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
              double a = (x1[i] - x1_old)/(x3[i] - x3_old);
              double r1 = i_n*dr;
              double delta_z1 = (r1 - x1_old)/a;
              double z1 = x3_old + delta_z1;
              double z2 = k_n*dz;
              double delta_r2 = (z2-x3_old)*a;
              double r2 = x1_old+ delta_r2;
              if (z1<k_n*dz)
              {
                simple_j_weighting(time1, j1, r1, z1 ,x1_old, x3_old, i_n-1, k_n-1,i);
                simple_j_weighting(time1, j1, r2, z2, r1, z1, i_n, k_n-1,i);
                simple_j_weighting(time1, j1, x1[i], x3[i], r2, z2, i_n, k_n,i);
              }
              else if (z1>k_n*dz)
              {
                simple_j_weighting(time1, j1, r2, z2 ,x1_old, x3_old, i_n-1, k_n-1,i);
                simple_j_weighting(time1, j1, r1, z1, r2, z2,i_n-1, k_n,i);
                simple_j_weighting(time1, j1, x1[i], x3[i], r1, z1, i_n, k_n,i);
              }
            }
            // case, when particle move from [i-1][k+1] -> [i][k] cell
            else
            {
              double a = (x1[i] - x1_old)/(x3[i] - x3_old);
              double r1 = i_n*dr;
              double delta_z1 = (r1 - x1_old)/a;
              double z1 = x3_old + delta_z1;

              double z2 = (k_n+1)*dz;
              double delta_r2 = -(x3_old-z2)*a;
              double r2 = x1_old+ delta_r2;
              if (z1>(k_n+1)*dz)
              {
                simple_j_weighting(time1, j1, r1, z1 ,x1_old, x3_old, i_n-1, k_n+1,i);
                simple_j_weighting(time1, j1, r2, z2, r1, z1, i_n, k_n+1,i);
                simple_j_weighting(time1, j1, x1[i], x3[i], r2, z2, i_n, k_n,i);
              }
              else if (z1<(k_n+1)*dz)
              {
                simple_j_weighting(time1, j1, r2, z2 ,x1_old, x3_old, i_n-1, k_n+1,i);
                simple_j_weighting(time1, j1, r1, z1, r2, z2,i_n-1, k_n,i);
                simple_j_weighting(time1, j1, x1[i], x3[i], r1, z1, i_n, k_n,i);
              }
            }
          }
          // case, when particle move from [i+1] cell to [i] cell
          else if (i_o>i_n)
          {
            // case, when particle move from [i+1][k-1] -> [i][k] cell
            if(k_o<k_n)
            {
              double a = (x1[i] - x1_old)/(x3[i] - x3_old);
              double r1 = (i_n+1)*dr;
              double delta_z1 = -(x1_old-r1)/a;
              double z1 = x3_old + delta_z1;

              double z2 = k_n*dz;
              double delta_r2 = -(z2-x3_old)*a;
              double r2 = x1_old- delta_r2;

              if (z1<(k_n)*dz)
              {
                simple_j_weighting(time1, j1, r1, z1 ,x1_old, x3_old, i_n+1, k_n-1,i);
                simple_j_weighting(time1,j1, r2, z2, r1, z1, i_n, k_n-1,i);
                simple_j_weighting(time1,j1, x1[i], x3[i], r2, z2, i_n, k_n,i);
              }
              else if (z1>(k_n)*dz)
              {
                simple_j_weighting(time1, j1, r2, z2 ,x1_old, x3_old, i_n+1, k_n-1,i);
                simple_j_weighting(time1, j1, r1, z1, r2, z2,i_n+1, k_n,i);
                simple_j_weighting(time1, j1, x1[i], x3[i], r1, z1, i_n, k_n,i);
              }

            }
            // case, when particle move from [i+1][k+1] -> [i][k] cell
            else if (k_o>k_n)
            {
              double a = (x1_old-x1[i])/(x3_old-x3[i]);
              double r1 = (i_n+1)*dr;
              double delta_z1 = (r1-x1[i])/a;
              double z1 = x3[i] + delta_z1;

              double z2 = (k_n+1)*dz;
              double delta_r2 = (z2-x3[i])*a;
              double r2 = x1[i] + delta_r2;

              if (z1>(k_n+1)*dz)
              {
                simple_j_weighting(time1, j1, r1, z1 ,x1_old,x3_old, i_n+1, k_n+1,i);
                simple_j_weighting(time1, j1, r2, z2, r1, z1, i_n, k_n+1,i);
                simple_j_weighting(time1, j1, x1[i], x3[i], r2, z2, i_n, k_n,i);
              }
              else if (z1<(k_n+1)*dz)
              {
                simple_j_weighting(time1, j1, r2, z2 ,x1_old, x3_old, i_n+1, k_n+1,i);
                simple_j_weighting(time1, j1, r1, z1, r2, z2,i_n+1, k_n,i);
                simple_j_weighting(time1, j1, x1[i], x3[i], r1, z1, i_n, k_n,i);
              }
            }
          }
        }
        break;
        }
      }
    }

}
void Particles::azimuthal_j_weighting(Time* time1, current *this_j)
{

  int r_i=0;  // number of particle i cell
  int z_k=0;  // number of particle k cell

  double dr = geom1->dr;
  double dz = geom1->dz;
  double r1, r2, r3; // temp variables for calculation
  double dz1, dz2;   // temp var.: width of k and k+1 cell

  double ro_v =0; // charge density Q/V, V - volume of particle
  double v_1 =0; // volume of [i][k] cell
  double v_2= 0; // volume of [i+1][k] cell
  //double ro_v_2=0; // charge density in i+1 cell

  double rho =0; //charge density in cell
  double current; // j_phi in cell
  // double **temp = this_j->get_j2();

  for(int i=0;i<number;i++)
    if (is_alive[i])
    {
      // finding number of i and k cell. example: dr = 0.5; r = 0.4; i =0
      r_i = (int)ceil((x1[i])/dr)-1;
      z_k = (int)ceil((x3[i])/dz)-1;
      // TODO: workaround: sometimes it gives -1.
      // Just get 0 cell if it happence
      if (r_i < 0) { r_i = 0; }
      if (z_k < 0) { z_k = 0; }

      // in first cell other alg. of ro_v calc
      if(x1[i]>dr)
      {
        r1 = x1[i] - 0.5*dr;
        r2 = (r_i+0.5)*dr;
        r3 = x1[i] + 0.5*dr;
        ro_v = charge_array[i]/(2.0*PI*dz*dr*x1[i]);
        v_1 = PI*dz*dr*dr*2.0*(r_i);
        v_2 = PI*dz*dr*dr*2.0*(r_i+1);
        dz1 = (z_k+1)*dz-x3[i];
        dz2 = x3[i] - z_k*dz;

        // weighting in j[i][k] cell
        rho = ro_v*PI*dz1*(r2*r2-r1*r1)/v_1;
        current = rho*v2[i];
        this_j->set_j2(r_i, z_k, current);

        // weighting in j[i+1][k] cell
        rho = ro_v*PI*dz1*(r3*r3-r2*r2)/v_2;
        current = rho*v2[i];
        this_j->set_j2(r_i+1,z_k, current);

        // weighting in j[i][k+1] cell
        rho = ro_v*PI*dz2*(r2*r2-r1*r1)/v_1;
        current = rho*v2[i];
        this_j->set_j2(r_i, z_k+1,  current);

        // weighting in j[i+1][k+1] cell
        rho = ro_v*PI*dz2*(r3*r3-r2*r2)/v_2;
        current = rho*v2[i];
        this_j->set_j2(r_i+1, z_k+1, current);

      }
      else
      {
        r1 = x1[i] - 0.5*dr;
        r2 = (r_i+0.5)*dr;
        r3 = x1[i]+0.5*dr;
        dz1 = (z_k+1)*dz-x3[i];
        dz2 = x3[i] - z_k*dz;
        ro_v = charge_array[i]/(2.0*PI*dz*dr*x1[i]);
        v_1 = PI*dz*dr*dr/4.0;
        v_2 = PI*dz*dr*dr*2.0*(r_i+1);

        // weighting in j[i][k] cell
        rho = ro_v*PI*dz1*(r2*r2-r1*r1)/v_1;
        current = rho*v2[i];
        this_j->set_j2(r_i, z_k,  current);

        // weighting in j[i+1][k] cell
        rho = ro_v*PI*dz1*(r3*r3-r2*r2)/v_2;
        current = rho*v2[i];
        this_j->set_j2(r_i+1,z_k,    current);

        // weighting in j[i][k+1] cell
        rho = ro_v*PI*dz2*(r2*r2-r1*r1)/v_1;
        current = rho*v2[i];
        this_j->set_j2(r_i, z_k+1,  current);

        // weighting in j[i+1][k+1] cell
        rho = ro_v*PI*dz2*(r3*r3-r2*r2)/v_2;
        current = rho*v2[i];
        this_j->set_j2(r_i+1, z_k+1,  current);
      }
    }
}

void Particles::strict_motion_weighting(Time *time1,
                                        current *this_j,
                                        double x1_new,
                                        double x3_new,
                                        double x1_old,
                                        double x3_old,
                                        int p_number)
{
  double dr = geom1->dr;
  double dz = geom1->dz;

  // defining number of cell
  int i_n = (int)ceil((x1_new)/dr)-1;
  int k_n =(int)ceil((x3_new)/dz)-1;
  int i_o = (int)ceil((x1_old)/dr)-1;
  int k_o =(int)ceil((x3_old)/dz)-1;
  // TODO: workaround: sometimes it gives -1.
  // Just get 0 cell if it happence
  if (i_n < 0) { i_n = 0; }
  if (k_n < 0) { k_n = 0; }
  if (i_o < 0) { i_o = 0; }
  if (k_o < 0) { k_o = 0; }

  if ((abs(x1_new-x1_old)<1e-15)&&(abs(x3_new-x3_old)<1e-15))
  {
    cout<<"WARNING! zero velocity!\n";
    return;
  }
  // stirct axis motion
  if (abs(x1_new-x1_old)<1e-15)
  {
    double r1=0, r2=0,r3=0;
    double delta_z = 0.0;
    double  value_part = 2.0*PI*x1_new*dr*dz;
    double wj_lower =0;
    r1 = x1_new-0.5*dr;
    r2 = (i_n+0.5)*dr;
    r3 = x1_new+0.5*dr;
    if (i_n==0)
    {
      wj_lower = charge_array[p_number]/(time1->delta_t*PI*dr*dr/4.0) * PI*(r2*r2-r1*r1)/value_part;
    }
    else
    {
      wj_lower = charge_array[p_number]/(time1->delta_t*2.0*PI*i_n*dr*dr) * PI*(r2*r2-r1*r1)/value_part;
    }
    double wj_upper =    charge_array[p_number]/(time1->delta_t*2*PI*(i_n+1)*dr*dr) *PI*(r3*r3-r2*r2)/value_part;
    double wj=0;
    this_j->set_j1(i_n,k_n,0.0);
    this_j->set_j1(i_n,k_n+1,0.0);
    int res_k = k_n-k_o;
    switch(res_k)
    {
    case 0:
    {
      delta_z = x3_new - x3_old;
      wj = wj_lower*delta_z;
      this_j->set_j3(i_n,k_n,wj);
      wj = wj_upper*delta_z;
      this_j->set_j3(i_n+1,k_n,wj);
    }
    break;

    case 1:
    {
      delta_z = k_n*dz - x3_old;
      wj = wj_lower*delta_z;
      this_j->set_j3(i_n,k_n-1,wj);
      wj = wj_upper*delta_z;
      this_j->set_j3(i_n+1,k_n-1,wj);

      delta_z = x3_new - k_n*dz;
      wj = wj_lower*delta_z;
      this_j->set_j3(i_n,k_n,wj);
      wj = wj_upper*delta_z;
      this_j->set_j3(i_n+1,k_n,wj);
    }
    break;

    case -1:
    {

      delta_z = (k_n+1)*dz - x3_old;
      wj = wj_lower*delta_z;
      this_j->set_j3(i_n,k_n+1,wj);
      wj = wj_upper*delta_z;
      this_j->set_j3(i_n+1,k_n+1,wj);

      delta_z = x3_new - (k_n+1)*dz;
      wj = wj_lower*delta_z;
      this_j->set_j3(i_n,k_n,wj);
      wj = wj_upper*delta_z;
      this_j->set_j3(i_n+1,k_n,wj);

    }
    break;
    }
  }

  // stirct radial motion
  else if (abs(x3_new-x3_old)<1e-15)
  {
    double r0    =(i_n+0.5)*dr;
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
      wj = charge_array[p_number]/(PI*4.0*r0*dz*dz*dr*time1->delta_t)*(delta_r - r0*r0/(x1_old+delta_r) +r0*r0/x1_old + dr*dr/(4.0*(x1_old+delta_r)) - dr*dr/(4.0*x1_old));
      res_j = wj*left_delta_z;
      this_j->set_j1(i_n,k_n,res_j);
      res_j = wj*right_delta_z;
      this_j->set_j1(i_n,k_n+1,res_j);

    }
    break;
    case 1:
    {
      delta_r = (i_n)*dr- x1_old;
      left_delta_z = (k_n+1)*dz-x3_new;
      right_delta_z = x3_new - k_n*dz;
      r0 = (i_n-0.5)*dr;
      wj = charge_array[p_number]/(PI*4.0*r0*dz*dz*dr*time1->delta_t)*(delta_r - r0*r0/(x1_old+delta_r) +r0*r0/x1_old + dr*dr/(4.0*(x1_old+delta_r)) - dr*dr/(4.0*x1_old));
      res_j = wj*left_delta_z;
      this_j->set_j1(i_n-1,k_n,res_j);
      res_j = wj*right_delta_z;
      this_j->set_j1(i_n-1,k_n+1,res_j);

      delta_r = x1_new - i_n*dr;
      r0 = (i_n+0.5)*dr;
      wj = charge_array[p_number]/(PI*4.0*r0*dz*dz*dr*time1->delta_t)*(delta_r - r0*r0/(i_n*dr+delta_r) +r0*r0/i_n*dr + dr*dr/(4.0*(i_n*dr+delta_r)) - dr*dr/(4.0*i_n*dr));
      res_j = wj*left_delta_z;
      this_j->set_j1(i_n,k_n,res_j);
      res_j = wj*right_delta_z;
      this_j->set_j1(i_n,k_n+1,res_j);
    }
    break;
    case -1:
    {
      delta_r = (i_n+1)*dr - x1_old ;
      left_delta_z = (k_n+1)*dz-x3_new;
      right_delta_z = x3_new - k_n*dz;
      r0 = (i_n+1.5)*dr;
      wj = charge_array[p_number]/(PI*4.0*r0*dz*dz*dr*time1->delta_t)*(delta_r - r0*r0/(x1_old+delta_r) +r0*r0/x1_old + dr*dr/(4.0*(x1_old+delta_r)) - dr*dr/(4.0*x1_old));
      res_j = wj*left_delta_z;
      this_j->set_j1(i_n+1,k_n,res_j);
      res_j = wj*right_delta_z;
      this_j->set_j1(i_n+1,k_n+1,res_j);

      delta_r = x1_new - (i_n+1)*dr;
      r0 = (i_n+0.5)*dr;
      wj = charge_array[p_number]/(PI*4.0*r0*dz*dz*dr*time1->delta_t)*(delta_r - r0*r0/((i_n+1)*dr+delta_r) +r0*r0/(i_n+1)*dr + dr*dr/(4.0*((i_n+1)*dr+delta_r)) - dr*dr/(4.0*(i_n+1)*dr));
      res_j = wj*left_delta_z;
      this_j->set_j1(i_n,k_n,res_j);
      res_j = wj*right_delta_z;
      this_j->set_j1(i_n,k_n+1,res_j);
    }
    break;
    }
  }
}

bool continuity_equation(Time *input_time,
                         Geometry *input_geometry,
                         current *input_J,
                         charge_density *rho_old,
                         charge_density *rho_new)
{
  double **rho_old_array = rho_old->get_ro() ;
  double **rho_new_array = rho_new->get_ro() ;
  double **J1 = input_J->get_j1() ;
  double **J3 = input_J->get_j3() ;
  // double delta_rho = 1.0/(input_geometry->dz*4.0*PI*input_geometry->dr*input_geometry->dr) ;
  bool ok = true;
  double res, tolerance = 1e-3 ;
  for (int i=1;i<input_geometry->n_grid_1-1;i++)

    for (int k=1;k<input_geometry->n_grid_2-1;k++)
    {
      res = (rho_new_array[i][k] - rho_old_array[i][k])/input_time->delta_t +(J3[i][k] - J3[i][k-1])/input_geometry->dz     +(J1[i][k] - J1[i-1][k])/input_geometry->dr + (J1[i][k] + J1[i-1][k])/(2.0*i*input_geometry->dr);
      if (res > tolerance)
      {
        ok = false;
        // i = input_geometry->n_grid_1;

      }
    }
  return ok;
}
//int* Particles::get_cell_numbers_jr(double x1_new,double x3_new, double x1_old, double x3_old)
/*
  int* Particles::get_cell_numbers_jr(double x1_new,double x3_new, double x1_old, double x3_old)
  {

  int i_new = (int)ceil((x1_new)/geom1->dr)-1;
  int k_new =(int)ceil((x3_new)/geom1->dz)-1;
  int i_old = (int)ceil((x1_old)/geom1->dr)-1;
  int k_old =(int)ceil((x3_old)/geom1->dz)-1;
  if (x1_old==(i_old+1)*geom1->dr)
  i_old=i_new;
  if(x3_old==(k_old+1)*geom1->dz)
  k_old=k_new;
  if (x1_new==(i_new+1)*geom1->dr)
  i_new=i_old;
  if(x3_new==(k_new+1)*geom1->dz)
  k_new=k_old;


  /// 4 cells

  if (i_new==i_old&&k_new==k_old)
  {
  int * res_cells = new int[4];

  res_cells[0]=i_new;
  res_cells[1]=k_new;
  res_cells[2]=i_new;
  res_cells[3]=k_new+1;

  }
  /// 7 cells
  if (i_new!=i_old||k_new!=k_old)
  {
  if(i_new!=i_old)
  {
  int i_new_max=0;
  if (i_new<i_old)
  int i_new_max=i_old;
  else
  i_new_max=i_new;
  int * res_cells = new int[8];
  res_cells[0]= i_new_max;
  res_cells[1]= k_new;
  res_cells[2]=i_new_max;
  res_cells[3]=k_new+1;
  res_cells[4]=i_new_max-1;
  res_cells[5]=k_new;
  res_cells[6]=i_new_max-1;
  res_cells[7]=k_new+1;

  }

  if(k_new!=k_old)
  {
  int k_new_max=0;
  if (k_new<k_old)
  int k_new_max=k_old;
  else
  k_new_max=k_new;

  int * res_cells = new int[6];
  res_cells[0]=i_new;
  res_cells[1]=k_new_max-1;
  res_cells[2]=i_new;
  res_cells[3]=k_new_max;
  res_cells[4]=i_new;
  res_cells[5]=k_new_max+1;

  }

  }
  // 10 cells

  {
  if((i_new>i_old)&&(k_new>k_old)||(i_new<i_old)&&(k_new<k_old))
  {
  int ik_temp=0;
  double x_temp= 0;
  if (i_new<i_old)
  {
  x_temp=x1_new;
  ik_temp=i_new;
  i_new=i_old;
  x1_new=x1_old;
  i_old=ik_temp;
  x1_old=x_temp;

  }


  if (k_new<k_old)
  {
  x_temp=x3_new;
  ik_temp=k_new;
  k_new=k_old;
  x3_new=x3_old;
  k_old=ik_temp;
  x3_old=x_temp;
  }
  double a = (x1_new - x1_old)/(x3_new - x3_old);
  double r1 = i_new*geom1->dr;
  double delta_z1 = (r1 - x1_old)/a;
  double z1 = x3_old + delta_z1;
  double z2 = k_new*geom1->dz;
  double delta_r2 = (z2-x3_old)*a;
  double r2 = x1_old+ delta_r2;
  if (z1<k_new*geom1->dz)
  {
  int * res_cells = new int[10];
  res_cells[0]=i_new-1;
  res_cells[1]=k_new-1;
  res_cells[2]=i_new-1;
  res_cells[3]=k_new;

  res_cells[4]=i_new;
  res_cells[5]=k_new-1;
  res_cells[6]=i_new;
  res_cells[7]=k_new;

  res_cells[8]=i_new;
  res_cells[9]=k_new+1;



  }
  else if (z1>k_new*geom1->dz)
  {
  int * res_cells = new int[10];
  res_cells[0]=i_new-1;
  res_cells[1]=k_new-1;

  res_cells[2]=i_new-1;
  res_cells[3]=k_new;


  res_cells[6]=i_new-1;
  res_cells[7]=k_new+1;

  res_cells[8]=i_new;
  res_cells[9]=k_new+1;

  }




  }

  }

  }
  void Particles::get_cell_numbers_jr_1(double x1_new,double x3_new, double x1_old, double x3_old,int *i_return, int* k_return,int* accur)
  {

  int i_new = (int)ceil((x1_new)/geom1->dr)-1;
  int k_new =(int)ceil((x3_new)/geom1->dz)-1;
  int i_old = (int)ceil((x1_old)/geom1->dr)-1;
  int k_old =(int)ceil((x3_old)/geom1->dz)-1;
  if (x1_old==(i_old+1)*geom1->dr)
  i_old=i_new;
  if(x3_old==(k_old+1)*geom1->dz)
  k_old=k_new;
  if (x1_new==(i_new+1)*geom1->dr)
  i_new=i_old;
  if(x3_new==(k_new+1)*geom1->dz)
  k_new=k_old;

  int cell_value = abs(i_new-i_old)+abs(k_new-k_old);

  /// 4 cells

  if (cell_value==0)
  {
  *i_return=i_new;
  *k_return=k_new;
  *accur=1;
  }
  /// 7 cells
  if (cell_value==1)
  {
  *i_return=i_new;
  *k_return=k_new;
  *accur=2;

  }
  // 10 cells

  {

  *i_return=i_new;
  *k_return=k_new;
  *accur=3;


  }

  }

  void Particles::get_cell_numbers_jr_2(double x1_new,double x3_new, double x1_old, double x3_old,int **cell_arr_jr, int **cell_arr_jz, int* number)
  {

  int i_new = (int)ceil((x1_new)/geom1->dr)-1;
  int k_new =(int)ceil((x3_new)/geom1->dz)-1;
  int i_old = (int)ceil((x1_old)/geom1->dr)-1;
  int k_old =(int)ceil((x3_old)/geom1->dz)-1;
  if (x1_old==(i_old+1)*geom1->dr)
  i_old=i_new;
  if(x3_old==(k_old+1)*geom1->dz)
  k_old=k_new;
  if (x1_new==(i_new+1)*geom1->dr)
  i_new=i_old;
  if(x3_new==(k_new+1)*geom1->dz)
  k_new=k_old;

  int cell_value = abs(i_new-i_old)+abs(k_new-k_old);

  /// 4 cells

  if (cell_value==0)
  {

  cell_arr_jr[0][0]=i_new;
  cell_arr_jr[0][1]=k_new;
  cell_arr_jr[1][0]=i_new;
  cell_arr_jr[1][1]=k_new+1;

  cell_arr_jz[0][0]=i_new;
  cell_arr_jz[0][1]=k_new;
  cell_arr_jz[1][0]=i_new+1;
  cell_arr_jz[1][1]=k_new;



  }
  /// 7 cells
  if (cell_value==1)
  {
  set_simple_cell(cell_arr_jr, cell_arr_jz, 0,i_old,k_old);
  set_simple_cell(cell_arr_jr, cell_arr_jz, 2,i_new,k_new);

  }
  // 10 cells

  {

  if (((i_new-i_old)+(k_new-k_old))!=0)
  {
  set_simple_cell(cell_arr_jr, cell_arr_jz, 0,i_old,k_old);
  set_simple_cell(cell_arr_jr, cell_arr_jz, 2,i_new,k_new);
  int i_temp= i_new;
  int k_temp = k_new;
  if (i_new <i_old)
  i_temp=i_old;
  if (k_new <k_old)
  k_temp=k_old;
  set_simple_cell(cell_arr_jr, cell_arr_jz, 4,i_temp,k_temp-1);
  set_simple_cell(cell_arr_jr, cell_arr_jz, 6,i_temp-1,k_temp);

  }
  else if (((i_new-i_old)+(k_new-k_old))==0)
  {
  set_simple_cell(cell_arr_jr, cell_arr_jz, 0,i_old,k_old);
  set_simple_cell(cell_arr_jr, cell_arr_jz, 2,i_new,k_new);
  int i_temp= i_new;
  int k_temp = k_new;
  if (i_new <i_old)
  i_temp=i_old;
  if (k_new >k_old)
  k_temp=k_old;
  set_simple_cell(cell_arr_jr, cell_arr_jz, 4,i_temp-1,k_temp);
  set_simple_cell(cell_arr_jr, cell_arr_jz, 6,i_temp,k_temp+1);

  }


  }

  }
  void Particles:: set_simple_cell(int** cell_arr_jr, int** cell_arr_jz,int start_number, int i_new, int k_new)
  {
  cell_arr_jr[0+start_number][0]=i_new;
  cell_arr_jr[0+start_number][1]=k_new;
  cell_arr_jr[1+start_number][0]=i_new;
  cell_arr_jr[1+start_number][1]=k_new+1;

  cell_arr_jz[0+start_number][0]=i_new;
  cell_arr_jz[0+start_number][1]=k_new;
  cell_arr_jz[1+start_number][0]=i_new+1;
  cell_arr_jz[1+start_number][1]=k_new;

  }

  void Particles::get_cell_numbers_new(double x1_new,double x3_new, double x1_old, double x3_old, int*i_min, int* i_max,int* k_min, int* k_max)
  {

  int i_new = (int)ceil((x1_new)/geom1->dr)-1;
  int k_new =(int)ceil((x3_new)/geom1->dz)-1;
  int i_old = (int)ceil((x1_old)/geom1->dr)-1;
  int k_old =(int)ceil((x3_old)/geom1->dz)-1;
  if (x1_old==(i_old+1)*geom1->dr)
  i_old=i_new;
  if(x3_old==(k_old+1)*geom1->dz)
  k_old=k_new;
  if (x1_new==(i_new+1)*geom1->dr)
  i_new=i_old;
  if(x3_new==(k_new+1)*geom1->dz)
  k_new=k_old;

  int cell_value = abs(i_new-i_old)+abs(k_new-k_old);

  /// 4 cells

  if (cell_value==0)
  {
  *i_min = i_new;
  *i_max  =i_new+1;
  *k_min = k_new;
  *k_min = k_new+1;


  }
  /// 7 cells
  if (cell_value==1)
  {
  if (i_new!=i_old)
  {
  int i_m = i_new;
  if (i_new<i_old)
  {
  i_m=i_old;
  }
  *i_min = i_m-1;
  *i_max  =i_m+1;
  *k_min = k_new;
  *k_max = k_new+1;

  }
  else {

  int k_m = k_new;
  if (k_new<k_old)
  {
  k_m=k_old;
  }
  *i_min = i_new;
  *i_max  =i_new+1;
  *k_min = k_m-1;
  *k_max = k_m+1;



  }


  }
  // 10 cells

  {

  if (((i_new-i_old)+(k_new-k_old))!=0)
  {
  *i_max=i_new+1;
  *i_min = i_new-1;
  if (i_new<i_old)
  {
  *i_max=i_old+1;
  * i_min = i_old-1;
  }

  *k_max=k_new+1;
  *k_min=k_new-1;
  if (i_new<i_old)
  {
  *k_max=k_old+1;
  *k_min = k_old-1;
  }

  }

  }
*/
