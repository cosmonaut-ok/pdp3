#include "particles.h"
#include "eField.h"
#include "hField.h"
#include "pdp3Time.h"
#include "triple.h"
#include <math.h>
#include <string.h>
#include <cstdlib>
#include "constant.h"

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
                     Geometry* geom):geom1(geom)
{
  name = (char *) p_name;
  charge = (double)p_charge*(double)EL_CHARGE;
  init_const_mass = p_mass;
  mass = (double)p_mass*(double)EL_MASS;
  number = p_number;

  // allocate memory for coordinates and velocities of particles
  mass_array = new double[number];
  charge_array = new double[number];
  x1 = new double[number];
  x3 = new double[number];
  v1 = new double[number];
  v2 = new double[number];
  v3 = new double[number];
  is_alive = new int[number];

#pragma omp parallel for
  for (int i = 0; i < number; i++)
    is_alive[i] = 1;
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

#pragma omp parallel for
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
#pragma omp parallel for
  for(int i=0; i<number; i++)
  {
    v1[i] = 0;
    v2[i] = 1e5; // TODO: WHY?
    v3[i] = 0;
  }
}

void Particles::set_x_0()
{
#pragma omp parallel for
  for(int i=0; i<number; i++)
  {
    x1[i] = 0.5; // TODO: WHY?
    x3[i] = 0.5; // TODO: WHY?
  }
}

// calculate Lorentz factor
double Particles::get_gamma(int i)
{
  double gamma, beta;

  beta = (pow(v1[i], 2) + pow(v2[i], 2) + pow(v3[i], 2)) / LIGHT_SPEED_POW_2;

  if (beta > 1) // it's VERY BAD! Beta should not be more, than 1
  {
    cerr << "CRITICAL!(get_gamma): Lorentz factor aka gamma is comples. Can not continue." << endl
         << "\tUsually it happens, when <Time> -> <delta_t> value is too big." << endl
         << "\tv1[" << i << "] = " << v1[i]
         << "\tv2[" << i << "] = " << v2[i]
         << "\tv3[" << i << "] = " << v3[i];
    exit(1);
  }

  gamma = pow(1.0 - beta, -0.5);

  if (isinf(gamma) == 1) { // avoid infinity values
    cerr << "WARNING(get_gamma): gamma (Lorenz factor) girects to infinity for [v1, v2, v3, i]: ["
         << v1[i] << ", " << v2[i] << ", " << v3[i] << ", " << i << "]" << endl;
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
  double gamma, beta;

  beta = (pow(v1[i], 2) + pow(v2[i], 2) + pow(v3[i], 2)) / LIGHT_SPEED_POW_2;

  if (beta > 1) // it's VERY BAD! Beta should not be more, than 1
  {
    cerr << "CRITICAL!(get_gamma_inv): Lorentz factor aka gamma is comples. Can not continue." << endl
         << "\tUsually it happens, when <Time> -> <delta_t> value is too big." << endl
         << "\tv1[" << i << "] = " << v1[i]
         << "\tv2[" << i << "] = " << v2[i]
         << "\tv3[" << i << "] = " << v3[i];
    exit(1);
  }

  gamma = pow(1.0 + beta, 0.5);

  if (isinf(gamma) == 1) { // avoid infinity values
    cerr << "WARNING(get_gamma_inv): reciprocal gamma (Lorenz factor) directs to infinity for [v1, v2, v3, i]: ["
         << v1[i] << ", " << v2[i] << ", " << v3[i] << ", " << i << "]" << endl;
    return 1e100; // just return some very big value
  }
  else
  {
    return gamma;
  }
}

void Particles::step_v(EField *e_fld, HField *h_fld, Time* t)
{

#pragma omp parallel for shared(e_fld, h_fld, t)
  for(int i=0; i<number; i++)
    if (is_alive[i])
    {
      // define vars directly in cycle, because of multithreading
      double gamma, b1, b2, b3, e1, e2, e3, vv1, vv2, vv3, const1, const2;
      Triple E_compon(0.0, 0.0, 0.0), B_compon(0.0, 0.0, 0.0);

      // check if x1 and x3 are correct
      if (isnan(x1[i]) || isinf(x1[i]) != 0 || isnan(x3[i]) || isinf(x3[i]) != 0)
      {
        cerr << "ERROR(step_v): x1[" << i << "] or x3[" << i << "] is not valid number. Can not continue." << endl;
        exit(1);
      }
      // q*t/2*m TODO: what constant is it?
      const1 = charge_array[i]*t->delta_t/2.0/mass_array[i];
      E_compon = e_fld->get_field(x1[i],x3[i]);
      B_compon = h_fld->get_field(x1[i],x3[i]);
      e1 = E_compon.first*const1;
      e2 = E_compon.second*const1;
      e3 = E_compon.third*const1;
      b1 = B_compon.first*MAGN_CONST*const1;
      b2 = B_compon.second*MAGN_CONST*const1;
      b3 = B_compon.third*MAGN_CONST*const1;

      // round very small velicities to avoid exceptions
      vv1 = (abs(v1[i]) < 1e-15) ? 0 : v1[i];
      vv2 = (abs(v2[i]) < 1e-15) ? 0 : v2[i];
      vv3 = (abs(v3[i]) < 1e-15) ? 0 : v3[i];

      // 1. Multiplication by relativistic factor
      // u(n-1/2) = gamma(n-1/2)*v(n-1/2)
      gamma = get_gamma(i);
      //
      v1[i] = gamma*vv1;
      v2[i] = gamma*vv2;
      v3[i] = gamma*vv3;

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
// #pragma omp parallel for shared(dr, dz, x1_wall, x3_wall, half_dr, half_dz, x1_wallX2, x3_wallX2, half_dt)
#pragma omp parallel for
  for(int i=0; i<number; i++)
    if (is_alive[i])
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
void Particles::charge_weighting(ChargeDensity* ro1)
{
  int r_i = 0;  // number of particle i cell
  int z_k = 0;  // number of particle k cell

  double dr = geom1->dr;
  double dz = geom1->dz;
  double r1, r2, r3; // temp variables for calculation
  double dz1, dz2; // temp var.: width of k and k+1 cell

  double ro_v = 0; // charge density Q/V, V - volume of particle
  double v_1 = 0; // volume of [i][k] cell
  double v_2 = 0; // volume of [i+1][k] cell
  //double ro_v_2 = 0; // charge density in i+1 cell

  double value = 0;
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

void Particles::velocity_distribution(double tempr_ev)
{
  double therm_vel = sqrt(tempr_ev*2.0*EL_CHARGE/
                          (this->init_const_mass*EL_MASS));
  // double R =0; // number from [0;1]
	// TODO: why 5e6?
  double dv = therm_vel/0.5e7; // velocity step in calculation integral
	// TODO: why 9.0?
  double cutoff_vel = 9.0*therm_vel; // cutoff velocity
  int lenght_arr = (int)cutoff_vel/dv;
  double s = 0;

  double* integ_array = new double [lenght_arr];

  double const1 = 2*therm_vel*therm_vel;
  // part of numerical integral calculation
  for (int i=0; i<lenght_arr; i++)
  {
		double ds = 0;
    ds = exp(-dv*i*dv*i/const1)*dv;
    s = s+ds;
    integ_array[i] = s;
  }

#pragma omp parallel for shared(integ_array, lenght_arr, const1, dv)
  for(int i_n=0; i_n<number; i_n++)
  {
    double Rr = random_reverse(i_n,3);
    double Rfi = random_reverse(i_n,5);
    double Rz = random_reverse(i_n,7);
    double t_z = sqrt(PI/2.0)*therm_vel;
    // R = rand()/(double)32768;
    double f_vr = Rr*t_z;
    double f_vfi = Rfi*t_z;
    double f_vz = Rz*t_z;
    int sign = 1;
    if (i_n%2 == 1)
      sign =-1;

    // binary search
    int i = 0;
    int j = lenght_arr;
    int k = 0;
    while(i <= j)
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

    i = 0;
    j = lenght_arr;
    k = 0;
    while(i <= j)
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
  }
  delete []integ_array;
}

double Particles::random_reverse(double vel, int power)
{
  int int_vel =(int) floor(vel);
  double ost = 0;
  double r = 0;
  int order = 1;
  while(int_vel >= 1)
  {
    ost = int_vel % power;

    r = r + ost*pow((double)power,(-order));

    int_vel = (int_vel - ost)/power;
    order = order+1;
  }
  return r;
}

void Particles::load_spatial_distribution(double n1, double n2, double left_plasma_boundary, int type)
{
  // calculate number of electrons in a big particle
  double rand_r;
  double rand_z;
  double dr = geom1->dr*1.00000001;// TODO: why?
  double dz = geom1->dz*1.00000001;// TODO: why?
  double dn = n2 - n1;
  switch (type)
  {
  case 0:
  {
    double n_in_big = (PI*geom1->first_size*geom1->first_size*geom1->second_size/number*(n2+n1)/2.0);
    charge *= n_in_big;
    mass *= n_in_big;
#pragma omp parallel for shared (dr, dz, dn, left_plasma_boundary, n1)
    for(int n = 0; n < number; n++)
    {
      rand_r = random_reverse(n, 13); // TODO: why 11 and 13?
      rand_z = random_reverse(number - 1 - n, 11);
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
#pragma omp parallel for shared (sigma, R_sq, dz)
    for (int i = 0; i<number;i++)
    {
      rand_r = random_reverse(i,13);
      double int_rd =   exp(-dr*dr/(8.0*sigma*sigma));
      x1[i]=sigma*sqrt(-2.0*log(int_rd - rand_r*(int_rd-exp(-R_sq/(2.0*sigma*sigma)))));

      //x1[i] = (geom1->first_size - dr)*(rand_r)*rand_r + dr/2.0;
      // double tt = exp(-R_sq/(2.0*sigma*sigma));
      rand_z = random_reverse(number - 1 - i, 11);
      x3[i] = (geom1->second_size - dz/2.0)*rand_z + dz/2.0;
      //x3[i] = (geom1->second_size - left_plasma_boundary - dz)/dn*(sqrt(n1*n1 + rand_z*(2*n1*dn + dn*dn)) - n1) + left_plasma_boundary + dz/2.0;
    }
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
  double dr = geom1->dr*1.00000001; // TODO: why?
  double dz = geom1->dz*1.00000001; // TODO: why?
  double dn = n2 - n1;
  switch (type)
  {
  case 0:
  {

    double N_big_for_cell=(double) number/( (double) geom1->n_grid_1*geom1->n_grid_2);
    double N_real_i = 8.0*PI*(n2+n1)/2.0*dr*dz;
    double n_in_big =0;

#pragma omp parallel for shared(dr, dz, dn, left_plasma_boundary, n1, n2) private(rand_r, rand_z, n_in_big)
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

void Particles::simple_j_weighting(Time* time1,
                                   Current *j1,
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
    // equation y = k*x+b;//
    // finding k & b//
    double k = delta_r/delta_z;
    double b = x1_old;
    //calculate current jz in [i,k] cell//
    wj = charge_array[p_number]/(2*dr*dz*delta_t*2*PI*i_n*dr*dr) *
      (dr*delta_z - k*delta_z*delta_z/2.0 - delta_z*b + dr*dr/k *
       ((i_n+0.5)*(i_n+0.5)-0.25)*log((k*delta_z+b)/b));
    // set new weighting current value
    j1->set_j3(i_n,k_n, wj);

    //calculate current in [i+1,k] cell//
    wj = charge_array[p_number]/(2*dr*dz*delta_t*2*PI*(i_n+1)*dr*dr) *
      (k*delta_z*delta_z/2.0+delta_z*b + delta_z*dr + dr*dr/k *
       (0.25-(i_n+0.5)*(i_n+0.5)) * log((k*delta_z+b)/b));
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
    wj = charge_array[p_number]/(2*PI*r0*dz*dz*dr*delta_t) *
      (r0*k*delta_r+k/2.0 * delta_r*(x1_old+delta_r/2.0)+0.5*delta_r*
       (b-k*(2*r0+r1))+ delta_r*(b-k*r1)*(4*r0*r0-dr*dr)/(8*x1_old*(x1_old+delta_r)) +
       (k*(r0*r0/2.0-dr*dr/8.0))*log((x1_old+delta_r)/x1_old));
    j1->set_j1(i_n,k_n, wj);

    b= x3_old- k_n*dz;;
    //weighting jr in [i][k+1] cell
    wj = charge_array[p_number]/(2*PI*r0*dz*dz*dr*delta_t) *
      (-r0*k*delta_r - k/2.0*delta_r*(x1_old+delta_r/2.0)+0.5*delta_r*
       (b+k*(2*r0+r1))+ delta_r*(b+k*r1)*(4*r0*r0-dr*dr)/(8*x1_old*(x1_old+delta_r)) -
       (k*(r0*r0/2.0-dr*dr/8.0))*log((x1_old+delta_r)/x1_old));
    j1->set_j1(i_n, k_n+1, wj);
  }
  // if i cell is equal 0
  else
  {
    // equation y = k*x+b;//
    // finding k & b//
    double k = delta_r/delta_z;
    double b = x1_old;
    //calculate current jz in [i,k] cell//
    wj = charge_array[p_number]/(2.0*dr*dz*delta_t*PI*dr*dr/4.0) *
      (dr*delta_z - k*delta_z*delta_z/2.0 - delta_z*b );
    // set new weighting current value
    j1->set_j3(i_n,k_n, wj);

    //calculate current in [i+1,k] cell//
    wj = charge_array[p_number]/(2.0*dr*dz*delta_t*2.0*PI*dr*dr) *
      (k*delta_z*delta_z/2.0 + delta_z*dr +delta_z*b);
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
    wj = charge_array[p_number]/(2*PI*r0*dz*dz*dr*delta_t) *
      (r0*k*delta_r+k/2.0 * delta_r*(x1_old+delta_r/2.0)+0.5*delta_r*(b-k*(2*r0+r1))+
       delta_r*(b-k*r1)*(4*r0*r0-dr*dr)/(8*x1_old*(x1_old+delta_r)) +
       (k*(r0*r0/2.0-dr*dr/8.0))*log((x1_old+delta_r)/x1_old));
    j1->set_j1(i_n,k_n, wj);

    b= x3_old- k_n*dz;;
    //weighting jr in [i][k+1] cell
    wj = charge_array[p_number]/(2*PI*r0*dz*dz*dr*delta_t) *
      (-r0*k*delta_r - k/2.0*delta_r*(x1_old+delta_r/2.0)+0.5*delta_r*(b+k*(2*r0+r1))+
       delta_r*(b+k*r1)*(4*r0*r0-dr*dr)/(8*x1_old*(x1_old+delta_r)) -
       (k*(r0*r0/2.0-dr*dr/8.0))*log((x1_old+delta_r)/x1_old));
    j1->set_j1(i_n, k_n+1, wj);
  }

  //}
}
void Particles::simple_constrho_j_weighting(Time* time1,
                                            Current *j1,
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

  // if i cell is not equal 0
  if (i_n>=1)
  {
    j1->set_j3(i_n,k_n, wj);

    //calculate current in [i+1,k] cell//
    j1->set_j3(i_n+1,k_n, wj);

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
    // calculate current jz in [i,k] cell
    // wj = charge/(2.0*dr*dz*delta_t*PI*dr*dr/4.0) * (dr*delta_z - k*delta_z*delta_z/2.0 - delta_z*b );
    // set new weighting current value
    j1->set_j3(i_n,k_n, wj);

    // calculate current in [i+1,k] cell
    // wj = charge/(2.0*dr*dz*delta_t*2.0*PI*dr*dr) * (k*delta_z*delta_z/2.0 + delta_z*dr +delta_z*b);
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
    wj = charge_array[p_number]/(2*PI*r0*dz*dz*dr*delta_t) *
      (r0*k*delta_r+k/2.0 * delta_r*(x1_old+delta_r/2.0)+0.5*delta_r*(b-k*(2*r0+r1))+
       delta_r*(b-k*r1)*(4*r0*r0-dr*dr)/(8*x1_old*(x1_old+delta_r)) +
       (k*(r0*r0/2.0-dr*dr/8.0))*log((x1_old+delta_r)/x1_old));
    j1->set_j1(i_n,k_n, wj);

    b = x3_old - k_n*dz;
    //weighting jr in [i][k+1] cell
    wj = charge_array[p_number]/(2*PI*r0*dz*dz*dr*delta_t) *
      (-r0*k*delta_r - k/2.0*delta_r*(x1_old+delta_r/2.0)+0.5*delta_r*(b+k*(2*r0+r1))+
       delta_r*(b+k*r1)*(4*r0*r0-dr*dr)/(8*x1_old*(x1_old+delta_r)) -
       (k*(r0*r0/2.0-dr*dr/8.0))*log((x1_old+delta_r)/x1_old));
    j1->set_j1(i_n, k_n+1, wj);
  }
}

void Particles::j_weighting(Time* time1, Current *j1, double* x1_o,double* x3_o)
{
  double dr = geom1->dr;
  double dz = geom1->dz;

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

        ///////// 3) charge in 10 cells /////////
        case 2:
        {
          // case, when particle move from [i-1] cell to [i] cell
          if (i_o<i_n)
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
void Particles::azimuthal_j_weighting(Time* time1, Current *this_j)
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
                                        Current *this_j,
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
    cerr<<"WARNING! zero velocity!" << endl;
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
      wj_lower = charge_array[p_number]/(time1->delta_t*PI*dr*dr/4.0) * PI*(r2*r2-r1*r1)/value_part;
    else
      wj_lower = charge_array[p_number]/(time1->delta_t*2.0*PI*i_n*dr*dr) * PI*(r2*r2-r1*r1)/value_part;
    double wj_upper = charge_array[p_number]/(time1->delta_t*2*PI*(i_n+1)*dr*dr) * PI*(r3*r3-r2*r2)/value_part;
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
      wj = charge_array[p_number]/(PI*4.0*r0*dz*dz*dr*time1->delta_t)*
        (delta_r - r0*r0/(x1_old+delta_r) + r0*r0/x1_old + dr*dr/(4.0*(x1_old+delta_r)) -
         dr*dr/(4.0*x1_old));
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
      wj = charge_array[p_number]/(PI*4.0*r0*dz*dz*dr*time1->delta_t) *
        (delta_r - r0*r0/(x1_old+delta_r) +r0*r0/x1_old + dr*dr/(4.0*(x1_old+delta_r)) -
         dr*dr/(4.0*x1_old));
      res_j = wj*left_delta_z;
      this_j->set_j1(i_n-1,k_n,res_j);
      res_j = wj*right_delta_z;
      this_j->set_j1(i_n-1,k_n+1,res_j);

      delta_r = x1_new - i_n*dr;
      r0 = (i_n+0.5)*dr;
      wj = charge_array[p_number]/(PI*4.0*r0*dz*dz*dr*time1->delta_t) *
        (delta_r - r0*r0/(i_n*dr+delta_r) +r0*r0/i_n*dr + dr*dr/(4.0*(i_n*dr+delta_r)) -
         dr*dr/(4.0*i_n*dr));
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
      wj = charge_array[p_number]/(PI*4.0*r0*dz*dz*dr*time1->delta_t) *
        (delta_r - r0*r0/(x1_old+delta_r) +r0*r0/x1_old + dr*dr/(4.0*(x1_old+delta_r)) -
         dr*dr/(4.0*x1_old));
      res_j = wj*left_delta_z;
      this_j->set_j1(i_n+1,k_n,res_j);
      res_j = wj*right_delta_z;
      this_j->set_j1(i_n+1,k_n+1,res_j);

      delta_r = x1_new - (i_n+1)*dr;
      r0 = (i_n+0.5)*dr;
      wj = charge_array[p_number]/(PI*4.0*r0*dz*dz*dr*time1->delta_t) *
        (delta_r - r0*r0/((i_n+1)*dr+delta_r) +r0*r0/(i_n+1)*dr + dr*dr/(4.0*((i_n+1)*dr+delta_r)) -
         dr*dr/(4.0*(i_n+1)*dr));
      res_j = wj*left_delta_z;
      this_j->set_j1(i_n,k_n,res_j);
      res_j = wj*right_delta_z;
      this_j->set_j1(i_n,k_n+1,res_j);
    }
    break;
    }
  }
}
