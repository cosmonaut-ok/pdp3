#include "particles.h"

using namespace std;
using namespace constant;
using namespace math;
using namespace math::random;

// use classical calculations, if velocity lower, than minimal
#ifndef REL_LIMIT
#define REL_LIMIT 5e7
#endif
// define REL_LIMIT^2 to decrease number of operations
const double REL_LIMIT_POW_2 = pow (REL_LIMIT, 2);

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
                     Geometry *geom):geom1(geom)
{
  name = (char*) p_name;
  charge = p_charge * EL_CHARGE;
  init_const_mass = p_mass;
  mass = p_mass * EL_MASS;
  number = p_number;

  // allocate memory for position and velocity of particles
  mass_array = new double[number];
  charge_array = new double[number];
  vel = new double *[number];
  pos = new double *[number];
  pos_old = new double *[number];

  is_alive = new bool[number];

  sin_theta_r = new double[number];
  cos_theta_r = new double[number];

#pragma omp parallel for
  for (unsigned int i = 0; i < number; i++)
  {
    is_alive[i] = true;
    vel[i] = new double[3];
    pos[i] = new double[3];
    pos_old[i] = new double[3];
  }
}

// Destructor
Particles::~Particles()
{
  for (unsigned int i = 0; i < number; i++)
  {
    delete [] pos[i];
    delete [] pos_old[i];
    delete [] vel[i];
  }

  delete [] pos;
  delete [] pos_old;
  delete [] vel;
  delete [] is_alive;

  delete [] mass_array;
  delete [] charge_array;
  delete [] sin_theta_r;
  delete [] cos_theta_r;
}

void Particles::half_step_pos(Time *t)
{
  double half_dt = t->delta_t / 2.;

#pragma omp parallel for shared (half_dt)
  for(unsigned int i = 0; i < number; i++)
    if (is_alive[i])
    {
      // check if radius and longitude are correct
      if (isnan(pos[i][0]) || isinf(pos[i][0]) != 0 || isnan(pos[i][2]) || isinf(pos[i][2]) != 0)
      {
        cerr << "ERROR(half_step_pos): radius[" << i << "] or longitude[" << i << "] is not valid number. Can not continue." << endl;
        exit(1);
      }
      //
      pos[i][0] = pos[i][0] + vel[i][0] * half_dt;
      //! we use "fake" rotation component to correct position from xy to rz pane
      pos[i][1] = pos[i][1] + vel[i][1] * half_dt;
      pos[i][2] = pos[i][2] + vel[i][2] * half_dt;
    }
}

// function for charge density weighting
void Particles::charge_weighting(ChargeDensity *ro1)
{
  int r_i = 0; // number of particle i cell
  int z_k = 0; // number of particle k cell

  double dr = geom1->dr;
  double dz = geom1->dz;
  double r1, r2, r3; // temp variables for calculation
  double dz1, dz2; // temp var.: width of k and k + 1 cell

  double ro_v = 0; // charge density Q/V, V - volume of particle
  double v_1 = 0; // volume of [i][k] cell
  double v_2 = 0; // volume of [i + 1][k] cell
  //double ro_v_2 = 0; // charge density in i + 1 cell

  double value = 0;
  // double **temp = ro1->get_rho();

  for(unsigned int i = 0; i < number; i++)
    if (is_alive[i])
    {
      // check if radius and longitude are correct
      if (isnan(pos[i][0]) || isinf(pos[i][0]) != 0 || isnan(pos[i][2]) || isinf(pos[i][2]) != 0)
      {
        cerr << "ERROR(charge_weighting): radius[" << i << "] or longitude[" << i << "] is not valid number. Can not continue." << endl;
        exit(1);
      }

      // finding number of i and k cell. example: dr = 0.5; r = 0.4; i =0
      r_i = CELL_NUMBER(pos[i][0], dr);
      z_k = CELL_NUMBER(pos[i][2], dz);
      // TODO: workaround: sometimes it gives -1.
      // Just get 0 cell if it happence
      if (r_i < 0) r_i = 0;
      if (z_k < 0) z_k = 0;

      // in first cell other alg. of ro_v calc
      if(pos[i][0]>dr)
      {
        r1 =  pos[i][0] - 0.5 * dr;
        r2 = (r_i + 0.5) * dr;
        r3 = pos[i][0] + 0.5 * dr;
        ro_v = charge_array[i] / (2. * PI * dz * dr * pos[i][0]);
        v_1 = CELL_VOLUME(r_i, dr, dz);
        v_2 = CELL_VOLUME(r_i + 1, dr, dz);
	dz1 = (z_k + 0.5) * dz - (pos[i][2] - 0.5 * dz);
	dz2 = (pos[i][2] + 0.5 * dz) - (z_k + 0.5) * dz;

        // weighting in ro[i][k] cell
        value = ro_v * CYL_RNG_VOL(dz1, r1, r2) / v_1;
        ro1->inc_rho(r_i, z_k, value);

        // weighting in ro[i + 1][k] cell
        value = ro_v * CYL_RNG_VOL(dz1, r2, r3) / v_2;
        ro1->inc_rho(r_i + 1,z_k, value);

        // weighting in ro[i][k + 1] cell
        value = ro_v * CYL_RNG_VOL(dz2, r1, r2) / v_1;
        ro1->inc_rho(r_i, z_k + 1, value);

        // weighting in ro[i + 1][k + 1] cell
        value = ro_v * CYL_RNG_VOL(dz2, r2, r3) / v_2;
        ro1->inc_rho(r_i + 1, z_k + 1, value);

      }
      else if (pos[i][0] <= dr / 2.)
      {
        r_i = 0;
        r1 =  0.;
        r2 = (r_i + 0.5) * dr;
        r3 = pos[i][0] + 0.5 * dr;
	dz1 = (z_k + 0.5) * dz - (pos[i][2] - 0.5 * dz);
	dz2 = (pos[i][2] + 0.5 * dz) - (z_k + 0.5) * dz;
        ro_v = charge_array[i] / (PI * dz * (2. * pos[i][0] * pos[i][0] + dr * dr / 2.));
        v_1 = CYL_VOL(dz, dr);
        v_2 = CELL_VOLUME(r_i + 1, dr, dz);
        ////////////////////////// /

        // weighting in ro[i][k] cell
        value = ro_v * PI * dz1 * (dr * dr / 2.-pos[i][0] * dr + pos[i][0] * pos[i][0]) / v_1;
        ro1->inc_rho(r_i, z_k, value);

        // weighting in ro[i + 1][k] cell
        value = ro_v * CYL_RNG_VOL(dz1, r2, r3) / v_2;
        ro1->inc_rho(r_i + 1,z_k, value);

        // weighting in ro[i][k + 1] cell
        value = ro_v*PI * dz2 * (dr * dr / 2.-pos[i][0] * dr + pos[i][0] * pos[i][0]) / v_1;
        ro1->inc_rho(r_i, z_k + 1, value);

        // weighting in ro[i + 1][k + 1] cell
        value = ro_v * CYL_RNG_VOL(dz2, r2, r3) / v_2;
        ro1->inc_rho(r_i + 1, z_k + 1, value);

      }
      else
      {
        ////////////////////////// /
        r1 = pos[i][0] - 0.5 * dr;
        r2 = (r_i + 0.5) * dr;
        r3 = pos[i][0] + 0.5 * dr;
	dz1 = (z_k + 0.5) * dz - (pos[i][2] - 0.5 * dz);
	dz2 = (pos[i][2] + 0.5 * dz) - (z_k + 0.5) * dz;
        ro_v = charge_array[i] / (2. * PI * dz * dr * pos[i][0]);
        v_1 = CYL_VOL(dz, dr);
        v_2 = CELL_VOLUME(r_i + 1, dr, dz);
        ////////////////////////// /

        // weighting in ro[i][k] cell
        value = ro_v * CYL_RNG_VOL(dz1, r1, r2) / v_1;
        ro1->inc_rho(r_i, z_k, value);

        // weighting in ro[i + 1][k] cell
        value = ro_v * CYL_RNG_VOL(dz1, r2, r3) / v_2;
        ro1->inc_rho(r_i + 1,z_k, value);

        // weighting in ro[i][k + 1] cell
        value = ro_v * CYL_RNG_VOL(dz2, r1, r2) / v_1;
        ro1->inc_rho(r_i, z_k + 1, value);

        // weighting in ro[i + 1][k + 1] cell
        value = ro_v * CYL_RNG_VOL(dz2, r2, r3) / v_2;
        ro1->inc_rho(r_i + 1, z_k + 1, value);
      }
    }
}

void Particles::velocity_distribution(double tempr_ev)
{
  double therm_vel = lib::sq_rt(2. * EL_CHARGE * tempr_ev / mass);

  // Sample the energies in the MJ distribution
  vector<double> energies = math::maxwell_juttner::maxwellJuttner(number, tempr_ev);

  for (unsigned int p=0; p < number; p++)
  {
    double therm_vel_el = lib::sq_rt(2 * EL_CHARGE * energies[p] / mass);

    // velocity for components. Require to normalize values
    double norm_coeff = 1; // 0.7071067811865475;
    therm_vel_el = therm_vel_el * norm_coeff;
      
    double rnd_0, rnd_1, rnd_2;

#ifdef TESTMODE
    rnd_0 = random_reverse(p, 3);
    rnd_1 = random_reverse(p, 5);
    rnd_2 = random_reverse(p, 7);

    // emulate negative numbers
    double k0 = random_reverse(p, 4);
    double k1 = random_reverse(p, 6);
    double k2 = random_reverse(p, 8);

    if (k0 > 0.5) rnd_0 = -rnd_0;
    if (k1 > 0.5) rnd_1 = -rnd_1;
    if (k2 > 0.5) rnd_2 = -rnd_2;
#else
    rnd_0 = uniform2();
    rnd_1 = uniform2();
    rnd_2 = uniform2();
#endif

    vel[p][0] = rnd_0 * therm_vel_el;
    vel[p][1] = rnd_1 * therm_vel_el;
    vel[p][2] = rnd_2 * therm_vel_el;

    // take into account relativistic factor
    // maxwellJuttner procedure returns relativistic momentums
    vel[p][0] *= lib::get_gamma_inv(vel[p][0] * vel[p][0]);
    vel[p][1] *= lib::get_gamma_inv(vel[p][1] * vel[p][1]);
    vel[p][2] *= lib::get_gamma_inv(vel[p][2] * vel[p][2]);
  }
}

void Particles::load_cylindrical_spatial_distribution(double n1, double n2, double left_plasma_boundary)
// ! calculate number of charged physical particles in macroparticle
{
  double dr = geom1->dr;
  double dz = geom1->dz;
  double dn = n2 - n1;

  // summary volume of all macroparticles
  double v_sum = 0;

  for (unsigned int n = 0; n < number; n++)
  {
    // very flat distribution, instead of random rectangle distribution
    double rand_r = random_reverse(n, 13);
    double rand_z = random_reverse(number - 1 - n, 11);
    // double rand_r = uniform();
    // double rand_z = uniform();

    pos[n][0] = (geom1->r_size - dr) * rand_r + dr / 2;
    pos[n][1] = 0;
    pos[n][2] = (geom1->z_size - left_plasma_boundary - dz)
      / dn * (lib::sq_rt(pow(n1, 2) + rand_z * (2 * n1 * dn + pow(dn, 2))) - n1)
      + left_plasma_boundary + dz / 2;

    v_sum += 2 * PI * pos[n][0] * dr * dz;
  }

  // average volume of single macroparticle
  double v_avg = v_sum / number;
  double N_total = (n1 + n2) / 2 * PI * geom1->r_size * geom1->r_size * geom1->z_size;

  double n_per_macro_avg = N_total / number;

  for (unsigned int n = 0; n < number; n++)
  {
    // coefitient of normalization
    double norm = 2 * PI * pos[n][0] * dr * dz / v_avg;

    // number of real particles per macroparticle
    double n_per_macro = n_per_macro_avg * norm;

    // set charge and mass of macroparticle
    charge_array[n] = charge * n_per_macro;
    mass_array[n] = mass * n_per_macro;
  }
}

void Particles::simple_current_distribution(Time *time1,
                                   Current *el_current,
                                   double radius_new,
                                   double longitude_new,
                                   double radius_old,
                                   double longitude_old,
                                   int i_n,
                                   int k_n,
                                   int p_number)
{
  double dr = geom1->dr;
  double dz = geom1->dz;
  double wj = 0;
  double delta_t = time1->delta_t;

  // distance of particle moving //
  double delta_r = radius_new - radius_old;
  double delta_z = longitude_new - longitude_old;

  double some_shit_density;

  if ((abs(delta_r) < MNZL) || (abs(delta_z) < MNZL)) // MNZL see constant.h
    return;
  // if i cell is not equal 0
  if (i_n>=1)
  {
    // equation y = k*x + b; //
    // finding k & b //
    double k = delta_r / delta_z;
    double b = radius_old;

    // calculate current jz in [i,k] cell
    some_shit_density = SOME_SHIT_DENSITY_R(charge_array[p_number], i_n * dr, dr, dz, delta_t);

    wj = some_shit_density
      * (dr * delta_z - k * delta_z * delta_z / 2. - delta_z * b + dr * dr / k
         * ((i_n + 0.5) * (i_n + 0.5) - 0.25) * log((k * delta_z + b) / b));
    // set new weighting current value
    el_current->inc_j_z(i_n, k_n, wj);

    some_shit_density = SOME_SHIT_DENSITY_R(charge_array[p_number], (i_n + 1) * dr, dr, dz, delta_t);

    // calculate current in [i+1,k] cell //
    wj = some_shit_density
      * (k * delta_z * delta_z / 2. + delta_z * b + delta_z * dr + dr * dr / k *
         (0.25-(i_n + 0.5) * (i_n + 0.5)) * log((k * delta_z + b) / b));
    // set new weighting current value
    el_current->inc_j_z(i_n+1, k_n, wj);

    ////////////////////////////////// /
    // calculate current jr in [i,k] cell //
    // equation y = k * x + b; //
    // finding k & b //
    k = -delta_z / delta_r;
    double r0 = (i_n + 0.5) * dr;
    double r1 =  radius_old;
    b = (k_n + 1.) * dz - longitude_old;

    // weighting jr in [i][k] cell
    some_shit_density = SOME_SHIT_DENSITY_Z(charge_array[p_number], r0, dr, dz, delta_t);

    wj = some_shit_density
      * (r0 * k * delta_r + k / 2. * delta_r * (radius_old + delta_r / 2.)
         + 0.5 * delta_r * (b - k * (2 * r0 + r1))
         + delta_r * (b-k * r1) * (4 * r0 * r0-dr * dr)
         / (8 * radius_old * (radius_old + delta_r)) +
         (k * (r0 * r0 / 2.-dr * dr / 8.))
         * log((radius_old + delta_r) / radius_old));
    el_current->inc_j_r(i_n,k_n, wj);

    b = longitude_old - k_n * dz;
    // weighting jr in [i][k+1] cell
    wj = some_shit_density
      * (-r0 * k * delta_r - k / 2. * delta_r * (radius_old + delta_r / 2.)
         + 0.5 * delta_r * (b + k * (2 * r0 + r1))
         + delta_r * (b + k * r1) * (4 * r0 * r0 - dr * dr)
         / (8 * radius_old * (radius_old + delta_r))
         - (k * (r0 * r0 / 2.-dr * dr / 8.))
         * log((radius_old + delta_r) / radius_old));
    el_current->inc_j_r(i_n, k_n+1, wj);
  }
  // if i cell is equal 0
  else
  {
    // equation y = k * x + b; //
    // finding k & b //
    double k = delta_r / delta_z;
    double b = radius_old;
    // calculate current jz in [i,k] cell //
    some_shit_density = charge_array[p_number]
      / (2. * PI * dr / 4. * dr * dz
         * delta_t
         * dr);

    wj = some_shit_density
      * (dr * delta_z - k * delta_z * delta_z / 2. - delta_z * b );
    // set new weighting current value
    el_current->inc_j_z(i_n,k_n, wj);

    // calculate current in [i + 1,k] cell //
    some_shit_density = SOME_SHIT_DENSITY_R(charge_array[p_number], dr, dr, dz, delta_t);

    wj = some_shit_density
      * (k * delta_z * delta_z / 2. + delta_z * dr + delta_z * b);
    // set new weighting current value
    el_current->inc_j_z(i_n + 1,k_n, wj);

    ////////////////////////////////// /
    // calculate current jr in [i,k] cell //
    // equation y = k * x + b; //
    // finding k & b //
    k = -delta_z / delta_r;
    double r0 = (i_n + 0.5) * dr;
    double r1 = radius_old;
    b= (k_n + 1.) * dz - longitude_old;

    // weighting jr in [i][k] cell
    some_shit_density = SOME_SHIT_DENSITY_Z(charge_array[p_number], r0, dr, dz, delta_t);

    wj = some_shit_density
      * (r0 * k * delta_r + k / 2. * delta_r * (radius_old + delta_r / 2.)
         + 0.5 * delta_r * (b-k * (2 * r0 + r1))
         + delta_r * (b-k * r1) * (4 * r0 * r0-dr * dr)
         / (8 * radius_old * (radius_old + delta_r))
         + (k * (r0 * r0 / 2.-dr * dr / 8.))
         * log((radius_old + delta_r) / radius_old));
    el_current->inc_j_r(i_n,k_n, wj);

    b = longitude_old- k_n * dz;;
    // weighting jr in [i][k + 1] cell
    // some_shit_density = SOME_SHIT_DENSITY_Z(charge_array[p_number], r0, dr, dz, delta_t);

    wj = some_shit_density
      * (-r0 * k * delta_r - k / 2. * delta_r * (radius_old + delta_r / 2.)
         + 0.5 * delta_r * (b + k * (2 * r0 + r1))
         + delta_r * (b + k * r1) * (4 * r0 * r0-dr * dr)
         / (8 * radius_old * (radius_old + delta_r))
         - (k * (r0 * r0 / 2.-dr * dr / 8.)) *
         log((radius_old + delta_r) / radius_old));
    el_current->inc_j_r(i_n, k_n + 1, wj);
  }
}

void Particles::current_distribution(Time * time1, Current * el_current)
{
  double dr = geom1->dr;
  double dz = geom1->dz;

  for (unsigned int i = 0;i < number;i++)
    if (is_alive[i])
    {
      // finding number new and old cells
      int i_n = CELL_NUMBER(pos[i][0], dr);
      int k_n = CELL_NUMBER(pos[i][2], dz);
      int i_o = CELL_NUMBER(pos_old[i][0], dr);
      int k_o = CELL_NUMBER(pos_old[i][2], dz);
      // TODO: workaround: sometimes it gives -1.
      // Just get 0 cell if it happence
      if (i_n < 0) i_n = 0;
      if (k_n < 0) k_n = 0;
      if (i_o < 0) i_o = 0;
      if (k_o < 0) k_o = 0;

      if (pos_old[i][0] == (i_o + 1) * dr) i_o = i_n;
      if (pos_old[i][2] == (k_o + 1) * dz) k_o = k_n;
      if (pos[i][0] == (i_n + 1) * dr) i_n = i_o;
      if (pos[i][2] == (k_n + 1) * dz) k_n = k_o;

      int res_cell = abs(i_n - i_o) + abs(k_n - k_o);

      if ((abs(pos[i][0] - pos_old[i][0]) < MNZL)
          || (abs(pos[i][2] - pos_old[i][2]) < MNZL))
        strict_motion_weighting(time1, el_current, pos[i][0], pos[i][2],
                                pos_old[i][0], pos_old[i][2], i);
      else
      {
        switch (res_cell)
        {
          // 1) charge in four nodes
        case 0: simple_current_distribution(time1, el_current, pos[i][0], pos[i][2],
                                   pos_old[i][0], pos_old[i][2], i_n, k_n, i);
          break;
          // 2) charge in 7 nodes
        case 1:
        {
          // charge in 7 nodes. Moving on r-axis (i_new != i_old)
          if ((i_n != i_o) && (k_n == k_o))
          {
            // moving to center from outer to innter cell
            if (pos_old[i][0] > (i_n + 1) * dr)
            {
              double a = (pos_old[i][0] - pos[i][0]) / (pos_old[i][2] - pos[i][2]);
              double r_boundary = (i_n + 1) * dr;
              double delta_r = r_boundary - pos[i][0];
              double z_boundary = pos[i][2] + delta_r / a;

              simple_current_distribution(time1, el_current, r_boundary, z_boundary,
                                 pos_old[i][0], pos_old[i][2], i_n + 1, k_n, i);
              simple_current_distribution(time1,el_current, pos[i][0], pos[i][2],
                                 r_boundary, z_boundary, i_n, k_n, i);
            }
            // moving to wall
            else
            {
              double a = (pos[i][0] - pos_old[i][0]) / (pos[i][2] - pos_old[i][2]);
              double r_boundary = (i_n) * dr;
              double delta_r = r_boundary - pos_old[i][0];
              double z_boundary = pos_old[i][2] + delta_r / a;

              simple_current_distribution(time1, el_current, r_boundary, z_boundary, pos_old[i][0], pos_old[i][2], i_n-1, k_n, i);
              simple_current_distribution(time1, el_current, pos[i][0], pos[i][2], r_boundary, z_boundary, i_n, k_n, i);
            }
          }
          // charge in seven cells. Moving on z-axis (k_new != k_old)
          else if ((i_n == i_o) && (k_n != k_o))
          {
            // moving forward from N to N + 1 cell
            if (pos_old[i][2] < k_n * dz)
            {
              double z_boundary = k_n * dz;
              double delta_z = z_boundary - pos_old[i][2];
              double a = (pos[i][0] - pos_old[i][0]) / (pos[i][2] - pos_old[i][2]);
              double r_boundary = pos_old[i][0] + a * delta_z;
              simple_current_distribution(time1, el_current, r_boundary, z_boundary ,pos_old[i][0], pos_old[i][2], i_n, k_n-1, i);
              simple_current_distribution(time1, el_current, pos[i][0], pos[i][2], r_boundary, z_boundary, i_n, k_n, i);
            }
            // moving backward
            else
            {
              double z_boundary = (k_n + 1) * dz;
              double delta_z = z_boundary - pos[i][2];
              double a = (pos_old[i][0] - pos[i][0]) / (pos_old[i][2] - pos[i][2]);
              double r_boundary = pos[i][0] + a * delta_z;
              simple_current_distribution(time1,el_current, r_boundary, z_boundary, pos_old[i][0], pos_old[i][2], i_n, k_n + 1, i);
              simple_current_distribution(time1,el_current, pos[i][0], pos[i][2], r_boundary, z_boundary, i_n, k_n, i);
            }
          }
        }
        break;
        // 3) charge in 10 nodes
        case 2:
        {
          // moving forward
          if (i_o < i_n)
          {
            // case, when particle move from [i-1][k-1] -> [i][k] cell
            if(k_o < k_n)
            {
              double a = (pos[i][0] - pos_old[i][0]) / (pos[i][2] - pos_old[i][2]);
              double r1 = i_n * dr;
              double delta_z1 = (r1 - pos_old[i][0]) / a;
              double z1 = pos_old[i][2] + delta_z1;
              double z2 = k_n * dz;
              double delta_r2 = (z2 - pos_old[i][2]) * a;
              double r2 = pos_old[i][0] + delta_r2;
              if (z1 < k_n * dz)
              {
                simple_current_distribution(time1, el_current, r1, z1 ,pos_old[i][0], pos_old[i][2], i_n-1, k_n-1,i);
                simple_current_distribution(time1, el_current, r2, z2, r1, z1, i_n, k_n-1,i);
                simple_current_distribution(time1, el_current, pos[i][0], pos[i][2], r2, z2, i_n, k_n,i);
              }
              else if (z1>k_n * dz)
              {
                simple_current_distribution(time1, el_current, r2, z2 ,pos_old[i][0], pos_old[i][2], i_n-1, k_n-1,i);
                simple_current_distribution(time1, el_current, r1, z1, r2, z2,i_n-1, k_n,i);
                simple_current_distribution(time1, el_current, pos[i][0], pos[i][2], r1, z1, i_n, k_n,i);
              }
            }
            // case, when particle move from [i-1][k + 1] -> [i][k] cell
            else
            {
              double a = (pos[i][0] - pos_old[i][0]) / (pos[i][2] - pos_old[i][2]);
              double r1 = i_n * dr;
              double delta_z1 = (r1 - pos_old[i][0]) / a;
              double z1 = pos_old[i][2] + delta_z1;

              double z2 = (k_n + 1) * dz;
              double delta_r2 = -(pos_old[i][2]-z2) * a;
              double r2 = pos_old[i][0] + delta_r2;
              if (z1>(k_n + 1) * dz)
              {
                simple_current_distribution(time1, el_current, r1, z1 ,pos_old[i][0], pos_old[i][2], i_n-1, k_n + 1,i);
                simple_current_distribution(time1, el_current, r2, z2, r1, z1, i_n, k_n + 1,i);
                simple_current_distribution(time1, el_current, pos[i][0], pos[i][2], r2, z2, i_n, k_n,i);
              }
              else if (z1<(k_n + 1) * dz)
              {
                simple_current_distribution(time1, el_current, r2, z2 ,pos_old[i][0], pos_old[i][2], i_n-1, k_n + 1,i);
                simple_current_distribution(time1, el_current, r1, z1, r2, z2,i_n-1, k_n,i);
                simple_current_distribution(time1, el_current, pos[i][0], pos[i][2], r1, z1, i_n, k_n,i);
              }
            }
          }
          // case, when particle move from [i + 1] cell to [i] cell
          else if (i_o > i_n)
          {
            // case, when particle move from [i + 1][k-1] -> [i][k] cell
            if(k_o<k_n)
            {
              double a = (pos[i][0] - pos_old[i][0]) / (pos[i][2] - pos_old[i][2]);
              double r1 = (i_n + 1) * dr;
              double delta_z1 = -(pos_old[i][0]-r1) / a;
              double z1 = pos_old[i][2] + delta_z1;

              double z2 = k_n * dz;
              double delta_r2 = -(z2-pos_old[i][2]) * a;
              double r2 = pos_old[i][0]- delta_r2;

              if (z1<(k_n) * dz)
              {
                simple_current_distribution(time1, el_current, r1, z1 ,pos_old[i][0], pos_old[i][2], i_n + 1, k_n-1,i);
                simple_current_distribution(time1,el_current, r2, z2, r1, z1, i_n, k_n-1,i);
                simple_current_distribution(time1,el_current, pos[i][0], pos[i][2], r2, z2, i_n, k_n,i);
              }
              else if (z1>(k_n) * dz)
              {
                simple_current_distribution(time1, el_current, r2, z2 ,pos_old[i][0], pos_old[i][2], i_n + 1, k_n-1,i);
                simple_current_distribution(time1, el_current, r1, z1, r2, z2,i_n + 1, k_n,i);
                simple_current_distribution(time1, el_current, pos[i][0], pos[i][2], r1, z1, i_n, k_n,i);
              }

            }
            // case, when particle move from [i + 1][k + 1] -> [i][k] cell
            else if (k_o>k_n)
            {
              double a = (pos_old[i][0]-pos[i][0]) / (pos_old[i][2]-pos[i][2]);
              double r1 = (i_n + 1) * dr;
              double delta_z1 = (r1-pos[i][0]) / a;
              double z1 = pos[i][2] + delta_z1;

              double z2 = (k_n + 1) * dz;
              double delta_r2 = (z2-pos[i][2]) * a;
              double r2 = pos[i][0] + delta_r2;

              if (z1>(k_n + 1) * dz)
              {
                simple_current_distribution(time1, el_current, r1, z1, pos_old[i][0],pos_old[i][2], i_n + 1, k_n + 1,i);
                simple_current_distribution(time1, el_current, r2, z2, r1, z1, i_n, k_n + 1,i);
                simple_current_distribution(time1, el_current, pos[i][0], pos[i][2], r2, z2, i_n, k_n,i);
              }
              else if (z1<(k_n + 1) * dz)
              {
                simple_current_distribution(time1, el_current, r2, z2, pos_old[i][0], pos_old[i][2], i_n + 1, k_n + 1,i);
                simple_current_distribution(time1, el_current, r1, z1, r2, z2,i_n + 1, k_n,i);
                simple_current_distribution(time1, el_current, pos[i][0], pos[i][2], r1, z1, i_n, k_n,i);
              }
            }
          }
        }
        break;
        }
      }
    }
}

void Particles::azimuthal_current_distribution(Current * this_j)
{
  double dr = geom1->dr;
  double dz = geom1->dz;

  for(unsigned int i = 0;i < number;i++)
    if (is_alive[i])
    {
      int r_i = 0;  // number of particle i cell
      int z_k = 0;  // number of particle k cell

      double r1, r2, r3; // temp variables for calculation
      double dz1, dz2;  // temp var.: width of k and k + 1 cell

      double ro_v = 0; // charge density Q / V, V - volume of particle
      double v_1 = 0; // volume of [i][k] cell
      double v_2 = 0; // volume of [i + 1][k] cell
      // double ro_v_2=0; // charge density in i + 1 cell

      double rho = 0; // charge density in cell
      double current; // j_phi in cell
      // double * * temp = this_j->get_j_phi();

      // finding number of i and k cell. example: dr = 0.5; r = 0.4; i =0
      r_i = CELL_NUMBER(pos[i][0], dr);
      z_k = CELL_NUMBER(pos[i][2], dz);
      // TODO: workaround: sometimes it gives -1.
      // Just get 0 cell if it happence
      if (r_i < 0) r_i = 0;
      if (z_k < 0) z_k = 0;

      // in first cell other alg. of ro_v calc
      if(pos[i][0] > dr)
      {
        r1 = pos[i][0] - 0.5 * dr;
        r2 = (r_i + 0.5) * dr;
        r3 = pos[i][0] + 0.5 * dr;
        ro_v = charge_array[i] / (2. * PI * dz * dr * pos[i][0]);
        v_1 = CELL_VOLUME(r_i, dr, dz);
        v_2 = CELL_VOLUME(r_i + 1, dr, dz);
	dz1 = (z_k + 0.5) * dz - (pos[i][2] - 0.5 * dz);
	dz2 = (pos[i][2] + 0.5 * dz) - (z_k + 0.5) * dz;
	
        // weighting in j[i][k] cell
        rho = ro_v * CYL_RNG_VOL(dz1, r1, r2) / v_1;
        current = rho * vel[i][1];
        this_j->inc_j_phi(r_i, z_k, current);

        // weighting in j[i + 1][k] cell
        rho = ro_v * CYL_RNG_VOL(dz1, r2, r3) / v_2;
        current = rho * vel[i][1];
        this_j->inc_j_phi(r_i + 1,z_k, current);

        // weighting in j[i][k + 1] cell
        rho = ro_v * CYL_RNG_VOL(dz2, r1, r2) / v_1;
        current = rho * vel[i][1];
        this_j->inc_j_phi(r_i, z_k + 1, current);

        // weighting in j[i + 1][k + 1] cell
        rho = ro_v * CYL_RNG_VOL(dz2, r2, r3) / v_2;
        current = rho * vel[i][1];
        this_j->inc_j_phi(r_i + 1, z_k + 1, current);

      }
      else
      {
        r1 = pos[i][0] - 0.5 * dr;
        r2 = (r_i + 0.5) * dr;
        r3 = pos[i][0] + 0.5 * dr;
	dz1 = (z_k + 0.5) * dz - (pos[i][2] - 0.5 * dz);
	dz2 = (pos[i][2] + 0.5 * dz) - (z_k + 0.5) * dz;
        ro_v = charge_array[i] / (2. * PI * dz * dr * pos[i][0]);
        v_1 = CYL_VOL(dz, dr);
        v_2 = CELL_VOLUME(r_i + 1, dr, dz);

        // weighting in j[i][k] cell
        rho = ro_v * CYL_RNG_VOL(dz1, r1, r2) / v_1;
        current = rho * vel[i][1];
        this_j->inc_j_phi(r_i, z_k, current);

        // weighting in j[i + 1][k] cell
        rho = ro_v * CYL_RNG_VOL(dz1, r2, r3) / v_2;
        current = rho * vel[i][1];
        this_j->inc_j_phi(r_i + 1,z_k, current);

        // weighting in j[i][k + 1] cell
        rho = ro_v * CYL_RNG_VOL(dz2, r1, r2) / v_1;
        current = rho * vel[i][1];
        this_j->inc_j_phi(r_i, z_k + 1, current);

        // weighting in j[i + 1][k + 1] cell
        rho = ro_v * CYL_RNG_VOL(dz2, r2, r3) / v_2;
        current = rho * vel[i][1];
        this_j->inc_j_phi(r_i + 1, z_k + 1, current);
      }
    }
}

void Particles::strict_motion_weighting(Time * time1,
                                        Current * this_j,
                                        double radius_new,
                                        double longitude_new,
                                        double radius_old,
                                        double longitude_old,
                                        int p_number)
{
  double dr = geom1->dr;
  double dz = geom1->dz;

  // defining number of cell
  int i_n = CELL_NUMBER(radius_new, dr);
  int k_n = CELL_NUMBER(longitude_new, dz);
  int i_o = CELL_NUMBER(radius_old, dr);;
  int k_o = CELL_NUMBER(longitude_old, dz);;
  // TODO: workaround: sometimes it gives -1.
  // Just get 0 cell if it happence
  if (i_n < 0) i_n = 0;
  if (k_n < 0) k_n = 0;
  if (i_o < 0) i_o = 0;
  if (k_o < 0) k_o = 0;

  if ((abs(radius_new - radius_old) < MNZL)
      && (abs(longitude_new - longitude_old) < MNZL))
    return;

  // stirct axis motion
  if (abs(radius_new - radius_old) < MNZL)
  {
    double delta_z = 0.;
    double value_part = 2. * PI * radius_new * dr * dz;
    double wj, wj_lower, wj_upper;
    double r1 = radius_new - 0.5 * dr;
    double r2 = (i_n + 0.5) * dr;
    double r3 = radius_new + 0.5 * dr;
    double delta_t = time1->delta_t;

    if (i_n == 0)
      wj_lower = charge_array[p_number]
        / (delta_t * PI * dr * dr / 4.)
        * PI * (r2 * r2 - r1 * r1) / value_part;
    else
      wj_lower = charge_array[p_number]
        / (delta_t * 2. * PI * i_n * dr * dr)
        * PI * (r2 * r2-r1 * r1) / value_part;
    wj_upper = charge_array[p_number]
      / (delta_t * 2 * PI * (i_n + 1) * dr * dr)
      * PI * (r3 * r3 - r2 * r2) / value_part;

    this_j->inc_j_r(i_n, k_n, 0.);
    this_j->inc_j_r(i_n, k_n + 1,0.);
    int res_k = k_n - k_o;

    switch(res_k)
    {
    case 0:
    {
      delta_z = longitude_new - longitude_old;
      wj = wj_lower * delta_z;
      this_j->inc_j_z(i_n, k_n, wj);
      wj = wj_upper * delta_z;
      this_j->inc_j_z(i_n+1, k_n, wj);
    }
    break;

    case 1:
    {
      delta_z = k_n * dz - longitude_old;
      wj = wj_lower * delta_z;
      this_j->inc_j_z(i_n, k_n-1, wj);
      wj = wj_upper * delta_z;
      this_j->inc_j_z(i_n+1, k_n-1, wj);

      delta_z = longitude_new - k_n * dz;
      wj = wj_lower * delta_z;
      this_j->inc_j_z(i_n, k_n, wj);
      wj = wj_upper * delta_z;
      this_j->inc_j_z(i_n+1, k_n, wj);
    }
    break;

    case -1:
    {
      delta_z = (k_n + 1) * dz - longitude_old;
      wj = wj_lower * delta_z;
      this_j->inc_j_z(i_n, k_n+1, wj);
      wj = wj_upper * delta_z;
      this_j->inc_j_z(i_n+1, k_n+1, wj);

      delta_z = longitude_new - (k_n+1) * dz;
      wj = wj_lower * delta_z;
      this_j->inc_j_z(i_n, k_n, wj);
      wj = wj_upper * delta_z;
      this_j->inc_j_z(i_n+1, k_n, wj);
    }
    break;
    }
  }

  // stirct radial motion
  else if (abs(longitude_new-longitude_old)<MNZL)
  {
    double wj, delta_r, left_delta_z, right_delta_z, res_j, r0, some_shit_density;
    int res_i = i_n - i_o;
    double delta_t = time1->delta_t;

    switch(res_i)
    {
    case 0:
    {
      delta_r = radius_new - radius_old;
      left_delta_z = (k_n + 1) * dz-longitude_new;
      right_delta_z = longitude_new - k_n * dz;
      r0 = (i_n + 0.5) * dr;
      some_shit_density = SOME_SHIT_DENSITY_STRICT(charge_array[p_number], r0, dr, dz, delta_t);

      wj = some_shit_density
        * (delta_r - r0 * r0 / (radius_old + delta_r)
           + r0 * r0 / radius_old
           + dr * dr / (4. * (radius_old + delta_r))
           - dr * dr
           / (4. * radius_old));

      res_j = wj * left_delta_z;
      this_j->inc_j_r(i_n,k_n,res_j);
      res_j = wj * right_delta_z;
      this_j->inc_j_r(i_n,k_n + 1,res_j);
    }
    break;
    case 1:
    {
      delta_r = (i_n) * dr - radius_old;
      left_delta_z = (k_n + 1) * dz-longitude_new;
      right_delta_z = longitude_new - k_n * dz;
      r0 = (i_n-0.5) * dr;
      some_shit_density = SOME_SHIT_DENSITY_STRICT(charge_array[p_number], r0, dr, dz, delta_t);

      wj = some_shit_density
        * (delta_r - r0 * r0 / (radius_old + delta_r)
           + r0 * r0 / radius_old
           + dr * dr / (4. * (radius_old + delta_r))
           - dr * dr / (4. * radius_old));

      res_j = wj * left_delta_z;
      this_j->inc_j_r(i_n-1,k_n,res_j);
      res_j = wj * right_delta_z;
      this_j->inc_j_r(i_n-1,k_n + 1,res_j);

      delta_r = radius_new - i_n * dr;
      r0 = (i_n + 0.5) * dr;
      some_shit_density = SOME_SHIT_DENSITY_STRICT(charge_array[p_number], r0, dr, dz, delta_t);

      wj = some_shit_density
        * (delta_r - r0 * r0 / (i_n * dr + delta_r)
           + r0 * r0 / i_n * dr
           + dr * dr / (4. * (i_n * dr + delta_r))
           - dr * dr / (4. * i_n * dr));

      res_j = wj * left_delta_z;
      this_j->inc_j_r(i_n,k_n,res_j);
      res_j = wj * right_delta_z;
      this_j->inc_j_r(i_n,k_n + 1,res_j);
    }
    break;
    case -1:
    {
      delta_r = (i_n + 1) * dr - radius_old ;
      left_delta_z = (k_n + 1) * dz-longitude_new;
      right_delta_z = longitude_new - k_n * dz;
      r0 = (i_n + 1.5) * dr;
      some_shit_density = SOME_SHIT_DENSITY_STRICT(charge_array[p_number], r0, dr, dz, delta_t);

      wj = some_shit_density
        * (delta_r - r0 * r0 / (radius_old + delta_r)
           + r0 * r0 / radius_old
           + dr * dr / (4. * (radius_old + delta_r))
           - dr * dr / (4. * radius_old));

      res_j = wj * left_delta_z;
      this_j->inc_j_r(i_n+1, k_n, res_j);
      res_j = wj * right_delta_z;
      this_j->inc_j_r(i_n+1, k_n+1, res_j);

      delta_r = radius_new - (i_n + 1) * dr;
      r0 = (i_n + 0.5) * dr;
      some_shit_density = SOME_SHIT_DENSITY_STRICT(charge_array[p_number], r0, dr, dz, delta_t);

      wj = some_shit_density
        * (delta_r - r0 * r0 / ((i_n + 1) * dr + delta_r)
           + r0 * r0 / (i_n + 1) * dr
           + dr * dr / (4. * ((i_n + 1) * dr + delta_r))
           - dr * dr / (4. * (i_n + 1) * dr));

      res_j = wj * left_delta_z;
      this_j->inc_j_r(i_n, k_n, res_j);
      res_j = wj * right_delta_z;
      this_j->inc_j_r(i_n, k_n+1, res_j);
    }
    break;
    }
  }
}

void Particles::back_position_to_rz()
{
  // ! implementation of backing coodrinates to rz pane
  // ! taken from https: // www.particleincell.com / 2015 / rz-pic /
#pragma omp parallel for
  for (unsigned int i = 0; i < number; i++)
    if (is_alive[i])
    {
      double r = lib::sq_rt(pos[i][0] * pos[i][0] + pos[i][1] * pos[i][1]);
      sin_theta_r[i] = pos[i][1] / r;
      pos[i][0] = r;
      pos[i][1] = 0;
    }
}

void Particles::back_velocity_to_rz_single(unsigned int i)
{
  // ! implementation of backing velocities to rz pane
  // ! taken from https: // www.particleincell.com / 2015 / rz-pic /
  // rotate velocity
  cos_theta_r[i] = lib::sq_rt( 1 - sin_theta_r[i] * sin_theta_r[i] );
  double u_2 = cos_theta_r[i] * vel[i][0] - sin_theta_r[i] * vel[i][1];
  double v_2 = sin_theta_r[i] * vel[i][0] + cos_theta_r[i] * vel[i][1];
  vel[i][0] = u_2;
  vel[i][1] = v_2;
}

void Particles::back_velocity_to_rz()
{
#pragma omp parallel for
  for (unsigned int i = 0; i < number; i++)
    if (is_alive[i])
      back_velocity_to_rz_single(i);
}

void Particles::boris_pusher_single(EField * e_fld, HField * h_fld,
                                    Time * t, unsigned int i)
{
  // !
  // ! boris pusher single
  // !
  // define vars directly in loop, because of multithreading
  double const1, const2, sq_velocity;
#ifdef PUSHER_BORIS_ADAPTIVE
  bool use_rel; // use relativistic calculations
#endif
#ifndef PUSHER_BORIS_CLASSIC // do not use gamma in non-relativistic boris pusher
  double gamma = 1;
#endif

  double * e;
  double * b;
  double vtmp[3];
  double velocity[3];

  // check if radius and longitude are correct
  if (isnan(pos[i][0]) ||
      isinf(pos[i][0]) != 0 ||
      isnan(pos[i][2]) ||
      isinf(pos[i][2]) != 0)
  {
    cerr << "ERROR(boris_pusher): radius[" << i << "] or longitude[" << i << "] is not valid number. Can not continue." << endl;
    exit(1);
  }

  e = e_fld->get_field(pos[i][0], pos[i][2]);
  b = h_fld->get_field(pos[i][0], pos[i][2]);

  // ! \f$ const1 = \frac{q t}{2 m} \f$
  // ! where \f$ q, m \f$ - particle charge and mass, \f$ t = \frac{\Delta t_{step}}{2} \f$
  const1 = charge_array[i] * t->delta_t / (2 * mass_array[i]);

  tinyvec3d::tv_product(e, const1);
  tinyvec3d::tv_product(b, MAGN_CONST * const1);

  // set velocity vector components and
  // round very small velicities to avoid exceptions
  velocity[0] = (abs(vel[i][0]) < MNZL) ? 0 : vel[i][0];
  velocity[1] = (abs(vel[i][1]) < MNZL) ? 0 : vel[i][1];
  velocity[2] = (abs(vel[i][2]) < MNZL) ? 0 : vel[i][2];

  // ! 0. check, if we should use classical calculations.
  // ! Required to increase modeling speed
#ifdef PUSHER_BORIS_ADAPTIVE
  if (pow(velocity[0], 2) + pow(velocity[1], 2) + pow(velocity[2], 2) > REL_LIMIT_POW_2)
    use_rel = true;
#endif

  // ! 1. Multiplication by relativistic factor (only for relativistic case)
  // ! \f$ u_{n-\frac{1}{2}} = \gamma_{n-\frac{1}{2}} * v_{n-\frac{1}{2}} \f$
#ifdef PUSHER_BORIS_ADAPTIVE
  if (use_rel)
#endif
#if defined (PUSHER_BORIS_RELATIVISTIC) || defined (PUSHER_BORIS_ADAPTIVE)
  {
    sq_velocity = tinyvec3d::tv_squared_sum(velocity);
    gamma = lib::get_gamma(sq_velocity);
    tinyvec3d::tv_product(velocity, gamma);
  }
#endif

  // ! 2. Half acceleration in the electric field
  // ! \f$ u'_n = u_{n-\frac{1}{2}} + \frac{q dt}{2 m E(n)} \f$
  // ! \f$ u'_n = u_{n-1 / 2} + \frac{q dt}{2 m E(n)} \f$
  tinyvec3d::tv_add(velocity, e);

  // ! 3. Rotation in the magnetic field
  // ! \f$ u" = u' + \frac{2}{1 + B'^2}  [(u' + [u' \times B'(n)] ) \times B'(n)] \f$,
  // ! \f$ B'(n) = \frac{B(n) q dt}{2 m * \gamma_n} \f$
#ifdef PUSHER_BORIS_ADAPTIVE
  if (use_rel)
#endif
#if defined (PUSHER_BORIS_RELATIVISTIC) || defined (PUSHER_BORIS_ADAPTIVE)
  {
    sq_velocity = tinyvec3d::tv_squared_sum(velocity);
    gamma = lib::get_gamma_inv(sq_velocity);
    tinyvec3d::tv_product(b, gamma);
  }
#endif
  // ! \f$ const2 = \frac{2}{1 + b_1^2 + b_2^2 + b_3^2} \f$
  const2 = 2. / (1. + tinyvec3d::tv_squared_sum(b));

  // set temporary velocity as old values
  // to calculate magnetic rotation
  tinyvec3d::tv_copy_components(velocity, vtmp);

  velocity[0] = vtmp[0] + const2 * (
    (vtmp[1] - vtmp[0] * b[2] + vtmp[2] * b[0]) * b[2]
    - (vtmp[2] + vtmp[0] * b[1] - vtmp[1] * b[0]) * b[1]
    );
  velocity[1] = vtmp[1] + const2 * (
    -(vtmp[0] + vtmp[1] * b[2] - vtmp[2] * b[1]) * b[2]
    + (vtmp[2] + vtmp[0] * b[1] - vtmp[1] * b[0]) * b[0]
    );
  velocity[2] = vtmp[2] + const2 * (
    (vtmp[0] + vtmp[1] * b[2] - vtmp[2] * b[1]) * b[1]
    - (vtmp[1] - vtmp[0] * b[2] + vtmp[2] * b[0]) * b[0]
    );

  // ! 4. Half acceleration in the electric field
  // ! \f$ u_{n + \frac{1}{2}} = u_n + \frac{q dt}{2 m E(n)} \f$
  tinyvec3d::tv_add(velocity, e);

  // ! 5. Division by relativistic factor
#ifdef PUSHER_BORIS_ADAPTIVE
  if (use_rel)
#endif
#if defined (PUSHER_BORIS_RELATIVISTIC) || defined (PUSHER_BORIS_ADAPTIVE)
  {
    sq_velocity = tinyvec3d::tv_squared_sum(velocity);
    gamma = lib::get_gamma_inv(sq_velocity);
    tinyvec3d::tv_product(velocity, gamma);
  }
#endif
  vel[i][0] = velocity[0];
  vel[i][1] = velocity[1];
  vel[i][2] = velocity[2];

  delete [] e;
  delete [] b;
}

void Particles::vay_pusher_single(EField * e_fld, HField * h_fld,
                                  Time * t, unsigned int i)
{
  // !
  // ! Vay pusher single
  // !

  double gamma;

  double const1, const2, sq_velocity;
  double * e;
  double * b;

  double momentum[3];
  double velocity[3];
  double vtmp[3];

  // check if radius and longitude are correct
  if (isnan(pos[i][0]) ||
      isinf(pos[i][0]) != 0 ||
      isnan(pos[i][2]) ||
      isinf(pos[i][2]) != 0)
  {
    cerr << "ERROR(vay_pusher): radius[" << i << "] or longitude[" << i << "] is not valid number. Can not continue." << endl;
    exit(1);
  }

  e = e_fld->get_field(pos[i][0], pos[i][2]);
  b = h_fld->get_field(pos[i][0], pos[i][2]);

  // ! \f$ const1 = \frac{q t}{2 m} \f$
  // ! where \f$ q, m \f$ - particle charge and mass, \f$ t = \frac{\Delta t_{step}}{2} \f$
  const1 = charge_array[i] * t->delta_t / (2 * mass_array[i]);

  tinyvec3d::tv_product(e, 2 * const1);
  tinyvec3d::tv_product(b, MAGN_CONST * const1);

  // set velocity vector components and
  // round very small velicities to avoid exceptions
  velocity[0] = (abs(vel[i][0]) < MNZL) ? 0 : vel[i][0];
  velocity[1] = (abs(vel[i][1]) < MNZL) ? 0 : vel[i][1];
  velocity[2] = (abs(vel[i][2]) < MNZL) ? 0 : vel[i][2];

  // get normalized relativistic momentum \f$ \gamma v \f$
  sq_velocity = tinyvec3d::tv_squared_sum(velocity);
  gamma = lib::get_gamma(sq_velocity);
  tinyvec3d::tv_product(velocity, gamma);

  // ____________________________________________
  // Part I: Computation of uprime

  tinyvec3d::tv_copy_components(velocity, vtmp);

  // Add Electric field
  tinyvec3d::tv_product(b, 2 * const1);
  tinyvec3d::tv_add(velocity, e);

  // Add magnetic field
  tinyvec3d::tv_product(b, const1);

  sq_velocity = tinyvec3d::tv_squared_sum(vtmp);
  gamma = lib::get_gamma_inv(sq_velocity);

  velocity[0] += gamma * (vtmp[1] * b[2] - vtmp[2] * b[1]);
  velocity[1] += gamma * (vtmp[2] * b[0] - vtmp[0] * b[2]);
  velocity[2] += gamma * (vtmp[0] * b[1] - vtmp[1] * b[0]);

  // alpha is 1 / gamma^2
  sq_velocity = tinyvec3d::tv_squared_sum(velocity);
  gamma = lib::get_gamma_inv(sq_velocity);

  double alpha = 1. / gamma / gamma;
  double B2 = tinyvec3d::tv_squared_sum(b);

  // ___________________________________________
  // Part II: Computation of Gamma^{i + 1}

  // s is sigma
  double s = alpha - B2;
  double us2 = pow(tinyvec3d::tv_dot(velocity, b), 2);

  // alpha becomes 1 / gamma^{i + 1}
  alpha = 1. / lib::sq_rt(0.5 * (s + lib::sq_rt(s * s + 4. * ( B2 + us2 ))));

  tinyvec3d::tv_product(b, alpha);

  s = 1. / (1. + tinyvec3d::tv_squared_sum(b));
  alpha = tinyvec3d::tv_dot(velocity, b);

  vtmp[0] = s * (velocity[0] + alpha * b[0] + b[2] * velocity[1] - b[1] * velocity[2]);
  vtmp[1] = s * (velocity[1] + alpha * b[1] + b[0] * velocity[2] - b[2] * velocity[0]);
  vtmp[2] = s * (velocity[2] + alpha * b[2] + b[1] * velocity[0] - b[0] * velocity[1]);

  // Inverse Gamma factor
  sq_velocity = tinyvec3d::tv_squared_sum(vtmp);
  gamma = lib::get_gamma_inv(sq_velocity);
  tinyvec3d::tv_product(vtmp, gamma);

  tinyvec3d::tv_copy_components(vtmp, vel[i]);

  delete [] e;
  delete [] b;
}

void Particles::dump_position_to_old_single(unsigned int i)
{
  pos_old[i][0] = pos[i][0];
  pos_old[i][2] = pos[i][2];
}

void Particles::dump_position_to_old()
{
#pragma omp parallel for
  for (unsigned int i = 0; i < number; i++)
    if (is_alive[i])
      dump_position_to_old_single(i);
}

void Particles::half_step_pos_single(Time * t, unsigned int i)
{
  double half_dt = t->delta_t / 2.;
  // check if radius and longitude are correct
  if (isnan(pos[i][0]) || isinf(pos[i][0]) != 0 || isnan(pos[i][2]) || isinf(pos[i][2]) != 0)
  {
    cerr << "ERROR(half_step_pos): radius[" << i << "] or longitude[" << i << "] is not valid number. Can not continue." << endl;
    exit(1);
  }
  //
  pos[i][0] = pos[i][0] + vel[i][0] * half_dt;
  // ! we use "fake" rotation component to correct position from xy to rz pane
  pos[i][1] = pos[i][1] + vel[i][1] * half_dt;
  pos[i][2] = pos[i][2] + vel[i][2] * half_dt;
}

void Particles::back_position_to_rz_single(unsigned int i)
{
  // ! implementation of backing coodrinates to rz pane
  // ! taken from https: // www.particleincell.com / 2015 / rz-pic /
  double r = lib::sq_rt(pos[i][0] * pos[i][0] + pos[i][1] * pos[i][1]);
  sin_theta_r[i] = pos[i][1] / r;
  pos[i][0] = r;
  pos[i][1] = 0;
}

void Particles::step_v(EField * e_fld, HField * h_fld, Time * t)
{
#pragma omp parallel for shared(e_fld, h_fld, t)
  for(unsigned int i = 0; i < number; i++)
  {
    if (is_alive[i])
    {
#if defined (PUSHER_BORIS_RELATIVISTIC) || defined (PUSHER_BORIS_ADAPTIVE) || (PUSHER_BORIS_CLASSIC)
      boris_pusher_single(e_fld, h_fld, t, i);
#elif defined (PUSHER_VAY)
      vay_pusher_single(e_fld, h_fld, t, i);
#endif
      // ! dump position to old
      dump_position_to_old_single(i);
      // ! half step position
      half_step_pos_single(t, i);
      // ! back position to rz
      back_position_to_rz_single(i);
    }
  }
}

void Particles::reflection_single(unsigned int i)
{
  double dr = geom1->dr;
  double dz = geom1->dz;
  double radius_wall = geom1->r_size - dr / 2.;
  double longitude_wall = geom1->z_size - dz / 2.;
  double half_dr = dr / 2.;
  double half_dz = dz / 2.;
  double radius_wallX2 = radius_wall * 2.;
  double longitude_wallX2 = longitude_wall * 2.;

  // ! FIXME: fix wall reflections for r-position
  if (pos[i][0] > radius_wall)
  {
    pos[i][0] = radius_wallX2 - pos[i][0];
    vel[i][0] = -vel[i][0];
  }

  if (pos[i][2] > longitude_wall)
  {
    pos[i][2] = longitude_wallX2 - pos[i][2];
    vel[i][2] = -vel[i][2];
  }

  if (pos[i][0] < half_dr)
  {
    pos[i][0] = dr - pos[i][0];
    vel[i][0] = -vel[i][0];
  }

  if (pos[i][2] < half_dz)
  {
    pos[i][2] = dz - pos[i][2];
    vel[i][2] = -vel[i][2];
  }
}

void Particles::move_half_reflect(Time *t)
{
#pragma omp parallel for shared(t)
  for (unsigned int i = 0; i < number; i++)
    if (is_alive[i])
    {
      half_step_pos_single(t, i);
      reflection_single(i);
      back_position_to_rz_single(i);
    }
}

void Particles::full_current_distribution(Current * current, Time *t)
{
  azimuthal_current_distribution(current);
  move_half_reflect(t);
  current_distribution(t, current);
  back_velocity_to_rz();
}
