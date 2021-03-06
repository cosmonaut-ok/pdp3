#include "particlesDensity.h"

//// constructor
ParticlesDensity::ParticlesDensity(Geometry *geom1):geom(geom1)
{
  density = new double *[geom->n_grid_r - 1];
  // t_sum = new double *[geom->n_grid_r];
  // t_vec_r = new double *[geom->n_grid_r];
  // t_vec_phi = new double *[geom->n_grid_r];
  // t_vec_z = new double *[geom->n_grid_r];
  // count = new double *[geom->n_grid_r];
  // count_sum = new double *[geom->n_grid_r];

  // filling second demension
  // for potential
  for (int i = 0; i < (geom->n_grid_r); i++)
  {
    density[i] = new double[geom->n_grid_z - 1];
    // t_sum[i] = new double[geom->n_grid_z];
    // t_vec_r[i] = new double[geom->n_grid_z];
    // t_vec_phi[i] = new double[geom->n_grid_z];
    // t_vec_z[i] = new double[geom->n_grid_z];
    // count[i] = new double[geom->n_grid_z];
    // count_sum[i] = new double[geom->n_grid_z];
  }
}

void ParticlesDensity::reset(void)
{
// #pragma omp parallel for
  for (unsigned int i = 0; i < geom->n_grid_r; i++)
    for (unsigned int j = 0; j < geom->n_grid_z; j++)
    {
      density[i][j] = 0;
      // t_sum[i][j] = 0;
      // t_vec_r[i][j] = 0;
      // t_vec_phi[i][j] = 0;
      // t_vec_z[i][j] = 0;
      // count[i][j] = 0;
      // count_sum[i][j] = 0;
    }
}

void ParticlesDensity::calc_density(Particles *prtls)
{
  reset();

  double dr = geom->dr;
  double dz = geom->dz;
  double r1, r2, r3; // temp variables for calculation
  double dz1, dz2; // temp var.: width of k and k + 1 cell

  double ro_v = 0; // charge density Q/V, V - volume of particle
  double ro_v_r = 0;
  double ro_v_phi = 0;
  double ro_v_z = 0;
  double ro_p = 0; // charge density Q/V, V - volume of particle
  double v_0 = 0; // charge density Q/V, V - volume of particle
  double v_1 = 0; // volume of [i][k] cell
  double v_2 = 0; // volume of [i + 1][k] cell
  double vel = 0;
  double vel_r = 0;
  double vel_phi = 0;
  double vel_z = 0;

  double value = 0;
  // double **temp = ro1->get_rho();

  for(unsigned int i = 0; i < prtls->number; i++)
  {
    // finding number of i and k cell. example: dr = 0.5; r = 0.4; i =0
    unsigned int r_i = CELL_NUMBER(prtls->pos[i][0], geom->dr);
    unsigned int z_k = CELL_NUMBER(prtls->pos[i][2], geom->dz);

    vel_r = prtls->vel[i][0];
    vel_phi = prtls->vel[i][1];
    vel_z = prtls->vel[i][2];

    vel = lib::sq_rt(vel_r * vel_r + vel_phi * vel_phi + vel_z * vel_z);

    // in first cell other alg. of ro_v calc
    if (prtls->pos[i][0] > dr)
    {
      r1 =  prtls->pos[i][0] - 0.5 * dr;
      r2 = (r_i + 0.5) * dr;
      r3 = prtls->pos[i][0] + 0.5 * dr;
      v_0 = 2. * PI * dz * dr * prtls->pos[i][0];

      // ro_v = prtls->mass_array[i] * vel / v_0;
      // ro_v_r = prtls->mass_array[i] * vel_r / v_0;
      // ro_v_phi = prtls->mass_array[i] * vel_phi / v_0;
      // ro_v_z = prtls->mass_array[i] * vel_z / v_0;
      ro_p = prtls->mass_array[i] / EL_MASS / v_0;

      v_1 = CELL_VOLUME(r_i, dr, dz);
      v_2 = CELL_VOLUME(r_i + 1, dr, dz);
      dz1 = (z_k + 0.5) * dz - (prtls->pos[i][2] - 0.5 * dz);
      dz2 = (prtls->pos[i][2] + 0.5 * dz) - (z_k + 0.5) * dz;

      // weighting in ro[i][k] cell
      value = CYL_RNG_VOL(dz1, r1, r2) / v_1;
      // inc_tmpr(r_i, z_k, ro_v * value);
      // inc_vec_r(r_i, z_k, ro_v_r * value);
      // inc_vec_phi(r_i, z_k, ro_v_phi * value);
      // inc_vec_z(r_i, z_k, ro_v_z * value);
      inc_count(r_i, z_k, ro_p * value);

      // weighting in ro[i + 1][k] cell
      value = CYL_RNG_VOL(dz1, r2, r3) / v_2;
      // inc_tmpr(r_i+1, z_k, ro_v * value);
      // inc_vec_r(r_i+1, z_k, ro_v_r * value);
      // inc_vec_phi(r_i+1, z_k, ro_v_phi * value);
      // inc_vec_z(r_i+1, z_k, ro_v_z * value);
      inc_count(r_i+1, z_k, ro_p * value);

      // weighting in ro[i][k + 1] cell
      value = CYL_RNG_VOL(dz2, r1, r2) / v_1;
      // inc_tmpr(r_i, z_k+1, ro_v * value);
      // inc_vec_r(r_i, z_k+1, ro_v_r * value);
      // inc_vec_phi(r_i, z_k+1, ro_v_phi * value);
      // inc_vec_z(r_i, z_k+1, ro_v_z * value);
      inc_count(r_i, z_k+1, ro_p * value);

      // weighting in ro[i + 1][k + 1] cell
      value = CYL_RNG_VOL(dz2, r2, r3) / v_2;
      // inc_tmpr(r_i + 1, z_k + 1, ro_v * value);
      // inc_vec_r(r_i + 1, z_k + 1, ro_v_r * value);
      // inc_vec_phi(r_i + 1, z_k + 1, ro_v_phi * value);
      // inc_vec_z(r_i + 1, z_k + 1, ro_v_z * value);
      inc_count(r_i + 1, z_k + 1, ro_p * value);
    }
    else if (prtls->pos[i][0] <= dr / 2.)
    {
      r_i = 0;
      r1 =  0.;
      r2 = (r_i + 0.5) * dr;
      r3 = prtls->pos[i][0] + 0.5 * dr;
      dz1 = (z_k + 0.5) * dz - (prtls->pos[i][2] - 0.5 * dz);
      dz2 = (prtls->pos[i][2] + 0.5 * dz) - (z_k + 0.5) * dz;
      v_0 = PI * dz * (2. * prtls->pos[i][0] * prtls->pos[i][0] + dr * dr / 2.);

      // ro_v = prtls->mass_array[i] * vel / v_0;
      // ro_v_r = prtls->mass_array[i] * vel_r / v_0;
      // ro_v_phi = prtls->mass_array[i] * vel_phi / v_0;
      // ro_v_z = prtls->mass_array[i] * vel_z / v_0;

      ro_p = prtls->mass_array[i] / EL_MASS / v_0;
      v_1 = CYL_VOL(dz, dr);
      v_2 = CELL_VOLUME(r_i + 1, dr, dz);
      ////////////////////////// /

      // weighting in ro[i][k] cell
      value = PI * dz1 * (dr * dr / 2. - prtls->pos[i][0] * dr + prtls->pos[i][0] * prtls->pos[i][0]) / v_1;
      // inc_tmpr(r_i, z_k, ro_v * value);
      // inc_vec_r(r_i, z_k, ro_v_r * value);
      // inc_vec_phi(r_i, z_k, ro_v_phi * value);
      // inc_vec_z(r_i, z_k, ro_v_z * value);
      inc_count(r_i, z_k, ro_p * value);

      // weighting in ro[i + 1][k] cell
      value = CYL_RNG_VOL(dz1, r2, r3) / v_2;
      // inc_tmpr(r_i + 1,z_k, ro_v * value);
      // inc_vec_r(r_i + 1,z_k, ro_v_r * value);
      // inc_vec_phi(r_i + 1,z_k, ro_v_phi * value);
      // inc_vec_z(r_i + 1,z_k, ro_v_z * value);
      inc_count(r_i + 1,z_k, ro_p * value);

      // weighting in ro[i][k + 1] cell
      value = PI * dz2 * (dr * dr / 2. - prtls->pos[i][0] * dr + prtls->pos[i][0] * prtls->pos[i][0]) / v_1;
      // inc_tmpr(r_i, z_k + 1, ro_v* value);
      // inc_vec_r(r_i, z_k + 1, ro_v_r* value);
      // inc_vec_phi(r_i, z_k + 1, ro_v_phi* value);
      // inc_vec_z(r_i, z_k + 1, ro_v_z* value);
      inc_count(r_i, z_k + 1, ro_p* value);

      // weighting in ro[i + 1][k + 1] cell
      value = CYL_RNG_VOL(dz2, r2, r3) / v_2;
      // inc_tmpr(r_i + 1, z_k + 1, ro_v * value);
      // inc_vec_r(r_i + 1, z_k + 1, ro_v_r * value);
      // inc_vec_phi(r_i + 1, z_k + 1, ro_v_phi * value);
      // inc_vec_z(r_i + 1, z_k + 1, ro_v_z * value);
      inc_count(r_i + 1, z_k + 1, ro_p * value);
    }
    else
    {
      ////////////////////////// /
      r1 = prtls->pos[i][0] - 0.5 * dr;
      r2 = (r_i + 0.5) * dr;
      r3 = prtls->pos[i][0] + 0.5 * dr;
      dz1 = (z_k + 0.5) * dz - (prtls->pos[i][2] - 0.5 * dz);
      dz2 = (prtls->pos[i][2] + 0.5 * dz) - (z_k + 0.5) * dz;
      v_0 = 2. * PI * dz * dr * prtls->pos[i][0];

      // ro_v = prtls->mass_array[i] * vel / v_0;
      // ro_v_r = prtls->mass_array[i] * vel_r / v_0;
      // ro_v_phi = prtls->mass_array[i] * vel_phi / v_0;
      // ro_v_z = prtls->mass_array[i] * vel_z / v_0;

      ro_p = prtls->mass_array[i] / EL_MASS / v_0;
      v_1 = CYL_VOL(dz, dr);
      v_2 = CELL_VOLUME(r_i + 1, dr, dz);
      ////////////////////////// /

      // weighting in ro[i][k] cell
      value = CYL_RNG_VOL(dz1, r1, r2) / v_1;
      // inc_tmpr(r_i, z_k, ro_v * value);
      // inc_vec_r(r_i, z_k, ro_v_r * value);
      // inc_vec_phi(r_i, z_k, ro_v_phi * value);
      // inc_vec_z(r_i, z_k, ro_v_z * value);
      inc_count(r_i, z_k, ro_p * value);

      // weighting in ro[i + 1][k] cell
      value = CYL_RNG_VOL(dz1, r2, r3) / v_2;
      // inc_tmpr(r_i + 1,z_k, ro_v * value);
      // inc_vec_r(r_i + 1,z_k, ro_v_r * value);
      // inc_vec_phi(r_i + 1,z_k, ro_v_phi * value);
      // inc_vec_z(r_i + 1,z_k, ro_v_z * value);
      inc_count(r_i + 1,z_k, ro_p * value);

      // weighting in ro[i][k + 1] cell
      value = CYL_RNG_VOL(dz2, r1, r2) / v_1;
      // inc_tmpr(r_i, z_k + 1, ro_v * value);
      // inc_vec_r(r_i, z_k + 1, ro_v_r * value);
      // inc_vec_phi(r_i, z_k + 1, ro_v_phi * value);
      // inc_vec_z(r_i, z_k + 1, ro_v_z * value);
      inc_count(r_i, z_k + 1, ro_p * value);

      // weighting in ro[i + 1][k + 1] cell
      value = CYL_RNG_VOL(dz2, r2, r3) / v_2;
      // inc_tmpr(r_i + 1, z_k + 1, ro_v * value);
      // inc_vec_r(r_i + 1, z_k + 1, ro_v_r * value);
      // inc_vec_phi(r_i + 1, z_k + 1, ro_v_phi * value);
      // inc_vec_z(r_i + 1, z_k + 1, ro_v_z * value);
      inc_count(r_i + 1, z_k + 1, ro_p * value);
    }
  }
}

// void ParticlesDensity::inc_tmpr(unsigned int r, unsigned int z, double value)
// {
//   t_sum[r][z] += value;
// }

// void ParticlesDensity::inc_vec_r(unsigned int r, unsigned int z, double value)
// {
//   t_vec_r[r][z] += value;
// }

// void ParticlesDensity::inc_vec_phi(unsigned int r, unsigned int z, double value)
// {
//   t_vec_phi[r][z] += value;
// }

// void ParticlesDensity::inc_vec_z(unsigned int r, unsigned int z, double value)
// {
//   t_vec_z[r][z] += value;
// }

void ParticlesDensity::inc_count(unsigned int r, unsigned int z, double value)
{
  density[r][z] += value;
}

// void ParticlesDensity::normalize(Particles *prtls)
// {
//   // lib::bilinear_interpolation(t_sum, t, geom->n_grid_r, geom->n_grid_z);
//   // lib::bilinear_interpolation(count_sum, count, geom->n_grid_r, geom->n_grid_z);

//   // lib::bicubic_interpolation(t_sum, t, geom->n_grid_r, geom->n_grid_z);
//   // lib::bicubic_interpolation(count_sum, count, geom->n_grid_r, geom->n_grid_z);

//   for (unsigned int r = 0; r < geom->n_grid_r - 1; r++)
//     for (unsigned int z = 0; z < geom->n_grid_z - 1; z++)
//     {
//       double p_sum2 = t_vec_r[r][z] * t_vec_r[r][z] + t_vec_phi[r][z] * t_vec_phi[r][z] + t_vec_z[r][z] * t_vec_z[r][z];
//       t[r][z] = (t_sum[r][z] * t_sum[r][z] - p_sum2) / (2 * EL_MASS * count_sum[r][z] * count_sum[r][z] * EL_CHARGE);
//     }
// }
