#include "temperature.h"

//// constructor
Temperature::Temperature(Geometry *geom1):geom(geom1)
{
  t = new double *[geom->n_grid_r];
  t_sum = new double *[geom->n_grid_r];
  count = new double *[geom->n_grid_r];

  // filling second demension
  // for potential
  for (int i = 0; i < (geom->n_grid_r); i++)
  {
    t[i] = new double[geom->n_grid_z];
    t_sum[i] = new double[geom->n_grid_z];
    count[i] = new double[geom->n_grid_z];
  }
}

void Temperature::reset(void)
{
// #pragma omp parallel for
  for (unsigned int i = 0; i < geom->n_grid_r; i++)
    for (unsigned int j = 0; j < geom->n_grid_z; j++)
    {
      t[i][j] = 0;
      t_sum[i][j] = 0;
      count[i][j] = 0;
    }
}

void Temperature::calc_t_r(Particles *prtls)
{
  reset();

  for (unsigned int i = 0; i < prtls->number; i++)
  {
    unsigned int i_r = CELL_NUMBER(prtls->pos[i][0], geom->dr);
    unsigned int k_z = CELL_NUMBER(prtls->pos[i][2], geom->dz);

    t_sum[i_r][k_z] = t_sum[i_r][k_z] + VEL_TO_TEMPR(abs(prtls->vel[i][0]), prtls->mass);
    count[i_r][k_z]++;
  }

// #pragma omp parallel for
  for (unsigned int i = 0; i < geom->n_grid_r; i++)
    for (unsigned int j = 0; j < geom->n_grid_z; j++)
      {
	if (count[i][j] == 0) count[i][j] = 1;
	t[i][j] = t_sum[i][j] / count[i][j];
      }
}

void Temperature::calc_t_phi(Particles *prtls)
{
  reset();

  for (unsigned int i = 0; i < prtls->number; i++)
  {
    unsigned int i_r = CELL_NUMBER(prtls->pos[i][0], geom->dr);
    unsigned int k_z = CELL_NUMBER(prtls->pos[i][2], geom->dz);

    t_sum[i_r][k_z] = t_sum[i_r][k_z] + VEL_TO_TEMPR(abs(prtls->vel[i][1]), prtls->mass);
    count[i_r][k_z]++;
  }

// #pragma omp parallel for
  for (unsigned int i = 0; i < geom->n_grid_r; i++)
    for (unsigned int j = 0; j < geom->n_grid_z; j++)
      {
	if (count[i][j] == 0) count[i][j] = 1;
	t[i][j] = t_sum[i][j] / count[i][j];
      }
}

void Temperature::calc_t_z(Particles *prtls)
{
  reset();

  for (unsigned int i = 0; i < prtls->number; i++)
  {
    unsigned int i_r = CELL_NUMBER(prtls->pos[i][0], geom->dr);
    unsigned int k_z = CELL_NUMBER(prtls->pos[i][2], geom->dz);

    t_sum[i_r][k_z] = t_sum[i_r][k_z] + VEL_TO_TEMPR(abs(prtls->vel[i][2]), prtls->mass);
    count[i_r][k_z]++;
  }

// #pragma omp parallel for
  for (unsigned int i = 0; i < geom->n_grid_r; i++)
    for (unsigned int j = 0; j < geom->n_grid_z; j++)
      {
	if (count[i][j] == 0) count[i][j] = 1;
	t[i][j] = t_sum[i][j] / count[i][j];
      }
}

void Temperature::calc_t(Particles *prtls)
{
  reset();

  for (unsigned int i = 0; i < prtls->number; i++)
  {
    unsigned int i_r = CELL_NUMBER(prtls->pos[i][0], geom->dr);
    unsigned int k_z = CELL_NUMBER(prtls->pos[i][2], geom->dz);

    double t_r_cmp = VEL_TO_TEMPR(prtls->vel[i][0], prtls->mass);
    double t_phi_cmp = VEL_TO_TEMPR(prtls->vel[i][1], prtls->mass);
    double t_z_cmp = VEL_TO_TEMPR(prtls->vel[i][2], prtls->mass);

    t_sum[i_r][k_z] = t_sum[i_r][k_z] + lib::sq_rt(t_r_cmp * t_r_cmp + t_phi_cmp * t_phi_cmp + t_z_cmp * t_z_cmp);
    count[i_r][k_z]++;
  }

// #pragma omp parallel for
  for (unsigned int i = 0; i < geom->n_grid_r; i++)
    for (unsigned int j = 0; j < geom->n_grid_z; j++)
    {
      if (count[i][j] == 0) count[i][j] = 1;
      t[i][j] = t_sum[i][j] / count[i][j];
    }
}

void Temperature::calc_t_r(Particles *prtls)
{
  reset();

  double dr = geom->dr;
  double dz = geom->dz;
  double r1, r2, r3; // temp variables for calculation
  double dz1, dz2; // temp var.: width of k and k + 1 cell

  double ro_v = 0; // charge density Q/V, V - volume of particle
  double v_1 = 0; // volume of [i][k] cell
  double v_2 = 0; // volume of [i + 1][k] cell
  double vel2 = 0;

  double value = 0;
  // double **temp = ro1->get_rho();

  for(unsigned int i = 0; i < prtls->number; i++)
  {
    // finding number of i and k cell. example: dr = 0.5; r = 0.4; i =0
    unsigned int r_i = CELL_NUMBER(prtls->pos[i][0], geom->dr);
    unsigned int z_k = CELL_NUMBER(prtls->pos[i][2], geom->dz);

    vel2 = prtls->vel[i][0] * prtls->vel[i][0];

    // in first cell other alg. of ro_v calc
    if (prtls->pos[i][0] > dr)
    {
      r1 =  prtls->pos[i][0] - 0.5 * dr;
      r2 = (r_i + 0.5) * dr;
      r3 = prtls->pos[i][0] + 0.5 * dr;
      ro_v = prtls->mass_array[i] * vel2 / (4. * PI * dz * dr * prtls->pos[i][0]);


      v_1 = CELL_VOLUME(r_i, dr, dz);
      v_2 = CELL_VOLUME(r_i + 1, dr, dz);
      dz1 = (z_k + 1) * dz - prtls->pos[i][2];
      dz2 = prtls->pos[i][2] - z_k * dz;

      // weighting in ro[i][k] cell
      value = ro_v * CYL_RNG_VOL(dz1, r1, r2) / v_1;
      inc_tmpr(r_i, z_k, value);

      // weighting in ro[i + 1][k] cell
      value = ro_v * CYL_RNG_VOL(dz1, r2, r3) / v_2;
      inc_tmpr(r_i+1, z_k, value);

      // weighting in ro[i][k + 1] cell
      value = ro_v * CYL_RNG_VOL(dz2, r1, r2) / v_1;
      inc_tmpr(r_i, z_k+1, value);

      // weighting in ro[i + 1][k + 1] cell
      value = ro_v * CYL_RNG_VOL(dz2, r2, r3) / v_2;
      inc_tmpr(r_i + 1, z_k + 1, value);
    }
    else if (prtls->pos[i][0] <= dr / 2.)
    {
      r_i = 0;
      r1 =  0.;
      r2 = (r_i + 0.5) * dr;
      r3 = prtls->pos[i][0] + 0.5 * dr;
      dz1 = (z_k + 1) * dz - prtls->pos[i][2];
      dz2 = prtls->pos[i][2] - z_k * dz;
      ro_v = prtls->mass_array[i] * vel2 / (2 * PI * dz * (2. * prtls->pos[i][0] * prtls->pos[i][0] + dr * dr / 2.));
      v_1 = CYL_VOL(dz, dr);
      v_2 = CELL_VOLUME(r_i + 1, dr, dz);
      ////////////////////////// /

      // weighting in ro[i][k] cell
      value = ro_v * PI * dz1 * (dr * dr / 2. - prtls->pos[i][0] * dr + prtls->pos[i][0] * prtls->pos[i][0]) / v_1;
      inc_tmpr(r_i, z_k, value);

      // weighting in ro[i + 1][k] cell
      value = ro_v * CYL_RNG_VOL(dz1, r2, r3) / v_2;
      inc_tmpr(r_i + 1,z_k, value);

      // weighting in ro[i][k + 1] cell
      value = ro_v*PI * dz2 * (dr * dr / 2. - prtls->pos[i][0] * dr + prtls->pos[i][0] * prtls->pos[i][0]) / v_1;
      inc_tmpr(r_i, z_k + 1, value);

      // weighting in ro[i + 1][k + 1] cell
      value = ro_v * CYL_RNG_VOL(dz2, r2, r3) / v_2;
      inc_tmpr(r_i + 1, z_k + 1, value);

    }
    else
    {
      ////////////////////////// /
      r1 = prtls->pos[i][0] - 0.5 * dr;
      r2 = (r_i + 0.5) * dr;
      r3 = prtls->pos[i][0] + 0.5 * dr;
      dz1 = (z_k + 1) * dz - prtls->pos[i][2];
      dz2 = prtls->pos[i][2] - z_k * dz;
      ro_v = prtls->mass_array[i] * vel2 / (4. * PI * dz * dr * prtls->pos[i][0]);
      v_1 = CYL_VOL(dz, dr);
      v_2 = CELL_VOLUME(r_i + 1, dr, dz);
      ////////////////////////// /

      // weighting in ro[i][k] cell
      value = ro_v * CYL_RNG_VOL(dz1, r1, r2) / v_1;
      inc_tmpr(r_i, z_k, value);

      // weighting in ro[i + 1][k] cell
      value = ro_v * CYL_RNG_VOL(dz1, r2, r3) / v_2;
      inc_tmpr(r_i + 1,z_k, value);

      // weighting in ro[i][k + 1] cell
      value = ro_v * CYL_RNG_VOL(dz2, r1, r2) / v_1;
      inc_tmpr(r_i, z_k + 1, value);

      // weighting in ro[i + 1][k + 1] cell
      value = ro_v * CYL_RNG_VOL(dz2, r2, r3) / v_2;
      inc_tmpr(r_i + 1, z_k + 1, value);
    }
  }

  for (unsigned int i = 0; i < geom->n_grid_r; i++)
    for (unsigned int j = 0; j < geom->n_grid_z; j++)
    {
      if (count[i][j] == 0) count[i][j] = 1;
      t[i][j] = t_sum[i][j] / count[i][j];
    }
}





void Temperature::calc_t_phi(Particles *prtls)
{
  reset();

  double dr = geom->dr;
  double dz = geom->dz;
  double r1, r2, r3; // temp variables for calculation
  double dz1, dz2; // temp var.: width of k and k + 1 cell

  double ro_v = 0; // charge density Q/V, V - volume of particle
  double v_1 = 0; // volume of [i][k] cell
  double v_2 = 0; // volume of [i + 1][k] cell
  double vel2 = 0;

  double value = 0;
  // double **temp = ro1->get_rho();

  for(unsigned int i = 0; i < prtls->number; i++)
  {
    // finding number of i and k cell. example: dr = 0.5; r = 0.4; i =0
    unsigned int r_i = CELL_NUMBER(prtls->pos[i][0], geom->dr);
    unsigned int z_k = CELL_NUMBER(prtls->pos[i][2], geom->dz);

    vel2 = prtls->vel[i][1] * prtls->vel[i][1];

    // in first cell other alg. of ro_v calc
    if (prtls->pos[i][0] > dr)
    {
      r1 =  prtls->pos[i][0] - 0.5 * dr;
      r2 = (r_i + 0.5) * dr;
      r3 = prtls->pos[i][0] + 0.5 * dr;
      ro_v = prtls->mass_array[i] * vel2 / (4. * PI * dz * dr * prtls->pos[i][0]);


      v_1 = CELL_VOLUME(r_i, dr, dz);
      v_2 = CELL_VOLUME(r_i + 1, dr, dz);
      dz1 = (z_k + 1) * dz - prtls->pos[i][2];
      dz2 = prtls->pos[i][2] - z_k * dz;

      // weighting in ro[i][k] cell
      value = ro_v * CYL_RNG_VOL(dz1, r1, r2) / v_1;
      inc_tmpr(r_i, z_k, value);

      // weighting in ro[i + 1][k] cell
      value = ro_v * CYL_RNG_VOL(dz1, r2, r3) / v_2;
      inc_tmpr(r_i+1, z_k, value);

      // weighting in ro[i][k + 1] cell
      value = ro_v * CYL_RNG_VOL(dz2, r1, r2) / v_1;
      inc_tmpr(r_i, z_k+1, value);

      // weighting in ro[i + 1][k + 1] cell
      value = ro_v * CYL_RNG_VOL(dz2, r2, r3) / v_2;
      inc_tmpr(r_i + 1, z_k + 1, value);
    }
    else if (prtls->pos[i][0] <= dr / 2.)
    {
      r_i = 0;
      r1 =  0.;
      r2 = (r_i + 0.5) * dr;
      r3 = prtls->pos[i][0] + 0.5 * dr;
      dz1 = (z_k + 1) * dz - prtls->pos[i][2];
      dz2 = prtls->pos[i][2] - z_k * dz;
      ro_v = prtls->mass_array[i] * vel2 / (2 * PI * dz * (2. * prtls->pos[i][0] * prtls->pos[i][0] + dr * dr / 2.));
      v_1 = CYL_VOL(dz, dr);
      v_2 = CELL_VOLUME(r_i + 1, dr, dz);
      ////////////////////////// /

      // weighting in ro[i][k] cell
      value = ro_v * PI * dz1 * (dr * dr / 2. - prtls->pos[i][0] * dr + prtls->pos[i][0] * prtls->pos[i][0]) / v_1;
      inc_tmpr(r_i, z_k, value);

      // weighting in ro[i + 1][k] cell
      value = ro_v * CYL_RNG_VOL(dz1, r2, r3) / v_2;
      inc_tmpr(r_i + 1,z_k, value);

      // weighting in ro[i][k + 1] cell
      value = ro_v*PI * dz2 * (dr * dr / 2. - prtls->pos[i][0] * dr + prtls->pos[i][0] * prtls->pos[i][0]) / v_1;
      inc_tmpr(r_i, z_k + 1, value);

      // weighting in ro[i + 1][k + 1] cell
      value = ro_v * CYL_RNG_VOL(dz2, r2, r3) / v_2;
      inc_tmpr(r_i + 1, z_k + 1, value);

    }
    else
    {
      ////////////////////////// /
      r1 = prtls->pos[i][0] - 0.5 * dr;
      r2 = (r_i + 0.5) * dr;
      r3 = prtls->pos[i][0] + 0.5 * dr;
      dz1 = (z_k + 1) * dz - prtls->pos[i][2];
      dz2 = prtls->pos[i][2] - z_k * dz;
      ro_v = prtls->mass_array[i] * vel2 / (4. * PI * dz * dr * prtls->pos[i][0]);
      v_1 = CYL_VOL(dz, dr);
      v_2 = CELL_VOLUME(r_i + 1, dr, dz);
      ////////////////////////// /

      // weighting in ro[i][k] cell
      value = ro_v * CYL_RNG_VOL(dz1, r1, r2) / v_1;
      inc_tmpr(r_i, z_k, value);

      // weighting in ro[i + 1][k] cell
      value = ro_v * CYL_RNG_VOL(dz1, r2, r3) / v_2;
      inc_tmpr(r_i + 1,z_k, value);

      // weighting in ro[i][k + 1] cell
      value = ro_v * CYL_RNG_VOL(dz2, r1, r2) / v_1;
      inc_tmpr(r_i, z_k + 1, value);

      // weighting in ro[i + 1][k + 1] cell
      value = ro_v * CYL_RNG_VOL(dz2, r2, r3) / v_2;
      inc_tmpr(r_i + 1, z_k + 1, value);
    }
  }

  for (unsigned int i = 0; i < geom->n_grid_r; i++)
    for (unsigned int j = 0; j < geom->n_grid_z; j++)
    {
      if (count[i][j] == 0) count[i][j] = 1;
      t[i][j] = t_sum[i][j] / count[i][j];
    }
}




void Temperature::calc_t_z(Particles *prtls)
{
  reset();

  double dr = geom->dr;
  double dz = geom->dz;
  double r1, r2, r3; // temp variables for calculation
  double dz1, dz2; // temp var.: width of k and k + 1 cell

  double ro_v = 0; // charge density Q/V, V - volume of particle
  double v_1 = 0; // volume of [i][k] cell
  double v_2 = 0; // volume of [i + 1][k] cell
  double vel2 = 0;

  double value = 0;
  // double **temp = ro1->get_rho();

  for(unsigned int i = 0; i < prtls->number; i++)
  {
    // finding number of i and k cell. example: dr = 0.5; r = 0.4; i =0
    unsigned int r_i = CELL_NUMBER(prtls->pos[i][0], geom->dr);
    unsigned int z_k = CELL_NUMBER(prtls->pos[i][2], geom->dz);

    vel2 = prtls->vel[i][2] * prtls->vel[i][2];

    // in first cell other alg. of ro_v calc
    if (prtls->pos[i][0] > dr)
    {
      r1 =  prtls->pos[i][0] - 0.5 * dr;
      r2 = (r_i + 0.5) * dr;
      r3 = prtls->pos[i][0] + 0.5 * dr;
      ro_v = prtls->mass_array[i] * vel2 / (4. * PI * dz * dr * prtls->pos[i][0]);


      v_1 = CELL_VOLUME(r_i, dr, dz);
      v_2 = CELL_VOLUME(r_i + 1, dr, dz);
      dz1 = (z_k + 1) * dz - prtls->pos[i][2];
      dz2 = prtls->pos[i][2] - z_k * dz;

      // weighting in ro[i][k] cell
      value = ro_v * CYL_RNG_VOL(dz1, r1, r2) / v_1;
      inc_tmpr(r_i, z_k, value);

      // weighting in ro[i + 1][k] cell
      value = ro_v * CYL_RNG_VOL(dz1, r2, r3) / v_2;
      inc_tmpr(r_i+1, z_k, value);

      // weighting in ro[i][k + 1] cell
      value = ro_v * CYL_RNG_VOL(dz2, r1, r2) / v_1;
      inc_tmpr(r_i, z_k+1, value);

      // weighting in ro[i + 1][k + 1] cell
      value = ro_v * CYL_RNG_VOL(dz2, r2, r3) / v_2;
      inc_tmpr(r_i + 1, z_k + 1, value);
    }
    else if (prtls->pos[i][0] <= dr / 2.)
    {
      r_i = 0;
      r1 =  0.;
      r2 = (r_i + 0.5) * dr;
      r3 = prtls->pos[i][0] + 0.5 * dr;
      dz1 = (z_k + 1) * dz - prtls->pos[i][2];
      dz2 = prtls->pos[i][2] - z_k * dz;
      ro_v = prtls->mass_array[i] * vel2 / (2 * PI * dz * (2. * prtls->pos[i][0] * prtls->pos[i][0] + dr * dr / 2.));
      v_1 = CYL_VOL(dz, dr);
      v_2 = CELL_VOLUME(r_i + 1, dr, dz);
      ////////////////////////// /

      // weighting in ro[i][k] cell
      value = ro_v * PI * dz1 * (dr * dr / 2. - prtls->pos[i][0] * dr + prtls->pos[i][0] * prtls->pos[i][0]) / v_1;
      inc_tmpr(r_i, z_k, value);

      // weighting in ro[i + 1][k] cell
      value = ro_v * CYL_RNG_VOL(dz1, r2, r3) / v_2;
      inc_tmpr(r_i + 1,z_k, value);

      // weighting in ro[i][k + 1] cell
      value = ro_v*PI * dz2 * (dr * dr / 2. - prtls->pos[i][0] * dr + prtls->pos[i][0] * prtls->pos[i][0]) / v_1;
      inc_tmpr(r_i, z_k + 1, value);

      // weighting in ro[i + 1][k + 1] cell
      value = ro_v * CYL_RNG_VOL(dz2, r2, r3) / v_2;
      inc_tmpr(r_i + 1, z_k + 1, value);

    }
    else
    {
      ////////////////////////// /
      r1 = prtls->pos[i][0] - 0.5 * dr;
      r2 = (r_i + 0.5) * dr;
      r3 = prtls->pos[i][0] + 0.5 * dr;
      dz1 = (z_k + 1) * dz - prtls->pos[i][2];
      dz2 = prtls->pos[i][2] - z_k * dz;
      ro_v = prtls->mass_array[i] * vel2 / (4. * PI * dz * dr * prtls->pos[i][0]);
      v_1 = CYL_VOL(dz, dr);
      v_2 = CELL_VOLUME(r_i + 1, dr, dz);
      ////////////////////////// /

      // weighting in ro[i][k] cell
      value = ro_v * CYL_RNG_VOL(dz1, r1, r2) / v_1;
      inc_tmpr(r_i, z_k, value);

      // weighting in ro[i + 1][k] cell
      value = ro_v * CYL_RNG_VOL(dz1, r2, r3) / v_2;
      inc_tmpr(r_i + 1,z_k, value);

      // weighting in ro[i][k + 1] cell
      value = ro_v * CYL_RNG_VOL(dz2, r1, r2) / v_1;
      inc_tmpr(r_i, z_k + 1, value);

      // weighting in ro[i + 1][k + 1] cell
      value = ro_v * CYL_RNG_VOL(dz2, r2, r3) / v_2;
      inc_tmpr(r_i + 1, z_k + 1, value);
    }
  }

  for (unsigned int i = 0; i < geom->n_grid_r; i++)
    for (unsigned int j = 0; j < geom->n_grid_z; j++)
    {
      if (count[i][j] == 0) count[i][j] = 1;
      t[i][j] = t_sum[i][j] / count[i][j];
    }
}



void Temperature::calc_t(Particles *prtls)
{
  reset();

  double dr = geom->dr;
  double dz = geom->dz;
  double r1, r2, r3; // temp variables for calculation
  double dz1, dz2; // temp var.: width of k and k + 1 cell

  double ro_v = 0; // charge density Q/V, V - volume of particle
  double v_1 = 0; // volume of [i][k] cell
  double v_2 = 0; // volume of [i + 1][k] cell
  double vel2 = 0;

  double value = 0;
  // double **temp = ro1->get_rho();

  for(unsigned int i = 0; i < prtls->number; i++)
  {
    // finding number of i and k cell. example: dr = 0.5; r = 0.4; i =0
    unsigned int r_i = CELL_NUMBER(prtls->pos[i][0], geom->dr);
    unsigned int z_k = CELL_NUMBER(prtls->pos[i][2], geom->dz);

    vel2 = prtls->vel[i][0] * prtls->vel[i][0]
      + prtls->vel[i][1] * prtls->vel[i][1]
      + prtls->vel[i][2] * prtls->vel[i][2];

    // in first cell other alg. of ro_v calc
    if (prtls->pos[i][0] > dr)
    {
      r1 =  prtls->pos[i][0] - 0.5 * dr;
      r2 = (r_i + 0.5) * dr;
      r3 = prtls->pos[i][0] + 0.5 * dr;
      ro_v = prtls->mass_array[i] * vel2 / (4. * PI * dz * dr * prtls->pos[i][0]);


      v_1 = CELL_VOLUME(r_i, dr, dz);
      v_2 = CELL_VOLUME(r_i + 1, dr, dz);
      dz1 = (z_k + 1) * dz - prtls->pos[i][2];
      dz2 = prtls->pos[i][2] - z_k * dz;

      // weighting in ro[i][k] cell
      value = ro_v * CYL_RNG_VOL(dz1, r1, r2) / v_1;
      inc_tmpr(r_i, z_k, value);

      // weighting in ro[i + 1][k] cell
      value = ro_v * CYL_RNG_VOL(dz1, r2, r3) / v_2;
      inc_tmpr(r_i+1, z_k, value);

      // weighting in ro[i][k + 1] cell
      value = ro_v * CYL_RNG_VOL(dz2, r1, r2) / v_1;
      inc_tmpr(r_i, z_k+1, value);

      // weighting in ro[i + 1][k + 1] cell
      value = ro_v * CYL_RNG_VOL(dz2, r2, r3) / v_2;
      inc_tmpr(r_i + 1, z_k + 1, value);
    }
    else if (prtls->pos[i][0] <= dr / 2.)
    {
      r_i = 0;
      r1 =  0.;
      r2 = (r_i + 0.5) * dr;
      r3 = prtls->pos[i][0] + 0.5 * dr;
      dz1 = (z_k + 1) * dz - prtls->pos[i][2];
      dz2 = prtls->pos[i][2] - z_k * dz;
      ro_v = prtls->mass_array[i] * vel2 / (2 * PI * dz * (2. * prtls->pos[i][0] * prtls->pos[i][0] + dr * dr / 2.));
      v_1 = CYL_VOL(dz, dr);
      v_2 = CELL_VOLUME(r_i + 1, dr, dz);
      ////////////////////////// /

      // weighting in ro[i][k] cell
      value = ro_v * PI * dz1 * (dr * dr / 2. - prtls->pos[i][0] * dr + prtls->pos[i][0] * prtls->pos[i][0]) / v_1;
      inc_tmpr(r_i, z_k, value);

      // weighting in ro[i + 1][k] cell
      value = ro_v * CYL_RNG_VOL(dz1, r2, r3) / v_2;
      inc_tmpr(r_i + 1,z_k, value);

      // weighting in ro[i][k + 1] cell
      value = ro_v*PI * dz2 * (dr * dr / 2. - prtls->pos[i][0] * dr + prtls->pos[i][0] * prtls->pos[i][0]) / v_1;
      inc_tmpr(r_i, z_k + 1, value);

      // weighting in ro[i + 1][k + 1] cell
      value = ro_v * CYL_RNG_VOL(dz2, r2, r3) / v_2;
      inc_tmpr(r_i + 1, z_k + 1, value);

    }
    else
    {
      ////////////////////////// /
      r1 = prtls->pos[i][0] - 0.5 * dr;
      r2 = (r_i + 0.5) * dr;
      r3 = prtls->pos[i][0] + 0.5 * dr;
      dz1 = (z_k + 1) * dz - prtls->pos[i][2];
      dz2 = prtls->pos[i][2] - z_k * dz;
      ro_v = prtls->mass_array[i] * vel2 / (4. * PI * dz * dr * prtls->pos[i][0]);
      v_1 = CYL_VOL(dz, dr);
      v_2 = CELL_VOLUME(r_i + 1, dr, dz);
      ////////////////////////// /

      // weighting in ro[i][k] cell
      value = ro_v * CYL_RNG_VOL(dz1, r1, r2) / v_1;
      inc_tmpr(r_i, z_k, value);

      // weighting in ro[i + 1][k] cell
      value = ro_v * CYL_RNG_VOL(dz1, r2, r3) / v_2;
      inc_tmpr(r_i + 1,z_k, value);

      // weighting in ro[i][k + 1] cell
      value = ro_v * CYL_RNG_VOL(dz2, r1, r2) / v_1;
      inc_tmpr(r_i, z_k + 1, value);

      // weighting in ro[i + 1][k + 1] cell
      value = ro_v * CYL_RNG_VOL(dz2, r2, r3) / v_2;
      inc_tmpr(r_i + 1, z_k + 1, value);
    }
  }

  for (unsigned int i = 0; i < geom->n_grid_r; i++)
    for (unsigned int j = 0; j < geom->n_grid_z; j++)
    {
      if (count[i][j] == 0) count[i][j] = 1;
      t[i][j] = t_sum[i][j] / count[i][j];
    }
}

void Temperature::inc_tmpr(unsigned int r, unsigned int z, double value)
{
  t_sum[r][z] = t_sum[r][z] + value;
  count[r][z]++;
}
