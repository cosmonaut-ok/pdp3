#include "limitationsChecker.h"

using namespace std;

LimitationsChecker::LimitationsChecker(Parameters *config)
{
  cfg = config;

  double electron_density = 0;
  double electron_temperature = 0;


  for (auto i = cfg->particle_species.begin(); i != cfg->particle_species.end(); ++i)
  {
    number_macro += (*i).number_macro; // calculate total number of macroparticles

    if (strcmp((*i).name, "Electrons") == 0
        || strcmp((*i).name, "electrons") == 0
        || strcmp((*i).name, "electron") == 0)
    {
      electron_density = ((*i).left_density + (*i).right_density) / 2;
      electron_temperature = (*i).temperature;
    }
  }

  if (electron_density == 0 || electron_temperature == 0)
  {
    cerr << "ERROR! there is no electrons present in system. Exiting.";
    exit(1);
  }

  plasma_freq = sqrt(electron_density * EL_CHARGE * EL_CHARGE / (EL_MASS * EPSILON0));
  debye_length = sqrt(EPSILON0 * BOLTZMANN * electron_temperature / (4 * PI * electron_density * EL_CHARGE * EL_CHARGE));
}


bool check()
{

}

bool LimitationsChecker::check_velocity_time_step()
{
  bool status;

  double multiplicator = 100;

  double max_time_c = min(cfg->geom->dr / LIGHT_SPEED, cfg->geom->dz / LIGHT_SPEED);
  double max_time_wp = 2 / plasma_freq / multiplicator;

  double max_time = min(max_time_c, max_time_wp);

  if (cfg->time->delta_t > max_time)
  {
    cerr << "ERROR! too big time step: ``"
         << cfg->time->delta_t
         << " s.''. Should be less, than ``"
         << max_time
         << " s.''. Exiting" << endl;
    exit(1);
  }
  else
    status = true;

  return status;
}

bool LimitationsChecker::check_grid_size()
// implement d{r,z} < \lambda debye
{
  bool status;
  int debye_multiplicator = 10;

  if (cfg->geom->dr < debye_length * debye_multiplicator
      || cfg->geom->dr < debye_length * debye_multiplicator)
  {
    cerr << "ERROR! too small grid size: ``"
         << cfg->geom->dr << " x " << cfg->geom->dz
         << " m.''. Should be more, than ``"
         << debye_length * debye_multiplicator
         << " m.''. Exiting" << endl;
    exit(1);
  }
  else
    status = true;

  return status;
}

bool LimitationsChecker::check_system_size()
// implement d{r,z} < \lambda debye
{
  bool status = true;
  int debye_multiplicator = 100;

  if (cfg->geom->r_size * debye_multiplicator < debye_length * number_macro
      || cfg->geom->z_size * debye_multiplicator < debye_length * number_macro)
  {
    cerr << "ERROR! too small system size: ``"
         << cfg->geom->r_size << " x " << cfg->geom->z_size
         << " m.''. Should be more, than ``"
         <<  debye_length * number_macro / debye_multiplicator
         << " m.''. Exiting" << endl;
    exit(1);
  }

  return status;
}

bool LimitationsChecker::check_macro_number()
{
  bool status = true;
  int grid_multiplicator = 10;

  if (cfg->geom->r_size * cfg->geom->r_size * grid_multiplicator > number_macro)
  {
    cerr << "WARNING! too small summary number of macroparticles: ``"
         << number_macro
         << "''. Should be more, than ``"
         <<  cfg->geom->r_size * cfg->geom->r_size * grid_multiplicator
         << "''. You could get not relevant results." << endl;

    status = false;
  }

  return status;
}

LimitationsChecker::~LimitationsChecker(void)
{
}
