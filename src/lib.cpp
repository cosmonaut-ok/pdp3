#include "lib.h"

using namespace std;
using namespace constant;

namespace lib
{
  // C^2 define c^2 to decrease number of operations
  const double LIGHT_SPEED_POW_2 = pow (LIGHT_SPEED, 2);

  bool to_bool(string str)
  {
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    std::istringstream is(str);
    bool b;
    is >> std::boolalpha >> b;
    return b;
  }

  double get_gamma (double sq_velocity)
  {
    //! get_gamma takes squared velocity,
    //! only because of some features of code
    //! and optimisation issues
    double gamma, beta;

    beta = sq_velocity / LIGHT_SPEED_POW_2;

    if (beta > 1) // it's VERY BAD! Beta should not be more, than 1
    {
      cerr << "CRITICAL!(get_gamma): Lorentz factor aka gamma is complex. Can not continue." << endl
           << "\tvelocity is: " << sqrt(sq_velocity) << endl;
      exit(1);
    }

    gamma = pow(1.0 - beta, -0.5);

    if (isinf(gamma) == 1) { // avoid infinity values
      cerr << "WARNING(get_gamma): gamma (Lorenz factor) girects to infinity for velocity" << sqrt(sq_velocity) << endl;
      return 1e100; // just return some very big value
    }
    else
    {
      return gamma;
    }
  }

  double get_gamma_inv (double sq_velocity)
  {
    //! get_gamma_inv takes squared velocity,
    //! only because of some features of code
    //! and optimisation issues
    double gamma = pow(1.0 + sq_velocity / LIGHT_SPEED_POW_2, 0.5);

    return gamma;
  }

  double random_reverse(double vel, int power)
  {
    int int_vel =(int) floor(vel);
    double ost = 0;
    double r = 0;
    int order = 1;
    while(int_vel >= 1)
    {
      ost = int_vel % power;

      r = r + ost * pow((double)power, (-order));

      int_vel = (int_vel - ost)/power;
      order = order+1;
    }
    return r;
  }

  char* get_simulation_time()
  // get the time, spent since simulation launched (in "d h m s" format)
  {
    double time_sec = (clock() - SIMULATION_START_TIME) / (double)CLOCKS_PER_SEC;

    int d, hr, min;
    double sec;
    char* the_time = new char;

    if (time_sec > 60 && time_sec < 3600)
    {
      min = (int)time_sec / 60;
      sec = time_sec - min * 60;
      sprintf(the_time, "%dm %.2fs", min, sec);
    }
    else if (time_sec > 3600 && time_sec < 86400)
    {
      hr = (int)time_sec / 60 / 60;
      min = time_sec - hr * 60;
      sec = time_sec - hr * 60 * 60 - min * 60;
      sprintf(the_time, "%dh %dm %.1fs", hr, min, sec);
    }
    else if (time_sec > 86400)
    {
      d = (int)time_sec / 60 / 60 / 24;
      hr = time_sec - d * 24;
      min = time_sec - d * 24 - hr * 60;
      sec = time_sec - d * 24 - hr * 60 * 60 - min * 60;
      sprintf(the_time, "%dd %dh %dm %.0fs", d, hr, min, sec);
    }
    else
      sprintf(the_time, "%.2fs", time_sec);

    return the_time;
  }
}
