#include "lib.h"

using namespace constant;

namespace lib
{
  // C^2 define c^2 to decrease number of operations
  const double LIGHT_SPEED_POW_2 = pow (LIGHT_SPEED, 2);

#if __cplusplus >= 201103L
  using std::isnan;
#endif

#ifdef __SSE__
  float sqrt_recip(float x)
  {
    //! 1/sqrt(x), using MMX and Newton-Raphson approximation
    return _mm_cvtss_f32( _mm_rsqrt_ss( _mm_set_ps1(x) ) ); //same as _mm_set1_ps
  }
#endif

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

#ifdef __SSE__
    gamma = sqrt_recip(1.0 - beta);
#else
    gamma = pow(1.0 - beta, -0.5);
#endif

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
    double time_sec = std::time(nullptr) - SIMULATION_START_TIME;

    int d, hr, min;
    double sec;
    char* the_time = new char;

    int sec_in_min = 60;
    int sec_in_hr = 3600;
    int sec_in_day = 86400;

    if (time_sec > sec_in_min && time_sec < sec_in_hr)
    {
      min = (int)time_sec / sec_in_min;
      sec = time_sec - min * sec_in_min;
      sprintf(the_time, "%dm %.0fs", min, sec);
    }
    else if (time_sec > sec_in_hr && time_sec < sec_in_day)
    {
      hr = (int)time_sec / sec_in_hr;
      min = (int)((time_sec - hr * sec_in_hr) / sec_in_min);
      sec = time_sec - hr * sec_in_hr - min * sec_in_min;
      sprintf(the_time, "%dh %dm %.0fs", hr, min, sec);
    }
    else if (time_sec > sec_in_day)
    {
      d = (int)time_sec / sec_in_day;
      hr = (int)((time_sec - d * sec_in_day) / sec_in_hr);
      min = (int)((time_sec - d * sec_in_day - hr * sec_in_hr) / sec_in_min);
      sec = time_sec - d * sec_in_day - hr * sec_in_hr - min * sec_in_min;
      sprintf(the_time, "%dd %dh %dm %.0fs", d, hr, min, sec);
    }
    else
      sprintf(the_time, "%.0fs", time_sec);

    return the_time;
  }
}
