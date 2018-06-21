#include <stdio.h>
#include <stdlib.h>

#include <iostream>

#include <algorithm>
#include <iomanip>
#include "constant.h"

#include <math.h>
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
    double gamma, beta;

    beta = sq_velocity / LIGHT_SPEED_POW_2;

    if (beta > 1) // it's VERY BAD! Beta should not be more, than 1
    {
      cerr << "CRITICAL!(get_gamma_inv): Lorentz factor aka gamma is complex. Can not continue." << endl
           << "\tvelocity is: " << sqrt(sq_velocity) << endl;
      exit(1);
    }

    // gamma = pow(1.0 - beta, -0.5);
    gamma = pow(1.0 + beta, 0.5);

    if (isinf(gamma) == 1) { // avoid infinity values
      cerr << "WARNING(get_gamma_inv): gamma (Lorenz factor) girects to infinity for velocity" << sqrt(sq_velocity) << endl;
      return 1e100; // just return some very big value
    }
    else
    {
      return gamma;
    }
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

}
