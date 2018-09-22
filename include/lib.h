#include <sstream>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <ctime>

#ifdef __SSE__
#include <pmmintrin.h>
#endif

#include "constant.h"

using namespace std;

namespace lib
{
  float sqrt_recip(float x);
  bool to_bool(string str);
  double get_gamma (double velocity);
  double get_gamma_inv (double velocity);
  double random_reverse(double vel, int power);
  char* get_simulation_time();
}
