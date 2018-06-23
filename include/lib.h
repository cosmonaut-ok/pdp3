#include <sstream>
#include <iostream>
#include <math.h>
#include <algorithm>

#include "constant.h"

using namespace std;

namespace lib
{
  bool to_bool(string str);
  double get_gamma (double velocity);
  double get_gamma_inv (double velocity);
  double random_reverse(double vel, int power);
}
