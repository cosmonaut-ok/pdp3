#pragma once

// #include <fstream.h>
#include <sstream>
#include <iostream>
#include <math.h>

using namespace std;

namespace lib
{
  bool to_bool(string str);
  double get_gamma (double velocity);
  double get_gamma_inv (double velocity);
}
