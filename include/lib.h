#include <sstream>
#include <iostream>
#include <string>
#include <sys/stat.h> // stat
#include <errno.h> // errno, ENOENT, EEXIST
#include <math.h>
#include <algorithm>
#include <ctime>

#ifdef __SSE__
#include <pmmintrin.h>
#endif

#ifdef _WIN32
#include <direct.h> // _mkdir
#endif

#include "constant.h"
#include "config.h"

using namespace std;

namespace lib
{
  float sqrt_recip(float x);
  bool to_bool(string str);
  double get_gamma (double velocity);
  double get_gamma_inv (double velocity);
  double random_reverse(double vel, int power);
  char* get_simulation_time();
  bool directoryExists(const std::string& path);
  bool makeDirectory(const std::string& path);
}
