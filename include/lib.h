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

#ifdef __AVX__
#include <immintrin.h>//AVX Intrinsic Functions
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
  char* get_simulation_duration();
  bool directory_exists(const std::string& path);
  bool make_directory(const std::string& path);
  double sq_rt(double x);
  char *get_cmd_option(char **begin, char **end, const std::string &option);
  bool cmd_option_exists(char **begin, char **end, const std::string &option);
}
