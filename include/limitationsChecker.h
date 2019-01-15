//! Class for reading configfile and calculate parameters,
//! requires for modeling

/* #include <vector> */

/* #include "tinyxml2.h" */

#include <string.h>

#include "constant.h"
#include "lib.h"
#include "geometry.h"
#include "parameters.h"

using namespace std;
using namespace constant;
// using namespace tinyxml2;

/* <probe type="frame" r_start="0" z_start="0" r_end="10" z_end="200" schedule="5" component="E_r" /> */

/* struct probe { */
/*   char *component; */
/*   unsigned int type; */
/*   unsigned int r_start; */
/*   unsigned int r_end; */
/*   unsigned int z_start; */
/*   unsigned int z_end; */
/*   unsigned int schedule; */
/* }; */

/* struct particle_specie { */
/*   char *name; */
/*   unsigned int mass; */
/*   int charge; */
/*   double number_macro; */
/*   double left_density; */
/*   double right_density; */
/*   double temperature; */
/* }; */

class LimitationsChecker
{
public:
  // LimitationsChecker(void);
  LimitationsChecker(Parameters *cfg);
  ~LimitationsChecker(void);

public:
  double plasma_freq;
  double debye_length;
  double plasma_density;
  double number_macro;

public:
  bool check_velocity_time_step();
  bool check_grid_size();
  bool check_system_size();

private:
  Parameters *cfg;
};
