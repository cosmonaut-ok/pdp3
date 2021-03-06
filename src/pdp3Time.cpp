#include "pdp3Time.h"

Time::Time(void)
{
}

Time::~Time(void)
{
}

Time::Time(double ct, double st, double rt, double et, double dt)
{
  start_time = st;
  relaxation_time = rt;
  current_time = ct;
  end_time = et;
  delta_t = dt;
}
