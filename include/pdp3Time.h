#pragma once
#include "particlesStruct.h"

class Time
{
public:
  double current_time;
  double start_time;
  double relaxation_time;
  double end_time;
  double delta_t;
public:
  Time(void);
  Time(double ct, double st, double rt, double et, double dt);

  ~Time(void);
};
