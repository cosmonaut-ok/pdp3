#include <gtest/gtest.h>
#include "pdp3Time.h"
#include "pdp3Time.cpp"

// using namespace std;

double current_time = 1e-8;
double start_time = 0;
double relaxation_time = 0;
double end_time = 5e-8;
double delta_t = 1e-12;

TEST(Time, time_init)
{
  Time* pdp_time = new Time(current_time, start_time, relaxation_time, end_time, delta_t);

  EXPECT_TRUE(current_time == pdp_time->current_time);
  EXPECT_TRUE(start_time == pdp_time->start_time);
  EXPECT_TRUE(relaxation_time == pdp_time->relaxation_time);
  EXPECT_TRUE(end_time == pdp_time->end_time);
  EXPECT_TRUE(delta_t == pdp_time->delta_t);
}
