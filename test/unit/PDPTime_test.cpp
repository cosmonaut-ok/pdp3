#include <gtest/gtest.h>
#include "pdp3_time.h"

// using namespace std;

double current_time = 1e-8;
double start_time = 0;
double relaxation_time = 0;
double end_time = 5e-8;
double delta_t = 1e-12;
  
Time* pdp_time = new Time(current_time, start_time, relaxation_time, end_time, delta_t);

TEST(Time, init)
{
  EXPECT_TRUE(2 == 2);
}

int main(int ac, char* av[])
{
  testing::InitGoogleTest(&ac, av);
  return RUN_ALL_TESTS();
}
