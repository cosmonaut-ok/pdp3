#pragma once

#include <sstream>
#include <iostream>
#include <math.h>

namespace math::triple_vector
{
  Triple cross(Triple a, Triple b);
  Triple sum(Triple a, Triple b);
  Triple subst(Triple a, Triple b);
  double dot(Triple a, Triple b);
  Triple product(Triple a, double b);
  Triple pow(Triple a, double exp);
  Triple abs(Triple a);
  Triple abs2(Triple a);
}
