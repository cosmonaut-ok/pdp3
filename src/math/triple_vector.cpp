// some triple vector math functions. Inspired by picongpu's libpmacc

#include <cmath>
#include "math/triple.h"
#include "math/triple_vector.h"

// using namespace std;

namespace math
{
  namespace triple_vector
  {

    Triple cross(Triple a, Triple b)
    {
      Triple res(0.0, 0.0, 0.0);

      res.first = a.second * b.third - a.third * b.second;
      res.second = a.third * b.first - a.first * b.third;
      res.third = a.first * b.second - a.second * b.first;

      return res;
    }

    Triple sum(Triple a, Triple b)
    {
      Triple res(0.0, 0.0, 0.0);

      res.first = a.first + b.first;
      res.second = a.second + b.second;
      res.third = a.third + b.third;

      return res;
    }

    Triple subst(Triple a, Triple b)
    {
      Triple res(0.0, 0.0, 0.0);

      res.first = a.first - b.first;
      res.second = a.second - b.second;
      res.third = a.third - b.third;

      return res;
    }

    double dot(Triple a, Triple b)
    {
      double res = 0.0;

      res = a.first * b.first + a.second * b.second + a.third * b.third;
    }

    Triple product(Triple a, double b)
    {
      Triple res(0.0, 0.0, 0.0);

      res.first = a.first * b;
      res.second = a.second * b;
      res.third = a.third * b;

      return res;
    }

    Triple pow(Triple a, double exp)
    {
      Triple res(0.0, 0.0, 0.0);

      res.first = std::pow(a.first, exp);
      res.second = std::pow(a.second, exp);
      res.third = std::pow(a.third, exp);

      return res;
    }

    Triple abs(Triple a)
    {
      Triple res(0.0, 0.0, 0.0);

      res.first = std::abs(a.first);
      res.second = std::abs(a.second);
      res.third = std::abs(a.third);

      return res;
    }

    Triple abs2(Triple a)
    {
      Triple res(0.0, 0.0, 0.0);

      res.first = a.first * a.first;
      res.second = a.second * a.second;
      res.third = a.third * a.third;

      return res;
    }
  }
}
