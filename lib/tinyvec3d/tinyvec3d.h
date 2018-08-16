#ifndef TINYVEC3D_H_INCLUDED
#define TINYVEC3D_H_INCLUDED
#endif

#include <sstream>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <cmath>

extern "C"
{
  namespace tinyvec3d
  {
    double* mkvector3d(double x0, double x1, double x2);
    double* copy(double* v_this);
    void copy_components(double* v_this, double* v_other);
    void cross(double* v_this, double* v_other);
    double dot(double* v_this, double* v_other);
    void add(double* v_this, double* v_other);
    void subst(double* v_this, double* v_other);
    void plus(double* v_this, double number);
    void minus(double* v_this, double number);
    void product(double* v_this, double number);
    void div(double* v_this, double number);
    void power(double* v_this, double number);
    void abs(double* v_this);
    void abs2(double* v_this);
    double squared_sum(double* v_this);
  }
}
