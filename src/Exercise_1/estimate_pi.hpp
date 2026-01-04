#ifndef ESTIMATE_PI_HPP
#define ESTIMATE_PI_HPP
#define PI 3.141592653589793238462643383279502884L
#include "../General/Statistics.hpp"

Statistics random_point_method(double M, double N);

Statistics buffons_method(double M, double N, double t, double l);

#endif
