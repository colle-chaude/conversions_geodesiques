#ifndef CONV_GEO_H_INCLUDED
#define CONV_GEO_H_INCLUDED

#include <cmath>

namespace geo_lib
{

#define LAMBERT_I_TYPE 0
#define LAMBERT_II_TYPE 1
#define LAMBERT_III_TYPE 2
#define LAMBERT_IV_TYPE 3
#define LAMBERT_93_TYPE 4
#define LAMBERT_IT_EXT_TYPE 5

#define PI_GEO 3.1415926535897931

double ALG0002(double _l_iso, double _e, double _tol);
void ALG0004(double _n, double _e, double _c, double _lamc, double _Xs, double _Ys, double _X, double _Y, double _tol, double* _lambda, double* _phi);
void lambert_to_geo(double _X, double _Y, double* _longitude, double* _latitude, int lambert_type);

}//namespace geo_lib
#endif // CONV_GEO_H_INCLUDED
