#include "conv_geo.h"

namespace geo_lib{


// Valeurs des constantes Lambert France :
//          |   Lambert I       |   Lambert II      |   Lambert III     |   Lambert IV      |   Lambert-93
// n        |  0.760 405 965 6  |  0,728 968 627 4  |  0,695 912 796 6  |  0,671 267 932 2  |  0,725 607 765 0
// c (m)    |  11 603 796,98    |  11 745 793,39    |  11 947 992,52    |  12 136 281,99    |  11 754 255,426
// Xs (m)   |  600 000,0        |  600 000,0        |  600 000,0        |  234,358          |  700 000,0
// Ys (m)   |  5 657 616,674    |  6 199 695,768    |  6 791 905,085    |  7 239 161,542    |  12 655 612,050
// λ0 = 2°20’14,025” E par rapport au méridien de Greenwich
// = 0 grades par rapport au méridien de Paris.
// e = 0,082 483 256 76 (première excentricité de l’ellipsoïde Clarke 1880
// français).
// Les constantes en usage pour le Lambert II étendu sont celles du
// Lambert II avec Ys = 8 199 695,768.

    //              Lambert I       Lambert II      Lambert III     Lambert IV      Lambert_93    Lambert II etendu
const double n[6]  ={0.7604059656,  0.7289686274,   0.6959127966,   0.6712679322,   0.7256077650,   0.7289686274};
const double c[6]  ={11603796.98,   11745793.39,    11947992.52,    12136281.99,    11754255.426,   11745793.39};
const double Xs[6] ={600000.0,      600000.0,       600000.0,       234.358,        700000.0,       600000.0};
const double Ys[6] ={5657616.674,   6199695.768,    6791905.085,    7239161.542,    12655612.050,   8199695.768};
const double lamc = 0.04079234433;
const double tol = pow(10, -11);
const double e = 0.08248325676 ;

double ALG0002(double _l_iso, double _e, double _tol)
{
    double phi0 = 2* atan(exp(_l_iso)) - PI_GEO/2;
    double phi1 = phi0 + 2;
    int i = 0;
    
    //( (phi0 - phi1) > _tol || (phi0 - phi1) < -_tol )
    while( fabs(phi0 - phi1) && i < 5000)
    {
        phi1 = phi0;
        phi0 = 2* atan(
                        pow(
                            (1+_e*sin(phi1))/(1-_e*sin(phi1))
                            ,_e/2)
                        *exp(_l_iso)
                        ) 
                    - PI_GEO/2;

    }
    return phi0;
}

void ALG0004(double _n, double _e, double _c, double _lamc, double _Xs, double _Ys, double _X, double _Y, double _tol, double* _lambda, double* _phi)
{
    double R = sqrt((_X-_Xs)*(_X-_Xs) + (_Y-_Ys)*(_Y-_Ys) );
    double gama= atan((_X-_Xs)/(_Ys-_Y));

    *_lambda = _lamc + gama/_n ;

    double l_iso = -1/_n * log(fabs(R/_c));

    *_phi = ALG0002(l_iso, _e, _tol);
}




void lambert_to_geo(double _X, double _Y, double* _longitude, double* _latitude, int lambert_type)
{
    ALG0004(n[lambert_type], e, c[lambert_type], lamc, Xs[lambert_type], Ys[lambert_type], _X, _Y, tol, _longitude, _latitude); 
}


}//namespace geo_lib