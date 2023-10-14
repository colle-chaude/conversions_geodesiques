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

/// @brief Calcul de la latitude isométrique sur un ellipsoïde de première excentricité eau point de latitude ϕ.
/// @param _phi latitude en radian
/// @param _e première excentricité de l’ellipsoïde
/// @return latitude isométrique
double ALG0001(double _phi, double _e)
{
    return log( tan(PI_GEO/4+_phi/2) * pow( (1-_e*sin(_phi)) / (1+_e*sin(_phi)), _e/2));
}

/// @brief Calcul de la latitude ϕ à partir de la latitude isométrique L.
/// @param _l_iso latitude isométrique
/// @param _e première excentricité de l’ellipsoïde
/// @param _tol tolérance de convergence
/// @return latitude en radian
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

/// @brief 
/// @param _lambda longitude par rapport au méridien origine (Greenwich)
/// @param _phi latitude
/// @param _n exposant de la projection
/// @param _e première excentricité de l’ellipsoïde
/// @param _c constante de la projection
/// @param _lamc longitude de l’origine par rapport au méridien origine (Greenwich)
/// @param _Xs coordonnée X en projection du pôle
/// @param _Ys coordonnée Y en projection du pôle
/// @param _X coordonnée X de lambert (sortie)
/// @param _Y coordonnée Y de lambert (sortie)
void ALG0003(double _lambda, double _phi, double _n, double _e, double _c, double _lamc, double _Xs, double _Ys, double* _X, double* _Y)
{
    double l_iso = ALG0001(_phi, _e);
    *_X = _Xs + _c * exp(-_n*l_iso)*sin(_n*(_lambda - _lamc));
    *_Y = _Ys - _c * exp(-_n*l_iso)*cos(_n*(_lambda - _lamc));
}

/// @brief 
/// @param _X coordonnée X de lambert
/// @param _Y coordonnée Y de lambert
/// @param _n exposant de la projection
/// @param _e première excentricité de l’ellipsoïde
/// @param _c constante de la projection
/// @param _lamc longitude de l’origine par rapport au méridien origine (Greenwich)
/// @param _Xs coordonnée X en projection du pôle
/// @param _Ys coordonnée Y en projection du pôle$
/// @param _tol 
/// @param _lambda longitude par rapport au méridien origine (Greenwich) (sortie)
/// @param _phi latitude (sortie)
void ALG0004(double _X, double _Y, double _n, double _e, double _c, double _lamc, double _Xs, double _Ys, double _tol, double* _lambda, double* _phi)
{
    double R = sqrt((_X-_Xs)*(_X-_Xs) + (_Y-_Ys)*(_Y-_Ys) );
    double gama= atan((_X-_Xs)/(_Ys-_Y));

    *_lambda = _lamc + gama/_n ;

    double l_iso = -1/_n * log(fabs(R/_c));

    *_phi = ALG0002(l_iso, _e, _tol);
}



/// @brief convertir coordonées Lambert en coordonées géographiques longitude/latitude WGS84
/// @param _X coordonée X lambert
/// @param _Y coordonée Y lambert
/// @param _longitude longitude en degrés decimaux (Sortie)
/// @param _latitude  Latitude en degrés décimaux (Sortie)
/// @param lambert_type projection lambert utilisée (I, II, III, IV, 93, II Etendue)
void lambert_to_geo(double _X, double _Y, double* _longitude, double* _latitude, int lambert_type)
{
    ALG0004(_X, _Y, n[lambert_type], e, c[lambert_type], lamc, Xs[lambert_type], Ys[lambert_type], tol, _longitude, _latitude); 
    *_longitude *= 180.0/PI_GEO;
    *_latitude *= 180.0/PI_GEO;
}

/// @brief convertir coordonées géographiques longitude/latitude WGS84 en coordonées Lambert
/// @param _longitude longitude en degrés decimaux
/// @param _latitude Latitude en degrés decimaux
/// @param _X coordonée X lambert (Sortie)
/// @param _Y coordonée Y lambert (Sortie)
/// @param lambert_type 
void geo_to_lambert(double _longitude, double _latitude, double* _X, double* _Y, int lambert_type)
{
    ALG0003(_longitude*(PI_GEO/180.0), _latitude*(PI_GEO/180.0), n[lambert_type], e, c[lambert_type], lamc, Xs[lambert_type], Ys[lambert_type], _X, _Y); 
}


}//namespace geo_lib