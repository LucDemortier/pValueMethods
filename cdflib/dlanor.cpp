# include <iomanip>
# include <cmath>
using namespace std;
# include "cdflib.hpp"

//****************************************************************************80

double dlanor ( double *x )

//****************************************************************************80
//
//  Purpose:
//
//    DLANOR evaluates the logarithm of the asymptotic Normal CDF.
//
//  Discussion:
//
//    This routine computes the logarithm of the cumulative normal distribution
//    from abs ( x ) to infinity for  5 <= abs ( X ).
//
//    The relative error at X = 5 is about 0.5D-5.
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.2.12.
//
//  Parameters:
//
//    Input, double *X, the value at which the Normal CDF is to be
//    evaluated.  It is assumed that 5 <= abs ( X ).
//
//    Output, double DLANOR, the logarithm of the asymptotic
//    Normal CDF.
//
{
# define dlsqpi 0.91893853320467274177e0

  static double coef[12] = {
    -1.0e0,3.0e0,-15.0e0,105.0e0,-945.0e0,10395.0e0,-135135.0e0,2027025.0e0,
    -34459425.0e0,654729075.0e0,-13749310575.e0,316234143225.0e0
  };
  static int K1 = 12;
  static double dlanor,approx,correc,xx,xx2,T2;

  xx = fabs(*x);
  if ( xx < 5.0e0 )
  {
    ftnstop(" Argument too small in DLANOR");
  }
  approx = -dlsqpi-0.5e0*xx*xx-log(xx);
  xx2 = xx*xx;
  T2 = 1.0e0/xx2;
  correc = eval_pol ( coef, &K1, &T2 ) / xx2;
  correc = alnrel ( &correc );
  dlanor = approx+correc;
  return dlanor;
# undef dlsqpi
}
