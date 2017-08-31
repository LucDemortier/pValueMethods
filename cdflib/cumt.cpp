# include <iomanip>
# include <cmath>
using namespace std;
# include "cdflib.hpp"

//****************************************************************************80

void cumt ( double *t, double *df, double *cum, double *ccum )

//****************************************************************************80
//
//  Purpose:
//
//    CUMT evaluates the cumulative T distribution.
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    Formula 26.5.27.
//
//  Parameters:
//
//    Input, double *T, the upper limit of integration.
//
//    Input, double *DF, the number of degrees of freedom of
//    the T distribution.
//
//    Output, double *CUM, *CCUM, the T distribution CDF and
//    complementary CDF.
//
{
  static double a;
  static double dfptt;
  static double K2 = 0.5e0;
  static double oma;
  static double T1;
  static double tt;
  static double xx;
  static double yy;

  tt = (*t) * (*t);
  dfptt = ( *df ) + tt;
  xx = *df / dfptt;
  yy = tt / dfptt;
  T1 = 0.5e0 * ( *df );
  cumbet ( &xx, &yy, &T1, &K2, &a, &oma );

  if ( *t <= 0.0e0 )
  {
    *cum = 0.5e0 * a;
    *ccum = oma + ( *cum );
  }
  else
  {
    *ccum = 0.5e0 * a;
    *cum = oma + ( *ccum );
  }
  return;
}
