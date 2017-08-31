# include <iomanip>
# include <cmath>
using namespace std;
# include "cdflib.hpp"

//****************************************************************************80

double gamma_log ( double *a )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_LOG evaluates ln ( Gamma ( A ) ) for positive A.
//
//  Author:
//
//    Alfred H Morris, Jr,
//    Naval Surface Weapons Center,
//    Dahlgren, Virginia.
//
//  Reference:
//
//    Armido DiDinato and Alfred Morris,
//    Algorithm 708:
//    Significant Digit Computation of the Incomplete Beta Function Ratios,
//    ACM Transactions on Mathematical Software,
//    Volume 18, 1993, pages 360-373.
//
//  Parameters:
//
//    Input, double *A, the argument of the function.
//    A should be positive.
//
//    Output, double GAMMA_LOG, the value of ln ( Gamma ( A ) ).
//
{
  static double c0 = .833333333333333e-01;
  static double c1 = -.277777777760991e-02;
  static double c2 = .793650666825390e-03;
  static double c3 = -.595202931351870e-03;
  static double c4 = .837308034031215e-03;
  static double c5 = -.165322962780713e-02;
  static double d = .418938533204673e0;
  static double gamln,t,w;
  static int i,n;
  static double T1;

    if(*a > 0.8e0) goto S10;
    gamln = gamma_ln1 ( a ) - log ( *a );
    return gamln;
S10:
    if(*a > 2.25e0) goto S20;
    t = *a-0.5e0-0.5e0;
    gamln = gamma_ln1 ( &t );
    return gamln;
S20:
    if(*a >= 10.0e0) goto S40;
    n = ( int ) ( *a - 1.25e0 );
    t = *a;
    w = 1.0e0;
    for ( i = 1; i <= n; i++ )
    {
        t -= 1.0e0;
        w = t*w;
    }
    T1 = t-1.0e0;
    gamln = gamma_ln1 ( &T1 ) + log ( w );
    return gamln;
S40:
    t = pow(1.0e0/ *a,2.0);
    w = (((((c5*t+c4)*t+c3)*t+c2)*t+c1)*t+c0)/ *a;
    gamln = d+w+(*a-0.5e0)*(log(*a)-1.0e0);
    return gamln;
}
