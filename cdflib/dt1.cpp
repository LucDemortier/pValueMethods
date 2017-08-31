# include <iomanip>
# include <cmath>
using namespace std;
# include "cdflib.hpp"

//****************************************************************************80

double dt1 ( double *p, double *q, double *df )

//****************************************************************************80
//
//  Purpose:
//
//    DT1 computes an approximate inverse of the cumulative T distribution.
//
//  Discussion:
//
//    Returns the inverse of the T distribution function, i.e.,
//    the integral from 0 to INVT of the T density is P. This is an
//    initial approximation.
//
//  Parameters:
//
//    Input, double *P, *Q, the value whose inverse from the
//    T distribution CDF is desired, and the value (1-P).
//
//    Input, double *DF, the number of degrees of freedom of the
//    T distribution.
//
//    Output, double DT1, the approximate value of X for which
//    the T density CDF with DF degrees of freedom has value P.
//
{
  static double coef[4][5] = {
    {1.0e0,1.0e0,0.0e0,0.0e0,0.0e0},{3.0e0,16.0e0,5.0e0,0.0e0,0.0e0},{-15.0e0,17.0e0,
    19.0e0,3.0e0,0.0e0},{-945.0e0,-1920.0e0,1482.0e0,776.0e0,79.0e0}
  };
  static double denom[4] = {
    4.0e0,96.0e0,384.0e0,92160.0e0
  };
  static int ideg[4] = {
    2,3,4,5
  };
  static double dt1,denpow,sum,term,x,xp,xx;
  static int i;

    x = fabs(dinvnr(p,q));
    xx = x*x;
    sum = x;
    denpow = 1.0e0;
    for ( i = 0; i < 4; i++ )
    {
        term = eval_pol ( &coef[i][0], &ideg[i], &xx ) * x;
        denpow *= *df;
        sum += (term/(denpow*denom[i]));
    }
    if(!(*p >= 0.5e0)) goto S20;
    xp = sum;
    goto S30;
S20:
    xp = -sum;
S30:
    dt1 = xp;
    return dt1;
}
