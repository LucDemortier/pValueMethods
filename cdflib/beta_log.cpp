# include <iomanip>
# include <cmath>
using namespace std;
# include "cdflib.hpp"

//****************************************************************************80

double beta_log ( double *a0, double *b0 )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_LOG evaluates the logarithm of the beta function.
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
//    Input, double *A0, *B0, the parameters of the function.
//    A0 and B0 should be nonnegative.
//
//    Output, double *BETA_LOG, the value of the logarithm
//    of the Beta function.
//
{
  static double e = .918938533204673e0;
  static double value,a,b,c,h,u,v,w,z;
  static int i,n;
  static double T1;

    a = fifdmin1(*a0,*b0);
    b = fifdmax1(*a0,*b0);
    if(a >= 8.0e0) goto S100;
    if(a >= 1.0e0) goto S20;
//
//  PROCEDURE WHEN A .LT. 1
//
    if(b >= 8.0e0) goto S10;
    T1 = a+b;
    value = gamma_log ( &a )+( gamma_log ( &b )- gamma_log ( &T1 ));
    return value;
S10:
    value = gamma_log ( &a )+algdiv(&a,&b);
    return value;
S20:
//
//  PROCEDURE WHEN 1 .LE. A .LT. 8
//
    if(a > 2.0e0) goto S40;
    if(b > 2.0e0) goto S30;
    value = gamma_log ( &a )+ gamma_log ( &b )-gsumln(&a,&b);
    return value;
S30:
    w = 0.0e0;
    if(b < 8.0e0) goto S60;
    value = gamma_log ( &a )+algdiv(&a,&b);
    return value;
S40:
//
//  REDUCTION OF A WHEN B .LE. 1000
//
    if(b > 1000.0e0) goto S80;
    n = ( int ) ( a - 1.0e0 );
    w = 1.0e0;
    for ( i = 1; i <= n; i++ )
    {
        a -= 1.0e0;
        h = a/b;
        w *= (h/(1.0e0+h));
    }
    w = log(w);
    if(b < 8.0e0) goto S60;
    value = w+ gamma_log ( &a )+algdiv(&a,&b);
    return value;
S60:
//
//  REDUCTION OF B WHEN B .LT. 8
//
    n = ( int ) ( b - 1.0e0 );
    z = 1.0e0;
    for ( i = 1; i <= n; i++ )
    {
        b -= 1.0e0;
        z *= (b/(a+b));
    }
    value = w+log(z)+( gamma_log ( &a )+( gamma_log ( &b )-gsumln(&a,&b)));
    return value;
S80:
//
//  REDUCTION OF A WHEN B .GT. 1000
//
    n = ( int ) ( a - 1.0e0 );
    w = 1.0e0;
    for ( i = 1; i <= n; i++ )
    {
        a -= 1.0e0;
        w *= (a/(1.0e0+a/b));
    }
    value = log(w)-(double)n*log(b)+( gamma_log ( &a )+algdiv(&a,&b));
    return value;
S100:
//
//  PROCEDURE WHEN A .GE. 8
//
    w = bcorr(&a,&b);
    h = a/b;
    c = h/(1.0e0+h);
    u = -((a-0.5e0)*log(c));
    v = b*alnrel(&h);
    if(u <= v) goto S110;
    value = -(0.5e0*log(b))+e+w-v-u;
    return value;
S110:
    value = -(0.5e0*log(b))+e+w-u-v;
    return value;
}
