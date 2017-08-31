# include <iomanip>
# include <cmath>
using namespace std;
# include "cdflib.hpp"

//****************************************************************************80

double gsumln ( double *a, double *b )

//****************************************************************************80
//
//  Purpose:
//
//    GSUMLN evaluates the function ln(Gamma(A + B)).
//
//  Discussion:
//
//    GSUMLN is used for 1 <= A <= 2 and 1 <= B <= 2
//
//  Parameters:
//
//    Input, double *A, *B, values whose sum is the argument of
//    the Gamma function.
//
//    Output, double GSUMLN, the value of ln(Gamma(A+B)).
//
{
  static double gsumln,x,T1,T2;

    x = *a+*b-2.e0;
    if(x > 0.25e0) goto S10;
    T1 = 1.0e0+x;
    gsumln = gamma_ln1 ( &T1 );
    return gsumln;
S10:
    if(x > 1.25e0) goto S20;
    gsumln = gamma_ln1 ( &x ) + alnrel ( &x );
    return gsumln;
S20:
    T2 = x-1.0e0;
    gsumln = gamma_ln1 ( &T2 ) + log ( x * ( 1.0e0 + x ) );
    return gsumln;
}
