# include <iomanip>
# include <cmath>
using namespace std;
# include "cdflib.hpp"

//****************************************************************************80

double fpser ( double *a, double *b, double *x, double *eps )

//****************************************************************************80
//
//  Purpose:
//
//    FPSER evaluates IX(A,B)(X) for very small B.
//
//  Discussion:
//
//    This routine is appropriate for use when
//
//      B < min ( EPS, EPS * A )
//
//    and
//
//      X <= 0.5.
//
//  Parameters:
//
//    Input, double *A, *B, parameters of the function.
//
//    Input, double *X, the point at which the function is to
//    be evaluated.
//
//    Input, double *EPS, a tolerance.
//
//    Output, double FPSER, the value of IX(A,B)(X).
//
{
  static int K1 = 1;
  static double fpser,an,c,s,t,tol;

    fpser = 1.0e0;
    if(*a <= 1.e-3**eps) goto S10;
    fpser = 0.0e0;
    t = *a*log(*x);
    if(t < exparg(&K1)) return fpser;
    fpser = exp(t);
S10:
//
//                NOTE THAT 1/B(A,B) = B
//
    fpser = *b/ *a*fpser;
    tol = *eps/ *a;
    an = *a+1.0e0;
    t = *x;
    s = t/an;
S20:
    an += 1.0e0;
    t = *x*t;
    c = t/an;
    s += c;
    if(fabs(c) > tol) goto S20;
    fpser *= (1.0e0+*a*s);
    return fpser;
}
