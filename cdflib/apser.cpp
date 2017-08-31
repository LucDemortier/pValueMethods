# include <iomanip>
# include <cmath>
using namespace std;
# include "cdflib.hpp"

//****************************************************************************80

double apser ( double *a, double *b, double *x, double *eps )

//****************************************************************************80
//
//  Purpose:
//
//    APSER computes the incomplete beta ratio I(SUB(1-X))(B,A).
//
//  Discussion:
//
//    APSER is used only for cases where
//
//      A <= min ( EPS, EPS * B ),
//      B * X <= 1, and
//      X <= 0.5.
//
//  Parameters:
//
//    Input, double *A, *B, *X, the parameters of teh
//    incomplete beta ratio.
//
//    Input, double *EPS, a tolerance.
//
//    Output, double APSER, the computed value of the
//    incomplete beta ratio.
//
{
  static double g = 0.577215664901533e0;
  static double apser,aj,bx,c,j,s,t,tol;

    bx = *b**x;
    t = *x-bx;
    if(*b**eps > 2.e-2) goto S10;
    c = log(*x)+psi(b)+g+t;
    goto S20;
S10:
    c = log(bx)+g+t;
S20:
    tol = 5.0e0**eps*fabs(c);
    j = 1.0e0;
    s = 0.0e0;
S30:
    j = j + 1.0e0;
    t = t * (*x-bx/j);
    aj = t/j;
    s = s + aj;
    if(fabs(aj) > tol) goto S30;
    apser = -(*a*(c+s));
    return apser;
}
