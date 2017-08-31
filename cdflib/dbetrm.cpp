# include <iomanip>
# include <cmath>
using namespace std;
# include "cdflib.hpp"

//****************************************************************************80

double dbetrm ( double *a, double *b )

//****************************************************************************80
//
//  Purpose:
//
//    DBETRM computes the Sterling remainder for the complete beta function.
//
//  Discussion:
//
//    Log(Beta(A,B)) = Lgamma(A) + Lgamma(B) - Lgamma(A+B)
//    where Lgamma is the log of the (complete) gamma function
//
//    Let ZZ be approximation obtained if each log gamma is approximated
//    by Sterling's formula, i.e.,
//    Sterling(Z) = LOG( SQRT( 2*PI ) ) + ( Z-0.5D+00 ) * LOG( Z ) - Z
//
//    The Sterling remainder is Log(Beta(A,B)) - ZZ.
//
//  Parameters:
//
//    Input, double *A, *B, the parameters of the Beta function.
//
//    Output, double DBETRM, the Sterling remainder.
//
{
  static double dbetrm,T1,T2,T3;
//
//     Try to sum from smallest to largest
//
    T1 = *a+*b;
    dbetrm = -dstrem(&T1);
    T2 = fifdmax1(*a,*b);
    dbetrm += dstrem(&T2);
    T3 = fifdmin1(*a,*b);
    dbetrm += dstrem(&T3);
    return dbetrm;
}
