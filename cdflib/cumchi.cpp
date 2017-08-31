# include <iomanip>
# include <cmath>
using namespace std;
# include "cdflib.hpp"

//****************************************************************************80

void cumchi ( double *x, double *df, double *cum, double *ccum )

//****************************************************************************80
//
//  Purpose:
//
//    CUMCHI evaluates the cumulative chi-square distribution.
//
//  Parameters:
//
//    Input, double *X, the upper limit of integration.
//
//    Input, double *DF, the degrees of freedom of the
//    chi-square distribution.
//
//    Output, double *CUM, the cumulative chi-square distribution.
//
//    Output, double *CCUM, the complement of the cumulative
//    chi-square distribution.
//
{
  static double a;
  static double xx;

  a = *df * 0.5;
  xx = *x * 0.5;
  cumgam ( &xx, &a, cum, ccum );
  return;
}
