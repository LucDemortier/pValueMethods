# include <iomanip>
# include <cmath>
using namespace std;
# include "cdflib.hpp"

//****************************************************************************80

void cumbet ( double *x, double *y, double *a, double *b, double *cum,
  double *ccum )

//****************************************************************************80
//
//  Purpose:
//
//    CUMBET evaluates the cumulative incomplete beta distribution.
//
//  Discussion:
//
//    This routine calculates the CDF to X of the incomplete beta distribution
//    with parameters A and B.  This is the integral from 0 to x
//    of (1/B(a,b))*f(t)) where f(t) = t**(a-1) * (1-t)**(b-1)
//
//  Modified:
//
//    14 March 2006
//
//  Reference:
//
//    A R Didonato and Alfred Morris,
//    Algorithm 708:
//    Significant Digit Computation of the Incomplete Beta Function Ratios.
//    ACM Transactions on Mathematical Software,
//    Volume 18, Number 3, September 1992, pages 360-373.
//
//  Parameters:
//
//    Input, double *X, the upper limit of integration.
//
//    Input, double *Y, the value of 1-X.
//
//    Input, double *A, *B, the parameters of the distribution.
//
//    Output, double *CUM, *CCUM, the values of the cumulative
//    density function and complementary cumulative density function.
//
{
  static int ierr;

  if ( *x <= 0.0 )
  {
    *cum = 0.0;
    *ccum = 1.0;
  }
  else if ( *y <= 0.0 )
  {
    *cum = 1.0;
    *ccum = 0.0;
  }
  else
  {
    beta_inc ( a, b, x, y, cum, ccum, &ierr );
  }
  return;
}
