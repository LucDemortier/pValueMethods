# include <iomanip>
# include <cmath>
using namespace std;
# include "cdflib.hpp"

//****************************************************************************80

void cumf ( double *f, double *dfn, double *dfd, double *cum, double *ccum )

//****************************************************************************80
//
//  Purpose:
//
//    CUMF evaluates the cumulative F distribution.
//
//  Discussion:
//
//    CUMF computes the integral from 0 to F of the F density with DFN
//    numerator and DFD denominator degrees of freedom.
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.5.28.
//
//  Parameters:
//
//    Input, double *F, the upper limit of integration.
//
//    Input, double *DFN, *DFD, the number of degrees of
//    freedom for the numerator and denominator.
//
//    Output, double *CUM, *CCUM, the value of the F CDF and
//    the complementary F CDF.
//
{
# define half 0.5e0
# define done 1.0e0

  static double dsum,prod,xx,yy;
  static int ierr;
  static double T1,T2;

  if(!(*f <= 0.0e0)) goto S10;
  *cum = 0.0e0;
  *ccum = 1.0e0;
  return;
S10:
  prod = *dfn**f;
//
//     XX is such that the incomplete beta with parameters
//     DFD/2 and DFN/2 evaluated at XX is 1 - CUM or CCUM
//     YY is 1 - XX
//     Calculate the smaller of XX and YY accurately
//
  dsum = *dfd+prod;
  xx = *dfd/dsum;

  if ( xx > half )
  {
    yy = prod/dsum;
    xx = done-yy;
  }
  else
  {
    yy = done-xx;
  }

  T1 = *dfd*half;
  T2 = *dfn*half;
  beta_inc ( &T1, &T2, &xx, &yy, ccum, cum, &ierr );
  return;
# undef half
# undef done
}
