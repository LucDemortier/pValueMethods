# include <iomanip>
# include <cmath>
using namespace std;
# include "cdflib.hpp"

//****************************************************************************80

void cumpoi ( double *s, double *xlam, double *cum, double *ccum )

//****************************************************************************80
//
//  Purpose:
//
//    CUMPOI evaluates the cumulative Poisson distribution.
//
//  Discussion:
//
//    CUMPOI returns the probability of S or fewer events in a Poisson
//    distribution with mean XLAM.
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    Formula 26.4.21.
//
//  Parameters:
//
//    Input, double *S, the upper limit of cumulation of the
//    Poisson density function.
//
//    Input, double *XLAM, the mean of the Poisson distribution.
//
//    Output, double *CUM, *CCUM, the Poisson density CDF and
//    complementary CDF.
//
{
  static double chi,df;

  df = 2.0e0*(*s+1.0e0);
  chi = 2.0e0**xlam;
  cumchi(&chi,&df,ccum,cum);
  return;
}
