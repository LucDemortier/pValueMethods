# include <iomanip>
# include <cmath>
using namespace std;
# include "cdflib.hpp"

//****************************************************************************80

void cumbin ( double *s, double *xn, double *pr, double *ompr,
  double *cum, double *ccum )

//****************************************************************************80
//
//  Purpose:
//
//    CUMBIN evaluates the cumulative binomial distribution.
//
//  Discussion:
//
//    This routine returns the probability of 0 to S successes in XN binomial
//    trials, each of which has a probability of success, PR.
//
//  Modified:
//
//    14 March 2006
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.5.24.
//
//  Parameters:
//
//    Input, double *S, the upper limit of summation.
//
//    Input, double *XN, the number of trials.
//
//    Input, double *PR, the probability of success in one trial.
//
//    Input, double *OMPR, equals ( 1 - PR ).
//
//    Output, double *CUM, the cumulative binomial distribution.
//
//    Output, double *CCUM, the complement of the cumulative
//    binomial distribution.
//
{
  static double T1,T2;

  if ( *s < *xn )
  {
    T1 = *s + 1.0;
    T2 = *xn - *s;
    cumbet ( pr, ompr, &T1, &T2, ccum, cum );
  }
  else
  {
    *cum = 1.0;
    *ccum = 0.0;
  }
  return;
}
