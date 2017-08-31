# include <iomanip>
# include <cmath>
using namespace std;
# include "cdflib.hpp"

//****************************************************************************80

void cumnbn ( double *s, double *xn, double *pr, double *ompr,
  double *cum, double *ccum )

//****************************************************************************80
//
//  Purpose:
//
//    CUMNBN evaluates the cumulative negative binomial distribution.
//
//  Discussion:
//
//    This routine returns the probability that there will be F or
//    fewer failures before there are S successes, with each binomial
//    trial having a probability of success PR.
//
//    Prob(# failures = F | S successes, PR)  =
//                        ( S + F - 1 )
//                        (            ) * PR^S * (1-PR)^F
//                        (      F     )
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.5.26.
//
//  Parameters:
//
//    Input, double *F, the number of failures.
//
//    Input, double *S, the number of successes.
//
//    Input, double *PR, *OMPR, the probability of success on
//    each binomial trial, and the value of (1-PR).
//
//    Output, double *CUM, *CCUM, the negative binomial CDF,
//    and the complementary CDF.
//
{
  static double T1;

  T1 = *s+1.e0;
  cumbet(pr,ompr,xn,&T1,cum,ccum);
  return;
}
