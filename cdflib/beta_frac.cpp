# include <iomanip>
# include <cmath>
using namespace std;
# include "cdflib.hpp"

//****************************************************************************80

double beta_frac ( double *a, double *b, double *x, double *y, double *lambda,
  double *eps )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_FRAC evaluates a continued fraction expansion for IX(A,B).
//
//  Parameters:
//
//    Input, double *A, *B, the parameters of the function.
//    A and B should be nonnegative.  It is assumed that both A and
//    B are greater than 1.
//
//    Input, double *X, *Y.  X is the argument of the
//    function, and should satisy 0 <= X <= 1.  Y should equal 1 - X.
//
//    Input, double *LAMBDA, the value of ( A + B ) * Y - B.
//
//    Input, double *EPS, a tolerance.
//
//    Output, double BETA_FRAC, the value of the continued
//    fraction approximation for IX(A,B).
//
{
  static double bfrac,alpha,an,anp1,beta,bn,bnp1,c,c0,c1,e,n,p,r,r0,s,t,w,yp1;

  bfrac = beta_rcomp ( a, b, x, y );

  if ( bfrac == 0.0e0 )
  {
    return bfrac;
  }

  c = 1.0e0+*lambda;
  c0 = *b/ *a;
  c1 = 1.0e0+1.0e0/ *a;
  yp1 = *y+1.0e0;
  n = 0.0e0;
  p = 1.0e0;
  s = *a+1.0e0;
  an = 0.0e0;
  bn = anp1 = 1.0e0;
  bnp1 = c/c1;
  r = c1/c;
//
//  CONTINUED FRACTION CALCULATION
//
S10:
  n = n + 1.0e0;
  t = n/ *a;
  w = n*(*b-n)**x;
  e = *a/s;
  alpha = p*(p+c0)*e*e*(w**x);
  e = (1.0e0+t)/(c1+t+t);
  beta = n+w/s+e*(c+n*yp1);
  p = 1.0e0+t;
  s += 2.0e0;
//
//  UPDATE AN, BN, ANP1, AND BNP1
//
  t = alpha*an+beta*anp1;
  an = anp1;
  anp1 = t;
  t = alpha*bn+beta*bnp1;
  bn = bnp1;
  bnp1 = t;
  r0 = r;
  r = anp1/bnp1;

  if ( fabs(r-r0) <= (*eps) * r )
  {
    goto S20;
  }
//
//  RESCALE AN, BN, ANP1, AND BNP1
//
  an /= bnp1;
  bn /= bnp1;
  anp1 = r;
  bnp1 = 1.0e0;
  goto S10;
//
//  TERMINATION
//
S20:
  bfrac = bfrac * r;
  return bfrac;
}
