# include <cmath>

//****************************************************************************80

double dexpm1 ( double *x )

//****************************************************************************80
//
//  Purpose:
//
//    DEXPM1 evaluates the function EXP(X) - 1.
//
//  Reference:
//
//    Armido DiDinato and Alfred Morris,
//    Algorithm 708:
//    Significant Digit Computation of the Incomplete Beta Function Ratios,
//    ACM Transactions on Mathematical Software,
//    Volume 18, 1993, pages 360-373.
//
//  Parameters:
//
//    Input, double *X, the value at which exp(X)-1 is desired.
//
//    Output, double DEXPM1, the value of exp(X)-1.
//
{
  static double p1 = .914041914819518e-09;
  static double p2 = .238082361044469e-01;
  static double q1 = -.499999999085958e+00;
  static double q2 = .107141568980644e+00;
  static double q3 = -.119041179760821e-01;
  static double q4 = .595130811860248e-03;
  static double dexpm1;
  double w;

  if ( fabs(*x) <= 0.15e0 )
  {
    dexpm1 =   *x * ( ( (
        p2   * *x
      + p1 ) * *x
      + 1.0e0 )
      /((((
        q4   * *x
      + q3 ) * *x
      + q2 ) * *x
      + q1 ) * *x
      + 1.0e0 ) );
  }
  else if ( *x <= 0.0e0 )
  {
    w = exp(*x);
    dexpm1 = w-0.5e0-0.5e0;
  }
  else
  {
    w = exp(*x);
    dexpm1 = w*(0.5e0+(0.5e0-1.0e0/w));
  }

  return dexpm1;
}
