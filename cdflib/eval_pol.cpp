# include <cmath>

//****************************************************************************80

double eval_pol ( double a[], int *n, double *x )

//****************************************************************************80
//
//  Purpose:
//
//    EVAL_POL evaluates a polynomial at X.
//
//  Discussion:
//
//    EVAL_POL = A(0) + A(1)*X + ... + A(N)*X**N
//
//  Modified:
//
//    15 December 1999
//
//  Parameters:
//
//    Input, double precision A(0:N), coefficients of the polynomial.
//
//    Input, int *N, length of A.
//
//    Input, double *X, the point at which the polynomial
//    is to be evaluated.
//
//    Output, double EVAL_POL, the value of the polynomial at X.
//
{
  static double devlpl,term;
  static int i;

  term = a[*n-1];
  for ( i = *n-1-1; i >= 0; i-- )
  {
    term = a[i]+term**x;
  }

  devlpl = term;
  return devlpl;
}
