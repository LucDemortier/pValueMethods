//****************************************************************************80

void gamma_inc_values ( int *n_data, double *a, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_INC_VALUES returns some values of the incomplete Gamma function.
//
//  Discussion:
//
//    The (normalized) incomplete Gamma function P(A,X) is defined as:
//
//      PN(A,X) = 1/GAMMA(A) * Integral ( 0 <= T <= X ) T**(A-1) * exp(-T) dT.
//
//    With this definition, for all A and X,
//
//      0 <= PN(A,X) <= 1
//
//    and
//
//      PN(A,INFINITY) = 1.0
//
//    Mathematica can compute this value as
//
//      1 - GammaRegularized[A,X]
//
//  Modified:
//
//    31 May 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double *A, the parameter of the function.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 20

  double a_vec[N_MAX] = {
    0.1E+00,  0.1E+00,  0.1E+00,  0.5E+00,
    0.5E+00,  0.5E+00,  1.0E+00,  1.0E+00,
    1.0E+00,  1.1E+00,  1.1E+00,  1.1E+00,
    2.0E+00,  2.0E+00,  2.0E+00,  6.0E+00,
    6.0E+00, 11.0E+00, 26.0E+00, 41.0E+00 };
  double fx_vec[N_MAX] = {
    0.7420263E+00, 0.9119753E+00, 0.9898955E+00, 0.2931279E+00,
    0.7656418E+00, 0.9921661E+00, 0.0951626E+00, 0.6321206E+00,
    0.9932621E+00, 0.0757471E+00, 0.6076457E+00, 0.9933425E+00,
    0.0091054E+00, 0.4130643E+00, 0.9931450E+00, 0.0387318E+00,
    0.9825937E+00, 0.9404267E+00, 0.4863866E+00, 0.7359709E+00 };
  double x_vec[N_MAX] = {
    3.1622777E-02, 3.1622777E-01, 1.5811388E+00, 7.0710678E-02,
    7.0710678E-01, 3.5355339E+00, 0.1000000E+00, 1.0000000E+00,
    5.0000000E+00, 1.0488088E-01, 1.0488088E+00, 5.2440442E+00,
    1.4142136E-01, 1.4142136E+00, 7.0710678E+00, 2.4494897E+00,
    1.2247449E+01, 1.6583124E+01, 2.5495098E+01, 4.4821870E+01 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *a = 0.0E+00;
    *x = 0.0E+00;
    *fx = 0.0E+00;
  }
  else
  {
    *a = a_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }
  return;
# undef N_MAX
}
