//****************************************************************************80

void gamma_values ( int *n_data, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_VALUES returns some values of the Gamma function.
//
//  Definition:
//
//    GAMMA(Z) = Integral ( 0 <= T < Infinity) T**(Z-1) EXP(-T) dT
//
//  Recursion:
//
//    GAMMA(X+1) = X*GAMMA(X)
//
//  Restrictions:
//
//    0 < X ( a software restriction).
//
//  Special values:
//
//    GAMMA(0.5) = sqrt(PI)
//
//    For N a positive integer, GAMMA(N+1) = N!, the standard factorial.
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
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 18

  double fx_vec[N_MAX] = {
    4.590845E+00,     2.218160E+00,     1.489192E+00,     1.164230E+00,
    1.0000000000E+00, 0.9513507699E+00, 0.9181687424E+00, 0.8974706963E+00,
    0.8872638175E+00, 0.8862269255E+00, 0.8935153493E+00, 0.9086387329E+00,
    0.9313837710E+00, 0.9617658319E+00, 1.0000000000E+00, 3.6288000E+05,
    1.2164510E+17,    8.8417620E+30 };
  double x_vec[N_MAX] = {
    0.2E+00,  0.4E+00,  0.6E+00,  0.8E+00,
    1.0E+00,  1.1E+00,  1.2E+00,  1.3E+00,
    1.4E+00,  1.5E+00,  1.6E+00,  1.7E+00,
    1.8E+00,  1.9E+00,  2.0E+00, 10.0E+00,
   20.0E+00, 30.0E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *x = 0.0E+00;
    *fx = 0.0E+00;
  }
  else
  {
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }
  return;
# undef N_MAX
}
