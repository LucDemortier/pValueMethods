//****************************************************************************80

void psi_values ( int *n_data, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    PSI_VALUES returns some values of the Psi or Digamma function.
//
//  Discussion:
//
//    PSI(X) = d LN ( Gamma ( X ) ) / d X = Gamma'(X) / Gamma(X)
//
//    PSI(1) = - Euler's constant.
//
//    PSI(X+1) = PSI(X) + 1 / X.
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
# define N_MAX 11

  double fx_vec[N_MAX] = {
    -0.5772156649E+00, -0.4237549404E+00, -0.2890398966E+00,
    -0.1691908889E+00, -0.0613845446E+00, -0.0364899740E+00,
     0.1260474528E+00,  0.2085478749E+00,  0.2849914333E+00,
     0.3561841612E+00,  0.4227843351E+00 };
  double x_vec[N_MAX] = {
    1.0E+00,  1.1E+00,  1.2E+00,
    1.3E+00,  1.4E+00,  1.5E+00,
    1.6E+00,  1.7E+00,  1.8E+00,
    1.9E+00,  2.0E+00 };

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
