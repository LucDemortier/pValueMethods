# include <cmath>

//****************************************************************************80

void binomial_cdf_values ( int *n_data, int *a, double *b, int *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    BINOMIAL_CDF_VALUES returns some values of the binomial CDF.
//
//  Discussion:
//
//    CDF(X)(A,B) is the probability of at most X successes in A trials,
//    given that the probability of success on a single trial is B.
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
//    Daniel Zwillinger,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition, CRC Press, 1996, pages 651-652.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int *A, double *B, the parameters of the function.
//
//    Output, int *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 17

  int a_vec[N_MAX] = {
     2,  2,  2,  2,
     2,  4,  4,  4,
     4, 10, 10, 10,
    10, 10, 10, 10,
    10 };
  double b_vec[N_MAX] = {
    0.05E+00, 0.05E+00, 0.05E+00, 0.50E+00,
    0.50E+00, 0.25E+00, 0.25E+00, 0.25E+00,
    0.25E+00, 0.05E+00, 0.10E+00, 0.15E+00,
    0.20E+00, 0.25E+00, 0.30E+00, 0.40E+00,
    0.50E+00 };
  double fx_vec[N_MAX] = {
    0.9025E+00, 0.9975E+00, 1.0000E+00, 0.2500E+00,
    0.7500E+00, 0.3164E+00, 0.7383E+00, 0.9492E+00,
    0.9961E+00, 0.9999E+00, 0.9984E+00, 0.9901E+00,
    0.9672E+00, 0.9219E+00, 0.8497E+00, 0.6331E+00,
    0.3770E+00 };
  int x_vec[N_MAX] = {
     0, 1, 2, 0,
     1, 0, 1, 2,
     3, 4, 4, 4,
     4, 4, 4, 4,
     4 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *a = 0;
    *b = 0.0E+00;
    *x = 0;
    *fx = 0.0E+00;
  }
  else
  {
    *a = a_vec[*n_data-1];
    *b = b_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }
  return;
# undef N_MAX
}
