# include <cmath>

//****************************************************************************80

void f_noncentral_cdf_values ( int *n_data, int *a, int *b, double *lambda,
  double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    F_NONCENTRAL_CDF_VALUES returns some values of the F CDF test function.
//
//  Discussion:
//
//    The value of NONCENTRAL_F_CDF ( DFN, DFD, LAMDA, X ) can be evaluated
//    in Mathematica by commands like:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      CDF[NoncentralFRatioDistribution[ DFN, DFD, LAMBDA ], X ]
//
//  Modified:
//
//    12 June 2004
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
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Wolfram Media / Cambridge University Press, 1999.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int *A, int *B, double *LAMBDA, the
//    parameters of the function.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 22

  int a_vec[N_MAX] = {
     1,  1,  1,  1,
     1,  1,  1,  1,
     1,  1,  2,  2,
     3,  3,  4,  4,
     5,  5,  6,  6,
     8, 16 };
  int b_vec[N_MAX] = {
     1,  5,  5,  5,
     5,  5,  5,  5,
     5,  5,  5, 10,
     5,  5,  5,  5,
     1,  5,  6, 12,
    16,  8 };
  double fx_vec[N_MAX] = {
    0.500000E+00, 0.636783E+00, 0.584092E+00, 0.323443E+00,
    0.450119E+00, 0.607888E+00, 0.705928E+00, 0.772178E+00,
    0.819105E+00, 0.317035E+00, 0.432722E+00, 0.450270E+00,
    0.426188E+00, 0.337744E+00, 0.422911E+00, 0.692767E+00,
    0.363217E+00, 0.421005E+00, 0.426667E+00, 0.446402E+00,
    0.844589E+00, 0.816368E+00 };
  double lambda_vec[N_MAX] = {
    0.00E+00,  0.000E+00, 0.25E+00,  1.00E+00,
    1.00E+00,  1.00E+00,  1.00E+00,  1.00E+00,
    1.00E+00,  2.00E+00,  1.00E+00,  1.00E+00,
    1.00E+00,  2.00E+00,  1.00E+00,  1.00E+00,
    0.00E+00,  1.00E+00,  1.00E+00,  1.00E+00,
    1.00E+00,  1.00E+00 };
  double x_vec[N_MAX] = {
    1.00E+00,  1.00E+00, 1.00E+00,  0.50E+00,
    1.00E+00,  2.00E+00, 3.00E+00,  4.00E+00,
    5.00E+00,  1.00E+00, 1.00E+00,  1.00E+00,
    1.00E+00,  1.00E+00, 1.00E+00,  2.00E+00,
    1.00E+00,  1.00E+00, 1.00E+00,  1.00E+00,
    2.00E+00,  2.00E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *a = 0;
    *b = 0;
    *lambda = 0.0E+00;
    *x = 0.0E+00;
    *fx = 0.0E+00;
  }
  else
  {
    *a = a_vec[*n_data-1];
    *b = b_vec[*n_data-1];
    *lambda = lambda_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
