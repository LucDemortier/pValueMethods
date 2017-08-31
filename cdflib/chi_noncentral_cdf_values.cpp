# include <cmath>

//****************************************************************************80

void chi_noncentral_cdf_values ( int *n_data, double *x, double *lambda,
  int *df, double *cdf )

//****************************************************************************80
//
//  Purpose:
//
//    CHI_NONCENTRAL_CDF_VALUES returns values of the noncentral chi CDF.
//
//  Discussion:
//
//    The CDF of the noncentral chi square distribution can be evaluated
//    within Mathematica by commands such as:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      CDF [ NoncentralChiSquareDistribution [ DF, LAMBDA ], X ]
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
//    Output, double *X, the argument of the function.
//
//    Output, double *LAMBDA, the noncentrality parameter.
//
//    Output, int *DF, the number of degrees of freedom.
//
//    Output, double *CDF, the noncentral chi CDF.
//
{
# define N_MAX 27

  double cdf_vec[N_MAX] = {
    0.839944E+00, 0.695906E+00, 0.535088E+00,
    0.764784E+00, 0.620644E+00, 0.469167E+00,
    0.307088E+00, 0.220382E+00, 0.150025E+00,
    0.307116E-02, 0.176398E-02, 0.981679E-03,
    0.165175E-01, 0.202342E-03, 0.498448E-06,
    0.151325E-01, 0.209041E-02, 0.246502E-03,
    0.263684E-01, 0.185798E-01, 0.130574E-01,
    0.583804E-01, 0.424978E-01, 0.308214E-01,
    0.105788E+00, 0.794084E-01, 0.593201E-01 };
  int df_vec[N_MAX] = {
      1,   2,   3,
      1,   2,   3,
      1,   2,   3,
      1,   2,   3,
     60,  80, 100,
      1,   2,   3,
     10,  10,  10,
     10,  10,  10,
     10,  10,  10 };
  double lambda_vec[N_MAX] = {
     0.5E+00,  0.5E+00,  0.5E+00,
     1.0E+00,  1.0E+00,  1.0E+00,
     5.0E+00,  5.0E+00,  5.0E+00,
    20.0E+00, 20.0E+00, 20.0E+00,
    30.0E+00, 30.0E+00, 30.0E+00,
     5.0E+00,  5.0E+00,  5.0E+00,
     2.0E+00,  3.0E+00,  4.0E+00,
     2.0E+00,  3.0E+00,  4.0E+00,
     2.0E+00,  3.0E+00,  4.0E+00 };
  double x_vec[N_MAX] = {
     3.000E+00,  3.000E+00,  3.000E+00,
     3.000E+00,  3.000E+00,  3.000E+00,
     3.000E+00,  3.000E+00,  3.000E+00,
     3.000E+00,  3.000E+00,  3.000E+00,
    60.000E+00, 60.000E+00, 60.000E+00,
     0.050E+00,  0.050E+00,  0.050E+00,
     4.000E+00,  4.000E+00,  4.000E+00,
     5.000E+00,  5.000E+00,  5.000E+00,
     6.000E+00,  6.000E+00,  6.000E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *x = 0.0E+00;
    *lambda = 0.0E+00;
    *df = 0;
    *cdf = 0.0E+00;
  }
  else
  {
    *x = x_vec[*n_data-1];
    *lambda = lambda_vec[*n_data-1];
    *df = df_vec[*n_data-1];
    *cdf = cdf_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
