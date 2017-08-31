# include <cmath>

//****************************************************************************80

void chi_square_cdf_values ( int *n_data, int *a, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    CHI_SQUARE_CDF_VALUES returns some values of the Chi-Square CDF.
//
//  Discussion:
//
//    The value of CHI_CDF ( DF, X ) can be evaluated in Mathematica by
//    commands like:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      CDF[ChiSquareDistribution[DF], X ]
//
//  Modified:
//
//    11 June 2004
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
//    Output, int *A, the parameter of the function.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 21

  int a_vec[N_MAX] = {
     1,  2,  1,  2,
     1,  2,  3,  4,
     1,  2,  3,  4,
     5,  3,  3,  3,
     3,  3, 10, 10,
    10 };
  double fx_vec[N_MAX] = {
    0.0796557E+00, 0.00498752E+00, 0.112463E+00,    0.00995017E+00,
    0.472911E+00,  0.181269E+00,   0.0597575E+00,   0.0175231E+00,
    0.682689E+00,  0.393469E+00,   0.198748E+00,    0.090204E+00,
    0.0374342E+00, 0.427593E+00,   0.608375E+00,    0.738536E+00,
    0.828203E+00,  0.88839E+00,    0.000172116E+00, 0.00365985E+00,
    0.0185759E+00 };
  double x_vec[N_MAX] = {
    0.01E+00, 0.01E+00, 0.02E+00, 0.02E+00,
    0.40E+00, 0.40E+00, 0.40E+00, 0.40E+00,
    1.00E+00, 1.00E+00, 1.00E+00, 1.00E+00,
    1.00E+00, 2.00E+00, 3.00E+00, 4.00E+00,
    5.00E+00, 6.00E+00, 1.00E+00, 2.00E+00,
    3.00E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *a = 0;
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
