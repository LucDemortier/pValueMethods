# include <cmath>

//****************************************************************************80

void beta_inc_values ( int *n_data, double *a, double *b, double *x,
  double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_INC_VALUES returns some values of the incomplete Beta function.
//
//  Discussion:
//
//    The incomplete Beta function may be written
//
//      BETA_INC(A,B,X) = Integral (0 to X) T**(A-1) * (1-T)**(B-1) dT
//                      / Integral (0 to 1) T**(A-1) * (1-T)**(B-1) dT
//
//    Thus,
//
//      BETA_INC(A,B,0.0) = 0.0
//      BETA_INC(A,B,1.0) = 1.0
//
//    Note that in Mathematica, the expressions:
//
//      BETA[A,B]   = Integral (0 to 1) T**(A-1) * (1-T)**(B-1) dT
//      BETA[X,A,B] = Integral (0 to X) T**(A-1) * (1-T)**(B-1) dT
//
//    and thus, to evaluate the incomplete Beta function requires:
//
//      BETA_INC(A,B,X) = BETA[X,A,B] / BETA[A,B]
//
//  Modified:
//
//    09 June 2004
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
//    Karl Pearson,
//    Tables of the Incomplete Beta Function,
//    Cambridge University Press, 1968.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double *A, *B, the parameters of the function.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 30

  double a_vec[N_MAX] = {
     0.5E+00,  0.5E+00,  0.5E+00,  1.0E+00,
     1.0E+00,  1.0E+00,  1.0E+00,  1.0E+00,
     2.0E+00,  2.0E+00,  2.0E+00,  2.0E+00,
     2.0E+00,  2.0E+00,  2.0E+00,  2.0E+00,
     2.0E+00,  5.5E+00, 10.0E+00, 10.0E+00,
    10.0E+00, 10.0E+00, 20.0E+00, 20.0E+00,
    20.0E+00, 20.0E+00, 20.0E+00, 30.0E+00,
    30.0E+00, 40.0E+00 };
  double b_vec[N_MAX] = {
     0.5E+00,  0.5E+00,  0.5E+00,  0.5E+00,
     0.5E+00,  0.5E+00,  0.5E+00,  1.0E+00,
     2.0E+00,  2.0E+00,  2.0E+00,  2.0E+00,
     2.0E+00,  2.0E+00,  2.0E+00,  2.0E+00,
     2.0E+00,  5.0E+00,  0.5E+00,  5.0E+00,
     5.0E+00, 10.0E+00,  5.0E+00, 10.0E+00,
    10.0E+00, 20.0E+00, 20.0E+00, 10.0E+00,
    10.0E+00, 20.0E+00 };
  double fx_vec[N_MAX] = {
    0.0637686E+00, 0.2048328E+00, 1.0000000E+00, 0.0E+00,
    0.0050126E+00, 0.0513167E+00, 0.2928932E+00, 0.5000000E+00,
    0.028E+00,     0.104E+00,     0.216E+00,     0.352E+00,
    0.500E+00,     0.648E+00,     0.784E+00,     0.896E+00,
    0.972E+00,     0.4361909E+00, 0.1516409E+00, 0.0897827E+00,
    1.0000000E+00, 0.5000000E+00, 0.4598773E+00, 0.2146816E+00,
    0.9507365E+00, 0.5000000E+00, 0.8979414E+00, 0.2241297E+00,
    0.7586405E+00, 0.7001783E+00 };
  double x_vec[N_MAX] = {
    0.01E+00, 0.10E+00, 1.00E+00, 0.0E+00,
    0.01E+00, 0.10E+00, 0.50E+00, 0.50E+00,
    0.1E+00,  0.2E+00,  0.3E+00,  0.4E+00,
    0.5E+00,  0.6E+00,  0.7E+00,  0.8E+00,
    0.9E+00,  0.50E+00, 0.90E+00, 0.50E+00,
    1.00E+00, 0.50E+00, 0.80E+00, 0.60E+00,
    0.80E+00, 0.50E+00, 0.60E+00, 0.70E+00,
    0.80E+00, 0.70E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *a = 0.0E+00;
    *b = 0.0E+00;
    *x = 0.0E+00;
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
