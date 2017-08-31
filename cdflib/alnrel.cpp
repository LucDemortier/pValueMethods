# include <cmath>

//****************************************************************************80

double alnrel ( double *a )

//****************************************************************************80
//
//  Purpose:
//
//    ALNREL evaluates the function ln ( 1 + A ).
//
//  Modified:
//
//    17 November 2006
//
//  Reference:
//
//    Armido DiDinato, Alfred Morris,
//    Algorithm 708:
//    Significant Digit Computation of the Incomplete Beta Function Ratios,
//    ACM Transactions on Mathematical Software,
//    Volume 18, 1993, pages 360-373.
//
//  Parameters:
//
//    Input, double *A, the argument.
//
//    Output, double ALNREL, the value of ln ( 1 + A ).
//
{
  double alnrel;
  static double p1 = -0.129418923021993e+01;
  static double p2 =  0.405303492862024e+00;
  static double p3 = -0.178874546012214e-01;
  static double q1 = -0.162752256355323e+01;
  static double q2 =  0.747811014037616e+00;
  static double q3 = -0.845104217945565e-01;
  double t;
  double t2;
  double w;
  double x;

  if ( fabs ( *a ) <= 0.375e0 )
  {
    t = *a / ( *a + 2.0e0 );
    t2 = t * t;
    w = (((p3*t2+p2)*t2+p1)*t2+1.0e0)
      / (((q3*t2+q2)*t2+q1)*t2+1.0e0);
    alnrel = 2.0e0 * t * w;
  }
  else
  {
    x = 1.0e0 + *a;
    alnrel = log ( x );
  }
  return alnrel;
}
