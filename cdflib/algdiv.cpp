# include <iomanip>
# include <cmath>
using namespace std;
# include "cdflib.hpp"

//****************************************************************************80

double algdiv ( double *a, double *b )

//****************************************************************************80
//
//  Purpose:
//
//    ALGDIV computes ln ( Gamma ( B ) / Gamma ( A + B ) ) when 8 <= B.
//
//  Discussion:
//
//    In this algorithm, DEL(X) is the function defined by
//
//      ln ( Gamma(X) ) = ( X - 0.5 ) * ln ( X ) - X + 0.5 * ln ( 2 * PI )
//                      + DEL(X).
//
//  Parameters:
//
//    Input, double *A, *B, define the arguments.
//
//    Output, double ALGDIV, the value of ln(Gamma(B)/Gamma(A+B)).
//
{
  static double algdiv;
  static double c;
  static double c0 =  0.833333333333333e-01;
  static double c1 = -0.277777777760991e-02;
  static double c2 =  0.793650666825390e-03;
  static double c3 = -0.595202931351870e-03;
  static double c4 =  0.837308034031215e-03;
  static double c5 = -0.165322962780713e-02;
  static double d;
  static double h;
  static double s11;
  static double s3;
  static double s5;
  static double s7;
  static double s9;
  static double t;
  static double T1;
  static double u;
  static double v;
  static double w;
  static double x;
  static double x2;

  if ( *b <= *a )
  {
    h = *b / *a;
    c = 1.0e0 / ( 1.0e0 + h );
    x = h / ( 1.0e0 + h );
    d = *a + ( *b - 0.5e0 );
  }
  else
  {
    h = *a / *b;
    c = h / ( 1.0e0 + h );
    x = 1.0e0 / ( 1.0e0 + h );
    d = *b + ( *a - 0.5e0 );
  }
//
//  SET SN = (1 - X**N)/(1 - X)
//
  x2 = x * x;
  s3 = 1.0e0 + ( x + x2 );
  s5 = 1.0e0 + ( x + x2 * s3 );
  s7 = 1.0e0 + ( x + x2 * s5 );
  s9 = 1.0e0 + ( x + x2 * s7 );
  s11 = 1.0e0 + ( x + x2 * s9 );
//
//  SET W = DEL(B) - DEL(A + B)
//
  t = pow ( 1.0e0 / *b, 2.0 );

  w = (((( c5 * s11  * t
         + c4 * s9 ) * t
         + c3 * s7 ) * t
         + c2 * s5 ) * t
         + c1 * s3 ) * t
         + c0;

  w *= ( c / *b );
//
//  Combine the results.
//
  T1 = *a / *b;
  u = d * alnrel ( &T1 );
  v = *a * ( log ( *b ) - 1.0e0 );

  if ( v < u )
  {
    algdiv = w - v - u;
  }
  else
  {
    algdiv = w - u - v;
  }
  return algdiv;
}
