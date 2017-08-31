# include <iomanip>
# include <cmath>
using namespace std;
# include "cdflib.hpp"

//****************************************************************************80

double beta ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    BETA evaluates the beta function.
//
//  Modified:
//
//    03 December 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the arguments of the beta function.
//
//    Output, double BETA, the value of the beta function.
//
{
  return ( exp ( beta_log ( &a, &b ) ) );
}
