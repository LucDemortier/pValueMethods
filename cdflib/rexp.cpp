# include <cmath>

//****************************************************************************80

double rexp ( double *x )

//****************************************************************************80
//
//  Purpose:
//
//    REXP evaluates the function EXP(X) - 1.
//
//  Modified:
//
//    09 December 1999
//
//  Parameters:
//
//    Input, double *X, the argument of the function.
//
//    Output, double REXP, the value of EXP(X)-1.
//
{
  static double p1 = .914041914819518e-09;
  static double p2 = .238082361044469e-01;
  static double q1 = -.499999999085958e+00;
  static double q2 = .107141568980644e+00;
  static double q3 = -.119041179760821e-01;
  static double q4 = .595130811860248e-03;
  static double rexp,w;

    if(fabs(*x) > 0.15e0) goto S10;
    rexp = *x*(((p2**x+p1)**x+1.0e0)/((((q4**x+q3)**x+q2)**x+q1)**x+1.0e0));
    return rexp;
S10:
    w = exp(*x);
    if(*x > 0.0e0) goto S20;
    rexp = w-0.5e0-0.5e0;
    return rexp;
S20:
    rexp = w*(0.5e0+(0.5e0-1.0e0/w));
    return rexp;
}
