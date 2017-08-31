# include <cmath>

//****************************************************************************80

double esum ( int *mu, double *x )

//****************************************************************************80
//
//  Purpose:
//
//    ESUM evaluates exp ( MU + X ).
//
//  Parameters:
//
//    Input, int *MU, part of the argument.
//
//    Input, double *X, part of the argument.
//
//    Output, double ESUM, the value of exp ( MU + X ).
//
{
  static double esum,w;

    if(*x > 0.0e0) goto S10;
    if(*mu < 0) goto S20;
    w = (double)*mu+*x;
    if(w > 0.0e0) goto S20;
    esum = exp(w);
    return esum;
S10:
    if(*mu > 0) goto S20;
    w = (double)*mu+*x;
    if(w < 0.0e0) goto S20;
    esum = exp(w);
    return esum;
S20:
    w = *mu;
    esum = exp(w)*exp(*x);
    return esum;
}
