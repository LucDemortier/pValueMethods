# include <iomanip>
# include <cmath>
using namespace std;
# include "cdflib.hpp"

//****************************************************************************80

double exparg ( int *l )

//****************************************************************************80
//
//  Purpose:
//
//    EXPARG returns the largest or smallest legal argument for EXP.
//
//  Discussion:
//
//    Only an approximate limit for the argument of EXP is desired.
//
//  Modified:
//
//    09 December 1999
//
//  Parameters:
//
//    Input, int *L, indicates which limit is desired.
//    If L = 0, then the largest positive argument for EXP is desired.
//    Otherwise, the largest negative argument for EXP for which the
//    result is nonzero is desired.
//
//    Output, double EXPARG, the desired value.
//
{
  static int K1 = 4;
  static int K2 = 9;
  static int K3 = 10;
  static double exparg,lnb;
  static int b,m;

    b = ipmpar(&K1);
    if(b != 2) goto S10;
    lnb = .69314718055995e0;
    goto S40;
S10:
    if(b != 8) goto S20;
    lnb = 2.0794415416798e0;
    goto S40;
S20:
    if(b != 16) goto S30;
    lnb = 2.7725887222398e0;
    goto S40;
S30:
    lnb = log((double)b);
S40:
    if(*l == 0) goto S50;
    m = ipmpar(&K2)-1;
    exparg = 0.99999e0*((double)m*lnb);
    return exparg;
S50:
    m = ipmpar(&K3);
    exparg = 0.99999e0*((double)m*lnb);
    return exparg;
}
