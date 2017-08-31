# include <iomanip>
# include <cmath>
using namespace std;
# include "cdflib.hpp"

//****************************************************************************80

void cdfchn ( int *which, double *p, double *q, double *x, double *df,
  double *pnonc, int *status, double *bound )

//****************************************************************************80
//
//  Purpose:
//
//    CDFCHN evaluates the CDF of the Noncentral Chi-Square.
//
//  Discussion:
//
//    This routine calculates any one parameter of the noncentral chi-square
//    distribution given values for the others.
//
//    The value P of the cumulative distribution function is calculated
//    directly.
//
//    Computation of the other parameters involves a seach for a value that
//    produces the desired value of P.  The search relies on the
//    monotonicity of P with respect to the other parameters.
//
//    The computation time required for this routine is proportional
//    to the noncentrality parameter (PNONC).  Very large values of
//    this parameter can consume immense computer resources.  This is
//    why the search range is bounded by 10,000.
//
//    The CDF of the noncentral chi square distribution can be evaluated
//    within Mathematica by commands such as:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      CDF[ NoncentralChiSquareDistribution [ DF, LAMBDA ], X ]
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.5.25.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Wolfram Media / Cambridge University Press, 1999.
//
//  Parameters:
//
//    Input, int *WHICH, indicates which argument is to be calculated
//    from the others.
//    1: Calculate P and Q from X, DF and PNONC;
//    2: Calculate X from P, DF and PNONC;
//    3: Calculate DF from P, X and PNONC;
//    4: Calculate PNONC from P, X and DF.
//
//    Input/output, double *P, the integral from 0 to X of
//    the noncentral chi-square distribution.  If this is an input
//    value, it should lie in the range: [0, 1.0-1.0D-16).
//
//    Input/output, double *Q, is generally not used by this
//    subroutine and is only included for similarity with other routines.
//    However, if P is to be computed, then a value will also be computed
//    for Q.
//
//    Input, double *X, the upper limit of integration of the
//    noncentral chi-square distribution.  If this is an input value, it
//    should lie in the range: [0, +infinity).  If it is an output value,
//    it will be sought in the range: [0,1.0D+300].
//
//    Input/output, double *DF, the number of degrees of freedom
//    of the noncentral chi-square distribution.  If this is an input value,
//    it should lie in the range: (0, +infinity).  If it is an output value,
//    it will be searched for in the range: [ 1.0D-300, 1.0D+300].
//
//    Input/output, double *PNONC, the noncentrality parameter of
//    the noncentral chi-square distribution.  If this is an input value, it
//    should lie in the range: [0, +infinity).  If it is an output value,
//    it will be searched for in the range: [0,1.0D+4]
//
//    Output, int *STATUS, reports on the calculation.
//    0, if calculation completed correctly;
//    -I, if input parameter number I is out of range;
//    1, if the answer appears to be lower than the lowest search bound;
//    2, if the answer appears to be higher than the greatest search bound.
//
//    Output, double *BOUND, is only defined if STATUS is nonzero.
//    If STATUS is negative, then this is the value exceeded by parameter I.
//    if STATUS is 1 or 2, this is the search bound that was exceeded.
//
{
# define tent4 1.0e4
# define tol (1.0e-8)
# define atol (1.0e-50)
# define zero (1.0e-300)
# define one (1.0e0-1.0e-16)
# define inf 1.0e300

  static double K1 = 0.0e0;
  static double K3 = 0.5e0;
  static double K4 = 5.0e0;
  static double fx,cum,ccum;
  static unsigned long qhi,qleft;
  static double T2,T5,T6,T7,T8,T9,T10,T11,T12,T13;

  *status = 0;
  *bound = 0.0;
//
//     Check arguments
//
    if(!(*which < 1 || *which > 4)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
S10:
    *bound = 4.0e0;
S20:
    *status = -1;
    return;
S30:
    if(*which == 1) goto S70;
//
//     P
//
    if(!(*p < 0.0e0 || *p > one)) goto S60;
    if(!(*p < 0.0e0)) goto S40;
    *bound = 0.0e0;
    goto S50;
S40:
    *bound = one;
S50:
    *status = -2;
    return;
S70:
S60:
    if(*which == 2) goto S90;
//
//     X
//
    if(!(*x < 0.0e0)) goto S80;
    *bound = 0.0e0;
    *status = -4;
    return;
S90:
S80:
    if(*which == 3) goto S110;
//
//     DF
//
    if(!(*df <= 0.0e0)) goto S100;
    *bound = 0.0e0;
    *status = -5;
    return;
S110:
S100:
    if(*which == 4) goto S130;
//
//     PNONC
//
    if(!(*pnonc < 0.0e0)) goto S120;
    *bound = 0.0e0;
    *status = -6;
    return;
S130:
S120:
//
//     Calculate ANSWERS
//
    if(1 == *which) {
//
//     Calculating P and Q
//
        cumchn(x,df,pnonc,p,q);
        *status = 0;
    }
    else if(2 == *which) {
//
//     Calculating X
//
        *x = 5.0e0;
        T2 = inf;
        T5 = atol;
        T6 = tol;
        dstinv(&K1,&T2,&K3,&K3,&K4,&T5,&T6);
        *status = 0;
        dinvr(status,x,&fx,&qleft,&qhi);
S140:
        if(!(*status == 1)) goto S150;
        cumchn(x,df,pnonc,&cum,&ccum);
        fx = cum-*p;
        dinvr(status,x,&fx,&qleft,&qhi);
        goto S140;
S150:
        if(!(*status == -1)) goto S180;
        if(!qleft) goto S160;
        *status = 1;
        *bound = 0.0e0;
        goto S170;
S160:
        *status = 2;
        *bound = inf;
S180:
S170:
        ;
    }
    else if(3 == *which) {
//
//     Calculating DF
//
        *df = 5.0e0;
        T7 = zero;
        T8 = inf;
        T9 = atol;
        T10 = tol;
        dstinv(&T7,&T8,&K3,&K3,&K4,&T9,&T10);
        *status = 0;
        dinvr(status,df,&fx,&qleft,&qhi);
S190:
        if(!(*status == 1)) goto S200;
        cumchn(x,df,pnonc,&cum,&ccum);
        fx = cum-*p;
        dinvr(status,df,&fx,&qleft,&qhi);
        goto S190;
S200:
        if(!(*status == -1)) goto S230;
        if(!qleft) goto S210;
        *status = 1;
        *bound = zero;
        goto S220;
S210:
        *status = 2;
        *bound = inf;
S230:
S220:
        ;
    }
    else if(4 == *which) {
//
//     Calculating PNONC
//
        *pnonc = 5.0e0;
        T11 = tent4;
        T12 = atol;
        T13 = tol;
        dstinv(&K1,&T11,&K3,&K3,&K4,&T12,&T13);
        *status = 0;
        dinvr(status,pnonc,&fx,&qleft,&qhi);
S240:
        if(!(*status == 1)) goto S250;
        cumchn(x,df,pnonc,&cum,&ccum);
        fx = cum-*p;
        dinvr(status,pnonc,&fx,&qleft,&qhi);
        goto S240;
S250:
        if(!(*status == -1)) goto S280;
        if(!qleft) goto S260;
        *status = 1;
        *bound = zero;
        goto S270;
S260:
        *status = 2;
        *bound = tent4;
S270:
        ;
    }
S280:
    return;
# undef tent4
# undef tol
# undef atol
# undef zero
# undef one
# undef inf
}
