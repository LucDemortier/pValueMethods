# include <iomanip>
# include <cmath>
using namespace std;
# include "cdflib.hpp"

//****************************************************************************80

void cdfpoi ( int *which, double *p, double *q, double *s, double *xlam,
  int *status, double *bound )

//****************************************************************************80
//
//  Purpose:
//
//    CDFPOI evaluates the CDF of the Poisson distribution.
//
//  Discussion:
//
//    This routine calculates any one parameter of the Poisson distribution
//    given the others.
//
//    The value P of the cumulative distribution function is calculated
//    directly.
//
//    Computation of other parameters involve a seach for a value that
//    produces the desired value of P.  The search relies on the
//    monotonicity of P with respect to the other parameters.
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.4.21.
//
//  Parameters:
//
//    Input, int *WHICH, indicates which argument is to be calculated
//    from the others.
//    1: Calculate P and Q from S and XLAM;
//    2: Calculate A from P, Q and XLAM;
//    3: Calculate XLAM from P, Q and S.
//
//    Input/output, double *P, the cumulation from 0 to S of the
//    Poisson density.  Whether this is an input or output value, it will
//    lie in the range [0,1].
//
//    Input/output, double *Q, equal to 1-P.  If Q is an input
//    value, it should lie in the range [0,1].  If Q is an output value,
//    it will lie in the range [0,1].
//
//    Input/output, double *S, the upper limit of cumulation of
//    the Poisson CDF.  If this is an input value, it should lie in
//    the range: [0, +infinity).  If it is an output value, it will be
//    searched for in the range: [0,1.0D+300].
//
//    Input/output, double *XLAM, the mean of the Poisson
//    distribution.  If this is an input value, it should lie in the range
//    [0, +infinity).  If it is an output value, it will be searched for
//    in the range: [0,1E300].
//
//    Output, int *STATUS, reports the status of the computation.
//     0, if the calculation completed correctly;
//    -I, if the input parameter number I is out of range;
//    +1, if the answer appears to be lower than lowest search bound;
//    +2, if the answer appears to be higher than greatest search bound;
//    +3, if P + Q /= 1.
//
//    Output, double *BOUND, is only defined if STATUS is nonzero.
//    If STATUS is negative, then this is the value exceeded by parameter I.
//    if STATUS is 1 or 2, this is the search bound that was exceeded.
//
{
# define tol (1.0e-8)
# define atol (1.0e-50)
# define inf 1.0e300

  static int K1 = 1;
  static double K2 = 0.0e0;
  static double K4 = 0.5e0;
  static double K5 = 5.0e0;
  static double fx,cum,ccum,pq;
  static unsigned long qhi,qleft,qporq;
  static double T3,T6,T7,T8,T9,T10;

  *status = 0;
  *bound = 0.0;
//
//     Check arguments
//
    if(!(*which < 1 || *which > 3)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
S10:
    *bound = 3.0e0;
S20:
    *status = -1;
    return;
S30:
    if(*which == 1) goto S70;
//
//     P
//
    if(!(*p < 0.0e0 || *p > 1.0e0)) goto S60;
    if(!(*p < 0.0e0)) goto S40;
    *bound = 0.0e0;
    goto S50;
S40:
    *bound = 1.0e0;
S50:
    *status = -2;
    return;
S70:
S60:
    if(*which == 1) goto S110;
//
//     Q
//
    if(!(*q <= 0.0e0 || *q > 1.0e0)) goto S100;
    if(!(*q <= 0.0e0)) goto S80;
    *bound = 0.0e0;
    goto S90;
S80:
    *bound = 1.0e0;
S90:
    *status = -3;
    return;
S110:
S100:
    if(*which == 2) goto S130;
//
//     S
//
    if(!(*s < 0.0e0)) goto S120;
    *bound = 0.0e0;
    *status = -4;
    return;
S130:
S120:
    if(*which == 3) goto S150;
//
//     XLAM
//
    if(!(*xlam < 0.0e0)) goto S140;
    *bound = 0.0e0;
    *status = -5;
    return;
S150:
S140:
    if(*which == 1) goto S190;
//
//     P + Q
//
    pq = *p+*q;
    if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0*dpmpar(&K1))) goto S180;
    if(!(pq < 0.0e0)) goto S160;
    *bound = 0.0e0;
    goto S170;
S160:
    *bound = 1.0e0;
S170:
    *status = 3;
    return;
S190:
S180:
    if(!(*which == 1)) qporq = *p <= *q;
//
//     Select the minimum of P or Q
//     Calculate ANSWERS
//
    if(1 == *which) {
//
//     Calculating P
//
        cumpoi(s,xlam,p,q);
        *status = 0;
    }
    else if(2 == *which) {
//
//     Calculating S
//
        *s = 5.0e0;
        T3 = inf;
        T6 = atol;
        T7 = tol;
        dstinv(&K2,&T3,&K4,&K4,&K5,&T6,&T7);
        *status = 0;
        dinvr(status,s,&fx,&qleft,&qhi);
S200:
        if(!(*status == 1)) goto S230;
        cumpoi(s,xlam,&cum,&ccum);
        if(!qporq) goto S210;
        fx = cum-*p;
        goto S220;
S210:
        fx = ccum-*q;
S220:
        dinvr(status,s,&fx,&qleft,&qhi);
        goto S200;
S230:
        if(!(*status == -1)) goto S260;
        if(!qleft) goto S240;
        *status = 1;
        *bound = 0.0e0;
        goto S250;
S240:
        *status = 2;
        *bound = inf;
S260:
S250:
        ;
    }
    else if(3 == *which) {
//
//     Calculating XLAM
//
        *xlam = 5.0e0;
        T8 = inf;
        T9 = atol;
        T10 = tol;
        dstinv(&K2,&T8,&K4,&K4,&K5,&T9,&T10);
        *status = 0;
        dinvr(status,xlam,&fx,&qleft,&qhi);
S270:
        if(!(*status == 1)) goto S300;
        cumpoi(s,xlam,&cum,&ccum);
        if(!qporq) goto S280;
        fx = cum-*p;
        goto S290;
S280:
        fx = ccum-*q;
S290:
        dinvr(status,xlam,&fx,&qleft,&qhi);
        goto S270;
S300:
        if(!(*status == -1)) goto S330;
        if(!qleft) goto S310;
        *status = 1;
        *bound = 0.0e0;
        goto S320;
S310:
        *status = 2;
        *bound = inf;
S320:
        ;
    }
S330:
    return;
# undef tol
# undef atol
# undef inf
}
