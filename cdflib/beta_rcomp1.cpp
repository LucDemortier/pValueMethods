# include <iomanip>
# include <cmath>
using namespace std;
# include "cdflib.hpp"

//****************************************************************************80

double beta_rcomp1 ( int *mu, double *a, double *b, double *x, double *y )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_RCOMP1 evaluates exp(MU) * X**A * Y**B / Beta(A,B).
//
//  Parameters:
//
//    Input, int MU, ?
//
//    Input, double A, B, the parameters of the Beta function.
//    A and B should be nonnegative.
//
//    Input, double X, Y, ?
//
//    Output, double BETA_RCOMP1, the value of
//    exp(MU) * X**A * Y**B / Beta(A,B).
//
{
  static double Const = .398942280401433e0;
  static double brcmp1,a0,apb,b0,c,e,h,lambda,lnx,lny,t,u,v,x0,y0,z;
  static int i,n;
//
//     CONST = 1/SQRT(2*PI)
//
  static double T1,T2,T3,T4;

    a0 = fifdmin1(*a,*b);
    if(a0 >= 8.0e0) goto S130;
    if(*x > 0.375e0) goto S10;
    lnx = log(*x);
    T1 = -*x;
    lny = alnrel(&T1);
    goto S30;
S10:
    if(*y > 0.375e0) goto S20;
    T2 = -*y;
    lnx = alnrel(&T2);
    lny = log(*y);
    goto S30;
S20:
    lnx = log(*x);
    lny = log(*y);
S30:
    z = *a*lnx+*b*lny;
    if(a0 < 1.0e0) goto S40;
    z -= beta_log(a,b);
    brcmp1 = esum(mu,&z);
    return brcmp1;
S40:
//
//   PROCEDURE FOR A .LT. 1 OR B .LT. 1
//
    b0 = fifdmax1(*a,*b);
    if(b0 >= 8.0e0) goto S120;
    if(b0 > 1.0e0) goto S70;
//
//  ALGORITHM FOR B0 .LE. 1
//
    brcmp1 = esum(mu,&z);
    if(brcmp1 == 0.0e0) return brcmp1;
    apb = *a+*b;
    if(apb > 1.0e0) goto S50;
    z = 1.0e0+gam1(&apb);
    goto S60;
S50:
    u = *a+*b-1.e0;
    z = (1.0e0+gam1(&u))/apb;
S60:
    c = (1.0e0+gam1(a))*(1.0e0+gam1(b))/z;
    brcmp1 = brcmp1*(a0*c)/(1.0e0+a0/b0);
    return brcmp1;
S70:
//
//  ALGORITHM FOR 1 .LT. B0 .LT. 8
//
    u = gamma_ln1 ( &a0 );
    n = ( int ) ( b0 - 1.0e0 );
    if(n < 1) goto S90;
    c = 1.0e0;
    for ( i = 1; i <= n; i++ )
    {
        b0 -= 1.0e0;
        c *= (b0/(a0+b0));
    }
    u = log(c)+u;
S90:
    z -= u;
    b0 -= 1.0e0;
    apb = a0+b0;
    if(apb > 1.0e0) goto S100;
    t = 1.0e0+gam1(&apb);
    goto S110;
S100:
    u = a0+b0-1.e0;
    t = (1.0e0+gam1(&u))/apb;
S110:
    brcmp1 = a0*esum(mu,&z)*(1.0e0+gam1(&b0))/t;
    return brcmp1;
S120:
//
//  ALGORITHM FOR B0 .GE. 8
//
    u = gamma_ln1 ( &a0 ) + algdiv ( &a0, &b0 );
    T3 = z-u;
    brcmp1 = a0*esum(mu,&T3);
    return brcmp1;
S130:
//
//    PROCEDURE FOR A .GE. 8 AND B .GE. 8
//
    if(*a > *b) goto S140;
    h = *a/ *b;
    x0 = h/(1.0e0+h);
    y0 = 1.0e0/(1.0e0+h);
    lambda = *a-(*a+*b)**x;
    goto S150;
S140:
    h = *b/ *a;
    x0 = 1.0e0/(1.0e0+h);
    y0 = h/(1.0e0+h);
    lambda = (*a+*b)**y-*b;
S150:
    e = -(lambda/ *a);
    if(fabs(e) > 0.6e0) goto S160;
    u = rlog1(&e);
    goto S170;
S160:
    u = e-log(*x/x0);
S170:
    e = lambda/ *b;
    if(fabs(e) > 0.6e0) goto S180;
    v = rlog1(&e);
    goto S190;
S180:
    v = e-log(*y/y0);
S190:
    T4 = -(*a*u+*b*v);
    z = esum(mu,&T4);
    brcmp1 = Const*sqrt(*b*x0)*z*exp(-bcorr(a,b));
    return brcmp1;
}
