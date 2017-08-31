# include <iomanip>
# include <cmath>
using namespace std;
# include "cdflib.hpp"

//****************************************************************************80

void beta_inc ( double *a, double *b, double *x, double *y, double *w,
  double *w1, int *ierr )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_INC evaluates the incomplete beta function IX(A,B).
//
//  Author:
//
//    Alfred H Morris, Jr,
//    Naval Surface Weapons Center,
//    Dahlgren, Virginia.
//
//  Parameters:
//
//    Input, double *A, *B, the parameters of the function.
//    A and B should be nonnegative.
//
//    Input, double *X, *Y.  X is the argument of the
//    function, and should satisy 0 <= X <= 1.  Y should equal 1 - X.
//
//    Output, double *W, *W1, the values of IX(A,B) and
//    1-IX(A,B).
//
//    Output, int *IERR, the error flag.
//    0, no error was detected.
//    1, A or B is negative;
//    2, A = B = 0;
//    3, X < 0 or 1 < X;
//    4, Y < 0 or 1 < Y;
//    5, X + Y /= 1;
//    6, X = A = 0;
//    7, Y = B = 0.
//
{
  static int K1 = 1;
  static double a0,b0,eps,lambda,t,x0,y0,z;
  static int ierr1,ind,n;
  static double T2,T3,T4,T5;
//
//  EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE SMALLEST
//  NUMBER FOR WHICH 1.0 + EPS .GT. 1.0
//
    eps = dpmpar ( &K1 );
    *w = *w1 = 0.0e0;
    if(*a < 0.0e0 || *b < 0.0e0) goto S270;
    if(*a == 0.0e0 && *b == 0.0e0) goto S280;
    if(*x < 0.0e0 || *x > 1.0e0) goto S290;
    if(*y < 0.0e0 || *y > 1.0e0) goto S300;
    z = *x+*y-0.5e0-0.5e0;
    if(fabs(z) > 3.0e0*eps) goto S310;
    *ierr = 0;
    if(*x == 0.0e0) goto S210;
    if(*y == 0.0e0) goto S230;
    if(*a == 0.0e0) goto S240;
    if(*b == 0.0e0) goto S220;
    eps = fifdmax1(eps,1.e-15);
    if(fifdmax1(*a,*b) < 1.e-3*eps) goto S260;
    ind = 0;
    a0 = *a;
    b0 = *b;
    x0 = *x;
    y0 = *y;
    if(fifdmin1(a0,b0) > 1.0e0) goto S40;
//
//  PROCEDURE FOR A0 .LE. 1 OR B0 .LE. 1
//
    if(*x <= 0.5e0) goto S10;
    ind = 1;
    a0 = *b;
    b0 = *a;
    x0 = *y;
    y0 = *x;
S10:
    if(b0 < fifdmin1(eps,eps*a0)) goto S90;
    if(a0 < fifdmin1(eps,eps*b0) && b0*x0 <= 1.0e0) goto S100;
    if(fifdmax1(a0,b0) > 1.0e0) goto S20;
    if(a0 >= fifdmin1(0.2e0,b0)) goto S110;
    if(pow(x0,a0) <= 0.9e0) goto S110;
    if(x0 >= 0.3e0) goto S120;
    n = 20;
    goto S140;
S20:
    if(b0 <= 1.0e0) goto S110;
    if(x0 >= 0.3e0) goto S120;
    if(x0 >= 0.1e0) goto S30;
    if(pow(x0*b0,a0) <= 0.7e0) goto S110;
S30:
    if(b0 > 15.0e0) goto S150;
    n = 20;
    goto S140;
S40:
//
//  PROCEDURE FOR A0 .GT. 1 AND B0 .GT. 1
//
    if(*a > *b) goto S50;
    lambda = *a-(*a+*b)**x;
    goto S60;
S50:
    lambda = (*a+*b)**y-*b;
S60:
    if(lambda >= 0.0e0) goto S70;
    ind = 1;
    a0 = *b;
    b0 = *a;
    x0 = *y;
    y0 = *x;
    lambda = fabs(lambda);
S70:
    if(b0 < 40.0e0 && b0*x0 <= 0.7e0) goto S110;
    if(b0 < 40.0e0) goto S160;
    if(a0 > b0) goto S80;
    if(a0 <= 100.0e0) goto S130;
    if(lambda > 0.03e0*a0) goto S130;
    goto S200;
S80:
    if(b0 <= 100.0e0) goto S130;
    if(lambda > 0.03e0*b0) goto S130;
    goto S200;
S90:
//
//  EVALUATION OF THE APPROPRIATE ALGORITHM
//
    *w = fpser(&a0,&b0,&x0,&eps);
    *w1 = 0.5e0+(0.5e0-*w);
    goto S250;
S100:
    *w1 = apser(&a0,&b0,&x0,&eps);
    *w = 0.5e0+(0.5e0-*w1);
    goto S250;
S110:
    *w = beta_pser(&a0,&b0,&x0,&eps);
    *w1 = 0.5e0+(0.5e0-*w);
    goto S250;
S120:
    *w1 = beta_pser(&b0,&a0,&y0,&eps);
    *w = 0.5e0+(0.5e0-*w1);
    goto S250;
S130:
    T2 = 15.0e0*eps;
    *w = beta_frac ( &a0,&b0,&x0,&y0,&lambda,&T2 );
    *w1 = 0.5e0+(0.5e0-*w);
    goto S250;
S140:
    *w1 = beta_up ( &b0, &a0, &y0, &x0, &n, &eps );
    b0 = b0 + (double)n;
S150:
    T3 = 15.0e0*eps;
    beta_grat (&b0,&a0,&y0,&x0,w1,&T3,&ierr1);
    *w = 0.5e0+(0.5e0-*w1);
    goto S250;
S160:
    n = ( int ) b0;
    b0 -= (double)n;
    if(b0 != 0.0e0) goto S170;
    n -= 1;
    b0 = 1.0e0;
S170:
    *w = beta_up ( &b0, &a0, &y0, &x0, &n, &eps );
    if(x0 > 0.7e0) goto S180;
    *w = *w + beta_pser(&a0,&b0,&x0,&eps);
    *w1 = 0.5e0+(0.5e0-*w);
    goto S250;
S180:
    if(a0 > 15.0e0) goto S190;
    n = 20;
    *w = *w + beta_up ( &a0, &b0, &x0, &y0, &n, &eps );
    a0 = a0 + (double)n;
S190:
    T4 = 15.0e0*eps;
    beta_grat ( &a0, &b0, &x0, &y0, w, &T4, &ierr1 );
    *w1 = 0.5e0+(0.5e0-*w);
    goto S250;
S200:
    T5 = 100.0e0*eps;
    *w = beta_asym ( &a0, &b0, &lambda, &T5 );
    *w1 = 0.5e0+(0.5e0-*w);
    goto S250;
S210:
//
//  TERMINATION OF THE PROCEDURE
//
    if(*a == 0.0e0) goto S320;
S220:
    *w = 0.0e0;
    *w1 = 1.0e0;
    return;
S230:
    if(*b == 0.0e0) goto S330;
S240:
    *w = 1.0e0;
    *w1 = 0.0e0;
    return;
S250:
    if(ind == 0) return;
    t = *w;
    *w = *w1;
    *w1 = t;
    return;
S260:
//
//  PROCEDURE FOR A AND B .LT. 1.E-3*EPS
//
    *w = *b/(*a+*b);
    *w1 = *a/(*a+*b);
    return;
S270:
//
//  ERROR RETURN
//
    *ierr = 1;
    return;
S280:
    *ierr = 2;
    return;
S290:
    *ierr = 3;
    return;
S300:
    *ierr = 4;
    return;
S310:
    *ierr = 5;
    return;
S320:
    *ierr = 6;
    return;
S330:
    *ierr = 7;
    return;
}
