# include <iomanip>
# include <cmath>
using namespace std;
# include "cdflib.hpp"

//****************************************************************************80

void beta_grat ( double *a, double *b, double *x, double *y, double *w,
  double *eps,int *ierr )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_GRAT evaluates an asymptotic expansion for IX(A,B).
//
//  Parameters:
//
//    Input, double *A, *B, the parameters of the function.
//    A and B should be nonnegative.  It is assumed that 15 <= A
//    and B <= 1, and that B is less than A.
//
//    Input, double *X, *Y.  X is the argument of the
//    function, and should satisy 0 <= X <= 1.  Y should equal 1 - X.
//
//    Input/output, double *W, a quantity to which the
//    result of the computation is to be added on output.
//
//    Input, double *EPS, a tolerance.
//
//    Output, int *IERR, an error flag, which is 0 if no error
//    was detected.
//
{
  static double bm1,bp2n,cn,coef,dj,j,l,lnx,n2,nu,p,q,r,s,sum,t,t2,u,v,z;
  static int i,n,nm1;
  static double c[30],d[30],T1;

    bm1 = *b-0.5e0-0.5e0;
    nu = *a+0.5e0*bm1;
    if(*y > 0.375e0) goto S10;
    T1 = -*y;
    lnx = alnrel(&T1);
    goto S20;
S10:
    lnx = log(*x);
S20:
    z = -(nu*lnx);
    if(*b*z == 0.0e0) goto S70;
//
//  COMPUTATION OF THE EXPANSION
//  SET R = EXP(-Z)*Z**B/GAMMA(B)
//
    r = *b*(1.0e0+gam1(b))*exp(*b*log(z));
    r *= (exp(*a*lnx)*exp(0.5e0*bm1*lnx));
    u = algdiv(b,a)+*b*log(nu);
    u = r*exp(-u);
    if(u == 0.0e0) goto S70;
    gamma_rat1 ( b, &z, &r, &p, &q, eps );
    v = 0.25e0*pow(1.0e0/nu,2.0);
    t2 = 0.25e0*lnx*lnx;
    l = *w/u;
    j = q/r;
    sum = j;
    t = cn = 1.0e0;
    n2 = 0.0e0;
    for ( n = 1; n <= 30; n++ )
    {
        bp2n = *b+n2;
        j = (bp2n*(bp2n+1.0e0)*j+(z+bp2n+1.0e0)*t)*v;
        n2 = n2 + 2.0e0;
        t *= t2;
        cn /= (n2*(n2+1.0e0));
        c[n-1] = cn;
        s = 0.0e0;
        if(n == 1) goto S40;
        nm1 = n-1;
        coef = *b-(double)n;
        for ( i = 1; i <= nm1; i++ )
        {
            s = s + (coef*c[i-1]*d[n-i-1]);
            coef = coef + *b;
        }
S40:
        d[n-1] = bm1*cn+s/(double)n;
        dj = d[n-1]*j;
        sum = sum + dj;
        if(sum <= 0.0e0) goto S70;
        if(fabs(dj) <= *eps*(sum+l)) goto S60;
    }
S60:
//
//  ADD THE RESULTS TO W
//
    *ierr = 0;
    *w = *w + (u*sum);
    return;
S70:
//
//  THE EXPANSION CANNOT BE COMPUTED
//
    *ierr = 1;
    return;
}
