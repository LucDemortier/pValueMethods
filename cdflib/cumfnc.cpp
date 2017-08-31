# include <iomanip>
# include <cmath>
using namespace std;
# include "cdflib.hpp"

//****************************************************************************80

void cumfnc ( double *f, double *dfn, double *dfd, double *pnonc,
  double *cum, double *ccum )

//****************************************************************************80
//
//  Purpose:
//
//    CUMFNC evaluates the cumulative noncentral F distribution.
//
//  Discussion:
//
//    This routine computes the noncentral F distribution with DFN and DFD
//    degrees of freedom and noncentrality parameter PNONC.
//
//    The series is calculated backward and forward from J = LAMBDA/2
//    (this is the term with the largest Poisson weight) until
//    the convergence criterion is met.
//
//    The sum continues until a succeeding term is less than EPS
//    times the sum (or the sum is less than 1.0e-20).  EPS is
//    set to 1.0e-4 in a data statement which can be changed.
//
//
//    The original version of this routine allowed the input values
//    of DFN and DFD to be negative (nonsensical) or zero (which
//    caused numerical overflow.)  I have forced both these values
//    to be at least 1.
//
//  Modified:
//
//    15 June 2004
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.5.16, 26.6.17, 26.6.18, 26.6.20.
//
//  Parameters:
//
//    Input, double *F, the upper limit of integration.
//
//    Input, double *DFN, *DFD, the number of degrees of freedom
//    in the numerator and denominator.  Both DFN and DFD must be positive,
//    and normally would be integers.  This routine requires that they
//    be no less than 1.
//
//    Input, double *PNONC, the noncentrality parameter.
//
//    Output, double *CUM, *CCUM, the noncentral F CDF and
//    complementary CDF.
//
{
# define qsmall(x) (int)(sum < 1.0e-20 || (x) < eps*sum)
# define half 0.5e0
# define done 1.0e0

  static double eps = 1.0e-4;
  static double dsum,dummy,prod,xx,yy,adn,aup,b,betdn,betup,centwt,dnterm,sum,
    upterm,xmult,xnonc;
  static int i,icent,ierr;
  static double T1,T2,T3,T4,T5,T6;

    if(!(*f <= 0.0e0)) goto S10;
    *cum = 0.0e0;
    *ccum = 1.0e0;
    return;
S10:
    if(!(*pnonc < 1.0e-10)) goto S20;
//
//  Handle case in which the non-centrality parameter is
//  (essentially) zero.
//
    cumf(f,dfn,dfd,cum,ccum);
    return;
S20:
    xnonc = *pnonc/2.0e0;
//
//  Calculate the central term of the poisson weighting factor.
//
    icent = ( int ) xnonc;
    if(icent == 0) icent = 1;
//
//  Compute central weight term
//
    T1 = (double)(icent+1);
    centwt = exp(-xnonc+(double)icent*log(xnonc)- gamma_log ( &T1 ) );
//
//  Compute central incomplete beta term
//  Assure that minimum of arg to beta and 1 - arg is computed
//  accurately.
//
    prod = *dfn**f;
    dsum = *dfd+prod;
    yy = *dfd/dsum;
    if(yy > half) {
        xx = prod/dsum;
        yy = done-xx;
    }
    else  xx = done-yy;
    T2 = *dfn*half+(double)icent;
    T3 = *dfd*half;
    beta_inc ( &T2, &T3, &xx, &yy, &betdn, &dummy, &ierr );
    adn = *dfn/2.0e0+(double)icent;
    aup = adn;
    b = *dfd/2.0e0;
    betup = betdn;
    sum = centwt*betdn;
//
//  Now sum terms backward from icent until convergence or all done
//
    xmult = centwt;
    i = icent;
    T4 = adn+b;
    T5 = adn+1.0e0;
    dnterm = exp( gamma_log ( &T4 ) - gamma_log ( &T5 )
      - gamma_log ( &b ) + adn * log ( xx ) + b * log(yy));
S30:
    if(qsmall(xmult*betdn) || i <= 0) goto S40;
    xmult *= ((double)i/xnonc);
    i -= 1;
    adn -= 1.0;
    dnterm = (adn+1.0)/((adn+b)*xx)*dnterm;
    betdn += dnterm;
    sum += (xmult*betdn);
    goto S30;
S40:
    i = icent+1;
//
//  Now sum forwards until convergence
//
    xmult = centwt;
    if(aup-1.0+b == 0) upterm = exp(-gamma_log ( &aup )
      - gamma_log ( &b ) + (aup-1.0)*log(xx)+
      b*log(yy));
    else  {
        T6 = aup-1.0+b;
        upterm = exp( gamma_log ( &T6 ) - gamma_log ( &aup )
          - gamma_log ( &b ) + (aup-1.0)*log(xx)+b*
          log(yy));
    }
    goto S60;
S50:
    if(qsmall(xmult*betup)) goto S70;
S60:
    xmult *= (xnonc/(double)i);
    i += 1;
    aup += 1.0;
    upterm = (aup+b-2.0e0)*xx/(aup-1.0)*upterm;
    betup -= upterm;
    sum += (xmult*betup);
    goto S50;
S70:
    *cum = sum;
    *ccum = 0.5e0+(0.5e0-*cum);
    return;
# undef qsmall
# undef half
# undef done
}
