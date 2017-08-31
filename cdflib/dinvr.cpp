# include <iomanip>
# include <cmath>
using namespace std;
# include "cdflib.hpp"

//****************************************************************************80

void dinvr ( int *status, double *x, double *fx,
  unsigned long *qleft, unsigned long *qhi )

//****************************************************************************80
//
//  Purpose:
//
//    DINVR bounds the zero of the function and invokes DZROR.
//
//  Discussion:
//
//    This routine seeks to find bounds on a root of the function and
//    invokes ZROR to perform the zero finding.  STINVR must have been
//    called before this routine in order to set its parameters.
//
//  Reference:
//
//    J C P Bus and T J Dekker,
//    Two Efficient Algorithms with Guaranteed Convergence for
//      Finding a Zero of a Function,
//    ACM Transactions on Mathematical Software,
//    Volume 1, Number 4, pages 330-345, 1975.
//
//  Parameters:
//
//    Input/output, integer STATUS.  At the beginning of a zero finding
//    problem, STATUS should be set to 0 and INVR invoked.  The value
//    of parameters other than X will be ignored on this call.
//    If INVR needs the function to be evaluated, it will set STATUS to 1
//    and return.  The value of the function should be set in FX and INVR
//    again called without changing any of its other parameters.
//    If INVR finishes without error, it returns with STATUS 0, and X an
//    approximate root of F(X).
//    If INVR cannot bound the function, it returns a negative STATUS and
//    sets QLEFT and QHI.
//
//    Output, double precision X, the value at which F(X) is to be evaluated.
//
//    Input, double precision FX, the value of F(X) calculated by the user
//    on the previous call, when INVR returned with STATUS = 1.
//
//    Output, logical QLEFT, is defined only if QMFINV returns FALSE.  In that
//    case, QLEFT is TRUE if the stepping search terminated unsucessfully
//    at SMALL, and FALSE if the search terminated unsucessfully at BIG.
//
//    Output, logical QHI, is defined only if QMFINV returns FALSE.  In that
//    case, it is TRUE if Y < F(X) at the termination of the search and FALSE
//    if F(X) < Y.
//
{
  E0000(0,status,x,fx,qleft,qhi,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
}
