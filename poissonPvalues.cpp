#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_integration.h>

double ppp_n_int(double x, void * p);
double ppp_nru_int(double x, void * p);
double ppp_logn_int(double x, void * p);
double fid_p_int(double x, void * p);
double api_pvalue(void * p);

using namespace std;

#include "cdflib/cdflib.hpp"

struct poiParams { double nObs; double poiMean; double poiUnc; double gauPoiRatio; double coeffOfVar; bool excess;};

int main()
{
    const double relError = 1.0e-08;
    const int workSize = 1000;
    gsl_integration_workspace * workPtr = gsl_integration_workspace_alloc(workSize);
    const int work2Size = 1000;
    gsl_integration_cquad_workspace * work2Ptr = gsl_integration_cquad_workspace_alloc(work2Size);
    struct poiParams par;
    int    status, which=2, acc=0;
    double bound, xMean=0.0, xStD=1.0, pAdjustment, n1Obs;

    cout << "Number of events observed: ";
    cin  >> par.nObs;
    cout << "Poisson mean: ";
    cin  >> par.poiMean;
    cout << "Uncertainty on mean: ";
    cin  >> par.poiUnc;
    cout << "P-value adjustment factor: ";
    cin  >> pAdjustment;

    n1Obs = par.nObs + 1;
    par.excess = (par.nObs >= par.poiMean);
    par.gauPoiRatio = 1.0;
    par.coeffOfVar = par.poiUnc/par.poiMean;

    cout << "\nPoisson mean: " << par.poiMean << " +/- " << par.poiUnc << ", observation: " << par.nObs << ", p-value adjustment: " << pAdjustment << endl;
    if (par.excess) {
        cout << "Computing the significance of an *excess*." << endl;
    } else {
        cout << "Computing the significance of a *deficit*." << endl;
    }
    cout << "\nP-Value      Nsigmas" << endl;
    cout << "--------------------" << endl;

// First ignore uncertainty on Poisson mean when computing p-value
    double pVal0, qVal0, nSig0;
    if (par.excess) {
        if(par.nObs > 0) {
            gamma_inc( &par.nObs, &par.poiMean, &pVal0, &qVal0, &acc );
        } else {
            pVal0  = 1;
        }
    } else {
        gamma_inc( &n1Obs, &par.poiMean, &qVal0, &pVal0, &acc );
    }
    pVal0 *= pAdjustment;
    qVal0 = 1.0 - pVal0;
    cdfnor( &which, &qVal0, &pVal0, &nSig0, &xMean, &xStD, &status, &bound );
    cout << setw(11) << left << pVal0 << "  " << setw(7) << left << nSig0 << "  (ignoring uncertainty on Poisson mean)" << endl;

    if (par.poiUnc != 0) {

// Try a truncated Gaussian prior for the Poisson mean
        double pVal1, qVal1, aErr1, rErr1, nSig1;
        if (!par.excess || (par.nObs > 0)) {
            gsl_function Fnor;
            Fnor.function = &ppp_n_int;
            Fnor.params   = &par;
            gsl_integration_qags(&Fnor, 0.0, 1.0, 0.0, relError, workSize, workPtr, &pVal1, &aErr1);
            rErr1 = aErr1/pVal1;
            pVal1 *= pAdjustment;
            qVal1 = 1.0 - pVal1;
            cdfnor( &which, &qVal1, &pVal1, &nSig1, &xMean, &xStD, &status, &bound );
        } else {
            pVal1 = 1.0;
            aErr1 = 0.0;
            rErr1 = 0.0;
            nSig1 = -INFINITY;
        }
        cout << setw(11) << left << pVal1 << "  " << setw(7) << left << nSig1 << "  (prior-pred., Gaussian prior)";
        if (rErr1 > relError) {cout << "; RP=" << rErr1;}
        cout << endl;

// Try a gamma prior for the Poisson mean
        double alpha, beta, betac;
        double pVal2, qVal2, nSig2;
        if (!par.excess || (par.nObs > 0)) {
            alpha = pow(par.poiMean/par.poiUnc, 2.0);
            betac = par.poiMean / (par.poiMean + (par.poiUnc*par.poiUnc));
            beta  = 1.0 - betac;
            if (par.excess) {
                cumbet( &beta, &betac, &par.nObs, &alpha, &pVal2, &qVal2 );
            } else {
                cumbet( &beta, &betac, &n1Obs, &alpha, &qVal2, &pVal2 );
            }
            pVal2 *= pAdjustment;
            qVal2 = 1.0 - pVal2;
            cdfnor( &which, &qVal2, &pVal2, &nSig2, &xMean, &xStD, &status, &bound );
        } else {
            pVal2 = 1.0;
            nSig2 = -INFINITY;
        }
        cout << setw(11) << left << pVal2 << "  " << setw(7) << left << nSig2 << "  (prior-pred., gamma prior)" << endl;

// Try a lognormal prior for the Poisson mean
        double pVal3, qVal3, aErr3, rErr3, nSig3;
        if (!par.excess || par.nObs > 0) {
            gsl_function Flognor;
            Flognor.function = &ppp_logn_int;
            Flognor.params = &par;
            gsl_integration_qags(&Flognor, 0.0, 1.0, 0.0, relError, workSize, workPtr, &pVal3, &aErr3);
            rErr3 = aErr3/pVal3;
            pVal3 *= pAdjustment;
            qVal3 = 1.0 - pVal3;
            cdfnor( &which, &qVal3, &pVal3, &nSig3, &xMean, &xStD, &status, &bound );
        } else {
            pVal3 = 1.0;
            aErr3 = 0.0;
            rErr3 = 0.0;
            nSig3 = -INFINITY;
        }
        cout << setw(11) << left << pVal3 << "  " << setw(7) << left << nSig3 << "  (prior-pred., lognormal prior)";
        if (rErr3 > relError) {cout << "; RP=" << rErr3;}
        cout << endl;

// Try a truncated Gaussian prior with *relative uncertainty* for the Poisson mean
        double pVal4, qVal4, aErr4, rErr4, nSig4;
        unsigned long nEvals;
        if (!par.excess || par.nObs > 0) {
            gsl_function Fnor_relunc;
            Fnor_relunc.function = &ppp_nru_int;
            Fnor_relunc.params = &par;
            gsl_integration_cquad (&Fnor_relunc, 0.0, 1.0, 0.0, relError, work2Ptr, &pVal4, &aErr4, &nEvals);
            rErr4 = aErr4/pVal4;
            pVal4 *= pAdjustment;
            qVal4 = 1.0 - pVal4;
            cdfnor( &which, &qVal4, &pVal4, &nSig4, &xMean, &xStD, &status, &bound );
        } else {
            pVal4 = 1.0;
            aErr4 = 0.0;
            rErr4 = 0.0;
            nSig4 = -INFINITY;
        }
        cout << setw(11) << left << pVal4 << "  " << setw(7) << left << nSig4 << "  (prior-pred., Gaussian prior on rel. unc.)";
        if (rErr4 > relError) {cout << "; RP=" << rErr4;}
        cout << endl;

// Try a fiducial p-value
        double pVal5, qVal5, aErr5, rErr5, nSig5;
//        if (par.nObs > 0) {
            gsl_function Ffid;
            Ffid.function = &fid_p_int;
            Ffid.params = &par;
            gsl_integration_qags(&Ffid, 0.0, 1.0, 0.0, relError, workSize, workPtr, &pVal5, &aErr5);
            rErr5 = aErr5/pVal5;
            if (!par.excess) {pVal5 = 1 - pVal5;}
            pVal5 *= pAdjustment;
            qVal5 = 1.0 - pVal5;
            cdfnor( &which, &qVal5, &pVal5, &nSig5, &xMean, &xStD, &status, &bound );
 //       } else {
 //           if (!par.excess) {
 //               pVal5 = gsl_cdf_ugaussian_Q(par.poiMean/par.poiUnc);
 //           } else {
 //               pVal5 = gsl_cdf_ugaussian_P(par.poiMean/par.poiUnc);
 //           }
 //           pVal5 *= pAdjustment;
 //           pErr5 = 0.0;
 //           nSig5 = gsl_cdf_ugaussian_Qinv(pVal5);
 //       }
        cout << setw(11) << left << pVal5 << "  " << setw(7) << left << nSig5 << "  (fiducial)";
        if (rErr5 > relError) {cout << "; RP=" << rErr5;}
        cout << endl;

// Try a plug-in p-value
        double pVal6, qVal6, nSig6;
        double dnu2 = pow(par.poiUnc, 2);
        double tmp = 0.5 * (par.poiMean - dnu2);
        double nuEst = tmp + sqrt(pow(tmp,2) + par.nObs*dnu2);
        if (par.excess) {
            if(par.nObs > 0) {
                gamma_inc( &par.nObs, &nuEst, &pVal6, &qVal6, &acc );
            } else {
                pVal6 = 1;
            }
        } else {
            gamma_inc( &n1Obs, &nuEst, &qVal6, &pVal6, &acc );
        }
        pVal6 *= pAdjustment;
        qVal6 = 1.0 - pVal6;
        cdfnor( &which, &qVal6, &pVal6, &nSig6, &xMean, &xStD, &status, &bound );
        cout << setw(11) << left << pVal6 << "  " << setw(7) << left << nSig6 << "  (plug-in)" << endl;

// Try an adjusted plug-in p-value
        double pVal7, qVal7, nSig7;
        pVal7 = api_pvalue(&par);
        pVal7 *= pAdjustment;
        qVal7 = 1.0 - pVal7;
        cdfnor( &which, &qVal7, &pVal7, &nSig7, &xMean, &xStD, &status, &bound );
        cout << setw(11) << left << pVal7 << "  " << setw(7) << left << nSig7 << "  (adjusted plug-in)" << endl;

    }

    gsl_integration_workspace_free(workPtr);
    gsl_integration_cquad_workspace_free(work2Ptr);
    return 0;
}

double ppp_n_int(double x, void * p) {
// Integrand of the prior-predictive p-value with truncated normal prior
    struct poiParams * params = (struct poiParams *)p;
    double nObs    = (params->nObs);
    double poiMean = (params->poiMean);
    double poiUnc  = (params->poiUnc);
    bool   excess  = (params->excess);

    int    which=1, status;
    double cval, y, tmp1, tmp2, tmp3, qval, bound, uLim, mean=0.0, sd=1.0;

    cval = max(1.0, nObs);
    y    = cval * (1.0-x)/x;
    uLim = poiMean / poiUnc;
    cdfnor( &which, &tmp3, &qval, &uLim, &mean, &sd, &status, &bound );
    uLim = (poiMean-y)/poiUnc;
    cdfnor( &which, &tmp2, &qval, &uLim, &mean, &sd, &status, &bound );
    if (excess) {
        tmp1 = (nObs-1)*log(y) - y - gsl_sf_lngamma(nObs);
    } else {
        tmp1 = nObs*log(y) - y - gsl_sf_lngamma(nObs+1);
        tmp2 = tmp3 - tmp2;
    }

    return (tmp2/tmp3)*exp(tmp1)*cval/pow(x,2);
}

double ppp_nru_int(double x, void * p) {
// Integrand of the prior-predictive p-value with truncated normal prior
// for the *relative uncertainty* on the Poisson mean. This version
// should be integrated from 0 to 1.
    struct poiParams * params = (struct poiParams *)p;
    double nObs        = (params->nObs);
    double poiMean     = (params->poiMean);
    double gauPoiRatio = (params->gauPoiRatio);
    double coeffOfVar  = (params->coeffOfVar);
    bool   excess      = (params->excess);

    double cval = max(1.0, nObs);
    double y    = cval * (1.0-x)/x;

    int    which, status;
    double tmp1, tmp2, tmp3, tmp4, qval, bound, ulim, mean, sd;
    which = 1;
    mean  = 0.0;
    sd    = 1.0;
    ulim  = 1.0/coeffOfVar;
    cdfnor( &which, &tmp3, &qval, &ulim, &mean, &sd, &status, &bound );
    ulim  = (gauPoiRatio*y-poiMean)/(gauPoiRatio*y*coeffOfVar);
    cdfnor( &which, &tmp4, &qval, &ulim, &mean, &sd, &status, &bound );
    if (excess) {
        tmp1 = (nObs-1)*log(y) - y - gsl_sf_lngamma(nObs);
        tmp2 = tmp3 - tmp4;
    } else {
        tmp1 = nObs*log(y) - y - gsl_sf_lngamma(nObs+1);
        tmp2 = tmp4;
    }

    return (tmp2/tmp3)*exp(tmp1)*cval/pow(x,2);
}

double ppp_logn_int(double x, void * p) {
// Integrand of the prior-predictive p-value with lognormal prior
    struct poiParams * params = (struct poiParams *)p;
    double nObs    = (params->nObs);
    double poiMean = (params->poiMean);
    double poiUnc  = (params->poiUnc);
    bool excess    = (params->excess);

    double relUnc2 = pow(poiUnc/poiMean, 2);
    double nu0     = poiMean/sqrt(1+relUnc2);
    double tau     = sqrt(log(1+relUnc2));
    double cval    = max(1.0, nObs);
    double y       = cval * (1.0-x)/x;
    double tmp1, tmp2;
    if (excess) {
        tmp1       = (nObs-1)*log(y) - y - gsl_sf_lngamma(nObs);
        tmp2       = gsl_cdf_ugaussian_P(-log(y/nu0)/tau);
    } else {
        tmp1       = nObs*log(y) - y - gsl_sf_lngamma(nObs+1);
        tmp2       = gsl_cdf_ugaussian_P(log(y/nu0)/tau);
    }

    return tmp2*exp(tmp1)*cval/pow(x,2);
}

double fid_p_int(double x, void * p) {
// Integrand of the fiducial p-value
    struct poiParams * params = (struct poiParams *)p;
    double nObs    = (params->nObs);
    double poiMean = (params->poiMean);
    double poiUnc  = (params->poiUnc);

    double cval    = max(1.0, nObs);
    double y       = cval * (1.0-x)/x;
//    double tmp1    = (nObs-1)*log(y) - y - gsl_sf_lngamma(nObs);
    double tmp1    = nObs*log(y) - y - gsl_sf_lngamma(nObs+1);
    double tmp2    = gsl_cdf_ugaussian_P((poiMean-y)/poiUnc);

    return tmp2*exp(tmp1)*cval/pow(x,2);
}

double api_pvalue(void * p) {
// Adjusted plug-in p-value
    struct poiParams * params = (struct poiParams *)p;
    double nObs    = (params->nObs);
    double poiMean = (params->poiMean);
    double poiUnc  = (params->poiUnc);
    bool excess    = (params->excess);

    const double epsi=1.0e-08;
    int ierror, acc=0;
    double pupi, qupi, xtld, tmp1, tmp2;
    double dnu2  = pow(poiUnc, 2);
    double tmp   = 0.5 * (poiMean - dnu2);
    double nuEst = tmp + sqrt(pow(tmp,2) + nObs*dnu2);
    double x0 = 0;
    double sum = 0;
    if (excess)
    {
        if (nObs > 0) {
            gamma_inc( &nObs, &nuEst, &pupi, &qupi, &acc );
            for (double nVal = 1, term = 1; (nVal <= 2*nuEst) || (term > epsi*sum); nVal++)
            {
                gamma_inc_inv( &nVal, &xtld, &x0, &pupi, &qupi, &ierror );
                if (ierror < 0) {
                    cout << "Error from gamma_inc_inv: " << ierror << endl;
                }
                xtld = xtld + (1-nVal/xtld)*dnu2;
                tmp1 = -nuEst + nVal*log(nuEst) - gsl_sf_lngamma(nVal+1);
                tmp2 = gsl_cdf_ugaussian_P((xtld-nuEst)/poiUnc);
                term = tmp2 * exp(tmp1);
                sum += term;
            }
        } else {
            sum = 1;
        }
    } else {
        double aVal=nObs+1;
        gamma_inc( &aVal, &nuEst, &qupi, &pupi, &acc );
        for (double nVal = 0, term = 1; (nVal <= 2*nuEst) || (term > epsi*sum); nVal++)
        {
            aVal = nVal + 1;
            gamma_inc_inv ( &aVal, &xtld, &x0, &qupi, &pupi, &ierror );
            if (ierror < 0) {
                cout << "Error from gamma_inc_inv: " << ierror << endl;
            }
            xtld = xtld + (1-nVal/xtld)*dnu2;
            tmp1 = -nuEst + nVal*log(nuEst) - gsl_sf_lngamma(aVal);
            tmp2 = gsl_cdf_ugaussian_P((nuEst-xtld)/poiUnc);
            term = tmp2 * exp(tmp1);
            sum += term;
        }
    }

    return sum;
}
