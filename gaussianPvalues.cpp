#include <iostream>
#include <iomanip>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cdf.h>

using namespace std;

int main()
{
    double Obs;
    double gauMean;
    double uncMean;
    double gauVar;
    double gauStD;
    double pAdjustment;
    bool excess;
    string bline(72, '-');

    cout << '\n' << bline << endl;
    cout << "Observation: ";
    cin  >> Obs;
    cout << "Gaussian mean: ";
    cin  >> gauMean;
    excess = Obs >= gauMean;
    cout << "         variance: ";
    cin  >> gauVar;
    gauStD = sqrt(gauVar);
    cout << "Uncertainty on Gaussian mean: ";
    cin  >> uncMean;
    cout << "P-value adjustment factor: ";
    cin  >> pAdjustment;

    cout << "\nGaussian mean: " << gauMean << " +/- " << uncMean << ", standard deviation: " << gauStD << ", observation: " << Obs << endl;
    cout << "P-value adjustment factor: " << pAdjustment << endl;
    if (excess) {
        cout << "Computing the significance of an *excess*." << endl;
    } else {
        cout << "Computing the significance of a *deficit*." << endl;
    }
    cout << "\nP-Value      Nsigmas" << endl;
    cout << "---------------------" << endl;

// First ignore uncertainty on Gaussian mean when computing p-value
    double pVal0, nSig0;
    if (excess) {
        pVal0 = gsl_cdf_ugaussian_Q((Obs-gauMean)/gauStD);
    } else {
        pVal0 = gsl_cdf_ugaussian_P((Obs-gauMean)/gauStD);
    }
    pVal0 *= pAdjustment;
    nSig0  = gsl_cdf_ugaussian_Qinv(pVal0);
    cout << setw(11) << left << pVal0 << "  " << setw(8) << left << nSig0 << "  (ignoring uncertainty on Gaussian mean)" << endl;

// Try a Gaussian prior for the Gaussian mean
    double combStD, pVal1, nSig1;
    combStD = sqrt( pow(uncMean,2) + pow(gauStD,2) );
    if (excess) {
        pVal1 = gsl_cdf_ugaussian_Q((Obs-gauMean)/combStD);
    } else {
        pVal1 = gsl_cdf_ugaussian_P((Obs-gauMean)/combStD);
    }
    pVal1 *= pAdjustment;
    nSig1  = gsl_cdf_ugaussian_Qinv(pVal1);
    cout << setw(11) << left << pVal1 << "  " << setw(8) << left << nSig1 << "  (prior-pred., Gaussian prior)" << endl;

    cout << bline << '\n' << endl;
    return 0;
}
