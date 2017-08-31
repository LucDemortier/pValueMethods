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
    double gauStD;
    bool excess;

    cout << "Observation: ";
    cin  >> Obs;
    cout << "Gaussian mean: ";
    cin  >> gauMean;
    excess = Obs >= gauMean;
    cout << "Gaussian standard deviation: ";
    cin  >> gauStD;

    cout << "\nGaussian mean: " << gauMean << " +/- " << gauStD << ", observation: " << Obs << endl;
    if (excess) {
        cout << "Computing the significance of an *excess*." << endl;
    } else {
        cout << "Computing the significance of a *deficit*." << endl;
    }
    cout << "\nP-Value      Nsigmas" << endl;
    cout << "--------------------" << endl;
    double pVal0, nSig0;
    if (excess) {
        pVal0 = gsl_cdf_ugaussian_Q((Obs-gauMean)/gauStD);
    } else {
        pVal0 = gsl_cdf_ugaussian_P((Obs-gauMean)/gauStD);
    }
    nSig0  = gsl_cdf_ugaussian_Qinv(pVal0);

    cout << setw(11) << left << pVal0 << "  " << setw(7) << left << nSig0 << endl;

    return 0;
}
