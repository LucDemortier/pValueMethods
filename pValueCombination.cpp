#include <iostream>
#include <iomanip>
#include <vector>
#include <sstream>
#include <math.h>
#include <gsl/gsl_cdf.h>

using namespace std;

int main()
{
    vector<double> pValues;
    string input;
    cout << "Enter p-values, one per line, end with an empty line:" << endl;
    while (getline(cin, input)) {
        if (input == "") {
            break;
        }
        double number;
        stringstream ss(input);
        while (ss >> number) {
            pValues.push_back(number);
        }
    }
    cout << "P-Value      Nsigmas" << endl;
    cout << "--------------------" << endl;

// Fisher's method
    double tStat1 = 0;
    for (vector<double>::const_iterator i = pValues.begin(); i != pValues.end(); ++i) {
        tStat1 += log(*i);
    }
    tStat1 *= -2;
    double nu = 2*pValues.size();
    double pComb1 = gsl_cdf_chisq_Q( tStat1, nu);
    double nSig1 = gsl_cdf_ugaussian_Qinv(pComb1);
    cout << setw(11) << left << pComb1 << "  " << setw(7) << left << nSig1 << "  (Fisher)" << endl;

// Fisher's method with chisquare distributions with different numbers of degrees of freedom
    double nDegF = 100;
    double sumDegF = 0;
    double tStat2 = 0;
    double chi2Val;
    for (vector<double>::const_iterator i = pValues.begin(); i != pValues.end(); ++i) {
        chi2Val = gsl_cdf_chisq_Qinv(*i, nDegF);
        tStat2 += chi2Val;
        sumDegF += nDegF;
    }
    double pComb2 = gsl_cdf_chisq_Q (tStat2, sumDegF);
    double nSig2 = gsl_cdf_ugaussian_Qinv(pComb2);
    cout << setw(11) << left << pComb2 << "  " << setw(7) << left << nSig2 << "  (Fisher with nDegF = " << nDegF << ")" << endl;

// Tippett's method
    double tStat3 = *min_element(pValues.begin(), pValues.end());
    double aval = 1.0;
    double bval = pValues.size();
    double pComb3 = gsl_cdf_beta_P(tStat3, aval, bval);
    double nSig3 = gsl_cdf_ugaussian_Qinv(pComb3);
    cout << setw(11) << left << pComb3 << "  " << setw(7) << left << nSig3 << "  (Tippett)" << endl;

// Stouffer's method
    double tStat4 = 0;
    for (vector<double>::const_iterator i = pValues.begin(); i != pValues.end(); ++i) {
        tStat4 += gsl_cdf_ugaussian_Qinv(*i);
    }
    tStat4 /= sqrt(pValues.size());
    double pComb4 = gsl_cdf_ugaussian_Q(tStat4);
    double nSig4 = gsl_cdf_ugaussian_Qinv(pComb4);
    cout << setw(11) << left << pComb4 << "  " << setw(7) << left << nSig4 << "  (Stouffer)" << endl;

    cout << endl;
    return 0;
}
