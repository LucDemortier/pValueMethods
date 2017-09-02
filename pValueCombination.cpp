#include <iostream>
#include <iomanip>
#include <vector>
#include <sstream>
#include <math.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_math.h>

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
    double numPvalues = pValues.size();
    cout << "Read in " << numPvalues << " p-values.\n" << endl;
    cout << "    Combination:    " << endl;
    cout << "P-Value      Nsigmas" << endl;
    cout << "---------------------" << endl;

// Fisher's method
    double tStat1 = 0;
    for (vector<double>::const_iterator i = pValues.begin(); i != pValues.end(); ++i) {
        tStat1 += log(*i);
    }
    tStat1 *= -2;
    double nu = 2*numPvalues;
    double pComb1 = gsl_cdf_chisq_Q( tStat1, nu);
    double nSig1 = gsl_cdf_ugaussian_Qinv(pComb1);
    cout << setw(11) << left << pComb1 << "  " << setw(8) << left << nSig1 << "  (Fisher)" << endl;

// Fisher's method with chisquare distributions with different numbers of degrees of freedom
    double nDegF2 = 100;
    double sumDegF = 0;
    double tStat2 = 0;
    double chi2Val;
    for (vector<double>::const_iterator i = pValues.begin(); i != pValues.end(); ++i) {
        chi2Val = gsl_cdf_chisq_Qinv(*i, nDegF2);
        tStat2 += chi2Val;
        sumDegF += nDegF2;
    }
    double pComb2 = gsl_cdf_chisq_Q (tStat2, sumDegF);
    double nSig2 = gsl_cdf_ugaussian_Qinv(pComb2);
    cout << setw(11) << left << pComb2 << "  " << setw(8) << left << nSig2 << "  (Fisher with nDegF = " << nDegF2 << ")" << endl;

// Tippett's method
    double tStat3 = *min_element(pValues.begin(), pValues.end());
    double aval = 1.0;
    double pComb3 = gsl_cdf_beta_P(tStat3, aval, numPvalues);
    double nSig3 = gsl_cdf_ugaussian_Qinv(pComb3);
    cout << setw(11) << left << pComb3 << "  " << setw(8) << left << nSig3 << "  (Tippett)" << endl;

// Stouffer's method
    double tStat4 = 0;
    for (vector<double>::const_iterator i = pValues.begin(); i != pValues.end(); ++i) {
        tStat4 += gsl_cdf_ugaussian_Qinv(*i);
    }
    tStat4 /= sqrt(numPvalues);
    double pComb4 = gsl_cdf_ugaussian_Q(tStat4);
    double nSig4 = gsl_cdf_ugaussian_Qinv(pComb4);
    cout << setw(11) << left << pComb4 << "  " << setw(8) << left << nSig4 << "  (Stouffer)" << endl;

// Logistic transformation
    double tStat5 = 0, nDegF5 = 5*numPvalues + 4;
    for (vector<double>::const_iterator i = pValues.begin(); i != pValues.end(); ++i) {
        tStat5 -= log(*i/(1-*i));
    }
    tStat5 /= M_PI * sqrt( numPvalues * (5*numPvalues+2) / (15*numPvalues+12) );
    double pComb5 = gsl_cdf_tdist_Q(tStat5, nDegF5);
    double nSig5 = gsl_cdf_ugaussian_Qinv(pComb5);
    cout << setw(11) << left << pComb5 << "  " << setw(8) << left << nSig5 << "  (Logit combination, approximate)" << endl;

    cout << endl;
    return 0;
}
