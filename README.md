# pValueMethods
This is a collection of methods for computing p-values and studying their properties. The following methods are currently available:

1. [**poissonPvalues:**](https://github.com/LucDemortier/pValueMethods/blob/master/poissonPvalues.cpp) computes the p-value corresponding to a Poisson observation, when the mean of the Poisson is uncertain. Several methods are used to quantify the uncertainty: prior-predictive with various forms for the prior, bootstrap (plug-in and adjusted plug-in), and fiducial.
2. [**gaussianPvalues:**](https://github.com/LucDemortier/pValueMethods/blob/master/gaussianPvalues.cpp) computes the p-value corresponding to a Gaussian observation, when the mean of the Gaussian is uncertain.
3. [**pValueCombination:**](https://github.com/LucDemortier/pValueMethods/blob/master/pValueCombination.cpp) combines an arbitrary number of p-values using a number of methods, for purposes of comparison. Includes the methods of Fisher, Tippett, and Stouffer.

This repository also includes [**cdflib**](https://github.com/LucDemortier/pValueMethods/tree/master/cdflib), a collection of routines for cumulative distribution functions, their inverses, and other parameters. This library was compiled and written by Barry W. Brown, James Lovato, and Kathy Russell. It was originally written in FORTRAN, but a C++ version is available on [John Burkardt's website](http://people.sc.fsu.edu/~jburkardt/cpp_src/cdflib/cdflib.html). I split the original file into individual routines and wrote a Makefile to create a static library. For the p-value project I initially intended to use the GNU Scientific Library (GSL) for all statistical computations, but GSL crashed on some calculations involving the gamma distribution (apparently this is a [known bug](https://lists.gnu.org/archive/html/bug-gsl/2011-10/msg00014.html)); cdflib appears to be more robust.

A document is being prepared that discusses p-value methods and properties in greater detail.
