# pValueMethods
This is a collection of methods for computing p-values and studying their properties. The following methods are currently available:

1. [**poissonPvalues:**](https://github.com/LucDemortier/pValueMethods/blob/master/poissonPvalues.cpp) computes the p-value corresponding to a Poisson observation, when the mean of the Poisson is uncertain. Several methods are used to incorporate this uncertainty into the p-value: prior-predictive (with truncated Gaussian, gamma, and log-normal priors); bootstrap (plug-in and adjusted plug-in); and fiducial.
2. [**gaussianPvalues:**](https://github.com/LucDemortier/pValueMethods/blob/master/gaussianPvalues.cpp) computes the p-value corresponding to a Gaussian observation, when the mean of the Gaussian is uncertain.
3. [**pValueCombination:**](https://github.com/LucDemortier/pValueMethods/blob/master/pValueCombination.cpp) combines an arbitrary number of *independent* p-values. Several combination methods are compared: Fisher, Tippett, Stouffer, the logit transform, Simes, Edgington, and Wilkinson.

This software uses the GNU Scientific Library (GSL) as well as  [**cdflib**](https://github.com/LucDemortier/pValueMethods/tree/master/cdflib), a collection of routines for cumulative distribution functions, their inverses, and other parameters, compiled and written by Barry W. Brown, James Lovato, and Kathy Russell.

A document is being prepared that discusses p-value methods and properties in greater detail.
