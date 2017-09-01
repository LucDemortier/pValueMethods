# pValueMethods
This is a collection of methods for computing p-values and studying their properties. The following methods are currently available:

1. [**poissonPvalues:**](https://github.com/LucDemortier/pValueMethods/blob/master/poissonPvalues.cpp) computes the p-value corresponding to a Poisson observation, when the mean of the Poisson is uncertain. Several methods are used to quantify the uncertainty: prior-predictive with various forms for the prior, bootstrap (plug-in and adjusted plug-in), and fiducial.
2. [**gaussianPvalues:**](https://github.com/LucDemortier/pValueMethods/blob/master/gaussianPvalues.cpp) computes the p-value corresponding to a Gaussian observation, when the mean of the Gaussian is uncertain.
3. [**pValueCombination:**](https://github.com/LucDemortier/pValueMethods/blob/master/pValueCombination.cpp) combines an arbitrary number of p-values using a number of methods, for purposes of comparison. Includes the methods of Fisher, Tippett, and Stouffer.

A document is being prepared that discusses these methods (and others) in greater detail.
