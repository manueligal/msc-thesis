# Estimating Causal Effects with Proxies and Domain Shifts

### Main Files

* `ConditionNumber.R`: Comparison of the point estimates from the causal and reduced parameterizations.
* `BaselineComparison.R`: Comparison of the point estimates of the two proposed estimators with three different baselines.
* `ReducedCI.R`: Representation of the coverage and interval length of the asymptotic confidence intervals.
* `Bootstrap.R`: Comparison of the asymptotic and bootstrap confidence intervals.
* `HotelsExample.R`: Application of the model and the reduced estimator to a real example.

### Auxiliary Files
The auxiliary files can be found in the `utils` folder.

* `misc.R`: General helper functions.
* `Lfunction.R`: Definition of the likelihood function.
* `FunctionsSampling.R`: Functions to generate the matrices that describe the model and samples from the entailed distribution.
* `FunctionsEstimation.R`: Functions to apply the estimators from the causal and reduced parameterizations.
* `FunctionsCI.R`: Function to calculate the asymptotic standard deviation in the reduced parameterization.
