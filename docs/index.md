### dynardl
Stata module to dynamically simulate autoregressive distributed lag (ARDL) models.

### Download 
You can download the most recent version of `dynardl` from the project site [here](https://github.com/andyphilips/dynardl/archive/master.zip). This program is part of a suite that also includes `pssbounds` (Jordan and Philips 2017), a Stata module to display the necessary critical values to conduct the Pesaran, Shin and Smith (2001) bounds test for cointegration.

### Version
Current version: 1.0.3. Note that a previous version of this program was called `dynpss`. `dynardl` provides a more flexible lag specification and adds additional plots.

### Table of Contents
 * [Description](#description)
 * [Reference](#reference)
 * [Authors](#authors)
 * [Citations](#citations)
 * [Examples](#examples)
 * [Example Papers](#example-papers)
 
### Description<a id="description"></a>
`dynardl` is a program to produce dynamic simulations of autoregressive distributed lag models (ARDL) of the sort recommended by Pesaran, Shin, and Smith (2001). See Philips (2017) for a discussion of this approach, and [Jordan and Philips](dynardl/pss Stata 2017.pdf) (2017) for an in-depth discussion of this program.

`dynardl` is designed to dynamically simulate the effects of a counterfactual change in one weakly exogenous regressor at a point in time, using stochastic simulation techniques. Since the ARDL procedure can produce models that are complicated to interpret, `dynardl` is designed to ease the burden of substantive interpretations through the creation of predicted (or expected) values of the dependent variable (along with associated confidence intervals), which can be plotted to show how a change in one variable "flows" through the model over time. `dynardl` takes 1000 (or however many simulations a user desires) draws of the set of parameters from a multivariate normal distribution, using the estimated parameters and the variance-covariance matrix from the linear regression. All covariates are set to certain values (typically means), which are used to create predicted Y-hat values plus stochastic uncertainty.

 
### Reference<a id="reference"></a>
If you use `dynardl`, please cite:

Jordan, Soren and Andrew Q. Philips. 2017 "[Cointegration testing and dynamic simulations of autoregressive distributed lag models](dynardl/pss Stata 2017.pdf)". Working Paper.

and

Philips, Andrew Q. 2017. "[Have your cake and eat it too? Cointegration and dynamic inference from autoregressive distributed lag models](http://dx.doi.org/10.1111/ajps.12318)." American Journal of Political Science.

### Authors<a id="authors"></a>

[Andrew Q. Philips](http://www.andyphilips.com), Department of Political Science, University of Colorado at Boulder. andrew.philips [AT] colorado.edu. @andyphilips

[Soren Jordan](http://sorenjordan.com), Department of Political Science, Auburn University. sorenjordanpols [AT] gmail.com.

### Citations<a id="citations"></a>

Pesaran, M Hashem, Yongcheol Shin and Richard J Smith. 2001. "Bounds testing approaches to the analysis of level relationships." Journal of Applied Econometrics 16(3):289-326.

### Examples<a id="examples"></a>

See the [working paper](dynardl/pss Stata 2017.pdf) for examples of `dynardl` in action.

### Example Papers<a id="example-papers"></a>
Use `dynardl` in one of your papers? Let me know (andrew.philips [AT] colorado.edu) and I will add it to the list below:
