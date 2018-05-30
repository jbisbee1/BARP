# BARP

Multilevel Regression and Poststratification (MRP) - Multilevel Regression + Bayesian Additive Regression Trees (BART) = BARP.

This package uses BARP to extrapolate a survey opinion to a finer level of geographic aggregation than the survey was originally intented to capture. It augments MRP (aka Mister P) by using a non-parametric function to estimate opinion. In so doing, it improves the predictive accuracy by relaxing the assumption that (1) the researcher has all necessary covariates and (2) the researcher knows the precise specification by which these covariates correlate with an opinion of interest. 

The user must provide both a survey dataset and a ``census'' dataset with identically named covariates. The census dataset can contain a column of proportions for each covariate bin (i.e., the number or percentage of the population that is a black woman between the ages of 30 and 40) but doesn't have to. The user must also identify which column contains the geographic unit at which she wants to extrapolate opinion and this column must also be included in the covariate list.

BARP produces estimates of an opinion of interest for each geographic unit. If the user chooses, she can also extract confidence intervals for these estimates. These confidence intervals should be used to account for uncertainty in the estimation of opinion when using the opinion estimates in conventional regression contexts. A detailed vignette is included to highlight features of the package.

The package can be installed by installing ```devtools``` and typing ```install_github('jbisbee1/BARP')```.