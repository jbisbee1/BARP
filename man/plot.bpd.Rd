% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.bpd.R
\name{plot.bpd}
\alias{plot.bpd}
\title{plot.bpd}
\usage{
\method{plot}{bpd}(
  bpd,
  var_names = NULL,
  var_labs = NULL,
  is_categorical = NULL,
  cols = c("#f2f0f7", "#cbc9e2", "#9e9ac8", "#6a51a3"),
  ...
)
}
\arguments{
\item{bpd}{A \code{\link{barp_partial_dependence}} partial dependence object of class "bpd".}

\item{var_names}{A vector of variable names. If \code{NULL} (the default), the variable names from the training data columns are used.}

\item{var_labs}{A list of variable labels corresponding to the levels at which the partial dependence was calculated. If \code{NULL} (the default), the raw values are used.}

\item{is_categorical}{A vector of logicals (T,F) corresponding to which variables are categorical and which are continuous.If \code{NULL} (the default), the function determines whether each variable is categorical or continuous based on the data class.}

\item{...}{Additional arguments to be passed to plot.}
}
\description{
This function plots the partial dependence results estimated returned by running barp_partial_dependence.
}
\examples{
data("gaymar")
barp.obj <- barp(y = "supp_gaymar",
             x = c("pvote","religcon","age","educ","gXr","stateid","region"),
             dat = svy,
             census = census06,
             geo.unit = "stateid",
             proportion = "n")
bpd <- barp_partial_dependence(BARP = barp.obj,
                               vars = c("age","educ"),
                               prop_data = .2,
                               levs = list(c(1:4),c(1:4),
                               credible_interval = c(0.025,0.975)))
plot(bpd,
     var_names = c("Age","Education"),
     var_labs = list(c("18-30","31-50","51-65","65+"),
                     c("LTHS","HS","Some Coll","Coll+")),
     is_categorical = c(T,T))
}
\seealso{
\code{\link{barp_partial_dependence}} which is a prerequisite to plotting.
}
\keyword{BARP}
\keyword{MRP}
\keyword{Mister}
\keyword{P}
\keyword{multilevel}
\keyword{poststratification}
\keyword{regression}
