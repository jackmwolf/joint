% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_sem.R
\name{estimate_effects}
\alias{estimate_effects}
\title{Estimate treatment effects from a structural equation model}
\usage{
estimate_effects(model, risk_difference = c(), sandwich = FALSE)
}
\arguments{
\item{model}{An object of class \link{joint_sem}}

\item{risk_difference}{Names of endpoints for which to calculate the risk
difference.}

\item{sandwich}{Logical indicator to use the sandwich variance estimator}
}
\value{
A list of effect estimates and variances
}
\description{
Estimate treatment effects from a structural equation model
}
\examples{
data(joint_example)
fit <- joint_sem(
  joint_example, endpoints = c("Y1", "Y2", "Y3_cat"))
fit
estimate_effects(fit, risk_difference = "Y3_cat")

}
