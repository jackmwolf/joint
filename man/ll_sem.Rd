% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_sem.R
\name{ll_sem}
\alias{ll_sem}
\title{Calculate the log likelihood for a joint model using SEM constraints}
\usage{
ll_sem(phi, data0, endpoints, categorical, treatment, .names = NULL)
}
\arguments{
\item{phi}{Named vector of parameter values}

\item{data0}{Data frame with treatment assignment and endpoints}

\item{endpoints}{Character vector of endpoint names as in \code{data0}}

\item{categorical}{Character vector of endpoint names to be treated as
categorical (binary or ordinal) variables}

\item{treatment}{Character name of the column name corresponding to binary
0-1 treatment indicator.}

\item{.names}{Optimal vector of names for phi}
}
\value{
The log likelihood for the input parameters and data
}
\description{
Calculate the log likelihood for a joint model using SEM constraints
}
