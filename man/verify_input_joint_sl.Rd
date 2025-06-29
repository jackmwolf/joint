% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/verify_input.R
\name{verify_input_joint_sl}
\alias{verify_input_joint_sl}
\title{Verify input to \code{\link[=joint_sl]{joint_sl()}}}
\usage{
verify_input_joint_sl(data0, endpoints, treatment, primary, n_boot, ci_level)
}
\arguments{
\item{data0}{Data frame with treatment assignment and endpoints}

\item{endpoints}{Character vector of endpoint names as in \code{data0}}

\item{treatment}{Character name of the column name corresponding to binary
0-1 treatment indicator.}

\item{primary}{Name of the primary endpoint}

\item{n_boot}{Number of bootstrap samples}

\item{ci_level}{Confidence level}
}
\value{
TRUE if all input is valid
}
\description{
Verify input to \code{\link[=joint_sl]{joint_sl()}}
}
