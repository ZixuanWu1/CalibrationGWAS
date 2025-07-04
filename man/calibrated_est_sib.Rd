% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sibling_calibration.R
\name{calibrated_est_sib}
\alias{calibrated_est_sib}
\title{Compute the calibrated estimator}
\usage{
calibrated_est_sib(X, data, alpha_ext, alpha_ext_var, N_ext, overlap_ratio = 0)
}
\arguments{
\item{X:}{a vector that contains transmitted allele values}

\item{data:}{output from format_data}

\item{alpha_ext:}{estimate of alpha from external data}

\item{alpha_ext_var:}{variance of alpha from external data}

\item{N_ext:}{number of samples in external data}

\item{overlap_ratio:}{proportion of internal data from external data}
}
\value{
a list of following elements
(1) beta_cal: calibrated estimate
(2) beta_cal_var: calibrated variance
(3) beta_int: within-family estimate
(4) beta_int_var: within-family variance
(5) alpha_int: marginal estimate
(6) alpha_int_var: marginal variance
}
\description{
Compute the calibrated estimator
}
