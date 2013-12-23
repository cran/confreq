\name{stirling_cfa}
\alias{stirling_cfa}
\title{Approximation to the binomial using Stirling's Formula}
\usage{
  stirling_cfa(observed, expected = NULL,
    n = sum(observed), p = NULL, cum = T, verb = T)
}
\arguments{
  \item{observed}{a integer vector with observed
  freqencies}

  \item{expected}{a vector giving the expected frequencies.
  expected can be set to \code{expected=NULL} if an vector
  of cell probabilities is given in argument \code{p}.}

  \item{n}{number of trials (scalar) default is \code{n =
  sum(observed) }.}

  \item{p}{a vector of cell probabilities. If p is not NULL
  the argument \code{expected} is ignored and this vector p
  of cell probabilities is used for calculatio instead of
  expected counts}

  \item{cum}{a logical - computation of cumulative density.
  If \code{cum=TRUE} (default) computes tail probability.
  If \code{cum=FALSE} computes prob. only for one cell
  (i.e. execute stircore only).}

  \item{verb}{logical - verbose results: If
  \code{verb=TRUE} (default) builds a results table. If
  \code{verb=FALSE} returns vector of cell p-values only.}
}
\description{
  Calculates the binomial aproximation using stirling's
  formula (Version of function: V 1.0 - November 2013)
}
\details{
  \itemize{ \item{Vector p must be of same length as
  observed _or_ p may be a scalar (e.g. in case of the
  zero-order CFA).} \item {The routine autoselects the
  upper or lower tail: \itemize{ \item{if obs > exp then
  sum obs:n} \item{ else sum 0:obs}}} \item {The stirling
  approximation cannot be evaluated if the observed
  frequency is 0 or n. Therefore, the proposal of A. von
  Eye (20xx) is adopted, taking the sum up to 1 or n-1,
  respectively.}}
}
\author{
  R.W. Alexandrowicz
}
\references{
  von Eye, A. (2002). \emph{Configural Frequency Analysis.
  Methods, Models, and Applications.} Mahwah, NJ, LEA.
}

