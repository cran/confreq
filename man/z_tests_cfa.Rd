% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/z_tests_cfa.R
\name{z_tests_cfa}
\alias{z_tests_cfa}
\title{Two z-Approximation Tests}
\usage{
z_tests_cfa(observed, expected, ccor = FALSE, ntotal = sum(observed))
}
\arguments{
\item{observed}{a vector giving the observed frequencies.}

\item{expected}{a vector giving the expected frequencies.}

\item{ccor}{either a logical (TRUE / FALSE) determining wether to apply a continuity correction or not to the Binomial Approximation of the z-Test. When set to \code{ccor=TRUE} continuity correction is applied for expected values 5 =< expected =< 10. For \code{ccor=FALSE} no continuity correction is applied. Another option is to set \code{ccor=c(x,y)} where x is the lower and y the upper bound for expected values where continuity correction is applied. So \code{ccor=c(5,10)} is equivalent to \code{ccor=TRUE}.}

\item{ntotal}{optional a numeric giving the total number of observations. By default ntotal is calculated as \code{ntotal=sum(observed)}.}
}
\value{
a list with z and p-values.
}
\description{
Calculates the Chi-square approximation to the z-test and the binomial approximation to the z-test.
}
\details{
An continuity correction can be applied to the binomial approximation -- see argument \code{ccor}.
}
\examples{
#######################################
# expected counts for LienertLSD data example.
designmatrix<-design_cfg_cfa(kat=c(2,2,2)) # generate an designmatrix (only main effects)
data(LienertLSD) # load example data
observed<-LienertLSD[,4] # extract observed counts
expected<-expected_cfa(des=designmatrix, observed=observed) # calculation of expected counts
z_tests_cfa(observed,expected)
####################################### 
}
\references{
No references in the moment
}
\keyword{misc}
