\name{binomial_test_cfa}
\alias{binomial_test_cfa}
\title{Binomial Test}
\usage{
  binomial_test_cfa(observed, expected,
    ntotal = sum(observed))
}
\arguments{
  \item{expected}{a vector giving the expected
  frequencies.}

  \item{observed}{a vector giving the observed
  frequencies.}

  \item{ntotal}{optional a numeric giving the total number
  of observations. By default ntotal is calculated as
  \code{ntotal=sum(observed)}.}
}
\value{
  a numeric giving the p-value.
}
\description{
  Calculates the (exact) binomial test based on obseved,
  expected frequencies an the total number of observations.
  #'
}
\details{
  No details
}
\examples{
#######################################
# first calculate expected counts for LienertLSD data example.
designmatrix<-design_cfg_cfa(kat=c(2,2,2)) # generate an designmatrix (only main effects)
data(LienertLSD) # load example data
observed<-LienertLSD[,4] # extract observed counts
expected<-expected_cfa(des=designmatrix, observed=observed) # calculation of expected counts
 binomial_test_cfa(observed,expected)
#######################################
}
\references{
  No references in the moment
}

