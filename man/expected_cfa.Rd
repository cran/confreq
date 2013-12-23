\name{expected_cfa}
\alias{expected_cfa}
\title{Expected frequencies with glm}
\usage{
  expected_cfa(des, observed, family = poisson(),
    intercept = FALSE, ...)
}
\arguments{
  \item{des}{a designmatrix (object of class
  \code{"matrix"}) as returned by function
  \code{design_cfg_cfa}.}

  \item{observed}{a integer vector with
  \code{lenght(observed) == dim(des)[1]}. WARNING: The
  observed frequencies counts must be in an order
  corresponding to the coding sheme in designmatix (see
  argument \code{des}).}

  \item{family}{argument passed to \code{\link{glm.fit}}
  with default set to \code{poisson()}}

  \item{intercept}{argument passed to \code{glm.fit} with
  default set to \code{FALSE}}

  \item{...}{aditional arguments optional passed to
  \code{\link{glm.fit}}}
}
\value{
  An vector object giving the expected counts.
}
\description{
  Calculates the expected frequencies of counts using log
  liniear model. #'
}
\details{
  No details
}
\examples{
#######################################
# expected counts for LienertLSD data example.
designmatrix<-design_cfg_cfa(kat=c(2,2,2)) # generate an designmatrix (only main effects)
data(LienertLSD) # load example data
observed<-LienertLSD[,4] # extract observed counts
expected_cfa(des=designmatrix, observed=observed) # calculation of expected counts
#######################################
}
\references{
  No references in the moment
}

