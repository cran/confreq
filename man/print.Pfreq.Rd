% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.Pfreq.R
\name{print.Pfreq}
\alias{print.Pfreq}
\title{S3 print for Pfreq}
\usage{
\method{print}{Pfreq}(
  x,
  ...,
  digits = NULL,
  quote = FALSE,
  right = TRUE,
  row.names = TRUE,
  max = NULL
)
}
\arguments{
\item{x}{object of class \code{"Pfreq"}}

\item{...}{further arguments passed to or from other methods.}

\item{digits}{minimal number of significant digits, see \code{\link{print.default}}.}

\item{quote}{logical, indicating whether or not strings should be printed with surrounding quotes.}

\item{right}{logical, indicating whether or not strings should be right aligned.}

\item{row.names}{logical (or character vector), indicating whether (or what) row names should be printed.}

\item{max}{numeric or \code{NULL}, specifying the maximal number of entries to be printed. By default, when \code{NULL}, \code{\link{getOption}("max.print")} used.}
}
\value{
output printed to the console
}
\description{
S3 print method for object of class\code{"Pfreq"}
}
\keyword{methods}
