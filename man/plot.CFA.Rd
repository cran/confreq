\name{plot.CFA}
\alias{plot.CFA}
\title{S3 plot for CFA}
\usage{
\method{plot}{CFA}(
  x,
  type = "z.pChi",
  holm = FALSE,
  fill=c("red", "blue", "grey"),
  adjalpha = "bonferroni",
  ...
)
}
\arguments{
\item{x}{object of class \code{"CFA"}}

\item{type}{character indicating which test to use for visualizing whether the observed pattern are 'Types', 'Antitypes' or not significant at all. Possible options for \code{type} are \code{"pChi"}, \code{"ex.bin.test"}, \code{"z.pChi"}, \code{"z.pBin"} and \code{"p.stir"}.}

\item{holm}{logical with default \code{holm = FALSE}. If set to \code{holm = TRUE}, significance testing is based on the holm procedure -- see references. This argument is deprecated (since version 1.5.6) and kept only for downward compatibility. Use argument \code{adjalpha} for any type of alpha adjustment.}

\item{fill}{a vector of (three) colors defining the coloring of significant 'Types' (default "red"), 'Antitypes' (default "blue") or not significant cells (default "grey") in the plot.}

\item{adjalpha}{character with default \code{adjalpha = "bonferroni"}. Selector for the type of alpha adjustment for multiple testing. Possible options are: \code{adjalpha = "none"}, for no adjustment; \code{adjalpha = "bonferroni"}, for bonferroni adjustment (default); \code{adjalpha = "holm"}, for alpha adjustment according to Holm (1979); other options to come ... .} 

\item{...}{other parameters passed trough}
}
\description{
S3 plot method for object of class\code{"CFA"}
}
\references{
Holm, S. (1979). A simple sequentially rejective multiple test procedure. \emph{Scandinavian Journal of Statistics, 6}(2), 65–70.
}
