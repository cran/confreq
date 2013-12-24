\name{lazar}
\alias{lazar}
\docType{data}
\title{The Data Example from Lazarsfeld and Henry}
\description{
data example by Lazarsfeld and Henry (1968) where \emph{N} = 1000 subjects need to solve questions or problems (i.e., A,B, and C). They either ‘1’ = solved or ‘2’ = did not solve the problems. The data is in pattern frequencies table representation (object of class c("data.frame", "Pfreq" )). 
}
\usage{
  data(lazar)
}

\format{
  A matrix with 4 columns and 8 rows. The last column gives the frequencies for the (response) pattern in column 1:3. 
}

\details{
No detail in the moment
}

\examples{
##############################################################
data(lazar)
dim(lazar)
##############################################################
}

\references{
Lazarsfeld, P. F., & Henry, N. W. (1968). \emph{Latent structure analysis}. Boston: Houghton Mifflin.
}
\keyword{datasets}
