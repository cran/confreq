\name{Lienert1978}
\alias{Lienert1978}
\docType{data}
\title{The Lienert (1978) Data}
\description{
Data used as an example for two-sample CFA in the textbook by Mark Stemmler (2014) taken from Linert (1978, p. 978). The data is in pattern frequencies table representation (object of class c("data.frame", "Pfreq" )). 
}
\usage{
  data(Lienert1978)
}

\format{
  A data frame (object of class \code{c("data.frame", "Pfreq")} ) with 4 columns and 12 rows. The last column gives the frequencies for the (response) pattern in column 1:2 of the respective 'Group' given in column 3. 
}

\details{
no details at the moment ...
}

\examples{
data(Lienert1978)
dim(Lienert1978)
##############################################################
colnames(Lienert1978) # show all variable names of Lienert1978
}

\references{
Lienert, G. A. (1978). \emph{Verteilungsfreie Methoden in der Biostatistik (Band II)} [Non-parametrical 168 methods in the field of biometrics (Vol. II)]. Meisenheim am Glan, Germany: Hain.

Stemmler, M. (2014). \emph{Person-Centered Methods -- Configural Frequency Analysis (CFA) and Other Methods for the Analysis of Contingency Tables}. Cham Heidelberg New York Dordrecht London: Springer.
}
\keyword{datasets}
