\name{LienertLSD}
\alias{LienertLSD}
\docType{data}
\title{The Linert LSD Data}
\description{
Data from the classical Linert LSD trial as an example for CFA. The data is in pattern frequencies table representation (object of class c("data.frame", "Pfreq" )). 
}
\usage{
  data(LienertLSD)
}

\format{
  A matrix with 4 columns and 8 rows. The last column gives the frequencies for the (response) pattern in column 1:3. 
}

\details{
The first three columns are named C, T and A which are abreviations for 
the observed symptoms after taking LSD:

C = narrowed consciousness

T = thought disturbance

A = affective disturbance
}

\examples{
data(LienertLSD)
dim(LienertLSD)
##############################################################
colnames(LienertLSD) # show all variable names of matrix LienertLSD
}

\references{
Lienert, G. A. (1971). Die Konfigurationsfrequenzanalyse: I. Ein neuer Weg zu Typen und Syndromen. \emph{Zeitschrift für Klinische Psychologie und Psychotherapie, 19}(2), 99-115.
}
\keyword{datasets}
