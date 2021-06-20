\name{LienertLSD}
\alias{LienertLSD}
\docType{data}
\title{The Lienert LSD Data}
\description{
Data from the classical Lienert LSD trial as an example for CFA (see Lienert, 1971, p. 103, ’Tabelle 1’). The data is in pattern frequencies table representation (object of class c("data.frame", "Pfreq" )). 
}
\usage{
  data(LienertLSD)
}

\format{
  A data frame (object of class \code{c("data.frame", "Pfreq")} ) with 4 columns and 8 rows. The last column gives the frequencies for the observed pattern of the psychotoxic basic syndrome (in column 1:3) due to the intake of lysergic acid diethylamide (LSD).
}

\details{
The first three columns are named C, T and A which are abbreviations for 
the observed symptoms after taking LSD:

C = narrowed consciousness

T = thought disturbance

A = affective disturbance

The coding of the observations is for all symptoms: present='+' and absent='-'
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
