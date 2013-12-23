\name{suicide}
\alias{suicide}
\docType{data}
\title{The Linert suicide Data}
\description{
Data from the Linert suicide example for CFA. The data is in data list representation (each row is one case). 
}
\usage{
  data(suicide)
}

\format{
  A \code{data.frame} with 3 columns (as factors). 
}

\details{
The three columns are named 'Geschlecht', 'Epoche' and 'Suizidart' which is 'gender', 'epoch' and 'type od suicide'. each of the variables are factors with the following levels:

Geschlecht: 'm' = 1 (male); 'w' = 2 (female) 

Epoche: '44' = 1 (the epoch 1944); '52' = 2 (the epoch 1952)

Suizidart: 'Eh' = 1(hang); 'Es' = 2 (shoot); 'Et' = 3(drown);  
'G' = 4(gas);  'H' = 5(crashing down);  'P' = 6(open vein);  'S' = 7(barbiturate);
}

\examples{
##############################################################
data(suicide) # to load the data.frame included in the package
class(suicide)
dim(suicide)
str(suicide)
}

\references{
Krauth, J., & Lienert, G. A. (1973). \emph{Die Konfigurationsfrequenzanalyse (KFA) und ihre Anwendung in Psychologie und Medizin: ein multivariates nichtparametrisches Verfahren zur Aufdeckung von Typen und Syndromen; mit 70 Tabellen}. Freiburg; München: Alber Karl.
}
\keyword{datasets}
