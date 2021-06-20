\name{newborns}
\alias{newborns}
\docType{data}
\title{The Data Example from Stemmler 2020}
\description{
data example by Stemmler (2020, table 4.1, p. 33) where \emph{n} = 56 newborns 'with seizures' = 1 or 'without seizures' = 2 (coded in the in first column named 'A') were tested with an intelligence test while they attended kindergarten. Children’s intelligence was divided into 'average or above' = 1 and 'below average' = 2 (coded in the in second column named 'B'). The third column gives the frequencies of the respective pattern.
}
\usage{
  data(newborns)
}

\format{
  A data.frame with 3 columns and 4 rows. The last column gives the frequencies for the observed pattern in column 1:2. The data is in pattern frequencies table representation (object of class c("data.frame", "Pfreq")) 
}

\details{
No detail in the moment
}

\examples{
##############################################################
data(newborns)
dim(newborns)
newborns
##############################################################
}

\references{
Stemmler, M. (2020). \emph{Person-Centered Methods -- Configural Frequency Analysis (CFA) and Other Methods for the Analysis of Contingency Tables}. Cham Heidelberg New York Dordrecht London: Springer.
}
\keyword{datasets}
