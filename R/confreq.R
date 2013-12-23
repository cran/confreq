#' @title Configural Frequencies Analysis Using Loglinear Modeling
#' @name confreq
#' @docType package
#' @description The package \code{confreq} offers some functions for Configural Frequencies Analysis (CFA) proposed by G.A. Lienert as an analysis of types and antitypes of response pattern. The core principle in the package \code{confreq} is to use the function \code{\link{glm}} to compute the expected counts based on a model (design) matrix.
#'
#' @details 
#' there are no further details at this point.
#' 
#' For further description see description of functions.
#' 
#' All the other things we already discussed in Klagenfurt at FGME 2013 with Rainer A., Mark S. \dots   
#'
#' @author \itemize{\item{Joerg-Henrik Heine <jhheine@@googlemail.com>}\item{R.W. Alexandrowicz}}
#' @references Krauth, J., & Lienert, G. A. (1973). \emph{Die Konfigurationsfrequenzanalyse (KFA) und ihre Anwendung in Psychologie und Medizin: ein multivariates nichtparametrisches Verfahren zur Aufdeckung von Typen und Syndromen; mit 70 Tabellen}. Freiburg; München: Alber Karl.
#' @references Lazarsfeld, P. F., & Henry, N. W. (1968). \emph{Latent structure analysis}. Boston: Houghton Mifflin.
#' @references Lienert, G. A. (1971). Die Konfigurationsfrequenzanalyse: I. Ein neuer Weg zu Typen und Syndromen. \emph{Zeitschrift für Klinische Psychologie und Psychotherapie, 19}(2), 99-115.
#' @references von Eye, A. (2002). \emph{Configural Frequency Analysis. Methods, Models, and Applications.} Mahwah, NJ, LEA.
#' @examples
#' #######################################
#' ######### some examples ########
#' data(LienertLSD)
#' LienertLSD
#' CFA(LienertLSD)
#' ## testing with (full) interactions
#' CFA(LienertLSD,form="~ C + T + A + C:T + C:A + T:A + C:T:A")
#' ## testing the null model
#' CFA(LienertLSD,form="null")
#' #######################
#' data(suicide)
#' suicide
#' # suicide data is in non tabulated data representation - so it must be tabulated !
#' CFA(dat2fre(suicide))
NULL
