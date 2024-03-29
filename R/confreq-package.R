#' @title Configural Frequencies Analysis Using Log-linear Modeling
#' @name confreq-package
#' @docType package
#' @importFrom gmp as.bigq
#' @importFrom gmp as.bigz
#' @importFrom gmp chooseZ
#' @importFrom gmp div.bigq 
#' @importFrom gmp prod.bigq
#' @importFrom methods is
#' @importFrom stats AIC
#' @importFrom stats BIC
#' @importFrom stats aggregate
#' @importFrom stats as.formula
#' @importFrom stats glm.fit
#' @importFrom stats logLik
#' @importFrom stats model.matrix
#' @importFrom stats na.omit 
#' @importFrom stats pchisq
#' @importFrom stats pnorm
#' @importFrom stats poisson
#' @importFrom stats quasipoisson
#' @importFrom stats coef
#' @importFrom vcd structable
#' @importFrom vcd strucplot
#' @import grid
#' @import vcd
#' @description The package \code{confreq} offers some functions for Configural Frequencies Analysis (CFA) proposed by G.A. Lienert as an analysis of types and antitypes of persons or objects grouped according to their characteristic (response) pattern. The core principle in the package \code{confreq} is to use the function \code{\link{glm}} to compute the expected counts based on a model (design) matrix. The main functions are \code{\link{CFA}} and \code{\link{S2CFA}} (see details).
#'
#' @details 
#' The simplest entry to the package \code{confreq} is to use the main function \code{\link{CFA}}, which will compute several coefficients of Configural Frequencies Analysis at once.
#' 
#' More sophisticated control can be achieved by using the several single functions like \code{\link{expected_cfa}}, \code{\link{design_cfg_cfa}}, \code{\link{chi_local_test_cfa}} , \code{\link{stirling_cfa}} , etc. \dots
#' 
#' Two-Sample-CFA, to detect discriminating pattern between two (sub-) samples, can be performed with the function \code{\link{S2CFA}}  
#' 
#' For further description see description of the respective functions.
#' 
#' A good introduction into the theory and applications of Configural Frequencies Analysis is given in the Textbook 'Person-Centered Methods' by Mark Stemmler (see references).
#' 
#' \emph{Additional Information}:
#' Some users running R on 'Linux like' OS distributions (like e.g. Ubuntu -- and in rare cases MAC OS) might report trouble during installation of \code{confreq} due to the package dependency \code{gmp}, which is used in \code{confreq} to perform the exact binomial test. 
#' This (miss-)behavior can usually traced back to a missing of 'the GNU Multiple Precision Arithmetic Library' in the respective OS installation. To fix this, users might consider to run the following Ubuntu Linux command in a terminal to install the latest GMP Library: 
#' 
#' '\code{sudo apt-get install libgmp3-dev}'.
#' 
#' For additional information see also \href{http://www.mathemagix.org/www/mmdoc/doc/html/external/gmp.en.html}{http://www.mathemagix.org/www/mmdoc/doc/html/external/gmp.en.html} and \href{https://gmplib.org/}{https://gmplib.org/}.
#' 
#' 
#' Annotation:
#' The foundations for this R-Package were established and discussed in Rothenberge (2011) and (finally) in Klagenfurt at FGME 2013 with Rainer Alexandrowicz and Mark Stemmler \dots   
#'
#' @author \itemize{\item{Joerg-Henrik Heine <jhheine@@googlemail.com>}\item{R.W. Alexandrowicz (function \code{stirling_cfa()})}}
#' @references von Eye, A. (2002). \emph{Configural Frequency Analysis. Methods, Models, and Applications.} Mahwah, NJ, LEA.
#' @references Krauth, J., & Lienert, G. A. (1973). \emph{Die Konfigurationsfrequenzanalyse (KFA) und ihre Anwendung in Psychologie und Medizin: ein multivariates nichtparametrisches Verfahren zur Aufdeckung von Typen und Syndromen; mit 70 Tabellen}. Freiburg; München: Alber Karl.
#' @references Lazarsfeld, P. F., & Henry, N. W. (1968). \emph{Latent structure analysis}. Boston: Houghton Mifflin.
#' @references Lienert, G. A. (1978). \emph{Verteilungsfreie Methoden in der Biostatistik (Band II)} [Non-parametrical methods in the field of biometrics (Vol. II)]. Meisenheim am Glan, Germany: Hain.
#' @references Lienert, G. A. (1971). Die Konfigurationsfrequenzanalyse: I. Ein neuer Weg zu Typen und Syndromen. \emph{Zeitschrift für Klinische Psychologie und Psychotherapie, 19}(2), 99-115.
#' @references Stemmler, M. (2020). \emph{Person-Centered Methods -- Configural Frequency Analysis (CFA) and Other Methods for the Analysis of Contingency Tables}. Cham Heidelberg New York Dordrecht London: Springer.
#' @references Stemmler, M., & Hammond, S. (1997). Configural frequency analysis of dependent samples for intra-patient treatment comparisons. \emph{Studia Psychologica, 39}, 167–175.
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
#' # suicide data is in non tabulated data representation 
#' # so it must be tabulated !
#' CFA(dat2fre(suicide))
NULL
