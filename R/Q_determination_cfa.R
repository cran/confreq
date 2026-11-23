#' @title Local Determination Coefficient (Q) for CFA
#' @keywords misc
#' @export Q_determination_cfa
#' @description Calculates Lienert's local determination coefficient (Q) for each cell of a contingency table, based on the observed and expected cell frequencies. 
#' The coefficient quantifies how strongly a single cell deviates from the independence model, relative to the largest deviation that is theoretically possible for that cell.
#' @details The determination coefficient for a cell \eqn{I} is defined (Betzin & Bollmann-Sdorra, 2003, eq. 1; going back to Lienert, 1971) as
#'   \deqn{Q_I = \frac{h_I - e_I}{\max\{e_I,\, N - e_I\}}}
#'   where \eqn{h_I} is the observed and \eqn{e_I} the expected (absolute) frequency of cell \eqn{I}, and \eqn{N} is the total sample size. 
#'   In relative frequencies (\eqn{f_I = h_I/N}, \eqn{p_I = e_I/N}) this is equivalent to
#'   \deqn{Q_I = \frac{f_I - p_I}{\max\{p_I,\, 1-p_I\}}.}
#'   \eqn{Q_I} lies between -1 and 1 for a single cell and coincides with Lienert's (1971) original formula \eqn{Q_I = (h_I-e_I)/e_I} whenever \eqn{2 e_I \geq N}.
#' 
#'   Note that values of \eqn{Q_I} from tables of different dimension are not directly comparable, even after normalizing by the bound -- see \code{\link{Q_determination_bound_cfa}}. 
#'
#' @param observed numeric vector of observed (absolute) cell frequencies \eqn{h_I}
#' @param expected numeric vector of expected (absolute) cell frequencies \eqn{e_I} of the same length as \code{observed}
#' @param n a numeric defining the sample size with default set to \code{sum(observed)}
#'
#' @return a numeric vector with the local determination coefficient \eqn{Q_I} for each cell, in the same order as \code{observed} / \code{expected}.
#'
#' @references Betzin, J., & Bollmann-Sdorra, P. (2003). On the determination coefficient in CFA. \emph{Psychology Science, 45}(2), 400-420.
#' @references Lienert, G. A. (1971). Die Konfigurationsfrequenzanalyse.\emph{Zeitschrift für Klinische Psychologie und Psychotherapie, 19}, 99-115.
#'
#' @seealso \code{\link{Q_determination_bound_cfa}}  
#'
#' @examples
#' #######################################
#' # first calculate expected counts for the LienertLSD data example.
#' # designmatrix with three main effects
#' designmatrix <- design_cfg_cfa(kat = c(2, 2, 2)) 
#' # load example data
#' data(LienertLSD)      
#' # extract observed counts                            
#' observed <- LienertLSD[, 4]                       
#' expected <- expected_cfa(des = designmatrix, observed = observed)
#' Q_determination_cfa(observed, expected)
#' #######################################
#' # an easy, hand-checkable 2x2 example:
#' # margins: row totals e1. = 60 and e2. = 40, 
#' # column totals e.1 = 50 and e.2 = 50, n = 100 -> e11 = 30
#' Q_determination_cfa(observed = 45, expected = 30, n = 100)
#' # (45-30) / max(30, 70) = 0.2142857
###############################################################
############### start of function definition ##################
Q_determination_cfa <- function(observed, expected, n = sum(observed)) { 
  ##### function programmed by Luciana Diaconescu

  if (length(observed) != length(expected)) {
    stop("'observed' and 'expected' must have the same length", call. = TRUE)
  } # die beiden Vektoren expected und observed sollten gleich lang sein

  # bedingung für den nenner (denominator) erstmal prüfen, wir wollen nicht durch 0 teilen
  if (any(expected <= 0) || any(expected >= n)) {  # || bedeutet OR, bei zwei Strichen wird nur ein Element des Vektors ausgespuckt --> eignet sich mehr bei if-loops für control flaw
                                                  # durch any() werden alle Elemente betrachtet --> braucht man wegen ||
    warning("Q_determination_cfa: at least one 'expected' value is <= 0 or >= n; the coefficient is undefined (division by zero) for such cells.", call. = FALSE)
  }

  # Formel (1) aus Betzin & Bollmann-Sdorra
  # neuer vektor 'denominator' mit allen nennern für jede zelle
  denominator <- pmax(expected, n - expected) # pmax() statt max(), weil für jede Zelle und nicht über gesamten vektor
  # hier ist kein loop notwendig, weil die operatoren in R anscheinend vektorisiert sind (habe ich mal selber ausprobiert, mal schauen, ob das für dieses beispiel tatsächlich klappt)
  Q <- abs(observed - expected) / denominator # abs() damit wir nur den Betrag bekommen, also keinen negativen Wert

  return(Q)
}
# End of function ------
