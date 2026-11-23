#' @title Theoretical Upper Bound for the Determination Coefficient
#' @keywords misc
#' @export Q_determination_bound_cfa
#' @description Calculates the theoretical upper bound of the local determination coefficient Q for a table with \code{T} variables (dimensions), following Theorem 1 in Betzin & Bollmann-Sdorra (2003).
#' @details For a \eqn{T}-dimensional table (i.e. \eqn{2^T} cells) the empirical value of the determination coefficient \code{\link{Q_determination_cfa}} is bounded by \deqn{Q_I < 2^{(T-1)/T} - 1} 
#' and this bound is sharp (attained in the limit, for any \eqn{\epsilon > 0} a table configuration exists that gets arbitrarily close to it). 
#' The bound increases with \code{T} but stays clearly below 1 even for moderately large \code{T} (e.g. ~0.682 for T = 4, ~0.866 for T = 10, ~0.986 for T = 100).
#' The authors explicitly warn against using this bound to *normalize* Q (see Section 4 / Example 1 (p. 409) of the paper, where normalizing even reverses the ranking of two Q values from tables of different dimension). 
#' Comparisons of Q across tables of different dimension therefore remain problematic even after normalization. 
#' 
#' @param D a numeric vector giving the number of variables (dimensions) of the table(s), \code{T >= 1}.
#' @return a numeric vector (same length as \code{T}) with the upper bound for Q.
#' 
#' @references Betzin, J., & Bollmann-Sdorra, P. (2003). On the determination coefficient in CFA. \emph{Psychology Science, 45}(2), 400-420.
#' 
#' @seealso \code{\link{Q_determination_cfa}}
#' 
#' @examples
#' #######################################
#' # reproduces Table 1 of Betzin & Bollmann-Sdorra (2003)
#' T <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 50, 75, 100)
#' round(Q_determination_bound_cfa(D = T), 3)
###############################################################
############### start of function definition ##################
Q_determination_bound_cfa <- function(D){
  ##### function programmed by Luciana Diaconescu

  # T als Notation genutzt, weil Betzin & Bollmann-Sdorra das so gemacht haben --> T entspricht der Anzahl an Dimensionen: changed JHH 12-11-2026 argument name is now D
  if (any(D < 1)) {
    stop("'D' (number of dimensions) must be >= 1", call. = TRUE) # 0 dimensionen oder negative ergeben keinen sinn
  }

  2^((D - 1) / D) - 1 # obere Grenze von Q für D (T) Dimensionen, nicht normalisiert, weil Autoren das nicht empfehlen!!
}
# End of function ------
