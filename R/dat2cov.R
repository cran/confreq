#' @title  Summarize Covariates by Configuration
#' @keywords utilities 
#' @export dat2cov
#' @description Compute one numeric summary per covariate and configuration, in the
#' order of the corresponding complete pattern-frequency table.
#' Given a dataset \code{x}, this function returns summary values for some (\code{numeric}) covariate variables in \code{x} for each pattern (configuration) defined by a set of factor variables in \code{x}. 
#' @details At least one configuration variable and one numeric covariate
#' are required. Explicit selection allows numeric or character configuration
#' variables without first converting them to factors. Factor levels can
#' specify unobserved categories. Configuration and covariate selections
#' must be disjoint and contain distinct existing column names.
#'
#' \code{notobs} applies only when no input record belongs to a configuration.
#' Missing or undefined summaries in observed groups remain \code{NA} or
#' \code{NaN}; their handling is controlled by \code{FUN} and its arguments.
#' A group containing records with zero weights is still observed, even if
#' its weighted mean is undefined. Missing configuration values are excluded
#' according to the default \code{dat2fre()} behavior.
#'
#' Names passed to \code{variables} and \code{covariates} must be written in
#' full because these arguments follow \code{...}. For use with \code{CFA()},
#' use the same configuration columns and ordering options in \code{dat2fre()}.
#'
#' 
#' @param x A data frame or matrix containing configuration variables
#'   and numeric covariates.
#' @param FUN A function or function name returning one numeric value per
#'   group and covariate. Defaults to \code{"mean"}. When \code{wgt} is
#'   supplied, only \code{mean} is supported and weighted means are computed.
#' @param ... Additional arguments passed to \code{FUN}, or to
#'   \code{\link[stats]{weighted.mean}} for weighted means, e.g.
#'   \code{na.rm = TRUE}. Not forwarded to \code{dat2fre()}.
#' @param notobs Numeric replacement for unobserved configurations.
#'   Defaults to zero. Supply one value for all covariates or one value per
#'   selected covariate, in the order of \code{covariates}. \code{NA} is allowed.
#' @param katorder Logical; see \code{\link{dat2fre}}. Use the same setting
#'   for the corresponding pattern-frequency table.
#' @param caseorder Logical; see \code{\link{dat2fre}}. Use the same setting
#'   for the corresponding pattern-frequency table.
#' @param wgt Optional vector of finite, nonnegative case weights, one per
#'   input row. Zero weights are allowed.
#' @param variables Optional character vector of configuration column names.
#'   By default, all factor columns of \code{x} are selected.
#' @param covariates Optional character vector of numeric covariate column
#'   names. By default, all numeric columns not selected in \code{variables}
#'   are used. Explicitly exclude identifiers and weight columns if necessary.   
#' @return A data frame of class \code{c("data.frame", "Pcov")}, containing
#'   the covariate summaries. The \code{configurations} attribute stores the
#'   corresponding configuration columns. Row names are readable labels made
#'   unique if necessary, not keys used for internal matching.
#' @seealso \code{\link{dat2fre}}, \code{\link{CFA}}
#' @examples
#' d <- data.frame(A = factor(c("yes", "yes", "no")),
#'                 B = factor(c("low", "low", "high")),
#'                 Age = c(20, 40, 60), Weight = c(1, 3, 2))
#' cov <- dat2cov(d, variables = c("A", "B"), covariates = "Age",
#'                wgt = d$Weight, notobs = NA_real_)
#' cov
#' attr(cov, "configurations")
############### start of function definition ##################
################ jhheine at googlemail.com ####################
# dat2cov: Explizite Spaltenauswahl ist optional. Ohne Auswahl gilt wie bisher:
# Faktoren definieren Konfigurationen, numerische Spalten die Kovariaten.
# Neue Auswahlargumente nehmen Spaltennamen entgegen.
dat2cov <- function(x, FUN = "mean", ..., notobs = 0, katorder = FALSE,
                    caseorder = TRUE, wgt = NULL, variables = NULL,
                    covariates = NULL) {
  d <- .cf_data(x)
  .cf_weights(wgt, nrow(d))
  if (is.null(variables)) variables <- names(d)[vapply(d, is.factor, logical(1))]
  if (is.null(covariates))
    covariates <- setdiff(names(d)[vapply(d, is.numeric, logical(1))], variables)
  valid_selection <- function(s) is.character(s) && length(s) > 0L &&
    !anyNA(s) && !anyDuplicated(s) && all(s %in% names(d))
  if (!valid_selection(variables) || !valid_selection(covariates) ||
      length(intersect(variables, covariates)))
    stop("Select distinct existing configuration and covariate columns by name.")
  if (!all(vapply(d[covariates], is.numeric, logical(1))))
    stop("Covariates must be numeric.")
  if (!is.numeric(notobs) || !length(notobs) ||
      !length(notobs) %in% c(1L, length(covariates)))
    stop("notobs must be numeric: one value or one value per covariate (NA allowed).")
  notobs <- rep(notobs, length.out = length(covariates))
  fun <- match.fun(FUN)
  if (!is.null(wgt) && !identical(fun, base::mean))
    stop("With wgt, only FUN = mean is supported; other functions need an explicit weighting definition.")
  pf <- dat2fre(d[variables], katorder = katorder, caseorder = caseorder)
  pattern <- pf[seq_len(ncol(pf) - 1L)]
  
  # Keine zusammengesetzten Textschluessel: pro Variable Kategorien vergleichen.
  # Fehlende Konfigurationswerte werden wie im dat2fre()-Standard ausgeschlossen.
  group <- rep(NA_integer_, nrow(d))
  for (i in seq_len(nrow(pattern))) {
    hit <- rep(TRUE, nrow(d))
    for (nm in names(pattern)) {
      a <- as.character(d[[nm]])
      hit <- hit & !is.na(a) & a == as.character(pattern[[nm]][i])
    }
    group[which(hit)] <- i
  }
  out <- as.data.frame(matrix(NA_real_, nrow(pattern), length(covariates)))
  names(out) <- covariates
  for (i in seq_len(nrow(pattern))) {
    take <- which(group == i)
    for (j in seq_along(covariates)) {
      if (!length(take)) {
        # notobs greift ausschliesslich bei unbeobachteten Konfigurationen.
        out[i, j] <- notobs[j]
      } else {
        values <- d[[covariates[j]]][take]
        value <- if (is.null(wgt)) fun(values, ...) else
          stats::weighted.mean(values, w = wgt[take], ...)
        if (!is.numeric(value) || length(value) != 1L)
          stop("FUN must return one numeric value per configuration and covariate.")
        out[i, j] <- value
      }
    }
  }
  # Lesbare Zeilennamen fuer bestehende Anwendungen; eindeutige Zuordnung
  # steht zusaetzlich als unverkuerzte Tabelle im Attribut configurations.
  rownames(out) <- make.unique(do.call(paste, c(lapply(pattern, as.character), sep = " ")))
  attr(out, "configurations") <- as.data.frame(pattern)
  class(out) <- c("data.frame", "Pcov")
  out
}





## old Stuff
# dat2cov <- function(x, FUN = "mean", ... , notobs = 0, katorder = FALSE, caseorder = TRUE, wgt=NULL ){
#   d <- x
#   #p_ind <- (sapply(d,class))=="factor"
#   p_ind <- (sapply(X = d,FUN = function(x){is.factor(x)}    ))
#   #cov_ind <- (sapply(d,class))=="numeric"
#   cov_ind <- (sapply(X = d,FUN = function(x){is.numeric(x)}    ))
#   nam <- names(d)[cov_ind]
#   
#   ta <- dat2fre(d[,p_ind],katorder = katorder,caseorder = caseorder,wgt = wgt)
#   if(is.null(wgt)){ 
#     fre <- ta[,ncol(ta)]
#   }
#   if(!is.null(wgt)){
#     fre <- attributes(ta)$WGT
#   }
#   pat <- apply(X =ta[,1:(ncol(ta)-1)],1,paste, collapse = " ")     
#   
#   dp <-  data.frame(pat = apply(X = d[,p_ind],1,paste, collapse = " "), d[,cov_ind])
#   da <- aggregate(x = dp[,-1], by=list(dp[,"pat"]),FUN = FUN, ...)
#   # db <- data.frame(pat=levels(da[,1]),da[,-1],stringsAsFactors =FALSE) # new AKA2021
#   db <- data.frame(pat=da[,1],da[,-1],stringsAsFactors =FALSE) 
#   names(db)[-1] <- nam
#   row.names(db) <-db[,1] 
#   # matrix(db[pat,-1],ncol=sum(cov_ind))
#   dc <- data.frame(db[pat,-1])
#   names(dc) <- nam
#   row.names(dc) <- pat
#   
#   if(!is.na(notobs)){ 
#     if(length(notobs)==1){notobs <- rep(notobs,times=sum(cov_ind))}
#     for(i in 1:sum(cov_ind)){
#       dc[,i][is.na(dc[,i])] <- notobs[i]  
#     }
#     }
#   
#   class(dc) <- c("data.frame", "Pcov")
#   return(dc)
# }