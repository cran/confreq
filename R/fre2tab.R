#' @title pattern frequency to table conversion 
#' @keywords utilities 
#' @export fre2tab
#' @description Given data as pattern frequencies (object of class class c("data.frame","Pfreq", see function \code{\link{dat2fre}}) this function returns a typical array representation (class "table" , see \code{\link{table}}) of it.
#' @details  This function was introduced in order to connect the typical confreq data representation in the objects of the class \code{c("data.frame","Pfreq")}, see function \code{\link{dat2fre}}, to the R-typical array representation as it exists in objects of the \code{"table"} class, see \code{\link{table}}. This array representation of multi-dimensional contingency tables is used more universally in R -- e.g. also in the R package \code{vcd}, see the examples section below.
#' 
#' It is assumed, that the last column of the object \code{patternfreq} represents the frequencies of the (response) pattern represented by the other columns in \code{patternfreq}.
#' 
#' Unused factor levels are retained. Duplicate configurations
#' are summed. Missing configuration categories are retained as table levels.
#' No \code{Pfreq} class assignment is needed for ordinary frequency data.
#' The choice of \code{weighted} affects only this conversion; it does not
#' change the handling of weights by \code{CFA()} or \code{S2CFA()}.
#' 
#' @return An object of class \code{"table"}, with named dimensions in
#'   the input order or the order specified by \code{form}.
#' @seealso \code{\link{dat2fre}}, \code{\link{fre2dat}},
#'   \code{\link[stats]{xtabs}}
#' 
#' @param patternfreq an object of class c("data.frame","Pfreq") 
#' @param form a formula object with possibly both left and right hand sides specifying the order of the variables in the resulting table. At default (\code{(formula=NULL)}) all variables in \code{(x)} are used in their respective order. With two sides,
#'   dimensions on the right precede those on the left, matching the previous
#'   \code{structable()} conversion. A left-hand side of \code{Freq} denotes
#'   frequencies rather than a configuration dimension, e.g. \code{Freq ~ A + B}.
#'   Omitting variables sums over those dimensions. Only variable names joined
#'   by \code{+} are supported; transformations and interactions are rejected.
#' @param weighted Logical. If \code{FALSE} (default), use the final frequency
#'   column. If \code{TRUE}, use the \code{wgtfreq} attribute and report an
#'   error if this attribute is absent or invalid.   
#' @return An object of class \code{"table"}, with named dimensions in
#'   the input order or the order specified by \code{form}; see \code{\link{table}}.
#' @seealso \code{\link{dat2fre}}, \code{\link{fre2dat}},
#'   \code{\link[stats]{xtabs}}
#' @references No references at the moment 
#' @examples #########################################
#' d <- data.frame(A = factor(c("yes", "no", "yes")),
#'                 B = factor(c("low", "low", "high")))
#' pf <- dat2fre(d, wgt = c(1, 2, 3))
#' fre2tab(pf)
#' fre2tab(pf, weighted = TRUE)
#' fre2tab(pf, form = ~ B + A)
#' fre2tab(pf, form = A ~ B)
#' fre2tab(pf, form = Freq ~ A + B)
#' fre2tab(pf, form = ~ A)
#' #########
#' data(LienertLSD)# loading example pattern frequencies table
#' fre2tab(LienertLSD)# coverting it into a table
#' 
#' ### examples using functions from package vcd
#' data(Lienert1978)# loading example pattern frequencies table
#' fre2tab(Lienert1978)# coverting it into a table
#' strucplot(fre2tab(Lienert1978))# plotting data with 'vcd'
#' structable(fre2tab(Lienert1978),direction = "v")# flatten table (vertical) with 'vcd'
#' 
#' # changing the vertical grouping when flattening the table by unsing a 'formula':
#' structable(fre2tab(Lienert1978, form=~Group + Student + Teacher),direction = "v")# flatten table
############### start of function definition ##################
########### pattern frequency to tabe conversion ##############
################ jhheine at googlemail.com ####################
########## new stuff 21-11-2026
fre2tab <- function(patternfreq, form = NULL, weighted = FALSE) {
  if (!is.logical(weighted) || length(weighted) != 1L || is.na(weighted))
    stop("weighted must be TRUE or FALSE.")
  z <- .cf_frequency(patternfreq, weighted = weighted)
  d <- z$data
  vars <- names(d)
  if (!is.null(form)) {
    if (!inherits(form, "formula") || !length(form) %in% c(2L, 3L))
      stop("form must be a formula, e.g. ~ A + B or A ~ B.")
    # Beide Formelseiten werden als Kategorienauswahl gelesen, nicht als
    # Regressionsmodell. Wie bei der bisherigen structable()-Umwandlung
    # kommen die Variablen der rechten Seite VOR denen der linken Seite.
    # Ein alleinstehendes Freq links ist dagegen der historische Sonderfall
    # fuer die Haeufigkeitsspalte und erzeugt keine zusaetzliche Dimension.
    form_names <- function(e) {
      if (is.symbol(e)) return(as.character(e))
      if (is.call(e) && identical(e[[1L]], as.name("+")) && length(e) == 3L)
        return(c(form_names(e[[2L]]), form_names(e[[3L]])))
      stop("form supports variable names joined by '+'.")
    }
    if (length(form) == 2L) {
      vars <- form_names(form[[2L]])
    } else {
      rhs <- form_names(form[[3L]])
      lhs <- if (identical(form[[2L]], as.name("Freq"))) character(0) else
        form_names(form[[2L]])
      vars <- c(rhs, lhs)
    }
    if (anyDuplicated(vars) || any(!vars %in% names(d)))
      stop("form must refer to distinct configuration columns.")
  }
  # Interne Namen vermeiden Formelprobleme bei Namen mit Leerzeichen.
  selected <- d[vars]
  names(selected) <- paste0(".v", seq_along(vars))
  selected$.freq <- z$freq
  f <- stats::reformulate(names(selected)[seq_along(vars)], response = ".freq")
  res <- stats::xtabs(f, data = selected, drop.unused.levels = FALSE,
                      na.action = stats::na.pass,
                      addNA = anyNA(selected[seq_along(vars)]))
  names(dimnames(res)) <- vars
  # Wie bisher eine normale table, ohne xtabs-spezifischen Aufruf.
  attr(res, "call") <- NULL
  class(res) <- "table"
  res
}

#### old stuff -----------------------------------------------------------------
# fre2tab <- function(patternfreq, form = NULL) {
#   if(any(class(patternfreq)=="Pfreq") != TRUE){stop("patternfreq must be an object of class 'Pfreq'","\n","see func. dat2fre()", call. = TRUE) }
#   if (is.null(form)){
#     form <- as.formula(paste("~", paste(colnames(patternfreq)[-ncol(patternfreq)],collapse=" + ")))
#   }
# res <- as.table(structable(formula = form, data = patternfreq))
# 
# return(res)
# }
