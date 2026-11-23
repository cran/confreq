# Neue interne Funktionen fuer confreq 1.7.0-5 -- Base R / stats.
# Stand: 21.11.2026. 

# Interne Hilfsfunktion: gemeinsame Eingabepruefung fuer Individualdaten.
.cf_data <- function(x) {
  if (!is.data.frame(x) && !is.matrix(x))
    stop("x must be a data.frame or matrix.")
  d <- as.data.frame(x, stringsAsFactors = FALSE)
  if (!ncol(d)) stop("At least one configuration variable is required.")
  if (anyDuplicated(names(d)) || anyNA(names(d)) || any(!nzchar(names(d))))
    stop("Column names must be nonempty and unique.")
  if (any(!vapply(d, function(z) is.factor(z) || is.character(z) ||
                 is.numeric(z) || is.logical(z), logical(1))))
    stop("Columns must be factors, character, numeric or logical vectors.")
  d
}

.cf_weights <- function(wgt, n) {
  if (!is.null(wgt) && (!is.numeric(wgt) || length(wgt) != n ||
      anyNA(wgt) || any(!is.finite(wgt)) || any(wgt < 0)))
    stop("wgt must contain one finite, nonnegative weight per row.")
}

# Gemeinsame Pruefung einer bereits tabellierten Eingabe.
# Es bleibt kompatibel: Die letzte Spalte enthaelt die Haeufigkeit.
.cf_frequency <- function(x, weighted = FALSE, integer = FALSE) {
  d <- .cf_data(x)
  if (ncol(d) < 2L) stop("Need configuration columns and a final frequency column.")
  f <- if (weighted) attr(x, "wgtfreq") else d[[ncol(d)]]
  if (is.null(f)) stop("No weighted frequencies available.")
  if (!is.numeric(f) || length(f) != nrow(d) || anyNA(f) ||
      any(!is.finite(f)) || any(f < 0))
    stop("Frequencies must be finite and nonnegative numeric values.")
  if (integer && any(f != floor(f)))
    stop("Individual records require integer frequencies; weighted sums cannot be rounded automatically.")
  list(data = d[seq_len(ncol(d) - 1L)], freq = f)
}

