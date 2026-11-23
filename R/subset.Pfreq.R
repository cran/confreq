#' @title Subset Pattern-Frequency Tables without Misaligning Weights
#' @keywords utilities
#' @description Keep cell-level weight attributes aligned when selecting or reordering rows of a \code{Pfreq} object.
#'
#' @param x A \code{Pfreq} object.
#' @param i,j Row and column selections, as for \code{[.data.frame}.
#' @param drop Passed to the data-frame method if supplied. Otherwise its
#'   ordinary default applies.
#'
#' @details The row attributes \code{wgtfreq}, \code{WGT}, \code{wgtfreq2}
#' and \code{wgtpositive} follow the selected rows, including named rows.
#' A row permutation retains the original individual-weight provenance
#' attribute \code{wgtdat}; selecting fewer configurations removes it because
#' the remaining individual weights cannot be reconstructed from cell totals.
#'
#' The \code{Pfreq} class is retained only for nonempty selections of unique,
#' existing rows that keep all configuration columns exactly once and the
#' original frequency column last. Other data-frame results lose the Pfreq
#' class and weight attributes. No automatic marginalization is performed
#' when configuration columns are removed. Vector results follow the ordinary
#' data-frame dropping rules. This method does not protect arbitrary editing
#' by \code{[<-}, \code{$<-}, \code{[[<-} or manual attribute replacement.
#'
#' An incomplete selection of unique configurations remains a frequency
#' table, but is not necessarily suitable for CFA, which needs a complete
#' grid. The corresponding analysis input checks should be respected.
#'
#' @return The selected object, with aligned weight attributes if it remains
#'   a Pfreq object; otherwise an ordinary data frame or dropped vector.
#' @seealso \code{\link{dat2fre}}, \code{\link{fre2tab}}
#' @method [ Pfreq
#' @export
#' @examples
#' d <- data.frame(A = factor(c("a", "b", "a")),
#'                 B = factor(c("x", "x", "y")))
#' pf <- dat2fre(d, wgt = c(1, 2, 3))
#' sorted <- pf[rev(seq_len(nrow(pf))), ]
#' attr(sorted, "wgtfreq")
#' 
################################################################################

`[.Pfreq` <- function(x, i, j, drop) {
  # Die Basismethode entscheidet unveraendert ueber Spaltenauswahl, drop,
  # Zeilennamen, negative/logische Indizes und Fehler bei ungueltigen Spalten.
  result <- NextMethod("[")
  if (!is.data.frame(result)) return(result)

  # x[i] ist Spaltenauswahl; x[i, ] ist Zeilenauswahl. Ein explizites drop
  # darf fuer diese Unterscheidung nicht als zusaetzlicher Index zaehlen.
  index_arguments <- nargs() - as.integer(!missing(drop))
  row_index <- seq_len(nrow(x))
  if (index_arguments >= 3L && !missing(i)) {
    # Eine kleine Index-Tabelle uebernimmt exakt die Zeilenauswahlregeln
    # von data.frame, auch bei Zeichenindizes und mehrfachen Zeilennamen.
    marker <- data.frame(.row = row_index, row.names = row.names(x))
    row_index <- marker[i, , drop = FALSE][[1L]]
  }
  weight_names <- c("wgtfreq", "WGT", "wgtfreq2", "wgtpositive")
  weights_present <- vapply(weight_names, function(nm)
    !is.null(attr(x, nm, exact = TRUE)), logical(1))

  # Ungueltige, zuvor bereits beschaedigte Attribute nicht still recyceln.
  for (nm in weight_names[weights_present]) {
    value <- attr(x, nm, exact = TRUE)
    if (!is.numeric(value) || length(value) != nrow(x))
      stop("Misaligned Pfreq weight attribute '", nm, "'; recreate the table with dat2fre().")
  }
  valid_rows <- length(row_index) > 0L && !anyNA(row_index) && !anyDuplicated(row_index)
  valid_columns <- ncol(result) == ncol(x) && !anyDuplicated(names(result)) &&
    setequal(names(result), names(x)) &&
    identical(tail(names(result), 1L), tail(names(x), 1L))
  if (valid_rows && valid_columns) {
    class(result) <- class(x)
    for (nm in weight_names[weights_present])
      attr(result, nm) <- attr(x, nm, exact = TRUE)[row_index]
    # wgtdat enthaelt Individual-, keine Zellgewichte: niemals mit Zellindizes
    # umsortieren. Nur bei unveraenderter Auswahl aller Zellen weitergeben.
    if (length(row_index) == nrow(x) && setequal(row_index, seq_len(nrow(x)))) {
      attr(result, "wgtdat") <- attr(x, "wgtdat", exact = TRUE)
    } else {
      attr(result, "wgtdat") <- NULL
    }
    comment(result) <- comment(x)
  } else {
    # Ein Ausschnitt ohne alle Merkmale ist keine automatisch aggregierte
    # Konfigurationstabelle. Ohne Pfreq-Klasse koennen CFA/S2CFA ihn nicht
    # versehentlich als vollwertige Eingabe interpretieren.
    class(result) <- setdiff(class(result), "Pfreq")
    for (nm in c(weight_names, "wgtdat")) attr(result, nm) <- NULL
    if (identical(comment(x), "using weighted data!")) comment(result) <- NULL
  }
  result
}
