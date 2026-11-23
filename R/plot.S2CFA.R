#' @exportS3Method plot S2CFA
#' @keywords methods
#' @method plot S2CFA
#' @title S3 plot for S2CFA
#' @description S3 plot method for object of class\code{"S2CFA"}
#' @param x object of class\code{"S2CFA"}
#' @param type Test column: \code{"ex.fisher.test"} (unweighted default),
#'   \code{"pChi"}, or \code{"p.wald"}. For survey-weighted results, an
#'   omitted type automatically selects \code{"p.wald"}; count-table tests
#'   cannot be requested for these weighted results. Conversely, \code{p.wald}
#'   is unavailable for unweighted results. Unavailable tests produce an
#'   informative error before adjustment or plotting starts.
#' @param fill a vector of (two) colors defining the coloring of discriminating 'Types' (default "red"), or not discriminating cells (default "grey") in the plot.
#' @param adjalpha Adjustment method: \code{"bonferroni"} (default),
#'   \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{"BH"},
#'   \code{"BY"}, \code{"fdr"} (alias for BH), or \code{"none"}.
#' @param ... other parameters passed trough.
#' @return Draws a plot and invisibly returns a list containing the table,
#'   fill array, significance decisions, configuration mapping, and
#'   additional entries \code{p.adj} and \code{adjalpha}.
#' @details Survey-weighted frequencies and Wald p-values are taken from
#' the weighted \code{S2CFA()} result; no integer rounding is performed.
#' The test assumes independent sampling units, without cluster or stratum
#' information. Non-computable p-values do not contribute to the Holm and
#' additional p.adjust test families. Summary and plot select the same test.
#' Additional adjustment methods are \code{"hochberg"},
#' \code{"hommel"}, \code{"BH"}, \code{"BY"} and \code{"fdr"} (an alias
#' for \code{"BH"}). These methods use \code{\link[stats]{p.adjust}} on
#' computable p-values only, following the same test-family rule as Holm.
#' Hochberg and Hommel require appropriate dependence assumptions; BH and BY
#' control the false discovery rate rather than the family-wise error rate.
#' See \code{p.adjust} for assumptions and statistical interpretation.
#' The historical Bonferroni denominator and the \code{"none"} behavior are
#' retained. In particular, the reduced test-family rule for Holm and the new
#' methods is not silently applied to the existing Bonferroni decisions.
#' With \code{adjalpha = "holm"}, decisions use
#' \code{stats::p.adjust(p, method = "holm")} on the non-missing p-values
#' only. Missing values (\code{NA}/\code{NaN}) do not contribute to the number
#' of comparisons. The historical strict rule, adjusted p-value less than
#' \code{alpha}, is retained. Decisions use full precision before rounding.
#' Only finite values in \code{[0, 1]} are computable p-values. Other numeric
#' values are excluded with a warning and treated as missing decisions.
#' Missing decisions are drawn using the neutral fill color.
#' @references Holm, S. (1979). A simple sequentially rejective multiple test procedure. \emph{Scandinavian Journal of Statistics, 6}(2), 65–70.
#' @references Bonferroni, C. E. (1935). Il calcolo delle assicurazioni su gruppi di teste. In S.O. Carboni (Ed.), \emph{Studi in Onore del Professore Salvatore Ortu Carboni} (S. 13–60). Roma, Tipografia del Senato: Bardi.
#' @examples
#' data(Lienert1978)
#' res <- S2CFA(Lienert1978)
#' plot(res, adjalpha = "BY")

########################### hier die plot method #class S2CFA #######################
#### JHH: komplett neu 21-11-2026 wegen kompletter Fehlerhaftigkeit der alten Vers. 
plot.S2CFA <- function(x, type = "ex.fisher.test", fill = c("red", "grey"),
                       adjalpha = "bonferroni", ...) {
  # 1. Einstellungen pruefen. Nicht berechenbare p-Werte werden unten
  #    neutral dargestellt, statt ungueltige Farben zu erzeugen.
  if (!inherits(x, "S2CFA")) stop("x must be an S2CFA result.")
  if (missing(type) && isTRUE(x$weighted)) type <- "p.wald"
  type <- match.arg(type, c("ex.fisher.test", "pChi", "p.wald"))
  # Testauswahl pruefen, bevor eine Adjustierung/Tabellenerstellung beginnt.
  # Ein nicht vorhandener Spaltenname darf keinen numeric(0)-Vektor erzeugen.
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("ex.fisher.test", "pChi", "p.wald"))
    stop("type must be 'ex.fisher.test', 'pChi' or 'p.wald'.", call. = FALSE)
  if (!isTRUE(x$weighted) && type == "p.wald")
    stop("type = 'p.wald' is only available for survey-weighted S2CFA results. For unweighted results use type = 'pChi' or type = 'ex.fisher.test'.", call. = FALSE)
  if (isTRUE(x$weighted) && type != "p.wald")
    stop("For survey-weighted results use type = 'p.wald'; count-table tests are not available.", call. = FALSE)
  if (!type %in% names(x$local.test))
    stop("The requested test column '", type, "' is missing from local.test; recompute S2CFA().", call. = FALSE)
  selected_p <- x$local.test[[type]]
  if (!is.numeric(selected_p) || length(selected_p) != nrow(x$local.test))
    stop("The selected test column must contain one numeric p-value per configuration; recompute S2CFA().", call. = FALSE)
  # Eine vorhandene Spalte voller NA ist zulaessig, z.B. bei alten gewichteten
  # Objekten ohne Varianzmomente. Nicht berechenbar ist nicht gleich fehlend.

  # Methodennamen pruefen; fdr ist der offizielle Alias von BH.
  adjalpha <- match.arg(adjalpha, stats::p.adjust.methods)
  if (adjalpha == "fdr") adjalpha <- "BH"
  if (length(fill) != 2L || anyNA(fill)) stop("fill requires two valid colors.")
  grDevices::col2rgb(fill) # Ungueltige Farbnamen bereits vor dem Plot melden.
  p <- x$local.test[[type]]
  extended <- adjalpha %in% c("holm", "hochberg", "hommel", "BH", "BY")
  if (!is.numeric(p) || (!extended && any(!is.na(p) & (!is.finite(p) | p < 0 | p > 1))))
    stop("Invalid local p-values.")

  # Neue Verfahren und Holm: nur berechenbare Werte in die Testfamilie.
  significant <- rep(FALSE, length(p))
  p_adjusted <- rep(NA_real_, length(p))
  usable <- !is.na(p)
  if (extended) {
    usable <- is.finite(p) & p >= 0 & p <= 1
    if (any(!is.na(p) & !usable))
      warning("Non-finite or out-of-range p-values excluded from adjustment.", call. = FALSE)
    p_adjusted[usable] <- stats::p.adjust(p[usable], method = adjalpha)
    significant[usable] <- p_adjusted[usable] < x$alpha
  }
  if (adjalpha == "bonferroni") {
    significant[usable] <- p[usable] < x$bonferroni.alpha
    p_adjusted[usable] <- pmin(1, p[usable] * x$alpha / x$bonferroni.alpha)
  }
  if (adjalpha == "none") {
    p_adjusted[usable] <- stats::p.adjust(p[usable], method = "none")
    significant[usable] <- p[usable] < x$alpha
  }

  input <- x$inputdata
  group <- x$grouping$var
  predictors <- x$pat
  if (length(group) != 1L || !group %in% names(input) || !length(predictors) ||
      any(!predictors %in% names(input)))
    stop("Missing configuration or grouping columns in S2CFA result.")
  
  # 3. WICHTIGE FEHLERKORREKTUR:
  # doubledecker() schiebt die abhaengige Variable an die LETZTE Dimension.
  # Ein bereits mehrdimensionales gp$fill wird dabei NICHT mit umgeordnet.
  # Deshalb stellen wir die Gruppierung selbst ans Ende, bevor Tabelle und
  # Farbfeld entstehen. doubledecker muss dann nichts mehr permutieren.
  # Die umgekehrte Reihenfolge der Praediktoren entspricht der bisherigen
  # Darstellung nach der internen Umordnung durch doubledecker().
  table_vars <- c(rev(predictors), group)
  # Grafik und Analyse muessen dieselbe Basis verwenden. data.frame-Subsetting
  # kann eigene Attribute entfernen; deshalb die unveraenderten Zeilengewichte
  # nach der reinen Spaltenumordnung ausdruecklich wieder anhaengen.
  plot_input <- input[, c(table_vars, names(input)[ncol(input)]), drop = FALSE]
  if (isTRUE(x$weighted))
    attr(plot_input, "wgtfreq") <- attr(input, "wgtfreq", exact = TRUE)
  tab <- fre2tab(plot_input, weighted = isTRUE(x$weighted))
  dims <- dim(tab)
  dn <- dimnames(tab)
  if (tail(dims, 1L) != 2L) stop("The grouping variable must have two categories.")
  
  # 4. Die Zeilen der lokalen Tests entsprechen in S2CFA den Eingabezeilen
  # der ersten Gruppe. Diese Konfigurationen liefern die Zuordnung zu p.
  # Es werden keine potentiell mehrdeutigen Textschluessel verkettet.
  first_group <- as.character(input[[group]]) == as.character(x$grouping$cats[1L])
  configurations <- input[which(first_group), rev(predictors), drop = FALSE]
  if (nrow(configurations) != length(p) || anyDuplicated(configurations))
    stop("Local tests cannot be matched uniquely to configurations.")
  
  # Gemischte Stellenwerte berechnen die Position einer Zelle im R-Array:
  # erste Dimension wechselt am schnellsten, letzte Dimension am langsamsten.
  stride <- c(1, cumprod(dims)[-length(dims)])
  cell <- rep(1, nrow(configurations))
  for (j in seq_along(configurations)) {
    position <- match(as.character(configurations[[j]]), dn[[j]])
    if (anyNA(position)) stop("Configuration category missing from plot table.")
    cell <- cell + (position - 1) * stride[j]
  }
  if (length(cell) != prod(head(dims, -1L)))
    stop("The plot requires a complete configuration table.")
  
  # Jede Konfiguration bekommt in beiden Gruppen dieselbe Farbe.
  # Explizite Zellzuordnung verhindert Fehler bei geaenderter Faktorordnung.
  colors <- array(fill[2L], dim = dims, dimnames = dn)
  configuration_colors <- ifelse(significant, fill[1L], fill[2L])
  for (g in seq_len(tail(dims, 1L)))
    colors[cell + (g - 1L) * tail(stride, 1L)] <- configuration_colors
  
  # 5. Farben und Tabelle haben dieselbe Dimensionsreihenfolge.
  # ... wird jetzt tatsaechlich weitergegeben, z.B. main = "Mein Titel".
  # Beide Gruppen haben je Konfiguration dieselbe Signifikanzfarbe.
  # Eine weisse Trennlinie macht ihre unterschiedlichen Anteile sichtbar;
  # mit col=NA wuerden die gleichfarbigen Gruppen optisch verschmelzen.
  gp <- grid::gpar(col = "white", lty = "solid", fill = colors)
  vcd::doubledecker(tab, depvar = length(dims), gp = gp, ...)
  cat("\nGrouping by variable:", group, ", with categories:", x$grouping$cats,
      "\npattern based on variables:", predictors, "\n-----------------------\n")
  
  # Unsichtbare Rueckgabe ermoeglicht die Pruefung der tatsaechlichen
  # Farbzuordnung, ohne die normale Anwendung plot(res) zu veraendern.
  invisible(list(table = tab, fill = colors, significant = significant,
                 configurations = configurations,
                 p.adj = p_adjusted, adjalpha = adjalpha))
}
