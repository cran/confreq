#' @exportS3Method summary S2CFA
#' @keywords methods
#' @method summary S2CFA
#' @title S3 Summary for S2CFA
#' @description S3 summary method for object of class\code{"S2CFA"}
#' @param object object of class\code{"S2CFA"}
#' @param digits integer rounds the values to the specified number of decimal places, default is \code{digits=3}.
#' @param type Test column: \code{"ex.fisher.test"} (unweighted default),
#'   \code{"pChi"}, or \code{"p.wald"}. For survey-weighted results, an
#'   omitted type automatically selects \code{"p.wald"}; count-table tests
#'   cannot be requested for these weighted results. Conversely, \code{p.wald}
#'   is unavailable for unweighted results. Unavailable tests produce an
#'   informative error before adjustment or plotting starts.
#' @param sorton sort results of local test by any column. By default the output is not sorted. Other options may be \code{"pat."}, \code{"disc.Type"}, etc. ... So all column names that can potentially appear in the result.
#' @param decreasing logical. Should the sort be increasing or decreasing? see \code{\link{order}}
#' @param showall logical with default \code{showall = TRUE}. To return only significant pattern (discriminating types) set it to \code{showall = FALSE}.
#' @param adjalpha Adjustment method: \code{"bonferroni"} (default),
#'   \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{"BH"},
#'   \code{"BY"}, \code{"fdr"} (alias for BH), or \code{"none"}.
#' @param ... other parameters passed trough.
#' @return A data frame of local results, with printed explanatory output.
#'   Includes a final \code{p.adj} column. For Holm adjustment, includes \code{Holm.crit} and \code{Holm.adj}.
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
#' The final column \code{p.adj} contains the selected adjusted p-values
#' (raw p-values for \code{"none"}); old columns keep their positions.
#' Existing \code{Holm.adj}, \code{Holm.crit} and, where applicable,
#' \code{Holland.adj} remain available. Classification precedes rounding.
#' With \code{adjalpha = "holm"}, decisions use
#' \code{stats::p.adjust(p, method = "holm")} on the non-missing p-values
#' only. Missing values (\code{NA}/\code{NaN}) do not contribute to the number
#' of comparisons. The historical strict rule, adjusted p-value less than
#' \code{alpha}, is retained. Decisions use full precision before rounding.
#' Only finite values in \code{[0, 1]} are computable p-values. Other numeric
#' values are excluded with a warning and treated as missing decisions.
#' The additional column \code{Holm.adj} contains adjusted p-values.
#' The existing \code{Holm.crit} column is retained as a diagnostic value:
#' \code{alpha / (m - r + 1)}, where \code{m} is the number of tested values
#' and \code{r} is their rank (minimum rank for ties). A row-wise comparison
#' with this column alone does not implement the sequential Holm procedure.
#' Missing decisions are displayed as \code{NA} and are omitted when
#' \code{showall = FALSE}. The existing Holm columns retain their positions.
#' @references Holm, S. (1979). A simple sequentially rejective multiple test procedure. \emph{Scandinavian Journal of Statistics, 6}(2), 65–70.
#' @references Bonferroni, C. E. (1935). Il calcolo delle assicurazioni su gruppi di teste. In S.O. Carboni (Ed.), \emph{Studi in Onore del Professore Salvatore Ortu Carboni} (S. 13–60). Roma, Tipografia del Senato: Bardi.
#' @examples
#' data(Lienert1978)
#' res <- S2CFA(Lienert1978)
#' summary(res, adjalpha = "BY")
########################### hier die summary method #class S2CFA #######################
summary.S2CFA<-function(object, digits=3, type="ex.fisher.test", sorton=NULL, decreasing=FALSE, showall=TRUE, adjalpha = "bonferroni", ...){
  #object, digits=3, type="z.pChi",sorton=NULL, decreasing=FALSE, showall=TRUE, holm=FALSE, wide=FALSE, adjalpha = "bonferroni", ...
  local.test <- object$local.test
  # Ein gewichtetes Ergebnis waehlt automatisch den Survey-Wald-Test.
  # Explizite ungeeignete Testwuensche werden nicht still umgedeutet.
  if (missing(type) && isTRUE(object$weighted)) type <- "p.wald"
  # Testauswahl pruefen, bevor eine Adjustierung/Tabellenerstellung beginnt.
  # Ein nicht vorhandener Spaltenname darf keinen numeric(0)-Vektor erzeugen.
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("ex.fisher.test", "pChi", "p.wald"))
    stop("type must be 'ex.fisher.test', 'pChi' or 'p.wald'.", call. = FALSE)
  if (!isTRUE(object$weighted) && type == "p.wald")
    stop("type = 'p.wald' is only available for survey-weighted S2CFA results. For unweighted results use type = 'pChi' or type = 'ex.fisher.test'.", call. = FALSE)
  if (isTRUE(object$weighted) && type != "p.wald")
    stop("For survey-weighted results use type = 'p.wald'; count-table tests are not available.", call. = FALSE)
  if (!type %in% names(object$local.test))
    stop("The requested test column '", type, "' is missing from local.test; recompute S2CFA().", call. = FALSE)
  selected_p <- object$local.test[[type]]
  if (!is.numeric(selected_p) || length(selected_p) != nrow(object$local.test))
    stop("The selected test column must contain one numeric p-value per configuration; recompute S2CFA().", call. = FALSE)
  # Eine vorhandene Spalte voller NA ist zulaessig, z.B. bei alten gewichteten
  # Objekten ohne Varianzmomente. Nicht berechenbar ist nicht gleich fehlend.

  stats_columns <- if (isTRUE(object$weighted))
    match(c("Wald", "df.wald", "p.wald", "prop.diff", "SE.diff"), names(local.test)) else 5:10
  if (anyNA(stats_columns))
    stop("The survey result is incomplete; recompute S2CFA().", call. = FALSE)

  global.test <- object$global.test
  # Methodennamen pruefen; fdr ist der offizielle Alias von BH.
  adjalpha <- match.arg(adjalpha, stats::p.adjust.methods)
  if (adjalpha == "fdr") adjalpha <- "BH"

  
  # object$bonferroni.alpha
  # disc.Type <- ifelse(local.test[,which(names(local.test)==type)] < object$bonferroni.alpha,yes="+", no=".")
  
  #### significant (anti)types cf. Holm ----------------------------------------
  ## added 29-06-2021 adapted from summary.CFA ........
  # Neue Verfahren: derselbe Rechenweg wie bei der korrigierten Holm-Methode.
  if (adjalpha %in% c("hochberg", "hommel", "BH", "BY")) {
    p_values <- local.test[[type]]
    if (!is.numeric(p_values)) stop("Local p-values must be numeric.")
    usable <- is.finite(p_values) & p_values >= 0 & p_values <= 1
    if (any(!is.na(p_values) & !usable))
      warning("Non-finite or out-of-range p-values excluded from adjustment.", call. = FALSE)
    p_adjusted <- rep(NA_real_, length(p_values))
    # Nur berechenbare Werte zaehlen zur Zahl der Vergleiche.
    p_adjusted[usable] <- stats::p.adjust(p_values[usable], method = adjalpha)
    significant <- p_adjusted < object$alpha
    disc.Type <- ifelse(significant, "+", ".")
    templocal <- data.frame(pat. = local.test[, 1L], disc.Type,
      round(local.test[, 2:3], digits = digits),
      round(local.test[,stats_columns,drop=FALSE], digits = digits), check.names = FALSE)
  }

  if(adjalpha=="holm") {
    # Nur berechenbare, tatsaechlich getestete Zellen bilden die Testfamilie.
    # NA/NaN bleiben im Ergebnis NA; funktionale CFA-Zellen werden entfernt.
    # Wichtig: p.adjust() auf den Teilvektor anwenden, damit fehlende Werte
    # nicht durch das voreingestellte n zur Zahl der Vergleiche gezaehlt werden.
    p_val <- local.test[[type]]
    if (!is.numeric(p_val)) stop("Local p-values must be numeric.")
    usable <- is.finite(p_val) & p_val >= 0 & p_val <= 1
    if (any(!is.na(p_val) & !usable))
      warning("Non-finite or out-of-range p-values excluded from Holm adjustment.", call. = FALSE)
    p_holm <- rep(NA_real_, length(p_val))
    p_holm[usable] <- stats::p.adjust(p_val[usable], method = "holm")
    # Die bisherige strikte Vergleichsregel '< alpha' bleibt erhalten.
    # p.adjust beruecksichtigt den notwendigen sequentiellen Holm-Abbruch.
    significant <- p_holm < object$alpha
    # Holm.crit bleibt fuer bestehende Auswerteskripte erhalten.
    # Bei Bindungen gilt fuer alle gleichen p-Werte derselbe kleinste Rang.
    # Dies ist nur ein rangbezogener Diagnosewert: ein isolierter Vergleich
    # p < Holm.crit ersetzt NICHT die sequentielle Holm-Entscheidung!
    alpha_bon <- rep(NA_real_, length(p_val))
    m <- sum(usable)
    if (m > 0L) {
      ranks <- rank(p_val[usable], ties.method = "min")
      alpha_bon[usable] <- object$alpha / (m - ranks + 1L)
    }
    disc.Type <- ifelse(significant, "+", ".")
    templocal <- data.frame(pat. = local.test[, 1L], disc.Type,
      Holm.crit = round(alpha_bon, digits = digits),
      Holm.adj = round(p_holm, digits = digits),
      round(local.test[, 2:3], digits = digits),
      round(local.test[,stats_columns,drop=FALSE], digits = digits), check.names = FALSE)
  }
  
  ### significant (anti)types cf. (classical) Bonferroni adjustment ------------
  if(adjalpha=="bonferroni"){ ## condition added 21-01-2019 # 30-04-2021 new type of alpha adjustment control
    disc.Type <- ifelse(test = local.test[,which(names(local.test)==type)] < object$bonferroni.alpha, yes="+",no=".")
    templocal <- data.frame(pat.=local.test[,1], disc.Type, round(local.test[,2:3],digits=digits), round(local.test[,stats_columns,drop=FALSE],digits=digits) , check.names = FALSE)
  }
  
  ### significant (anti)types without any adjustment for multiple testing ------
  if(adjalpha=="none"){ ## condition 30-04-2021
    disc.Type <- ifelse(test = local.test[,which(names(local.test)==type)] < object$alpha, yes="+",no=".")
    templocal <- data.frame(pat.=local.test[,1], disc.Type, round(local.test[,2:3],digits=digits), round(local.test[,stats_columns,drop=FALSE],digits=digits) , check.names = FALSE)
  }
  
  # Einheitliche Zusatzspalte am ENDE: bisherige Spaltenpositionen erhalten.
  # Holm.adj und Holm.crit bleiben als bestehende Spalten unveraendert verfuegbar.
  if (adjalpha %in% c("hochberg", "hommel", "BH", "BY")) {
    generic_adjusted <- p_adjusted
  } else if (adjalpha == "holm") {
    generic_adjusted <- p_holm
  } else {
    # Bonferroni/none behalten die bisherige Testfamilie und Entscheidung.
    # Insbesondere wird Bonferroni NICHT still auf weniger Vergleiche umgestellt.
    values <- local.test[[type]]
    usable <- is.finite(values) & values >= 0 & values <= 1
    generic_adjusted <- rep(NA_real_, length(values))
    if (adjalpha == "bonferroni") {
      historical_m <- object$alpha / object$bonferroni.alpha
      generic_adjusted[usable] <- pmin(1, values[usable] * historical_m)
    } else {
      generic_adjusted[usable] <- stats::p.adjust(values[usable], method = "none")
    }
  }
  templocal$p.adj <- generic_adjusted
  templocal$p.adj <- round(templocal$p.adj, digits = digits)

  ############## END of alpha adjustment methods -------------------------------
  
  # sort output by sorton ... --------------------------------------------------
  if (length(sorton)!=0 && any(names(templocal)%in%sorton) ){ #cond. changed 8-2-2020
    #   sorton <- "Holm.crit"
    #   sorton <- "p.stir"
    # any(names(templocal)%in%sorton)
    sorter <- order(templocal[,which(names(templocal)==sorton)],decreasing=decreasing)
    erg <- templocal[sorter,]
  }
  if (length(sorton)!=0 && !any(names(templocal)%in%sorton) ){ #cond. changed 8-2-2020
    erg <- templocal # dann wir das sortieren ignoriert wenns der user falsch eingibt -- ohne meldung
  }
  
  if (length(sorton)==0){
    erg <- templocal
    # muted the subsequent code on 23-06-2021
    # if (adjalpha=="holm"){
    #   sorton <- type
    #   sorter <- order(templocal[,which(names(templocal)==sorton)],decreasing=decreasing)
    #   erg <- templocal[sorter,]
    # }#added 'rev' 8-2-2020
  }
  if (showall==FALSE){
    erg <- erg[!is.na(erg$disc.Type) & erg$disc.Type != ".", , drop = FALSE]
  }
  
  
  # cat("results of global tests:","currently not implemented !","\n")
  cat("\n", "Grouping by variable:", object$grouping$var, ", with categories:", object$grouping$cats, "\n", "pattern based on variables:", object$pat,"\n", "-----------------------") 
  
  if (isTRUE(object$weighted)) {
    cat("\nSurvey-weighted frequencies; independent-unit linearized Wald test.\n")
  }
  # adapted from summary.CFA #29-06-2021 
  if(adjalpha=="bonferroni"){
    cat("\n","results of local tests:","\n", "-----------------------","\n","discriminating Type (+) / not discriminating Type (.) based on:", type,";","\n", "with Bonferroni adjusted alpha:", object$bonferroni.alpha,"\n")#changed 30-04-2021
  }
  if(adjalpha=="holm"){
    cat("\n","results of local tests:","\n", "-----------------------","\n","discriminating Type (+) / not discriminating Type (.) based on:", type, ";","\n","with Holm adjusted p-values - see \"Holm.adj\" \n")#added 21-01-2019 #changed 30-04-2021
  }
  if(adjalpha=="none"){
    cat("\n","results of local tests:","\n", "-----------------------","\n","discriminating Type (+) / not discriminating Type (.) based on:", type, ";","\n","with not adjusted alpha:", object$alpha,"\n")#added 30-04-2021
  }
  
  if (adjalpha %in% c("hochberg", "hommel", "BH", "BY")) {
    cat("\nAdjustment method:", adjalpha,
        "; compare p.adj with alpha =", object$alpha, "\n")
  }
  return(erg)
}
