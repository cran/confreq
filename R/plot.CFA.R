#' @exportS3Method plot CFA
#' @keywords methods
#' @method plot CFA
#' @title S3 plot for CFA
#' @description S3 plot method for object of class\code{"CFA"}
#' @param x object of class \code{"CFA"}
#' @param type character indicating which test to use for visualizing whether the observed pattern are 'Types', 'Antitypes' or not significant at all. Possible options for \code{type} are \code{"pChi"}, \code{"ex.bin.test"}, \code{"z.pChi"}, \code{"z.pBin"} and \code{"p.stir"}. 
#' @param fill a vector of (three) colors defining the coloring of significant 'Types' (default "red"), 'Antitypes' (default "blue") or not significant cells (default "grey") in the plot.
#' @param adjalpha Adjustment method: \code{"bonferroni"} (default),
#'   \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{"BH"},
#'   \code{"BY"}, \code{"fdr"} (alias for BH), or \code{"none"}.
#' @param ... other parameters passed trough
#' @return a plot visualizing the results.
#' @details Plot frequencies use the same raw or weighted frequency basis as
#' the CFA result. The existing CFA inferential calculations are unchanged;
#' this plotting correction does not supply a design-based survey variance.
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
#' Functional cells are excluded from the Holm test family and retain
#' their \code{"b"} label. Other pre-existing adjustment methods retain their decisions.
#' @references Holm, S. (1979). A simple sequentially rejective multiple test procedure. \emph{Scandinavian Journal of Statistics, 6}(2), 65–70.
#' @references Bonferroni, C. E. (1935). Il calcolo delle assicurazioni su gruppi di teste. In S.O. Carboni (Ed.), \emph{Studi in Onore del Professore Salvatore Ortu Carboni} (S. 13–60). Roma, Tipografia del Senato: Bardi.
#' @examples
#' data(LienertLSD)
#' res <- CFA(LienertLSD, bintest = FALSE)
#' plot(res, type = "pChi", adjalpha = "BH")

########################### hier die plot method #class CFA #######################
plot.CFA<-function(x, type="z.pChi", fill=c("red", "blue", "grey"), adjalpha = "bonferroni", ...){
  local.test <- x$local.test
  global.test <- x$global.test
  functional <- x$functional
  inputdata <- x$inputdata
  # Methodennamen pruefen; fdr ist der offizielle Alias von BH.
  adjalpha <- match.arg(adjalpha, stats::p.adjust.methods)
  if (adjalpha == "fdr") adjalpha <- "BH"

 
if(!is.null(functional)){   
  if(is.character(functional)){ ## added 22-01-2019 class(functional)=="character" changed 12-04-2022 ' inherits(x = class(functional),what = "character") ' again changed 19-08-2022 
    functional <- sapply(functional, function(y){which(y==as.character(x$local.test$pat.))})
 }
}
  
  
  # check for ex.bin.test present? # test added 23-06-2021
  if(all(is.na(x$local.test$ex.bin.test)) && type == "ex.bin.test"){stop("no results for 'ex.bin.test' available \n try to run CFA() again with argument 'bintest = TRUE' ")}
  
  out_wide <- c("Chi","df","pChi","ex.bin.test","z.Chi","z.pChi","z.Bin","z.pBin","cor.","p.stir","density.stir")

  ############## START of alpha adjustment methods -----------------------------
  

  #### significant (anti)types cf. Holm ----------------------------------------
  ## added 21-01-2019 #changed 30-04-2021 # new from summary 23-08-2022
  # Neue Verfahren: derselbe Rechenweg wie bei der korrigierten Holm-Methode.
  if (adjalpha %in% c("hochberg", "hommel", "BH", "BY")) {
    p_values <- local.test[[type]]
    if (!is.numeric(p_values)) stop("Local p-values must be numeric.")
    if (!is.null(functional)) p_values[functional] <- NA_real_
    usable <- is.finite(p_values) & p_values >= 0 & p_values <= 1
    if (any(!is.na(p_values) & !usable))
      warning("Non-finite or out-of-range p-values excluded from adjustment.", call. = FALSE)
    p_adjusted <- rep(NA_real_, length(p_values))
    # Nur berechenbare Werte zaehlen zur Zahl der Vergleiche.
    p_adjusted[usable] <- stats::p.adjust(p_values[usable], method = adjalpha)
    significant <- p_adjusted < x$alpha
    direction <- ifelse(local.test$obs. > local.test$exp., "+", "-")
    Type <- ifelse(significant, direction, ".")
    if (!is.null(functional)) Type[functional] <- "b"
  }

  if(adjalpha=="holm") {
    # Nur berechenbare, tatsaechlich getestete Zellen bilden die Testfamilie.
    # NA/NaN bleiben im Ergebnis NA; funktionale CFA-Zellen werden entfernt.
    # Wichtig: p.adjust() auf den Teilvektor anwenden, damit fehlende Werte
    # nicht durch das voreingestellte n zur Zahl der Vergleiche gezaehlt werden.
    p_val <- local.test[[type]]
    if (!is.numeric(p_val)) stop("Local p-values must be numeric.")
    if (!is.null(functional)) p_val[functional] <- NA_real_
    usable <- is.finite(p_val) & p_val >= 0 & p_val <= 1
    if (any(!is.na(p_val) & !usable))
      warning("Non-finite or out-of-range p-values excluded from Holm adjustment.", call. = FALSE)
    p_holm <- rep(NA_real_, length(p_val))
    p_holm[usable] <- stats::p.adjust(p_val[usable], method = "holm")
    # Die bisherige strikte Vergleichsregel '< alpha' bleibt erhalten.
    # p.adjust beruecksichtigt den notwendigen sequentiellen Holm-Abbruch.
    significant <- p_holm < x$alpha
    direction <- ifelse(local.test$obs. > local.test$exp., "+", "-")
    Type <- ifelse(significant, direction, ".")
    if (!is.null(functional)) Type[functional] <- "b"
    # Fehlende Entscheidungen bleiben in Type als NA erhalten;
    # der nachfolgende Zeichenteil verwendet dafuer die neutrale Farbe.
  }
  
  ### significant (anti)types cf. (classical) Bonferroni adjustment ------------
  if(adjalpha=="bonferroni"){ ## condition added 21-01-2019 # 30-04-2021 new type of alpha adjustment control
    temp1 <- local.test[,which(names(local.test)==type)] < x$bonferroni.alpha
    temp2 <- ifelse(test=local.test$obs. > local.test$exp., yes="+",no="-") 
    Type <- mapply(FUN=function(xx,y){ifelse(test=(xx==TRUE),yes=y, no="." )   },xx=temp1,y=temp2 )
    if(!is.null(functional)){Type[functional] <- "b"} # added 18-11-2016
    templocal <- data.frame(pat.=local.test[,"pat."],local.test[,c("obs.","exp.")],Type, local.test[,out_wide])
  }
  
  ### significant (anti)types without any adjustment for multiple testing ------
  if(adjalpha=="none"){ ## condition 30-04-2021
    temp1 <- local.test[,which(names(local.test)==type)] < x$alpha
    temp2 <- ifelse(test=local.test$obs. > local.test$exp., yes="+",no="-") 
    Type <- mapply(FUN=function(xx,y){ifelse(test=(xx==TRUE),yes=y, no="." )   },xx=temp1,y=temp2 )
    if(!is.null(functional)){Type[functional] <- "b"} # added 18-11-2016
    templocal <- data.frame(pat.=local.test[,"pat."],local.test[,c("obs.","exp.")],Type, local.test[,out_wide])    
  }
  
  if(any(Type=="b", na.rm = TRUE)){ cat("\n","Type (b): blanked out (functional CFA)","\n")}#(18-11-2016)}
  
  cat(attributes(x)$comment) # version 1.6
  
  inputdata_ <- inputdata

  ############## END of alpha adjustment methods -------------------------------
  
####### eigentlicher plotting teil ###################################
    # Die Grafik muss die tatsaechlich modellierten Haeufigkeiten zeigen.
    # Gewichtsattribute vor/nach der reinen Spaltenumordnung explizit sichern.
    use_weights <- !is.null(attr(inputdata_, "wgtfreq", exact = TRUE))
    plot_input <- inputdata_[, c(rev(colnames(inputdata_))[-1L], "Freq"), drop = FALSE]
    if (use_weights)
      attr(plot_input, "wgtfreq") <- attr(inputdata_, "wgtfreq", exact = TRUE)
    Preq_tab <- fre2tab(plot_input, weighted = use_weights)  
    fil <-sapply(1:length(Type), function(y){
      # Nicht berechenbare Tests werden neutral gezeichnet.
      if (is.na(Type[y])) return(fill[3])
      if(Type[y]=="+"){res <- fill[1]}
      if(Type[y]=="-"){res <- fill[2]}
      if(Type[y]=="."){res <- fill[3]}
      if(Type[y]=="b"){res <- fill[3]}
      res
    })
    # fil
    gp <- gpar(col = NA, lty = "solid", fill = array(fil, dim = dim(Preq_tab),dimnames = dimnames(Preq_tab)))
    strucplot(Preq_tab,gp = gp)
}
