#' @exportS3Method summary CFA
#' @keywords methods
#' @method summary CFA
#' @title S3 Summary for CFA
#' @description S3 summary method for object of class\code{"CFA"}
#' @param object object of class\code{"CFA"}
#' @param digits integer rounds the values to the specified number of decimal places, default is \code{digits=3}.
#' @param type character indicating which test to use for inference whether the observed pattern are 'Types', 'Antitypes' or not significant at all. Possible options for \code{type} are \code{"pChi"}, \code{"ex.bin.test"}, \code{"z.pChi"}, \code{"z.pBin"} and \code{"p.stir"}. 
#' @param sorton sort results of local test by any column. By default the output is not sorted. Other options may be \code{"pat."}, \code{"obs."}, \code{"exp."}, \code{"Type"}, \code{"Chi"}, etc. ... So all column names that can potentially appear in the result.
#' @param decreasing logical. Should the sort be increasing or decreasing? see \code{\link{order}}
#' @param showall logical with default \code{showall = TRUE}. To return only significant pattern ('Types' / 'Antitypes') set it to \code{showall = FALSE}.
#' @param holm logical with default \code{holm = FALSE}. If set to \code{holm = TRUE}, significance testing is based on the holm procedure -- see references. This argument is deprecated (since version 1.5.6) and kept only for downward compatibility. Use argument \code{adjalpha} for any type of alpha adjustment.  
#' @param wide logical with default \code{wide = FALSE}. If set to \code{wide = TRUE}, results for all significance tests are returned. 
#' @param adjalpha Adjustment method: \code{"bonferroni"} (default),
#'   \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{"BH"},
#'   \code{"BY"}, \code{"fdr"} (alias for BH), or \code{"none"}.
#'   \code{"holland"} remains available in this summary method.
#' @param ... other parameters passed trough.
#' @return A data frame of local results, with printed explanatory output.
#'   Includes a final \code{p.adj} column. For Holm adjustment, includes \code{Holm.crit} and \code{Holm.adj}.
#' @details Additional adjustment methods are \code{"hochberg"},
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
#' Functional cells are excluded from the Holm test family and retain
#' their \code{"b"} label. Other pre-existing adjustment methods retain their decisions.
#' @references Holm, S. (1979). A simple sequentially rejective multiple test procedure. \emph{Scandinavian Journal of Statistics, 6}(2), 65–70.
#' @references Bonferroni, C. E. (1935). Il calcolo delle assicurazioni su gruppi di teste. In S.O. Carboni (Ed.), \emph{Studi in Onore del Professore Salvatore Ortu Carboni} (S. 13–60). Roma, Tipografia del Senato: Bardi.
#' @references Holland, B. S., & Copenhaver, M. D. (1987). An Improved Sequentially Rejective Bonferroni Test Procedure. \emph{Biometrics, 43}(2), 417. https://doi.org/10.2307/2531823
#' @examples
#' data(LienertLSD)
#' res <- CFA(LienertLSD, bintest = FALSE)
#' summary(res, type = "pChi", adjalpha = "BH")


#sorted by the p-value selected in argument \code{type}
########################### hier die summary method #class CFA #######################
#digits=3; type="z.pChi";sorton=NULL; decreasing=FALSE; showall=TRUE; holm=TRUE; wide=FALSE; adjalpha = "bonferroni"

summary.CFA<-function(object, digits=3, type="z.pChi",sorton=NULL, decreasing=FALSE, showall=TRUE, holm=FALSE, wide=FALSE, adjalpha = "bonferroni", ...){
  local.test <- object$local.test
  global.test <- object$global.test
  functional <- object$functional
if(!is.null(functional)){   
 if(is.character(functional)){ ## added 22-01-2019 class(functional)=="character" changed 12-04-2022 ' inherits(x = class(functional),what = "character") ' again changed 19-08-2022 
   functional <- sapply(functional, function(x){which(x==as.character(object$local.test$pat.))})
 }
}
  # check for ex.bin.test present? # test added 23-06-2021
  if(all(is.na(object$local.test$ex.bin.test)) && type == "ex.bin.test"){stop("no results for 'ex.bin.test' available \n try to run CFA() again with argument 'bintest = TRUE' ")}

  
  # interpreting option "wide" # option added 8-2-2020 -------------------------
  if (wide==FALSE){
    if(type=="pChi"){out_wide <- c("Chi","df","pChi","Q.unst","Q.std")}# JHH added Q 13-11-2026
    if(type=="ex.bin.test"){out_wide <- c("df","ex.bin.test","Q.unst","Q.std")}# JHH added Q 13-11-2026
    if(type=="z.pChi"){out_wide <- c("df","z.Chi","z.pChi","Q.unst","Q.std")}# JHH added Q 13-11-2026
    if(type=="z.pBin"){out_wide <- c("df","z.Bin","z.pBin","cor.","Q.unst","Q.std")}# JHH added Q 13-11-2026
    if(type=="p.stir"){out_wide <- c("df","p.stir","density.stir","Q.unst","Q.std")}# JHH added Q 13-11-2026
    # if(type=="Q.unst"){out_wide <- c("df","z.Chi","z.pChi","Q.unst")} # JHH added 12-11-2026
    # if(type=="Q.std"){out_wide <- c("df","z.Chi","z.pChi","Q.std")} # JHH added 12-11-2026
  }else{out_wide <- c("Chi","df","pChi","ex.bin.test","z.Chi","z.pChi","z.Bin","z.pBin","cor.","p.stir","density.stir","Q.unst","Q.std")}  
  
  ############## START of alpha adjustment methods -----------------------------
  
  if(holm==TRUE){adjalpha <- "holm"} # 30-04-2021 for downward compatibility
  # Methodennamen pruefen; fdr ist der offizielle Alias von BH.
  adjalpha <- match.arg(adjalpha, c(stats::p.adjust.methods, "holland"))
  if (adjalpha == "fdr") adjalpha <- "BH"

  
  #### significant (anti)types cf. Holm ----------------------------------------
  ## added 21-01-2019 #changed 30-04-2021
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
    significant <- p_adjusted < object$alpha
    direction <- ifelse(local.test$obs. > local.test$exp., "+", "-")
    Type <- ifelse(significant, direction, ".")
    if (!is.null(functional)) Type[functional] <- "b"
    templocal <- data.frame(pat. = local.test$pat.,
      local.test[, c("obs.", "exp.")], Type, local.test[, out_wide],
      check.names = FALSE)
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
    direction <- ifelse(local.test$obs. > local.test$exp., "+", "-")
    Type <- ifelse(significant, direction, ".")
    if (!is.null(functional)) Type[functional] <- "b"
    templocal <- data.frame(pat. = local.test$pat.,
      local.test[, c("obs.", "exp.")], Type, Holm.crit = alpha_bon,
      Holm.adj = p_holm, local.test[, out_wide], check.names = FALSE)
  }
  
 
  ### significant (anti)types cf. (classical) Bonferroni adjustment ------------
  if(adjalpha=="bonferroni"){ ## condition added 21-01-2019 # 30-04-2021 new type of alpha adjustment control
  temp1 <- local.test[,which(names(local.test)==type)] < object$bonferroni.alpha
  temp2 <- ifelse(test=local.test$obs. > local.test$exp., yes="+",no="-") 
  Type <- mapply(FUN=function(xx,y){ifelse(test=(xx==TRUE),yes=y, no="." )   },xx=temp1,y=temp2 )
  if(!is.null(functional)){Type[functional] <- "b"} # added 18-11-2016
  templocal <- data.frame(pat.=local.test[,"pat."],local.test[,c("obs.","exp.")],Type, local.test[,out_wide])
  }
  
  ### significant (anti)types without any adjustment for multiple testing ------
  if(adjalpha=="none"){ ## condition 30-04-2021
    temp1 <- local.test[,which(names(local.test)==type)] < object$alpha
    temp2 <- ifelse(test=local.test$obs. > local.test$exp., yes="+",no="-") 
    Type <- mapply(FUN=function(xx,y){ifelse(test=(xx==TRUE),yes=y, no="." )   },xx=temp1,y=temp2 )
    if(!is.null(functional)){Type[functional] <- "b"} # added 18-11-2016
    templocal <- data.frame(pat.=local.test[,"pat."],local.test[,c("obs.","exp.")],Type, local.test[,out_wide])    
  }
  
  #### significant (anti)types as defined in Holland and Copenhaver (1987) ----------------------------------------
  ## added 24-11-2025 
  ### das WAR eine riesen B A U S T E L L E     ;-)))
  
  if(adjalpha=="holland"){ # 24-11-2025 new type of alpha adjustment control
    p_val <- local.test[,which(names(local.test)==type)]
    if(!is.null(functional)){p_val[functional] <- NA}# set all functional cells to NA
    
    
    ord <- order(p_val, na.last=NA)
    pvalues.sorted <- p_val[ord]
    k <- length(pvalues.sorted) + 1
    
    p.val.aux <- sapply(1:(k - 1), FUN=function(j, p.val, k){r <- 1 - (1 - p.val[j])^(k - j); return(r)},p.val=pvalues.sorted, k=k)
    
    p.adj.aux <- sapply(p.val.aux, FUN=function(x) {(min(1, x,na.rm = T))}) #     setUpperBound(p.val.aux, 1)
    
    p.adj.aux <- sapply(1:length(p.adj.aux), FUN=function(x) {(max(p.adj.aux[1:x],na.rm = T))} )# correctForMonotocity(p.adj.aux)
    
    p.adj <- rep(NA, length(p_val))
    suppressWarnings(expr = {
      p.adj[ord] <- p.adj.aux
    })
    
    # p.adj
    
    # ord_p <- order(p_val,decreasing = FALSE) # starting with the smallest 'p_val'
    # k <- length(p_val)
    # n <- length(na.omit(p_val))
    alpha <- object$alpha # 30-04-2021 according to new structure of CFA result object
    
    # alpha_bon_temp1 <-(alpha/1:n)[(order(rev(ord_p[1:(length(ord_p)-length(functional))])))]
    # alpha_bon_temp2 <- p_val
    # alpha_bon_temp2[!is.na(p_val)] <- alpha_bon_temp1
    # alpha_bon <- alpha_bon_temp2
    
    temp1 <- p.adj < alpha 
    temp2 <- ifelse(test=local.test$obs. > local.test$exp., yes="+",no="-") 
    Type <- mapply(FUN=function(xx,y){ifelse(test=(xx==TRUE),yes=y, no="." )   },xx=temp1,y=temp2 )
    if(!is.null(functional)){Type[functional] <- "b"} # added 18-11-2016
    
    templocal <- data.frame(pat.=local.test[,"pat."],local.test[,c("obs.","exp.")],Type, Holland.adj=p.adj ,local.test[,out_wide])  
  }
  
  # Einheitliche Zusatzspalte am ENDE: bisherige Spaltenpositionen erhalten.
  # Holm.adj und Holm.crit bleiben als bestehende Spalten unveraendert verfuegbar.
  if (adjalpha %in% c("hochberg", "hommel", "BH", "BY")) {
    generic_adjusted <- p_adjusted
  } else if (adjalpha == "holm") {
    generic_adjusted <- p_holm
  } else if (adjalpha == "holland") {
    generic_adjusted <- p.adj
  } else {
    # Bonferroni/none behalten die bisherige Testfamilie und Entscheidung.
    # Insbesondere wird Bonferroni NICHT still auf weniger Vergleiche umgestellt.
    values <- local.test[[type]]
    usable <- is.finite(values) & values >= 0 & values <= 1
    if (!is.null(functional)) usable[functional] <- FALSE
    generic_adjusted <- rep(NA_real_, length(values))
    if (adjalpha == "bonferroni") {
      historical_m <- object$alpha / object$bonferroni.alpha
      generic_adjusted[usable] <- pmin(1, values[usable] * historical_m)
    } else {
      generic_adjusted[usable] <- stats::p.adjust(values[usable], method = "none")
    }
  }
  templocal$p.adj <- generic_adjusted

  ############## END of alpha adjustment methods -------------------------------
  
  
  #print(templocal)
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
    erg <- erg[!is.na(erg$Type) & erg$Type != ".", , drop = FALSE]
  }
  
  #added (generic overall) rounding here 8-2-2020
  #round(erg,digits = digits)
  round_index <- sapply((names(erg)[-1]), function(x){is.numeric(erg[,x])},USE.NAMES = FALSE)
  round_index <- c(FALSE,round_index)
  erg[,round_index] <- round(erg[,round_index],digits = digits)
  erg
  
  cat("function Call:","\n","-------------","\n","Formula:",paste(object$used.formula,collapse = " "),"\n","Variables:", names(object$variables),"\n","Categories:", object$variables,"\n")
  
  cat("\n","results of global tests:",attributes(object)$comment,"\n", "-----------------------")
  cat("\n","pearson Chi-square test:",attributes(object)$comment,"\n")
  print(data.frame(global.test$pearson))
  
  cat("\n","likelihood ratio test:",attributes(object)$comment,"\n")
  print(data.frame(global.test$likelihood.ratio))
  
  cat("\n","Information Criteria:",attributes(object)$comment,"\n")
  print(data.frame(global.test$infocrit))
  
  if(adjalpha=="bonferroni"){
    cat("\n","results of local tests:",attributes(object)$comment,"\n", "-----------------------","\n","Type (+) / Antitype (-) based on:", type ,";","\n", "with Bonferroni adjusted alpha:", object$bonferroni.alpha,"\n")#changed 30-04-2021
  }
  if(adjalpha=="holm"){
    cat("\n","results of local tests:",attributes(object)$comment,"\n", "-----------------------","\n","Type (+) / Antitype (-) based on:", type, ";","\n","with Holm adjusted p-values - see \"Holm.adj\" \n")#added 21-01-2019 #changed 30-04-2021
  }
  if(adjalpha=="none"){
    cat("\n","results of local tests:",attributes(object)$comment,"\n", "-----------------------","\n","Type (+) / Antitype (-) based on:", type, ";","\n","with not adjusted alpha:", object$alpha,"\n")#added 30-04-2021
  }
  if(adjalpha=="holland"){
    cat("\n","results of local tests:",attributes(object)$comment,"\n", "-----------------------","\n","Type (+) / Antitype (-) based on:", "Holland and Copenhaver (1987) adjusted p-values from:", type, ";","\n","see column 'Holland.adj'","against not adjusted alpha =", object$alpha,"\n")#added 24-11-2025
  }
 
  cat("\n", "upper bound of determination coefficient Q is:", round(object$Qupper, digits),"\n" )#JHH added 12-11-2026
  #if(any(temp3==FALSE)){ cat("\n","Type (b): blanked out (functional CFA)","\n")}#(28-04-2015)}
  if(any(Type=="b", na.rm = TRUE)){ cat("\n","Type (b): blanked out (functional CFA)","\n")}#(18-11-2016)}
  
  if (adjalpha %in% c("hochberg", "hommel", "BH", "BY")) {
    cat("\nAdjustment method:", adjalpha,
        "; compare p.adj with alpha =", object$alpha, "\n")
  }
  return(erg)
}
