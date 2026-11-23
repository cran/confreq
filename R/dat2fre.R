#' @title Convert Individual Data to Pattern Frequencies
#' @keywords utilities 
#' @export dat2fre
#' @description Tabulate individual records as a complete table of configurations,
#' including configurations with zero observed frequency.
#' Given a dataset this function returns a (response) pattern frequencies table representation of it.
#' @details For data frames, factor levels define available categories and
#' their order, including unused levels. Other columns use observed categories
#' in the order determined by \code{table()}. Use factor levels to specify
#' additional categories. The arguments \code{kat} and \code{codes} are
#' rejected for data frames rather than silently ignored.
#'
#' Numeric matrices retain the historical assumption of codes starting at
#' one. Use a data frame for other coding schemes or character labels.
#' Missing matrix codes are rejected. For data frames, the default table
#' options exclude missing categories; \code{useNA = "ifany"} includes them.
#'
#' The \code{Freq} column always contains unweighted counts. If weights are
#' supplied, \code{wgtfreq} stores cell sums of weights, \code{WGT} stores
#' the mean weight per observed cell (one for unobserved cells), and
#' \code{wgtdat} stores the original weight vector. A comment marks the object
#' as weighted for compatibility with \code{CFA()}.
#' The additional attributes \code{wgtfreq2} (cell sums of squared weights)
#' and \code{wgtpositive} (cell counts of positive-weight records) support
#' the independent-unit survey variance calculation in \code{S2CFA()}.
#' Both attributes follow exactly the row order of the frequency table.
#' Use the default \code{caseorder = TRUE} for \code{CFA()} and use identical
#' variable selection and ordering options when preparing covariates.
#'
#' @param x A data frame containing factor, character, numeric or logical
#'   columns, or a numeric matrix of positive integer category codes.
#'   Column names must be unique; \code{Freq} is reserved for the output.
#' @param katorder Logical. Reorder variables by increasing number of
#'   categories? Defaults to \code{FALSE}.
#' @param caseorder Logical. If \code{TRUE} (default), order configurations
#'   by category order with the rightmost variable changing fastest.
#'   If \code{FALSE}, retain the order produced by \code{table()}.
#' @param kat For numeric matrices only: a vector specifying the number
#'   of categories per column. Codes must be in \code{1:kat[j]}.
#'   If omitted, column maxima determine the category counts.
#' @param codes For numeric matrices only: a list with one entry per column
#'   and one unique label per category. For named entries, names are used
#'   as labels; otherwise, entry values are used as labels.
#' @param wgt Optional numeric vector of finite, nonnegative case weights,
#'   one per input row. Zero weights are allowed.
#' @param ... Named options passed to \code{\link[base]{table}}:
#'   \code{exclude}, \code{useNA}, \code{dnn} and \code{deparse.level}.
#'   Other or unnamed arguments are rejected.
#' @return A data frame of class \code{c("Pfreq", "data.frame")}, with one
#'   factor column per configuration variable and a final \code{Freq} column.
#'   Weighted inputs additionally produce the attributes described above.
#' @seealso \code{\link{fre2dat}}, \code{\link{dat2cov}},
#'   \code{\link{fre2tab}}, \code{\link{CFA}}, \code{\link{S2CFA}}
#' @references No references in the moment 
#' @examples #######################################
#' d <- data.frame(A = factor(c("yes", "no", "yes")),
#'                 B = factor(c("low", "low", "high")))
#' pf <- dat2fre(d)
#' pf
#' weighted <- dat2fre(d, wgt = c(1, 2, 3))
#' attr(weighted, "wgtfreq")
#' ########################################
#' 
#' data(suicide)# loading data in data frame (702 cases) representation 
#' dat2fre(suicide) # converting it into a pattern frequencies table
#' 
#' ########### 
#' #######################################
#' data(LienertLSD)# loading example pattern frequencies table ..
#' test<-fre2dat(LienertLSD)# and coverting it into a simple (data) matrix
#' test<-test[sample(c(1:65),65),] # making a messy order 
#' ############
#' dat2fre(test) # making a proper ordered pattern frequencies table again
#' ##### try it with a data.frame too!
#' #######################################
############### start of function definition ##################
###########  data to pattern frequency conversion #############
################ jhheine at googlemail.com ####################
# dat2fre: Individualdaten -> vollstaendige Konfigurationshaeufigkeiten.
# Fuer numerische Matrizen bleibt die historische Codierung 1,...,kat erhalten.
# Bei Dataframes legen Faktorstufen die Kategorien fest; andere Spalten
# verwenden die sortierten beobachteten Werte. Labels werden nicht verkettet.
dat2fre <- function(x, katorder = FALSE, caseorder = TRUE, kat = NULL,
                    codes = NULL, wgt = NULL, ...) {
  matrix_input <- is.matrix(x)
  d <- .cf_data(x)
  .cf_weights(wgt, nrow(d))
  if ("Freq" %in% names(d)) stop("'Freq' is reserved for the frequency column.")
  if (!is.logical(katorder) || length(katorder) != 1L || is.na(katorder) ||
      !is.logical(caseorder) || length(caseorder) != 1L || is.na(caseorder))
    stop("katorder and caseorder must be TRUE or FALSE.")
  if (!matrix_input && (!is.null(kat) || !is.null(codes)))
    stop("kat and codes apply only to numeric matrices; use factor levels for dataframes.")
  
  if (matrix_input) {
    if (!is.numeric(x)) stop("Matrix input must be numeric; use a dataframe for labels.")
    if (anyNA(x) || any(!is.finite(x)) || any(x < 1 | x != floor(x)))
      stop("Numeric matrices require positive integer codes without missing values.")
    if (is.null(kat)) {
      if (!nrow(d)) stop("For empty matrices, supply kat.")
      kat <- vapply(d, max, numeric(1))
    }
    if (!is.numeric(kat) || length(kat) != ncol(d) || anyNA(kat) ||
        any(!is.finite(kat)) || any(kat < 1 | kat != floor(kat)))
      stop("kat must give one positive integer category count per column.")
    if (!is.null(codes) && (!is.list(codes) || length(codes) != ncol(d)))
      stop("codes must be a list with one entry per column.")
    for (j in seq_along(d)) {
      if (any(d[[j]] > kat[j])) stop("Observed code exceeds kat.")
      labels <- as.character(seq_len(kat[j]))
      if (!is.null(codes)) {
        # Historische benannte codes-Listen verwenden deren Namen als Labels.
        labels <- if (is.null(names(codes[[j]]))) as.character(codes[[j]]) else names(codes[[j]])
        if (length(labels) != kat[j] || anyNA(labels) || anyDuplicated(labels))
          stop("Each codes entry needs one unique label per category.")
      }
      d[[j]] <- factor(d[[j]], levels = seq_len(kat[j]), labels = labels)
    }
  }
  
  # table() uebernimmt Kategorien und den expliziten Umgang mit NA.
  # Unbekannte Argumente werden frueh gemeldet statt unbemerkt ignoriert.
  dots <- list(...)
  if (length(dots) && (is.null(names(dots)) || anyNA(names(dots)) ||
                       any(!nzchar(names(dots))) || anyDuplicated(names(dots)) ||
                       any(!names(dots) %in% c("exclude", "useNA", "dnn", "deparse.level"))))
    stop("... accepts named table options: exclude, useNA, dnn, deparse.level.")
  tab <- do.call(table, c(list(d), dots))
  if (any(dim(tab) == 0L)) stop("Every variable needs at least one category.")
  result <- as.data.frame(tab, stringsAsFactors = TRUE, responseName = "Freq")
  if (anyDuplicated(names(result)) || "Freq" %in% names(result)[-ncol(result)])
    stop("Table dimension names must be unique and must not be 'Freq'.")
  p <- ncol(result) - 1L
  
  # Die Zuordnung erfolgt spaltenweise durch match(), auch fuer NA-Kategorien.
  # Gemischte Stellenwerte ergeben eindeutige Zellnummern; Labels wie
  # (1,11), (11,1) oder Labels mit Leerzeichen koennen nicht kollidieren.
  dims <- dim(tab)
  stride <- c(1, cumprod(dims)[-length(dims)])
  cell <- rep(1, nrow(d))
  for (j in seq_len(p))
    cell <- cell + (match(as.character(d[[j]]), dimnames(tab)[[j]]) - 1) * stride[j]
  if (!is.null(wgt)) {
    sums <- numeric(nrow(result))
    # Survey-Varianzen brauchen neben Summe(w) auch Summe(w^2) je Zelle.
    # Die Zahl positiver Gewichtstraeger bestimmt die Stichprobenkorrektur.
    square_sums <- numeric(nrow(result))
    positive_counts <- numeric(nrow(result))
    keep <- !is.na(cell)
    if (any(keep)) {
      grouped <- rowsum(wgt[keep], group = cell[keep], reorder = FALSE)
      positions <- as.integer(rownames(grouped))
      sums[positions] <- grouped[, 1L]
      grouped_squares <- rowsum(wgt[keep]^2, group = cell[keep], reorder = FALSE)
      grouped_positive <- rowsum(as.numeric(wgt[keep] > 0), group = cell[keep], reorder = FALSE)
      square_sums[positions] <- grouped_squares[, 1L]
      positive_counts[positions] <- grouped_positive[, 1L]
    }
  }
  # Optional Variablen umordnen; standardmaessig bleibt die Eingabereihenfolge.
  if (katorder) {
    vars <- order(dims)
    result <- result[, c(vars, p + 1L), drop = FALSE]
  }
  # CFA erwartet: rechte Variable wechselt am schnellsten.
  # caseorder=FALSE behaelt die table()-Reihenfolge fuer bestehende Aufrufe.
  idx <- seq_len(nrow(result))
  if (caseorder) idx <- do.call(order, lapply(result[seq_len(p)], as.integer))
  result <- result[idx, , drop = FALSE]
  rownames(result) <- NULL
  class(result) <- c("Pfreq", "data.frame")
  if (!is.null(wgt)) {
    weighted <- sums[idx]
    mean_weight <- rep(1, nrow(result))
    observed <- result$Freq > 0
    mean_weight[observed] <- weighted[observed] / result$Freq[observed]
    comment(result) <- "using weighted data!"
    attr(result, "wgtdat") <- wgt
    attr(result, "wgtfreq") <- weighted
    attr(result, "wgtfreq2") <- square_sums[idx]
    attr(result, "wgtpositive") <- positive_counts[idx]
    attr(result, "WGT") <- mean_weight
  }
  result
}



### old Stuff
# dat2fre <- function(x, katorder=FALSE, caseorder=TRUE, kat=NULL, codes=NULL, wgt=NULL, ...) {
#   #### wenn x ein data.frame ist ... 
#   #### funktioniert bei: nur factors, nur integers, und gemischt
#   if( is(object = x,class2 = "data.frame") ){
#     # if(sum(is.na(x))!=0){cat("\n","x contains missing values (NA) implicitly treated with na.rm=TRUE","\n")}
#     unorderedresult<-as.data.frame(table(x, ...))# ungeordnetes tabulationsergebnis 
#     
#   if(katorder==TRUE){### sortieren der variablen nach deren kategoriezahl  
#     variablesorter<-order(sapply(  sapply(unorderedresult[,1:dim(unorderedresult)[2]-1],table,simplify=F)   ,length))
#     # variablesorter<-order(sapply(apply(unorderedresult[,1:dim(unorderedresult)[2]-1],2,table),length))
#     varorderedresult<-unorderedresult[,c(variablesorter,dim(unorderedresult)[2])]} else{varorderedresult<-unorderedresult} 
#     
#    if(caseorder==TRUE){ ### sortieren der 'fälle' bzw konfigurationen nach variablen kategorien aufsteigend
#     calllist<-lapply(varorderedresult,rank)
#     do.call("order", calllist)
#     result <- varorderedresult[do.call("order", calllist), ] 
#     rownames(result)<-1:dim(result)[1] } else{result <- varorderedresult }
#     
#     cat("Number of categories for each variable","\n","estimated from data are: " , "\n",sapply(lapply((result[,1:(dim(result)[2]-1)]),levels),length) ,"\n", "-->",dim(result)[1],"different configurations","\n")
#     
#     freq <- result$Freq # backup sample frequencies in case of weights
#         # version 1.6.0:
#     if(!is.null(wgt)){ 
#       all_pat <- apply(result[,1:(ncol(result)-1)],1, function(v){paste(v, collapse = "")})
#       tmp <- data.frame(pat.=factor(apply(x, 1, function(v){paste(v, collapse = "")}),levels=all_pat ), wgt=wgt)
#       Freq<-aggregate(x = tmp[,2],by = list(tmp$pat.),FUN = sum,drop = FALSE)[,2]
#       Freq[is.na(Freq)] <- 0
#       ### neuer gedanke: es werden gewichte für die tabulierten Häufigkeiten erstellt
#       # result$Freq <- Freq
#       frwgt <- Freq / freq ## frwgt sind die Gewicht für die (tabulierten) pattern Häufigkeiten
#       frwgt[is.na(frwgt)] <- 1 ## null mal beobachtete Patern bleiben bei einer gewichtung mit 1
#       # freq*frwgt # test OK
#     }
#     #end version 1.6.0  
#     }
#  #### ENDE von wenn x ein data.frame ist ...
#  
#  #### wenn x eine "matrix" ist ...     
#   if( is(object = x,class2 = "matrix") ){
#     if(is.numeric(x)!=TRUE){stop("x should be a numeric matrix !") }
#     X<-x
#     ### optional anzahl kategorien festlegen (falls nicht alle beobachtet wurden)
#   if(length(kat)==0){
#     
#     
#     kat<-apply(X,2,max)
#     cat("\n","Number of categories for each variable","\n","estimated from data are: " , "\n",names(kat),"\n", kat ,"\n", "-->",prod(kat),"different configurations","\n")
# 	}
#   
#     if(is.na(prod(kat))){stop("x contains missing values (NA) please recode or consider use na.rm=TRUE") }
#     
#       ### ENDE von optional anzahl kateg ....
#   konfig<-pos_cfg_cfa(kat)
#   if(length(colnames(x))!=0){colnames(konfig)<-colnames(X)} # variablennmen aus x uebertragen
#   if(length(colnames(x))==0){colnames(konfig)<-paste("V",1:length(kat),sep="")}
#   all_pat<-apply(konfig,1, function(v){paste(v, collapse = "")})
#   
#     Freq<-tabulate( factor(apply(X, 1, function(v){paste(v, collapse = "")}),levels=all_pat ),nbins=length(all_pat) )
#     freq <- Freq # backup sample frequencies in case of weights
#   # version 1.6.0:
#     if(!is.null(wgt)){
#     tmp <- data.frame(pat.=factor(apply(X, 1, function(v){paste(v, collapse = "")}),levels=all_pat ), wgt=wgt)
#     Freq<-aggregate(x = tmp[,2],by = list(tmp$pat.),FUN = sum,drop = FALSE)[,2]
#     Freq[is.na(Freq)] <- 0
#     ### neuer gedanke: es werden gewicht für die tabulierten Häufigkeiten erstellt
#     # result$Freq <- Freq
#     frwgt <- Freq / freq ## frwgt sind die Gewicht für die (tabulierten) pattern Häufigkeiten
#     frwgt[is.na(frwgt)] <- 1 ## null mal beobachtete Patern bleiben bei einer gewichtung mit 1
#     # freq*frwgt # test OK
#     }
#   #end version 1.6.0  
#   
#     konfig<-data.frame(apply(konfig,2,factor),stringsAsFactors=TRUE)# JHH 09-03-2021: ,stringsAsFactors=TRUE
#  
#        
#     result<-cbind(konfig,Freq) # the default return -- -6.8.2024 changed to Freq (instead of freq)
#     
#     
#     ### optionale umwandlung in data.frame mit codesx
#     if(length(codes)!=0){
#       konfigdata<-as.data.frame(matrix(NA,ncol=(dim(konfig)[2]),nrow=dim(konfig)[1]))
#       for (i in 1:dim(konfig)[2]){
#         konfigdata[,i]<-factor(konfig[,i],labels=names(codes[[i]])) 
#       }
#       colnames(konfigdata)<-colnames(x)  
#       result<-data.frame(konfigdata,Freq)    
#     } 
#   }
#   #### ENDE von wenn x eine "matrix" ist ...   
#   class(result)<-c("Pfreq","data.frame")  # was before 1.6 c("data.frame", "Pfreq")
#   
#   if(!is.null(wgt)){
#     comment(result) <- "using weighted data!"
#     attr(result,"wgtdat") <- wgt # save original data vector of weights
#     attr(result,"wgtfreq") <- Freq
#     attr(result,"WGT") <- frwgt
#     cat("using weighted data!", "\n")
#     }
#   
# return(result)
# }