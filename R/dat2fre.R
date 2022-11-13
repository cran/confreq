#' @title dataset to pattern frequency conversion 
#' @keywords utilities 
#' @export dat2fre
#' @description Given a dataset this function returns a (response) pattern frequencies table representation of it.
#' @details To use survey weights a vector of positive numeric values with length matching the number of rows in \code{x} must be assigned to the argument \code{wgt}. The individual case weights are then aggregated (respective sum of weights) for each pattern observed in the data (assigned to argument \code{x}).   
#' 
#' @param x an object of class "matrix" or "data.frame". If \code{x} is a "data.frame" each variable (column) must be an \code{integer} or a \code{factor}. If \code{x} is a "matrix" it is assumed that the categories for each variable in \code{x} start with \code{1} -- there is no check for that !!! 
#' @param katorder logical with default set to \code{katorder==FALSE}. When set to \code{katorder==TRUE} variables are ordered according to their number of categories (variable with most categories is the rightmost variable) in the resulting object. 
#' @param caseorder logical with default set to \code{caseorder==TRUE}. When set to \code{caseorder==FALSE} configurations are only ordered according to the categories of the rightmost variable in the resulting object.
#' @param kat ignored when \code{x} is a \code{data.frame}! If \code{x} is a \code{"matrix"} the optional argument \code{kat} must be an integer vector defining the number of categories for every variable in \code{x} (in the respective order). If left empty the (max) number of categories is estimated from the data given in \code{x}.
#' @param codes a list with character vectors containing coding for integers in matrix (if \code{x} is a numeric \code{matrix}). If \code{codes} is not empty (and the argument \code{x} is an object of class "matrix") the return object will be pattern frequencies table as \code{data.frame}.
#' @param wgt a numerical vector of survey weights to weight the cases (rows) in \code{x}   
#' @param ... other parameters passed through to \code{\link{table}} (in case of x  being a \code{data.frame}) or to \code{\link{tabulate}} (in case of x being a matrix).
#' @return An object of class c("data.frame","Pfreq") containing the (response) pattern frequencies table representation of the given dataset in the argument \code{x}. 
#' @references No references in the moment 
#' @examples #######################################
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
dat2fre <- function(x, katorder=FALSE, caseorder=TRUE, kat=NULL, codes=NULL, wgt=NULL, ...) {
  #### wenn x ein data.frame ist ... 
  #### funktioniert bei: nur factors, nur integers, und gemischt
  if( is(object = x,class2 = "data.frame") ){
    # if(sum(is.na(x))!=0){cat("\n","x contains missing values (NA) implicitly treated with na.rm=TRUE","\n")}
    unorderedresult<-as.data.frame(table(x, ...))# ungeordnetes tabulationsergebnis 
    
  if(katorder==TRUE){### sortieren der variablen nach deren kategoriezahl  
    variablesorter<-order(sapply(  sapply(unorderedresult[,1:dim(unorderedresult)[2]-1],table,simplify=F)   ,length))
    # variablesorter<-order(sapply(apply(unorderedresult[,1:dim(unorderedresult)[2]-1],2,table),length))
    varorderedresult<-unorderedresult[,c(variablesorter,dim(unorderedresult)[2])]} else{varorderedresult<-unorderedresult} 
    
   if(caseorder==TRUE){ ### sortieren der 'fälle' bzw konfigurationen nach variablen kategorien aufsteigend
    calllist<-lapply(varorderedresult,rank)
    do.call("order", calllist)
    result <- varorderedresult[do.call("order", calllist), ] 
    rownames(result)<-1:dim(result)[1] } else{result <- varorderedresult }
    
    cat("Number of categories for each variable","\n","estimated from data are: " , "\n",sapply(lapply((result[,1:(dim(result)[2]-1)]),levels),length) ,"\n", "-->",dim(result)[1],"different configurations","\n")
    
    freq <- result$Freq # backup sample frequencies in case of weights
        # version 1.6.0:
    if(!is.null(wgt)){ 
      all_pat <- apply(result[,1:(ncol(result)-1)],1, function(v){paste(v, collapse = "")})
      tmp <- data.frame(pat.=factor(apply(x, 1, function(v){paste(v, collapse = "")}),levels=all_pat ), wgt=wgt)
      Freq<-aggregate(x = tmp[,2],by = list(tmp$pat.),FUN = sum,drop = FALSE)[,2]
      Freq[is.na(Freq)] <- 0
      ### neuer gedanke: es werden gewichte für die tabulierten Häufigkeiten erstellt
      # result$Freq <- Freq
      frwgt <- Freq / freq ## frwgt sind die Gewicht für die (tabulierten) pattern Häufigkeiten
      frwgt[is.na(frwgt)] <- 1 ## null mal beobachtete Patern bleiben bei einer gewichtung mit 1
      # freq*frwgt # test OK
    }
    #end version 1.6.0  
    }
 #### ENDE von wenn x ein data.frame ist ...
 
 #### wenn x eine "matrix" ist ...     
  if( is(object = x,class2 = "matrix") ){
    if(is.numeric(x)!=TRUE){stop("x should be a numeric matrix !") }
    X<-x
    ### optional anzahl kategorien festlegen (falls nicht alle beobachtet wurden)
  if(length(kat)==0){
    
    
    kat<-apply(X,2,max, ...)
    cat("\n","Number of categories for each variable","\n","estimated from data are: " , "\n",names(kat),"\n", kat ,"\n", "-->",prod(kat),"different configurations","\n")
	}
  
    if(is.na(prod(kat))){stop("x contains missing values (NA) please recode or consider use na.rm=TRUE") }
    
      ### ENDE von optional anzahl kateg ....
  konfig<-pos_cfg_cfa(kat)
  if(length(colnames(x))!=0){colnames(konfig)<-colnames(X)} # variablennmen aus x uebertragen
  if(length(colnames(x))==0){colnames(konfig)<-paste("V",1:length(kat),sep="")}
  all_pat<-apply(konfig,1, function(v){paste(v, collapse = "")})
  
    Freq<-tabulate( factor(apply(X, 1, function(v){paste(v, collapse = "")}),levels=all_pat ),nbins=length(all_pat) )
    freq <- Freq # backup sample frequencies in case of weights
  # version 1.6.0:
    if(!is.null(wgt)){
    tmp <- data.frame(pat.=factor(apply(X, 1, function(v){paste(v, collapse = "")}),levels=all_pat ), wgt=wgt)
    Freq<-aggregate(x = tmp[,2],by = list(tmp$pat.),FUN = sum,drop = FALSE)[,2]
    Freq[is.na(Freq)] <- 0
    ### neuer gedanke: es werden gewicht für die tabulierten Häufigkeiten erstellt
    # result$Freq <- Freq
    frwgt <- Freq / freq ## frwgt sind die Gewicht für die (tabulierten) pattern Häufigkeiten
    frwgt[is.na(frwgt)] <- 1 ## null mal beobachtete Patern bleiben bei einer gewichtung mit 1
    # freq*frwgt # test OK
    }
  #end version 1.6.0  
  
    konfig<-data.frame(apply(konfig,2,factor),stringsAsFactors=TRUE)# JHH 09-03-2021: ,stringsAsFactors=TRUE
 
       
    result<-cbind(konfig,freq) # the default return
    
    
    ### optionale umwandlung in data.frame mit codesx
    if(length(codes)!=0){
      konfigdata<-as.data.frame(matrix(NA,ncol=(dim(konfig)[2]),nrow=dim(konfig)[1]))
      for (i in 1:dim(konfig)[2]){
        konfigdata[,i]<-factor(konfig[,i],labels=names(codes[[i]])) 
      }
      colnames(konfigdata)<-colnames(x)  
      result<-data.frame(konfigdata,Freq)    
    } 
  }
  #### ENDE von wenn x eine "matrix" ist ...   
  class(result)<-c("Pfreq","data.frame")  # was before 1.6 c("data.frame", "Pfreq")
  if(!is.null(wgt)){
    comment(result) <- "using weighted data!"
    attr(result,"wgtdat") <- wgt # save original data vector of weights
    attr(result,"wgtfreq") <- Freq
    attr(result,"WGT") <- frwgt
    cat("using weighted data!", "\n")
    }
  
return(result)
}