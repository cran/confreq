#' @export summary.CFA
#' @S3method summary CFA
#' @title S3 Summary for CFA
#' @description S3 summary method for object of class\code{"CFA"}
#' @param object object of class\code{"CFA"}
#' @param digits integer rounds the values to the specified number of decimal places, default is \code{digits=3}.
#' @param type character indicating which test to use for inference wether the observed pattern are 'Types', 'Antitypes' or not significant at all. Possible options for \code{type} are \code{"pChi"}, \code{"ex.bin.test"}, \code{"z.pChi"}, \code{"z.pBin"} and \code{"p.stir"}. 
#' @param sorton sort results of local test by any column. by default the output is not sorted. Other options may be \code{"pat."}, \code{"obs."}, \code{"exp."}, \code{"Type"}, \code{"Chi"}, etc. ...
#' @param decreasing logical. Should the sort be increasing or decreasing? see \code{\link{order}}
#' @param showall logical with default \code{showall = TRUE}. To return only significant pattern (types / antitypes) set it to \code{showall = FALSE}.
#' @param holm logical with default \code{holm = FALSE}. If set to \code{holm = TRUE}, significance testing is based on the holm procedure -- see references.  
#' @param ... other parameters passed trough
#' @references Holm, S. (1979). A simple sequentially rejective multiple test procedure. \emph{Scandinavian Journal of Statistics, 6}(2), 65–70.


#sorted by the p-value selected in argument \code{type}
########################### hier die summary method #class CFA #######################
summary.CFA<-function(object, digits=3, type="z.pChi",sorton=NULL, decreasing=FALSE, showall=TRUE, holm=FALSE, ...){
  local.test <- object$local.test
  global.test <- object$global.test
  functional <- object$functional
  
 if(class(functional)=="character"){ ## added 22-01-2019
   functional <- sapply(functional, function(x){which(x==as.character(object$local.test$pat.))})
 }
  
  # significant (anti)types cf. Holm 
  ## added 21-01-2019
  if(holm==TRUE){
  p_val <- local.test[,which(names(local.test)==type)]
    if(!is.null(functional)){p_val[functional] <- NA}
  
  ord_p <- order(p_val,decreasing = FALSE) # starting with the smallest
  k <- length(p_val)
  n <- length(na.omit(p_val))
  alpha <- object$global.test$pearson$alpha
  i=0
  holm_res <- vector(mode = "logical", length = k)
  if(i<k){
  while(p_val[ord_p[i+1]]  <  alpha/(n-(i))){
    holm_res[i+1] <- p_val[ord_p[i+1]]  <  alpha/(n-(i)) 
    i=i+1
    }
  }
  temp1 <- holm_res[order(ord_p)] # reorders it back
  temp2 <- ifelse(test=local.test$obs. > local.test$exp., yes="+",no="-") 
  Type <- mapply(FUN=function(x,y){ifelse(test=(x==TRUE),yes=y, no="." )   },x=temp1,y=temp2 )
  if(!is.null(functional)){Type[functional] <- "b"} # added 18-11-2016
  alpha_bon <-(c((alpha/1:n),rep(NA,times=sum(is.na(p_val)))))[order(ord_p)]
  # alpha_bon[is.na(alpha_bon)] <-0 
  }
  
  #object$bonferroni.alpha
  if(holm==FALSE){ ## condition added 21-01-2019
  temp1 <- local.test[,which(names(local.test)==type)] < object$bonferroni.alpha
  temp2 <- ifelse(test=local.test$obs. > local.test$exp., yes="+",no="-") 
  Type <- mapply(FUN=function(x,y){ifelse(test=(x==TRUE),yes=y, no="." )   },x=temp1,y=temp2 )
  
  if(!is.null(functional)){Type[functional] <- "b"} # added 18-11-2016
  }
  
  if(holm==FALSE){
  templocal <- data.frame(pat.=local.test[,1],round(local.test[,2:3],digits=digits),Type, round(local.test[,4:11],digits=digits)    ,cor.=local.test[,12],round(local.test[,13:14],digits=digits) )
  }
  
  if(holm==TRUE){
    templocal <- data.frame(pat.=local.test[,1],round(local.test[,2:3],digits=digits),Type, Holm.crit=round(alpha_bon,digits = digits),round(local.test[,4:11],digits=digits)    ,cor.=local.test[,12],round(local.test[,13:14],digits=digits) )
  }
  
  #print(templocal)
  # sort output by sorton ...
  if (length(sorton)!=0){
  sorter <- order(templocal[,which(names(templocal)==sorton)],decreasing=decreasing)
  erg <- templocal[sorter,]
  }
  if (length(sorton)==0){
    erg <- templocal
  }
  if (showall==FALSE){
    erg <- erg[erg$Type!=".",]
  }
  
  
  cat("function Call:","\n","-------------","\n","Formula:",paste(object$used.formula,collapse = " "),"\n","Variables:", names(object$variables),"\n","Categories:", object$variables,"\n")
  
  cat("\n","results of global tests:","\n", "-----------------------")
  cat("\n","pearson Chi-square test:","\n")
  print(data.frame(global.test$pearson))
  
  cat("\n","likelihood ratio test:","\n")
  print(data.frame(global.test$likelihood.ratio))
  
  cat("\n","Information Criteria:","\n")
  print(data.frame(global.test$infocrit))
  
  if(holm==FALSE){
    cat("\n","results of local tests:","\n", "-----------------------","\n","Type (+) / Antitype (-) based on:", type, "; Bonferoni adj. alpha:", object$bonferroni.alpha,"\n")
  }
  if(holm==TRUE){
    cat("\n","results of local tests:","\n", "-----------------------","\n","Type (+) / Antitype (-) based on:", type, "; with Bonferoni-Holm adjustment - see \"Holm.crit\" \n")#added 21-01-2019
  }
  
  #if(any(temp3==FALSE)){ cat("\n","Type (b): blanked out (functional CFA)","\n")}#(28-04-2015)}
  if(any(Type=="b")){ cat("\n","Type (b): blanked out (functional CFA)","\n")}#(18-11-2016)}
  
  return(erg)
}