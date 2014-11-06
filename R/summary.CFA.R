#' @export summary.CFA
#' @title S3 Summary for CFA
#' @description S3 summary method for object of class\code{"CFA"}
#' @param object object of class\code{"CFA"}
#' @param digits integer rounds the values to the specified number of decimal places, default is \code{digits=3}.
#' @param type character with default \code{type="z.pChi"}, to return wether the observed pattern are 'Types', 'Antitypes' or not significant at all. Possible options for \code{type} are \code{"pChi"}, \code{"ex.bin.test"}, \code{"z.pChi"}, \code{"z.pBin"} and \code{"p.stir"}.     
#' @param ... other parameters passed trough

########################### hier die summary method #class CFA #######################
summary.CFA<-function(object, digits=3, type="z.pChi",...){
  local.test <- object$local.test
  global.test <- object$global.test
    
  object$bonferroni.alpha
  temp1 <- local.test[,which(names(local.test)==type)] < object$bonferroni.alpha
  temp2 <- ifelse(test=local.test$obs. > local.test$exp., yes="+",no="-") 
  
  Type <- mapply(FUN=function(x,y){ifelse(test=(x==TRUE),yes=y, no="." )   },x=temp1,y=temp2 )
  
  templocal <- data.frame(pat.=local.test[,1],round(local.test[,2:3],digits=digits),Type, round(local.test[,4:11],digits=digits)    ,cor.=local.test[,12],round(local.test[,13:14],digits=digits) )
  #print(templocal)
  cat("results of global tests:", "\n")
  cat("pearson Chi-square test:","\n")
  print(data.frame(global.test$pearson))
  
  cat("likelihood ratio test:","\n")
  print(data.frame(global.test$likelihood.ratio))
  
  cat("Information Criteria:","\n")
  print(data.frame(global.test$infocrit))
  
  cat("results of local tests:","\n","Type (+) / Antitype (-) based on:", type, "; Bonferoni adj. alpha:", object$bonferroni.alpha,"\n")
  return(templocal)
}