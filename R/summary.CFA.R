#' @S3method summary CFA
#' @method summary CFA
#' @title S3 Summary for CFA
#' @description S3 summary method for object of class\code{"CFA"}
#' @param object object of class\code{"CFA"}
#' @param digits integer rounds the values to the specified number of decimal places, default is \code{digits=3}.
#' @param type character indicating which test to use for inference wether the observed pattern are 'Types', 'Antitypes' or not significant at all. Possible options for \code{type} are \code{"pChi"}, \code{"ex.bin.test"}, \code{"z.pChi"}, \code{"z.pBin"} and \code{"p.stir"}. 
#' @param sorton sort results of local test by any column. by default the output is not sorted. Other options may be \code{"pat."}, \code{"obs."}, \code{"exp."}, \code{"Type"}, \code{"Chi"}, etc. ...
#' @param decreasing logical. Should the sort be increasing or decreasing? see \code{\link{order}}
#' @param showall logical with default \code{showall = TRUE}. To return only significant pattern (types / antitypes) set it to \code{showall = FALSE}.
#' @param holm logical with default \code{holm = FALSE}. If set to \code{holm = TRUE}, significance testing is based on the holm procedure -- see references. 
#' @param wide logical with default \code{wide = FALSE}. If set to \code{wide = TRUE}, results for all significance tests are returned. 
#' @param ... other parameters passed trough
#' @references Holm, S. (1979). A simple sequentially rejective multiple test procedure. \emph{Scandinavian Journal of Statistics, 6}(2), 65–70.


#sorted by the p-value selected in argument \code{type}
########################### hier die summary method #class CFA #######################
#digits=4; type="z.pChi";sorton=NULL; decreasing=FALSE; showall=TRUE; holm=TRUE; wide=FALSE

summary.CFA<-function(object, digits=3, type="z.pChi",sorton=NULL, decreasing=FALSE, showall=TRUE, holm=FALSE, wide=FALSE, ...){
  local.test <- object$local.test
  global.test <- object$global.test
  functional <- object$functional
  
 if(class(functional)=="character"){ ## added 22-01-2019
   functional <- sapply(functional, function(x){which(x==as.character(object$local.test$pat.))})
 }
  
  # interpreting option "wide" # option added 8-2-2020
  if (wide==FALSE){
    if(type=="pChi"){out_wide <- c("Chi","df","pChi")}
    if(type=="ex.bin.test"){out_wide <- c("df","ex.bin.test")}
    if(type=="z.pChi"){out_wide <- c("df","z.Chi","z.pChi")}
    if(type=="z.pBin"){out_wide <- c("df","z.Bin","z.pBin","cor.")}
    if(type=="p.stir"){out_wide <- c("df","p.stir","density.stir")}
  }else{out_wide <- c("Chi","df","pChi","ex.bin.test","z.Chi","z.pChi","z.Bin","z.pBin","cor.","p.stir","density.stir")}  
  
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
  if(i<n){ # hier nur bis n (NAs werden nicht getestet)
  while(p_val[ord_p[i+1]]  <  alpha/(n-(i))){
    holm_res[i+1] <- p_val[ord_p[i+1]]  <  alpha/(n-(i)) 
    i=i+1
    }
  }
  temp1 <- holm_res[order(ord_p)] # reorders it back 
  temp2 <- ifelse(test=local.test$obs. > local.test$exp., yes="+",no="-") 
  Type <- mapply(FUN=function(x,y){ifelse(test=(x==TRUE),yes=y, no="." )   },x=temp1,y=temp2 )
  if(!is.null(functional)){Type[functional] <- "b"} # added 18-11-2016
  #subsequent added on 8-2-2020
  alpha_bon_temp1 <-(alpha/1:n)[(order(rev(ord_p[1:(length(ord_p)-length(functional))])))]
  alpha_bon_temp2 <- p_val
  alpha_bon_temp2[!is.na(p_val)] <- alpha_bon_temp1
  alpha_bon <- alpha_bon_temp2
  # round(alpha_bon_temp2,4)
  # round(p_val,4)
  # #alpha_bon <-((c((alpha/1:n),rep(NA,times=sum(is.na(p_val)))))[(order(rev(ord_p)))]) #added 'rev' 8-2-2020
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
  #templocal <- data.frame(pat.=local.test[,"pat."],round(local.test[,c("obs.","exp.")],digits=digits),Type, round(local.test[,out_wide],digits=digits), cor.=local.test[,"cor."])
  templocal <- data.frame(pat.=local.test[,"pat."],local.test[,c("obs.","exp.")],Type, local.test[,out_wide])
  }
  
  if(holm==TRUE){
  #templocal <- data.frame(pat.=local.test[,"pat."],round(local.test[,c("obs.","exp.")],digits=digits),Type, Holm.crit=round(alpha_bon,digits = digits) ,round(local.test[,out_wide],digits=digits), cor.=local.test[,"cor."])
  templocal <- data.frame(pat.=local.test[,"pat."],local.test[,c("obs.","exp.")],Type, Holm.crit=alpha_bon ,local.test[,out_wide])
  }
  
  #print(templocal)
  # sort output by sorton ...
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
    if (holm==TRUE){
      sorton <- type
      sorter <- order(templocal[,which(names(templocal)==sorton)],decreasing=decreasing)
      erg <- templocal[sorter,]
    }#added 'rev' 8-2-2020
  }
  if (showall==FALSE){
    erg <- erg[erg$Type!=".",]
  }
  
  #added (generic overall) rounding here 8-2-2020
  #round(erg,digits = digits)
  round_index <- sapply((names(erg)[-1]), function(x){is.numeric(erg[,x])},USE.NAMES = FALSE)
  round_index <- c(FALSE,round_index)
  erg[,round_index] <- round(erg[,round_index],digits = digits)
  erg
  
  cat("function Call:","\n","-------------","\n","Formula:",paste(object$used.formula,collapse = " "),"\n","Variables:", names(object$variables),"\n","Categories:", object$variables,"\n")
  
  cat("\n","results of global tests:","\n", "-----------------------")
  cat("\n","pearson Chi-square test:","\n")
  print(data.frame(global.test$pearson))
  
  cat("\n","likelihood ratio test:","\n")
  print(data.frame(global.test$likelihood.ratio))
  
  cat("\n","Information Criteria:","\n")
  print(data.frame(global.test$infocrit))
  
  if(holm==FALSE){
    cat("\n","results of local tests:","\n", "-----------------------","\n","Type (+) / Antitype (-) based on:", type ,";","\n", "Bonferoni adj. alpha:", object$bonferroni.alpha,"\n")
  }
  if(holm==TRUE){
    cat("\n","results of local tests:","\n", "-----------------------","\n","Type (+) / Antitype (-) based on:", type, ";","\n","with Holm-Bonferroni adjustment - see \"Holm.crit\" \n")#added 21-01-2019
  }
  
  #if(any(temp3==FALSE)){ cat("\n","Type (b): blanked out (functional CFA)","\n")}#(28-04-2015)}
  if(any(Type=="b")){ cat("\n","Type (b): blanked out (functional CFA)","\n")}#(18-11-2016)}
  
  return(erg)
}