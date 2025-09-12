#' @exportS3Method summary S2CFA
#' @keywords methods
#' @method summary S2CFA
#' @title S3 Summary for S2CFA
#' @description S3 summary method for object of class\code{"S2CFA"}
#' @param object object of class\code{"S2CFA"}
#' @param digits integer rounds the values to the specified number of decimal places, default is \code{digits=3}.
#' @param type character with default \code{type="ex.fisher.test"}, to return wether the observed pattern are 'discriminating Types' or not significant at all based on the respective p-value. Another option for \code{type} is \code{type="pChi"}.     
#' @param sorton sort results of local test by any column. By default the output is not sorted. Other options may be \code{"pat."}, \code{"disc.Type"}, etc. ... So all column names that can potentially appear in the result.
#' @param decreasing logical. Should the sort be increasing or decreasing? see \code{\link{order}}
#' @param showall logical with default \code{showall = TRUE}. To return only significant pattern (discriminating types) set it to \code{showall = FALSE}.
#' @param adjalpha character with default \code{adjalpha = "bonferroni"}. Selector for the type of alpha adjustment for multiple testing. Possible options are: \code{adjalpha = "none"}, for no adjustment; \code{adjalpha = "bonferroni"}, for bonferroni adjustment (default); \code{adjalpha = "holm"}, for alpha adjustment according to Holm (1979); other options to come ... . 
#' @param ... other parameters passed trough.
#' @return a summary of the results printed on the console.
#' @references Holm, S. (1979). A simple sequentially rejective multiple test procedure. \emph{Scandinavian Journal of Statistics, 6}(2), 65–70.
#' @references Bonferroni, C. E. (1935). Il calcolo delle assicurazioni su gruppi di teste. In S.O. Carboni (Ed.), \emph{Studi in Onore del Professore Salvatore Ortu Carboni} (S. 13–60). Roma, Tipografia del Senato: Bardi.
########################### hier die summary method #class S2CFA #######################
summary.S2CFA<-function(object, digits=3, type="ex.fisher.test", sorton=NULL, decreasing=FALSE, showall=TRUE, adjalpha = "bonferroni", ...){
  #object, digits=3, type="z.pChi",sorton=NULL, decreasing=FALSE, showall=TRUE, holm=FALSE, wide=FALSE, adjalpha = "bonferroni", ...
  local.test <- object$local.test
  global.test <- object$global.test
  
  # object$bonferroni.alpha
  # disc.Type <- ifelse(local.test[,which(names(local.test)==type)] < object$bonferroni.alpha,yes="+", no=".")
  
  #### significant (anti)types cf. Holm ----------------------------------------
  ## added 29-06-2021 adapted from summary.CFA ........
  if(adjalpha=="holm"){ # 30-04-2021 new type of alpha adjustment control
    p_val <- local.test[,which(names(local.test)==type)]
    # if(!is.null(functional)){p_val[functional] <- NA} # not implemented (yet)
    ord_p <- order(p_val,decreasing = FALSE) # starting with the smallest
    k <- length(p_val)
    n <- length(na.omit(p_val))
    alpha <- object$alpha # 30-04-2021 according to new structure of S2CFA result object
    i=0
    holm_res <- vector(mode = "logical", length = k)
    if(i<n){ # hier nur bis n (NAs werden nicht getestet)
      while(p_val[ord_p[i+1]]  <  alpha/(n-(i))){
        holm_res[i+1] <- p_val[ord_p[i+1]]  <  alpha/(n-(i)) 
        i=i+1
      }
    }
    temp1 <- holm_res[order(ord_p)] # reorders it back 
    disc.Type <- ifelse(test = local.test[,which(names(local.test)==type)] < temp1,yes="+", no=".")
    #subsequent added on 8-2-2020
    alpha_bon_temp1 <-(alpha/1:n)[(order(rev(ord_p[1:(length(ord_p))])))]
    alpha_bon_temp2 <- p_val
    alpha_bon_temp2[!is.na(p_val)] <- alpha_bon_temp1
    alpha_bon <- alpha_bon_temp2
    templocal <- data.frame(pat.=local.test[,1], disc.Type, Holm.crit=round(alpha_bon,digits=digits), round(local.test[,2:3],digits=digits), round(local.test[,5:10],digits=digits) , check.names = FALSE)
  }
  
  ### significant (anti)types cf. (classical) Bonferroni adjustment ------------
  if(adjalpha=="bonferroni"){ ## condition added 21-01-2019 # 30-04-2021 new type of alpha adjustment control
    disc.Type <- ifelse(test = local.test[,which(names(local.test)==type)] < object$bonferroni.alpha, yes="+",no=".")
    templocal <- data.frame(pat.=local.test[,1], disc.Type, round(local.test[,2:3],digits=digits), round(local.test[,5:10],digits=digits) , check.names = FALSE)
  }
  
  ### significant (anti)types without any adjustment for multiple testing ------
  if(adjalpha=="none"){ ## condition 30-04-2021
    disc.Type <- ifelse(test = local.test[,which(names(local.test)==type)] < object$alpha, yes="+",no=".")
    templocal <- data.frame(pat.=local.test[,1], disc.Type, round(local.test[,2:3],digits=digits), round(local.test[,5:10],digits=digits) , check.names = FALSE)
  }
  
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
    erg <- erg[erg$disc.Type!=".",]
  }
  
  
  # cat("results of global tests:","currently not implemented !","\n")
  cat("\n", "Grouping by variable:", object$grouping$var, ", with categories:", object$grouping$cats, "\n", "pattern based on variables:", object$pat,"\n", "-----------------------") 
  
  # adapted from summary.CFA #29-06-2021 
  if(adjalpha=="bonferroni"){
    cat("\n","results of local tests:","\n", "-----------------------","\n","discriminating Type (+) / not discriminating Type (.) based on:", type,";","\n", "with Bonferroni adjusted alpha:", object$bonferroni.alpha,"\n")#changed 30-04-2021
  }
  if(adjalpha=="holm"){
    cat("\n","results of local tests:","\n", "-----------------------","\n","discriminating Type (+) / not discriminating Type (.) based on:", type, ";","\n","with Holm adjusted alpha - see \"Holm.crit\" \n")#added 21-01-2019 #changed 30-04-2021
  }
  if(adjalpha=="none"){
    cat("\n","results of local tests:","\n", "-----------------------","\n","discriminating Type (+) / not discriminating Type (.) based on:", type, ";","\n","with not adjusted alpha:", object$alpha,"\n")#added 30-04-2021
  }
  
  return(erg)
}