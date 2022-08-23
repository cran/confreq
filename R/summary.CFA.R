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
#' @param adjalpha character with default \code{adjalpha = "bonferroni"}. Selector for the type of alpha adjustment for multiple testing. Possible options are: \code{adjalpha = "none"}, for no adjustment; \code{adjalpha = "bonferroni"}, for bonferroni adjustment (default); \code{adjalpha = "holm"}, for alpha adjustment according to Holm (1979); other options to come ... . 
#' @param ... other parameters passed trough.
#' @return a summary of the results printed on the console.
#' @references Holm, S. (1979). A simple sequentially rejective multiple test procedure. \emph{Scandinavian Journal of Statistics, 6}(2), 65–70.
#' @references Bonferroni, C. E. (1935). Il calcolo delle assicurazioni su gruppi di teste. In S.O. Carboni (Ed.), \emph{Studi in Onore del Professore Salvatore Ortu Carboni} (S. 13–60). Roma, Tipografia del Senato: Bardi.


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
    if(type=="pChi"){out_wide <- c("Chi","df","pChi")}
    if(type=="ex.bin.test"){out_wide <- c("df","ex.bin.test")}
    if(type=="z.pChi"){out_wide <- c("df","z.Chi","z.pChi")}
    if(type=="z.pBin"){out_wide <- c("df","z.Bin","z.pBin","cor.")}
    if(type=="p.stir"){out_wide <- c("df","p.stir","density.stir")}
  }else{out_wide <- c("Chi","df","pChi","ex.bin.test","z.Chi","z.pChi","z.Bin","z.pBin","cor.","p.stir","density.stir")}  
  
  ############## START of alpha adjustment methods -----------------------------
  
  if(holm==TRUE){adjalpha <- "holm"} # 30-04-2021 for downward compatibility
  
  #### significant (anti)types cf. Holm ----------------------------------------
  ## added 21-01-2019 #changed 30-04-2021
  if(adjalpha=="holm"){ # 30-04-2021 new type of alpha adjustment control
  p_val <- local.test[,which(names(local.test)==type)]
    if(!is.null(functional)){p_val[functional] <- NA}# set all functional cells to NA
  
  ord_p <- order(p_val,decreasing = FALSE) # starting with the smallest 'p_val'
  k <- length(p_val)
  n <- length(na.omit(p_val))
  alpha <- object$alpha # 30-04-2021 according to new structure of CFA result object
  # did this all new on 22-08-2022...........
  # i=0
  # holm_res <- vector(mode = "logical", length = k)
  # if(i<n){ # hier nur bis n (NAs werden nicht getestet)
  # while(p_val[ord_p[i+1]]  <  alpha/(n-(i))){
  #   holm_res[i+1] <- p_val[ord_p[i+1]]  <  alpha/(n-(i)) 
  #   i=i+1
  #   }
  # }
  # temp1 <- holm_res[order(ord_p)] # reorders it back 
  # temp2 <- ifelse(test=local.test$obs. > local.test$exp., yes="+",no="-") 
  # Type <- mapply(FUN=function(xx,y){ifelse(test=(xx==TRUE),yes=y, no="." )   },xx=temp1,y=temp2 )
  # if(!is.null(functional)){Type[functional] <- "b"} # added 18-11-2016
  #subsequent added on 8-2-2020
  alpha_bon_temp1 <-(alpha/1:n)[(order(rev(ord_p[1:(length(ord_p)-length(functional))])))]
  alpha_bon_temp2 <- p_val
  alpha_bon_temp2[!is.na(p_val)] <- alpha_bon_temp1
  alpha_bon <- alpha_bon_temp2
  # round(alpha_bon_temp2,4)
  # round(p_val,4)
  # #alpha_bon <-((c((alpha/1:n),rep(NA,times=sum(is.na(p_val)))))[(order(rev(ord_p)))]) #added 'rev' 8-2-2020
  # alpha_bon[is.na(alpha_bon)] <-0 
  # subsequent new added 22-08-2022 --------
   temp1 <- p_val < alpha_bon 
   temp2 <- ifelse(test=local.test$obs. > local.test$exp., yes="+",no="-") 
   Type <- mapply(FUN=function(xx,y){ifelse(test=(xx==TRUE),yes=y, no="." )   },xx=temp1,y=temp2 )
   if(!is.null(functional)){Type[functional] <- "b"} # added 18-11-2016
  # and here we go ...
  templocal <- data.frame(pat.=local.test[,"pat."],local.test[,c("obs.","exp.")],Type, Holm.crit=alpha_bon ,local.test[,out_wide])  
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
  
  ### significant (anti)types according to Holland and DiPonzio Copenhaver’s procedure any adjustment for multiple testing ------
### das ist eine riesen B A U S T E L L E
  
    # if(adjalpha=="holland"){ ## condition 30-04-2021
  #   
  #   setUpperBound <- function (vector, bound){
  #     # Simple function to limit the maximum value of a vector to a given value
  #     # Args:
  #     #   vector: Vector whose values will be bounded
  #     #   bound:  Maximum value in the result
  #     # Returns:
  #     #   A vector equal to 'vector' except for those values above 'bound', that are
  #     #   set to this upper bound
  #     res <- sapply(vector, 
  #                   FUN=function(x) {
  #                     return(min(bound, x))
  #                   })
  #     return(res)
  #   }
  #   correctForMonotonicity <- function (pvalues){
  #     # Function to ensure the monotonicity of the p-values after being corrected
  #     # Args:
  #     #   pvalues: Corrected p-values, in DECRESING ORDER ACCORDING TO THE RAW PVALUE
  #     # Returns:
  #     #   The corrected p-values such as there is not a p-value smaller than the any 
  #     #   of the previous ones.
  #     #
  #     pvalues <- sapply(1:length(pvalues), 
  #                       FUN=function(x) {
  #                         return(max(pvalues[1:x]))
  #                       })
  #     return(pvalues)
  #   }
  #   #################### end helper functions -------------
  #   p_val <- local.test[,which(names(local.test)==type)]
  #   names(p_val) <- local.test$pat. ## using names to control ordering
  #   if(!is.null(functional)){p_val[functional] <- NA}
  #   ord_p <- order(p_val, na.last=NA) # not sure about na.last option ... 'ord'
  #   pvalues.sorted <- p_val[ord_p]
  #   k <- length(pvalues.sorted) + 1
  #   
  #   p.val.aux <- sapply(1:(k - 1), 
  #                       FUN=function(j, p.val, k){
  #                         r <- 1 - (1 - p.val[j])^(k - j)
  #                         return(r)
  #                       },
  #                       p.val=pvalues.sorted, k=k)
  #   
  #   p.adj.aux <- setUpperBound(p.val.aux, 1)
  #   p.adj.aux <- correctForMonotonicity(p.adj.aux)
  #   
  #   p.adj <- rep(NA, length(p_val))
  #   suppressWarnings(expr = {
  #     p.adj[ord_p] <- p.adj.aux
  #   })
  # 
  #   
  # 
  #   }
  # 
  
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
    erg <- erg[erg$Type!=".",]
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
    cat("\n","results of local tests:",attributes(object)$comment,"\n", "-----------------------","\n","Type (+) / Antitype (-) based on:", type, ";","\n","with Holm adjusted alpha - see \"Holm.crit\" \n")#added 21-01-2019 #changed 30-04-2021
  }
  if(adjalpha=="none"){
    cat("\n","results of local tests:",attributes(object)$comment,"\n", "-----------------------","\n","Type (+) / Antitype (-) based on:", type, ";","\n","with not adjusted alpha:", object$alpha,"\n")#added 30-04-2021
  }
  
  #if(any(temp3==FALSE)){ cat("\n","Type (b): blanked out (functional CFA)","\n")}#(28-04-2015)}
  if(any(Type=="b")){ cat("\n","Type (b): blanked out (functional CFA)","\n")}#(18-11-2016)}
  
  return(erg)
}
