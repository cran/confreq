#' @S3method plot CFA
#' @method plot CFA
#' @title S3 plot for CFA
#' @description S3 plot method for object of class\code{"CFA"}
#' @param x object of class \code{"CFA"}
#' @param type character indicating which test to use for visualizing whether the observed pattern are 'Types', 'Antitypes' or not significant at all. Possible options for \code{type} are \code{"pChi"}, \code{"ex.bin.test"}, \code{"z.pChi"}, \code{"z.pBin"} and \code{"p.stir"}. 
#' @param holm logical with default \code{holm = FALSE}. If set to \code{holm = TRUE}, significance testing is based on the holm procedure -- see references. This argument is deprecated (since version 1.5.6) and kept only for downward compatibility. Use argument \code{adjalpha} for any type of alpha adjustment.  
#' @param fill a vector of (three) colors defining the coloring of significant 'Types' (default "red"), 'Antitypes' (default "blue") or not significant cells (default "grey") in the plot.
#' @param adjalpha character with default \code{adjalpha = "bonferroni"}. Selector for the type of alpha adjustment for multiple testing. Possible options are: \code{adjalpha = "none"}, for no adjustment; \code{adjalpha = "bonferroni"}, for bonferroni adjustment (default); \code{adjalpha = "holm"}, for alpha adjustment according to Holm (1979); other options to come ... . 
#' @param ... other parameters passed trough
#' @references Holm, S. (1979). A simple sequentially rejective multiple test procedure. \emph{Scandinavian Journal of Statistics, 6}(2), 65–70.

#sorted by the p-value selected in argument \code{type}
########################### hier die plot method #class CFA #######################
plot.CFA<-function(x, type="z.pChi", holm=FALSE, fill=c("red", "blue", "grey"), adjalpha = "bonferroni", ...){
  local.test <- x$local.test
  global.test <- x$global.test
  functional <- x$functional
  inputdata <- x$inputdata
  
 if(class(functional)=="character"){ ## added 22-01-2019
   functional <- sapply(functional, function(x){which(x==as.character(x$local.test$pat.))})
 }
  out_wide <- c("Chi","df","pChi","ex.bin.test","z.Chi","z.pChi","z.Bin","z.pBin","cor.","p.stir","density.stir")

  ############## START of alpha adjustment methods -----------------------------
  
  if(holm==TRUE){adjalpha <- "holm"} # 30-04-2021 for downward compatibility
  
  #### significant (anti)types cf. Holm ----------------------------------------
  ## added 21-01-2019 #changed 30-04-2021
  if(adjalpha=="holm"){ # 30-04-2021 new type of alpha adjustment control
    p_val <- local.test[,which(names(local.test)==type)]
    if(!is.null(functional)){p_val[functional] <- NA}
    
    ord_p <- order(p_val,decreasing = FALSE) # starting with the smallest
    k <- length(p_val)
    n <- length(na.omit(p_val))
    alpha <- x$alpha # 30-04-2021 according to new structure of CFA result object
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
    templocal <- data.frame(pat.=local.test[,"pat."],local.test[,c("obs.","exp.")],Type, Holm.crit=alpha_bon ,local.test[,out_wide])  
  }
  
  ### significant (anti)types cf. (classical) Bonferroni adjustment ------------
  if(adjalpha=="bonferroni"){ ## condition added 21-01-2019 # 30-04-2021 new type of alpha adjustment control
    temp1 <- local.test[,which(names(local.test)==type)] < x$bonferroni.alpha
    temp2 <- ifelse(test=local.test$obs. > local.test$exp., yes="+",no="-") 
    Type <- mapply(FUN=function(x,y){ifelse(test=(x==TRUE),yes=y, no="." )   },x=temp1,y=temp2 )
    if(!is.null(functional)){Type[functional] <- "b"} # added 18-11-2016
    templocal <- data.frame(pat.=local.test[,"pat."],local.test[,c("obs.","exp.")],Type, local.test[,out_wide])
  }
  
  ### significant (anti)types without any adjustment for multiple testing ------
  if(adjalpha=="none"){ ## condition 30-04-2021
    temp1 <- local.test[,which(names(local.test)==type)] < x$alpha
    temp2 <- ifelse(test=local.test$obs. > local.test$exp., yes="+",no="-") 
    Type <- mapply(FUN=function(x,y){ifelse(test=(x==TRUE),yes=y, no="." )   },x=temp1,y=temp2 )
    if(!is.null(functional)){Type[functional] <- "b"} # added 18-11-2016
    templocal <- data.frame(pat.=local.test[,"pat."],local.test[,c("obs.","exp.")],Type, local.test[,out_wide])    
  }
  
  if(any(Type=="b")){ cat("\n","Type (b): blanked out (functional CFA)","\n")}#(18-11-2016)}
  
  
  
  
  inputdata_ <- inputdata
  
  
  
  
  
  ############## END of alpha adjustment methods -------------------------------
  
  ############################# old code as comment (30-04-2021):
    # # significant (anti)types cf. Holm 
  # ## added 21-01-2019
  # if(holm==TRUE){
  # p_val <- local.test[,which(names(local.test)==type)]
  #   if(!is.null(functional)){p_val[functional] <- NA}
  # 
  # ord_p <- order(p_val,decreasing = FALSE) # starting with the smallest
  # k <- length(p_val)
  # n <- length(na.omit(p_val))
  # alpha <- x$global.test$pearson$alpha
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
  # Type <- mapply(FUN=function(x,y){ifelse(test=(x==TRUE),yes=y, no="." )   },x=temp1,y=temp2 )
  # if(!is.null(functional)){Type[functional] <- "b"} # added 18-11-2016
  # #subsequent added on 8-2-2020
  # alpha_bon_temp1 <-(alpha/1:n)[(order(rev(ord_p[1:(length(ord_p)-length(functional))])))]
  # alpha_bon_temp2 <- p_val
  # alpha_bon_temp2[!is.na(p_val)] <- alpha_bon_temp1
  # alpha_bon <- alpha_bon_temp2
  # }
  # 
  # if(holm==FALSE){ ## condition added 21-01-2019
  # temp1 <- local.test[,which(names(local.test)==type)] < x$bonferroni.alpha
  # temp2 <- ifelse(test=local.test$obs. > local.test$exp., yes="+",no="-") 
  # Type <- mapply(FUN=function(x,y){ifelse(test=(x==TRUE),yes=y, no="." )   },x=temp1,y=temp2 )
  # 
  # if(!is.null(functional)){Type[functional] <- "b"} # added 18-11-2016
  # }
  # 
  # if(holm==FALSE){
  #   templocal <- data.frame(pat.=local.test[,"pat."],local.test[,c("obs.","exp.")],Type, local.test[,out_wide])
  # }
  # 
  # if(holm==TRUE){
  #   templocal <- data.frame(pat.=local.test[,"pat."],local.test[,c("obs.","exp.")],Type, Holm.crit=alpha_bon ,local.test[,out_wide])
  # }
  #   erg <- templocal
  #   inputdata_ <- inputdata
  #   if (holm==TRUE){
  #     sorton <- type
  #     sorter <- order(templocal[,which(names(templocal)==sorton)],decreasing=TRUE)
  #     erg <- templocal[sorter,]
  #     inputdata_ <- inputdata[sorter,]
  #   }
  # if(any(Type=="b")){ cat("\n","Type (b): blanked out (functional CFA)","\n")}#(18-11-2016)}
  
    
####### eigentlicher plotting teil ###################################
    Preq_tab <- fre2tab(inputdata_[,c((rev(colnames(inputdata_)))[-1], "Freq")])  
    fil <-sapply(1:length(Type), function(x){
      if(Type[x]=="+"){res <- fill[1]}
      if(Type[x]=="-"){res <- fill[2]}
      if(Type[x]=="."){res <- fill[3]}
      res
    })
    # fil
    gp <- gpar(col = NA, lty = "solid", fill = array(fil, dim = dim(Preq_tab),dimnames = dimnames(Preq_tab)))
    strucplot(Preq_tab,gp = gp)
}