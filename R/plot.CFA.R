#' @exportS3Method plot CFA
#' @keywords methods
#' @method plot CFA
#' @title S3 plot for CFA
#' @description S3 plot method for object of class\code{"CFA"}
#' @param x object of class \code{"CFA"}
#' @param type character indicating which test to use for visualizing whether the observed pattern are 'Types', 'Antitypes' or not significant at all. Possible options for \code{type} are \code{"pChi"}, \code{"ex.bin.test"}, \code{"z.pChi"}, \code{"z.pBin"} and \code{"p.stir"}. 
#' @param fill a vector of (three) colors defining the coloring of significant 'Types' (default "red"), 'Antitypes' (default "blue") or not significant cells (default "grey") in the plot.
#' @param adjalpha character with default \code{adjalpha = "bonferroni"}. Selector for the type of alpha adjustment for multiple testing. Possible options are: \code{adjalpha = "none"}, for no adjustment; \code{adjalpha = "bonferroni"}, for bonferroni adjustment (default); \code{adjalpha = "holm"}, for alpha adjustment according to Holm (1979); other options to come ... . 
#' @param ... other parameters passed trough
#' @return a plot visualizing the results.
#' @references Holm, S. (1979). A simple sequentially rejective multiple test procedure. \emph{Scandinavian Journal of Statistics, 6}(2), 65–70.
#' @references Bonferroni, C. E. (1935). Il calcolo delle assicurazioni su gruppi di teste. In S.O. Carboni (Ed.), \emph{Studi in Onore del Professore Salvatore Ortu Carboni} (S. 13–60). Roma, Tipografia del Senato: Bardi.

########################### hier die plot method #class CFA #######################
plot.CFA<-function(x, type="z.pChi", fill=c("red", "blue", "grey"), adjalpha = "bonferroni", ...){
  local.test <- x$local.test
  global.test <- x$global.test
  functional <- x$functional
  inputdata <- x$inputdata
 
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
  if(adjalpha=="holm"){ # 30-04-2021 new type of alpha adjustment control
    p_val <- local.test[,which(names(local.test)==type)]
    if(!is.null(functional)){p_val[functional] <- NA}# set all functional cells to NA
    
    ord_p <- order(p_val,decreasing = FALSE) # starting with the smallest 'p_val'
    k <- length(p_val)
    n <- length(na.omit(p_val))
    alpha <- x$alpha # 30-04-2021 according to new structure of CFA result object
   
    alpha_bon_temp1 <-(alpha/1:n)[(order(rev(ord_p[1:(length(ord_p)-length(functional))])))]
    alpha_bon_temp2 <- p_val
    alpha_bon_temp2[!is.na(p_val)] <- alpha_bon_temp1
    alpha_bon <- alpha_bon_temp2
    
    temp1 <- p_val < alpha_bon 
    temp2 <- ifelse(test=local.test$obs. > local.test$exp., yes="+",no="-") 
    Type <- mapply(FUN=function(xx,y){ifelse(test=(xx==TRUE),yes=y, no="." )   },xx=temp1,y=temp2 )
    if(!is.null(functional)){Type[functional] <- "b"} # added 18-11-2016

    templocal <- data.frame(pat.=local.test[,"pat."],local.test[,c("obs.","exp.")],Type, Holm.crit=alpha_bon ,local.test[,out_wide])  
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
  
  if(any(Type=="b")){ cat("\n","Type (b): blanked out (functional CFA)","\n")}#(18-11-2016)}
  
  cat(attributes(x)$comment) # version 1.6
  
  inputdata_ <- inputdata

  ############## END of alpha adjustment methods -------------------------------
  
####### eigentlicher plotting teil ###################################
    Preq_tab <- fre2tab(inputdata_[,c((rev(colnames(inputdata_)))[-1], "Freq")])  
    fil <-sapply(1:length(Type), function(y){
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