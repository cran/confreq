#' @exportS3Method plot S2CFA
#' @keywords methods
#' @method plot S2CFA
#' @title S3 plot for S2CFA
#' @description S3 plot method for object of class\code{"S2CFA"}
#' @param x object of class\code{"S2CFA"}
#' @param type character with default \code{type="ex.fisher.test"}, to return wether the observed pattern are 'discriminating Types' or not significant at all based on the respective p-value. Another option for \code{type} is \code{type="pChi"}.     
#' @param fill a vector of (two) colors defining the coloring of discriminating 'Types' (default "red"), or not discriminating cells (default "grey") in the plot.
#' @param adjalpha character with default \code{adjalpha = "bonferroni"}. Selector for the type of alpha adjustment for multiple testing. Possible options are: \code{adjalpha = "none"}, for no adjustment; \code{adjalpha = "bonferroni"}, for bonferroni adjustment (default); \code{adjalpha = "holm"}, for alpha adjustment according to Holm (1979); other options to come ... . 
#' @param ... other parameters passed trough.
#' @return a plot visualizing the results.
#' @references Holm, S. (1979). A simple sequentially rejective multiple test procedure. \emph{Scandinavian Journal of Statistics, 6}(2), 65–70.
#' @references Bonferroni, C. E. (1935). Il calcolo delle assicurazioni su gruppi di teste. In S.O. Carboni (Ed.), \emph{Studi in Onore del Professore Salvatore Ortu Carboni} (S. 13–60). Roma, Tipografia del Senato: Bardi.

########################### hier die plot method #class S2CFA #######################
plot.S2CFA<-function(x, type="ex.fisher.test",  fill=c("red", "grey"), adjalpha = "bonferroni", ...){
  local.test <- x$local.test
  global.test <- x$global.test
  inputdata <- x$inputdata
  
  #### significant (anti)types cf. Holm ----------------------------------------
  ## added 29-06-2021 adapted from plot.CFA ........
  if(adjalpha=="holm"){ # 30-04-2021 new type of alpha adjustment control
    p_val <- local.test[,which(names(local.test)==type)]
    # if(!is.null(functional)){p_val[functional] <- NA} # not implemented (yet)
    ord_p <- order(p_val,decreasing = FALSE) # starting with the smallest
    k <- length(p_val)
    n <- length(na.omit(p_val))
    alpha <- x$alpha # 30-04-2021 according to new structure of S2CFA result object
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
    #templocal <- data.frame(pat.=local.test[,1], disc.Type, Holm.crit=round(alpha_bon,digits=digits), round(local.test[,2:3],digits=digits), round(local.test[,5:10],digits=digits) , check.names = FALSE)
  }
  
  ### significant (anti)types cf. (classical) Bonferroni adjustment ------------
  if(adjalpha=="bonferroni"){ ## condition added 21-01-2019 # 30-04-2021 new type of alpha adjustment control
    disc.Type <- ifelse(test = local.test[,which(names(local.test)==type)] < x$bonferroni.alpha, yes="+",no=".")
    #templocal <- data.frame(pat.=local.test[,1], disc.Type, round(local.test[,2:3],digits=digits), round(local.test[,5:10],digits=digits) , check.names = FALSE)
  }
  
  ### significant (anti)types without any adjustment for multiple testing ------
  if(adjalpha=="none"){ ## condition 30-04-2021
    disc.Type <- ifelse(test = local.test[,which(names(local.test)==type)] < x$alpha, yes="+",no=".")
    #templocal <- data.frame(pat.=local.test[,1], disc.Type, round(local.test[,2:3],digits=digits), round(local.test[,5:10],digits=digits) , check.names = FALSE)
  }
  
  inputdata_ <- inputdata
  Type <- disc.Type
  ############## END of alpha adjustment methods -------------------------------
  
  ####### eigentlicher plotting teil ###################################
  Preq_tab <- fre2tab(inputdata_[,c((rev(colnames(inputdata_)))[-1], "Freq")])  
  #Preq_tab <- fre2tab(inputdata_)  
  
  fil <-sapply(1:length(Type), function(y){
    if(Type[y]=="+"){res <- fill[1]}
    if(Type[y]=="."){res <- fill[2]}
    res
  })
  # fil
  gp <- gpar(col = NA, lty = "solid", fill = array(fil, dim = dim(Preq_tab),dimnames = dimnames(Preq_tab)))
  # strucplot(Preq_tab,gp = gp)
  doubledecker(Preq_tab,gp = gp,depvar = c((rev(colnames(inputdata_)))[-1], "Freq")[1])
  cat("\n", "Grouping by variable:", x$grouping$var, ", with categories:", x$grouping$cats, "\n", "pattern based on variables:", x$pat,"\n", "-----------------------") 
}