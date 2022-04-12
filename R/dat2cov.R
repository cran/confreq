#' @title conversion of a covariate dataset into summary covariate values
#' @keywords utilities 
#' @export dat2cov
#' @description Given a dataset \code{x}, this function returns summary values for some (\code{numeric}) covariate variables in \code{x} for each pattern (configuration) defined by a set of factor variables in \code{x}. 
#' @details No further details
#' 
#' @param x an object of class "data.frame" with at least 2 \code{factor} variables representing the pattern (configurations) and at least 1 \code{numeric} variable representing the covariate(s).
#' @param FUN a function to compute the summary statistics which can be applied to all covariate variables in \code{x}. See function \code{\link[stats]{aggregate}}.
#' @param ... further arguments passed to or used by methods un \code{FUN}.
#' @param notobs a numeric vector possibly with length equal to the number of \code{numeric} variables in \code{x}, defining the summary value for the respective covariate variable to use for unobserved pattern (configurations) defined by the \code{factor} variables in \code{x}. By default it is assumend that this value is 0. \code{notobs} is recycled if only one value is given.
#' @param katorder see \code{\link{dat2fre}}
#' @param caseorder see \code{\link{dat2fre}} 
#' @return An object of class \code{c("data.frame", "Pcov")} holding the summary statistics for the covariate variables corresponding to the pattern (configurations) of the given dataset in the argument \code{x}. 

############### start of function definition ##################
################ jhheine at googlemail.com ####################
dat2cov <- function(x, FUN = "mean", ... , notobs = 0, katorder = FALSE, caseorder = TRUE ){
  d <- x
  p_ind <- (sapply(d,class))=="factor"
  cov_ind <- (sapply(d,class))=="numeric"
  nam <- names(d)[cov_ind]
  
  ta <- dat2fre(d[,p_ind],katorder = katorder,caseorder = caseorder)
  fre <- ta[,ncol(ta)]
  pat <- apply(X =ta[,1:(ncol(ta)-1)],1,paste, collapse = " ")     
  
  dp <-  data.frame(pat = apply(X = d[,p_ind],1,paste, collapse = " "), d[,cov_ind])
  da <- aggregate(x = dp[,-1], by=list(dp[,"pat"]),FUN = FUN, ...)
  # db <- data.frame(pat=levels(da[,1]),da[,-1],stringsAsFactors =FALSE) # new AKA2021
  db <- data.frame(pat=da[,1],da[,-1],stringsAsFactors =FALSE) 
  names(db)[-1] <- nam
  row.names(db) <-db[,1] 
  # matrix(db[pat,-1],ncol=sum(cov_ind))
  dc <- data.frame(db[pat,-1])
  names(dc) <- nam
  row.names(dc) <- pat
  
  if(!is.na(notobs)){ 
    if(length(notobs)==1){notobs <- rep(notobs,times=sum(cov_ind))}
    for(i in 1:sum(cov_ind)){
      dc[,i][is.na(dc[,i])] <- notobs[i]  
    }
    }
  
  class(dc) <- c("data.frame", "Pcov")
  return(dc)
}