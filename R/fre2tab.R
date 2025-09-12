#' @title pattern frequency to table conversion 
#' @keywords utilities 
#' @export fre2tab
#' @description Given data as pattern frequencies (object of class class c("data.frame","Pfreq", see function \code{\link{dat2fre}}) this function returns a typical array representation (class "table" , see \code{\link{table}}) of it.
#' @details  This function was introduced in order to connect the typical confreq data representation in the objects of the class \code{c("data.frame","Pfreq")}, see function \code{\link{dat2fre}}, to the R-typical array representation as it exists in objects of the \code{"table"} class, see \code{\link{table}}. This array representation of multi-dimensional contingency tables is used more universally in R -- e.g. also in the R package \code{vcd}, see the examples section below.
#' 
#' It is assumed, that the last column of the object \code{patternfreq} represents the frequencies of the (response) pattern represented by the other columns in \code{patternfreq}.
#' 
#' @param patternfreq an object of class c("data.frame","Pfreq") 
#' @param form a formula object with possibly both left and right hand sides specifying the order of the variables in the resulting table. At default (\code{(formula=NULL)}) all variables in \code{(x)} are used in their respective order.
#' @return An object of class "table" see \code{\link{table}}. 
#' @references No references at the moment 
#' @examples #######################################
#' data(LienertLSD)# loading example pattern frequencies table
#' fre2tab(LienertLSD)# coverting it into a table
#' 
#' ### examples using functions from package vcd
#' data(Lienert1978)# loading example pattern frequencies table
#' fre2tab(Lienert1978)# coverting it into a table
#' strucplot(fre2tab(Lienert1978))# plotting data with 'vcd'
#' structable(fre2tab(Lienert1978),direction = "v")# flatten table (vertical) with 'vcd'
#' 
#' # changing the vertical grouping when flattening the table by unsing a 'formula':
#' structable(fre2tab(Lienert1978, form=~Group + Student + Teacher),direction = "v")# flatten table
############### start of function definition ##################
########### pattern frequency to tabe conversion ##############
################ jhheine at googlemail.com ####################
fre2tab <- function(patternfreq, form = NULL) {
  if(any(class(patternfreq)=="Pfreq") != TRUE){stop("patternfreq must be an object of class 'Pfreq'","\n","see func. dat2fre()", call. = TRUE) }
  if (is.null(form)){
    form <- as.formula(paste("~", paste(colnames(patternfreq)[-ncol(patternfreq)],collapse=" + ")))
  }
res <- as.table(structable(formula = form, data = patternfreq))

return(res)
}
