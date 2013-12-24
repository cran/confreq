#' @title Binomial Test 
#' @export binomial_test_cfa
#' @description Calculates the (exact) binomial test based on obseved, expected frequencies an the total number of observations. 
#' #' @details No details
#' 
#' @param expected a vector giving the expected frequencies.
#' @param observed a vector giving the observed frequencies.
#' @param ntotal optional a numeric giving the total number of observations. By default ntotal is calculated as \code{ntotal=sum(observed)}. 
#' @return a numeric giving the p-value.
#' @references No references in the moment
#' @examples #######################################
#' # first calculate expected counts for LienertLSD data example.
#' designmatrix<-design_cfg_cfa(kat=c(2,2,2)) # generate an designmatrix (only main effects)
#' data(LienertLSD) # load example data
#' observed<-LienertLSD[,4] # extract observed counts
#' expected<-expected_cfa(des=designmatrix, observed=observed) # calculation of expected counts
#'  binomial_test_cfa(observed,expected)
#' ####################################### 

############### start of function definition ##################
binomial_test_cfa<-function(observed,expected,ntotal=sum(observed)){
# func. by joerg-henrik heine jhheine(at)googlemail.com  
###############################################################
# helperfunction -- better than choose() !
bin.coeff <- function(n, k) exp(lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1))
# ENDE hilfsfunktion -- genauer als choose() !

exact<-function(obsi,expi,ntotal){ 
p_expi <- expi / ntotal   
sum.p <- 0

if(obsi>=expi){
 for (j in obsi:ntotal) {
    sum.p <- sum.p + ((p_expi^j) * (1 - p_expi)^(ntotal - j) * bin.coeff(ntotal, j))
  }
 exact = min(sum.p, 1-sum.p) 
}

if(obsi < expi){
  for (j in 0:obsi) {
    sum.p <- sum.p + ((p_expi^j) * (1 - p_expi)^(ntotal - j) * bin.coeff(ntotal, j))
  }
  exact = min(sum.p, 1-sum.p) 
}

return(exact)
}
#die eingegbenen argumente als matrix:
m<-cbind(observed,expected,rep(ntotal,length(observed)))
#rechnen des tests:
p.exact.binomial<-apply(m,1,function(x){exact(x[1],x[2],x[3]) })
####
#cat("p-value for exact binomial test:")
return(p.exact.binomial)
}
