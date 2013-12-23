#' @title Configural Frequencies Analysis Main Function 
#' @export CFA
#' @exportClass CFA
#' @description Calculates various coefficients for the Configural Frequencies Analysis (CFA) defining main- and (optionaly) interaction effects. The core principle is to use \code{\link{glm}} in package \code{stats} to calculate the expected counts considering a designmatrix, which is constructed based on an formular definition given in argument \code{form}. 
#' @details This is the main function of the package. It internaly calls several functions of the package \code{\link{confreq}} which are also available as single functions. 
#' 
#' @param patternfreq an object of class \code{"Pfreq"}, which is data in pattern frequencies representation - see function \code{\link{dat2fre}}.
#' 
#' @param alpha a numeric giving the alpha level for testing (default set to \code{alpha=.05})
#'    
#' @param form a character expression which can be coerced into a model formulae with the function \code{as.formula} in the package \code{stats}. If this argument is left empty (at default \code{form=NULL}) the function \code{design_cfg_cfa()} will return a designmatrix coding only main effects and no interactions -- for a designmatrix refering to  three variables (V1, V2, V3) for example, leaving the argument \code{form} empty will be equivalent to assigning the character \code{"~ V1 + V2 + V3"} to the argument (\code{form="~ V1 + V2 + V3"}).
#'  A special Case is to define a null-model or rather a cfa model of order zero. In such a model no (main) effects are considered. This can be achieved bei passing the character expression \code{"null"} to the argument \code{form} -- so: \code{form = "null"} -- not to be confound with the default setting of this argument \code{form=NULL}. 
#' 
#' @param ccor either a logical (TRUE / FALSE) determining wether to apply a continuity correction or not. When set to \code{ccor=TRUE} continuity correction is applied for expected values 5 =< expected =< 10. For \code{ccor=FALSE} no continuity correction is applied. Another option is to set \code{ccor=c(x,y)} where x is the lower and y the upper bound for expected values where continuity correction is applied. So \code{ccor=c(5,10)} is equivalent to \code{ccor=TRUE}.
#' 
#' @param family argument passed to \code{\link{glm.fit}} with default set to \code{poisson()}
#' 
#' @param intercept argument passed to \code{\link{glm.fit}} with default set to \code{FALSE} 
#'  
#' @param ... additional parameters passed through to other functions.
#' @return an object of class \code{CFA} with results.
#' @references Lienert, G. A. (1971). Die Konfigurationsfrequenzanalyse: I. Ein neuer Weg zu Typen und Syndromen. \emph{Zeitschrift für Klinische Psychologie und Psychotherapie, 19}(2), 99-115. 
#' 
#' @examples #######################################
#' ######### some examples ########
#' data(LienertLSD)
#' LienertLSD
#' CFA(LienertLSD)
#' ## testing with (full) interactions
#' CFA(LienertLSD,form="~ C + T + A + C:T + C:A + T:A + C:T:A")
#' #' ## testing the null model
#' CFA(LienertLSD,form="null")
#' #######################
#' data(suicide)
#' suicide
#' # suicide data is in non tabulated data representation - so it must be tabulated !
#' CFA(dat2fre(suicide))  

############### start of function definition ##################
CFA<-function(patternfreq, alpha=.05, form=NULL, ccor=FALSE, family=poisson(), intercept=FALSE, ...){
  
  if(length(form)==0){form<-paste("~", paste(names(patternfreq)[1:(length(patternfreq)-1)],collapse=" + "))}

  pattern <- do.call(paste, patternfreq[,1:(dim(patternfreq)[2]-1) ])
  
  observed <- patternfreq[,dim(patternfreq)[2]] 
  
  designmatrix <- design_cfg_cfa(kat=sapply(lapply(patternfreq[,1:(dim(patternfreq)[2]-1) ],levels),length) , form = form, ...)
  
  expected <- expected_cfa(des=designmatrix, observed=observed, family=family, intercept=intercept, ...)
  
  erg <- data.frame(pattern=pattern, observed=observed,expected=expected,do.call(cbind,chi_local_test_cfa(observed,expected)),ex.bin.test=round(binomial_test_cfa(observed,expected),5), z_tests_cfa(observed,expected,ccor=ccor),p.stirling=stirling_cfa(observed=observed, expected=expected,cum=TRUE,verb=FALSE),density.stirling=stirling_cfa(observed=observed, expected=expected,cum=FALSE,verb=FALSE))
  
  chi.square <- sum(erg$loc.chi.square)
  
  df <- df_des_cfa(designmatrix)
  
  chi.square.p <- (1-pchisq(chi.square,df))
  
  bonferroni <- alpha/length(expected)
  
  result <- list( local.test = erg, bonferroni.alpha=bonferroni, global.test = list(chi.square=chi.square,df=df,chi.square.p=chi.square.p,alpha=alpha)) 
  
  class(result)<-c("CFA","list")
 
  return(result)
}
# End of function ------ 