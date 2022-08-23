#' @title Configural Frequencies Analysis Main Function 
#' @keywords mainfunction
#' @export CFA
#' @description Calculates various coefficients for the Configural Frequencies Analysis (CFA) defining main- and (optional) interaction effects. The core principle is to use \code{\link{glm}} in package \code{stats} to calculate the expected counts considering a designmatrix, which is constructed based on an formula definition given in argument \code{form}. 
#' @details This is the main function of the package. It internal calls several functions of the package \code{\link{confreq-package}} which are also available as single functions. For classification of the observed patterns into 'Types' and 'Antitypes' according to Lienert  (1971), a S3 summary method for the resulting object of class \code{"CFA"} can be applied - see \code{\link{summary.CFA}}. An S3 plot method is useful for visualization of the contingency table and the 'Types' and 'Antitypes' -- see \code{\link{plot.CFA}}. Since version  1.6.0-1 of \code{confreq} survey weights are supported when tabluating a data set with function \code{\link{dat2fre}}. In case that for the resulting tabulated data in the object of class \code{c("data.frame","Pfreq")} survey weights were used the function \code{CFA} will take into account those weigts for estimation of the expected counts -- currently only when \code{method="log"}.  
#' 
#' @param patternfreq an object of class \code{"Pfreq"}, which is data in pattern frequencies representation - see function \code{\link{dat2fre}}.
#' 
#' @param alpha a numeric giving the alpha level for testing (default set to \code{alpha=.05})
#'    
#' @param form either a character expression which can be coerced into a model formula with the function \code{as.formula} in the package \code{stats}. If this argument is left empty (at default \code{form=NULL}) the (internal) function \code{design_cfg_cfa()} will return a designmatrix coding only main effects and no interactions -- for a designmatrix referring to  three variables (V1, V2, V3) for example, leaving the argument \code{form} empty will be equivalent to assigning the character \code{"~ V1 + V2 + V3"} to the argument (\code{form="~ V1 + V2 + V3"}).
#'  A special case is to define a null-model or rather a cfa model of order zero. In such a model no (main) effects are considered. This can be achieved bei passing the character expression \code{"null"} to the argument \code{form} -- so: \code{form = "null"} -- not to be confound with the default setting of this argument \code{form=NULL}. Another option is to define your own designmatrix and assign it to this argument (\code{form}) in this case the object assigned to \code{form} must be of class \code{"matrix"} and must logical match to the argument \code{patternfreq}, which is currently not checked! - but simply assumed.   
#' 
#' @param ccor either a logical (TRUE / FALSE) determining whether to apply a continuity correction or not for the Binomial Approximation to the z-Test. When set to \code{ccor=TRUE} continuity correction is applied for expected values 5 =< expected =< 10. For \code{ccor=FALSE} no continuity correction is applied. Another option is to set \code{ccor=c(x,y)} where x is the lower and y the upper bound for expected values where continuity correction is applied. So \code{ccor=c(5,10)} is equivalent to \code{ccor=TRUE}.
#' 
#' @param family argument passed to \code{\link{glm.fit}} with default set to \code{poisson()}
#' 
#' @param intercept argument passed to \code{\link{glm.fit}} with default set to \code{FALSE} 
#' 
#' @param method character defining the estimation method for expected frequencies with default set to \code{method="log"} to estimate the expected frequencies using \code{\link{glm}}. An other option is to set this argument to \code{method="margins"} which will result in expected frequencies calculated based on the margins of the multidimensional contingency table. Only main effects models are possible in this case and thus the arguments \code{form}, \code{family} \code{cova} and \code{intercept} are ignored.
#' 
#' @param blank can be used to indicate which pattern (configurations) are declared as structural cells (configurations) for functional CFA. Should be either (1) a character vector defining the pattern (with spaces between variable categories), which will be ignored for calculation of expected frequencies; or (2) a numeric vector defining the (row) positions of the pattern in an object of class \code{"Pfreq"} (see. argument \code{patternfreq}), which will be ignored for calculation of expected frequencies. At default (\code{blank=NULL}) all possible pattern, as listed in object of class \code{"Pfreq"}, are included for calculation of expected frequencies.  
#' 
#' @param cova a matrix (possibly with one or more columns) holding the covariate (mean) values for each pattern (configurations) see function \code{\link{dat2cov}}.
#' @param bintest a logical with default set to \code{bintest=TRUE}; if set to \code{bintest=FALSE} no calculations for the exact binomial test are performed, which can reduce processing time in some cases dramatically.
#' @param ... additional parameters passed through to other functions.
#' @return an object of class \code{CFA} with results.
#' @references Lienert, G. A. (1971). Die Konfigurationsfrequenzanalyse: I. Ein neuer Weg zu Typen und Syndromen. \emph{Zeitschrift für Klinische Psychologie und Psychotherapie, 19}(2), 99-115. 
#' @references Glück, J., & Von Eye, A. (2000). Including covariates in Configural Frequency Analysis. \emph{Psychologische Beitrage, 42}, 405–417.
#' @references Victor, N. (1989). An Alternativ Approach to Configural Frequency Analysis. \emph{Methodika, 3}, 61–73.
#' @references Stemmler, M. (2020). \emph{Person-Centered Methods}. Cham: Springer International Publishing.
#' @examples #######################################
#' ######### some examples ########
#' data(LienertLSD)
#' LienertLSD
#' res1 <- CFA(LienertLSD)
#' summary(res1)
#' ## testing with (full) interactions
#' res2 <- CFA(LienertLSD,form="~ C + T + A + C:T + C:A + T:A + C:T:A")
#' summary(res2)
#' #' ## testing the null model
#' res3 <- CFA(LienertLSD,form="null")
#' summary(res3)
#' #######################
#' data(suicide)
#' suicide
#' # suicide data is in non tabulated data representation - so it must be tabulated !
#' res4 <- CFA(dat2fre(suicide))  
#' summary(res4)
###############################################################
############### start of function definition ##################
CFA<-function(patternfreq, alpha=.05, form=NULL, ccor=FALSE, family=poisson(), intercept=FALSE, method="log", blank=NULL, cova=NULL,bintest=TRUE, ...){
  
if(any(class(patternfreq)=="Pfreq") != TRUE){stop("patternfreq must be an object of class 'Pfreq'","\n","see func. dat2fre()", call. = TRUE) }

# kategorie <- sapply(lapply(patternfreq[,1:(dim(patternfreq)[2]-1) ],levels),length)
# JHH: changed this to using 'unique' due to incident reported by M.Stemmler 15-10-2020 with R version 4.0.3  
kategorie <- sapply(lapply(patternfreq[,1:(dim(patternfreq)[2]-1) ],unique),length)

pattern <- do.call(paste, patternfreq[,1:(dim(patternfreq)[2]-1) ])
  
observed <- patternfreq[,dim(patternfreq)[2]] 

# version 1.6 weights
if(!is.null(attributes(patternfreq)$comment)){
if(attributes(patternfreq)$comment=="using weighted data!"){
 WGT <-  attributes(patternfreq)$WGT
}
}else{WGT <- NULL}
# condition added 22-06-2015 fixed class issue 2.2.2020
if(method=="log"){
  if(   !is(object = form,class2 = "matrix")   ){
    if(length(form)==0){form<-paste("~", paste(names(patternfreq)[1:(length(patternfreq)-1)],collapse=" + "))}
    
    designmatrix <- design_cfg_cfa(kat=kategorie , form = form, ...) 
    
    usedform <- form
    
    #option added 18-11-2016 
    if(!is.null(blank)){
      if(is.numeric(blank)){  #class(blank)!="character" , '!inherits(x = class(blank),what = "character")' again changed 19-08-2022
        blank_which <- blank
        Mi <- sapply(blank_which,function(x){b <- (rep(0,length.out=nrow(designmatrix))); b[x] <- 1; b})
        designmatrix <- cbind(designmatrix, Mi)
        usedform <- paste(paste(form),"; and functional pattern:", paste(pattern[blank],collapse = ", "))
      }
      if(is.character(blank)){ #class(blank)=="character", 'inherits(x = class(blank),what = "character"' again changed 19-08-2022
        blank_which <- sapply(blank, function(x){which(x==pattern)})
        Mi <- sapply(blank_which,function(x){b <- (rep(0,length.out=nrow(designmatrix))); b[x] <- 1; b})
        designmatrix <- cbind(designmatrix, Mi)
        usedform <- paste(paste(form),"; and functional pattern:", paste(blank,collapse = ", "))
      }  
    }else{blank_which <- blank}
    #option added 18-11-2016
    if(!is.null(cova)){
      designmatrix <- cbind(designmatrix, cova)
      usedform <- paste(usedform, "; and",ncol(cova) , "covariates") 
      }
  }
  
  if(is(object = form,class2 = "matrix")){
    designmatrix <- form
    usedform <- "designmatrix"
    #option added 08-01-2018 
    if(!is.null(blank)){
      if(!is.character(blank)){# '!inherits(x = class(blank),what = "character")' again changed 29-08-2022
        blank_which <- blank
        Mi <- sapply(blank_which,function(x){b <- (rep(0,length.out=nrow(designmatrix))); b[x] <- 1; b})
        designmatrix <- cbind(designmatrix, Mi)
        usedform <- paste(paste(usedform),"; and functional pattern:", paste(pattern[blank],collapse = ", "))
      } 
      if(is.character(blank)){# 'inherits(x = class(blank),what = "character")' again changed 29-08-2022
        blank_which <- sapply(blank, function(x){which(x==pattern)})
        Mi <- sapply(blank_which,function(x){b <- (rep(0,length.out=nrow(designmatrix))); b[x] <- 1; b})
        designmatrix <- cbind(designmatrix, Mi)
        usedform <- paste(paste(form),"; and functional pattern:", paste(blank,collapse = ", "))
      }  
    }else{blank_which <- blank}  
       
    #option added 11-05-2018
    if(!is.null(cova)){
      designmatrix <- cbind(designmatrix, cova)
      usedform <- paste(usedform, "; and",ncol(cova) , "covariates") 
    }
    } # !!! no further checks !!!
  
  # expected <- expected_cfa(des=designmatrix, observed=observed, family=family, intercept=intercept, ...) # padded out 24.10.2014
  glmfitres <- glm.fit(x=designmatrix, y=observed ,family=family, intercept = intercept,weights = WGT) #, ... added 24.10.2014
  expected <- glmfitres$fitted.value # added 24.10.2014
  #aic <- glmfitres$aic # added 24.10.2014
  class(glmfitres) <- c("glm", "lm" )# this is a trick!!! necessary for the following three lines of code - added 24.10.2014
  loglik <- logLik(glmfitres)# cheked against code below OK! added 24.10.2014
  aic <- AIC(glmfitres)# cheked against code below OK! added 24.10.2014
  bic <- BIC(glmfitres)# cheked against code below OK! added 24.10.2014
  # glmres <- glm(formula = paste("Freq",form) , data=patternfreq, , family = family) # added 24.10.2014
  # expected <- fitted(glmres) # OK
  # # resid(glmres)
  # # predict(glmres)
  # loglik <- logLik(glmres)
  # aic <- AIC(glmres)
  # bic <- BIC(glmres)
}

if(method=="margins"){
  expected <- expected_margin_cfa(Pfreq = patternfreq, blank = blank) # added 22-06-2015
  loglik <- NA # added 22-06-2015
  aic <- NA # added 22-06-2015
  bic <- NA # added 22-06-2015
  # nur für df
  form<-paste("~", paste(names(patternfreq)[1:(length(patternfreq)-1)],collapse=" + "))
  designmatrix <- design_cfg_cfa(kat=kategorie , form = form, ...) 
  df <- df_des_cfa(designmatrix)
  designmatrix <- NA # added 22-06-2015
  usedform <- "margins" # added 22-06-2015
}

if(bintest==TRUE){
  erg <- data.frame(pat.=pattern, obs.=observed,exp.=expected,do.call(cbind,chi_local_test_cfa(observed,expected)),ex.bin.test=binomial_test_cfa(observed,expected), z_tests_cfa(observed,expected,ccor=ccor),p.stir=stirling_cfa(observed=observed, expected=expected,cum=TRUE,verb=FALSE),density.stir=stirling_cfa(observed=observed, expected=expected,cum=FALSE,verb=FALSE))
}

if(bintest==FALSE){
  erg <- data.frame(pat.=pattern, obs.=observed,exp.=expected,do.call(cbind,chi_local_test_cfa(observed,expected)),ex.bin.test=as.numeric(rep(x = NA,length.out=length(observed))), z_tests_cfa(observed,expected,ccor=ccor),p.stir=stirling_cfa(observed=observed, expected=expected,cum=TRUE,verb=FALSE),density.stir=stirling_cfa(observed=observed, expected=expected,cum=FALSE,verb=FALSE))
}


chi.square <- sum(erg$Chi)
  
if(method=="log"){df <- df_des_cfa(designmatrix)} ## added 22-06-2015
  
chi.square.p <- (1-pchisq(chi.square,df))
  
lr.chi <- lr(observed,expected) ## added 20. October 2014 JHH
  
lr.p <- (1-pchisq(lr.chi,df)) ## added 20. October 2014 JHH
  
bonferroni <- alpha/length(expected)
  
result <- list( local.test = erg, bonferroni.alpha=bonferroni, global.test = list(pearson = list(Chi=chi.square,df=df,pChi=chi.square.p,alpha=alpha), likelihood.ratio = list(Chi=lr.chi,df=df,pChi=lr.p,alpha=alpha), infocrit=list(loglik=loglik, AIC=aic, BIC=bic)),designmatrix=designmatrix, variables=kategorie, used.formula=usedform, functional=blank, inputdata=patternfreq, alpha=alpha) 

if(!is.null(WGT)){
  comment(result) <- "using weighted data!"
}
class(result)<-c("CFA","list")
 
return(result)
}
# End of function ------ 