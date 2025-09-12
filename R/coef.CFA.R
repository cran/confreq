#' @exportS3Method coef CFA
#' @keywords methods
#' @aliases coefficients 
#' @method coef CFA
#' @title S3 coefficients for CFA
#' @description S3 coefficients method for object of class\code{"CFA"}.
#' @param object object of class \code{"CFA"}.
#' @param ... other parameters passed trough.
#' @return Coefficients extracted from the model object of class\code{"CFA"}.

########################### hier die plot method #class CFA #######################
coef.CFA<-function(object, ...){
  coef(object$glmfitres)
}