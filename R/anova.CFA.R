#' @exportS3Method anova CFA
#' @keywords methods
#' @aliases anova 
#' @method anova CFA
#' @title S3  analysis of variance (or deviance) tables for CFA objects
#' @description S3 anova method for object of class\code{"CFA"}.
#' @param object object of class \code{"CFA"}.
#' @param ... additional objects of the same type.
#' @param test Character string specifying the statistical test
#'   used for model comparison. See \code{\link[stats]{anova.glm}}
#'   for available options. Defaults to \code{"Chisq"}.
#' @return This (generic) function returns an object of class anova.
#' @details Model descriptions use the CFA specification stored in
#'   \code{used.formula}, including functional patterns and covariates.
#'   The model fits and comparison statistics are not changed. With a single
#'   object, the complete design matrix is treated as one GLM term; the table
#'   does not decompose the individual CFA main effects or interactions.

########################### hier die anova method #class CFA ###################
anova.CFA<-function(object, ..., test = "Chisq"){
  objects <- list(object, ...)
  if (!all(vapply(objects, inherits, logical(1), "CFA"))) {
    stop("All objects must inherit from class 'CFA'.", call. = FALSE)
  }
  # GLM-Ergebnisse aus allen CFA-Objekten extrahieren
  models <- lapply(objects, `[[`, "glmfitres")
  
  # Pruefen, ob gueltige glm-Objekte vorliegen
  valid <- vapply(models, inherits, logical(1), "glm")
  
  if (!all(valid)) {
    stop(
      "Model comparison using anova() is currently not possible.\n",
      "At least one CFA object does not contain a valid glm result.\n",
      "Expected frequencies may have been calculated directly\n",
      "from marginal frequencies.",
      call. = FALSE
    )
  }
  
  # ANOVA auf die glm-Objekte anwenden
  result <- do.call(stats::anova, c(unname(models), list(test = test)))
  
  # Nur die Darstellung anpassen, nicht die Formeln der GLM-Fits.
  descriptions <- vapply(objects, function(x) {
    description <- x$used.formula
    if (is.null(description) || length(description) == 0L) {
      return("CFA model (specification unavailable)")
    }
    if (inherits(description, "formula")) {
      description <- paste(deparse(description), collapse = " ")
    } else {
      description <- paste(as.character(description), collapse = " ")
    }
    if (identical(trimws(description), "null")) {
      description <- "~ 1 (order-zero CFA model)"
    }
    gsub("[\r\n]+", " ", description)
  }, character(1))
  
  if (length(models) > 1L) {
    attr(result, "heading") <- c(
      "Analysis of Deviance Table\n",
      paste0("Model ", seq_along(models), ": ", descriptions,
             collapse = "\n")
    )
  } else {
    heading <- attr(result, "heading")
    attr(result, "heading") <- c(
      heading,
      paste0("CFA specification: ", descriptions),
      "The complete CFA design is tested as one GLM term.\n"
    )
    rownames(result)[rownames(result) == "designmatrix"] <- "CFA design"
  }
  
  result
  
}