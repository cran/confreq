#' @exportS3Method print Pfreq
#' @keywords methods
#' @method print Pfreq
#' @title S3 print for Pfreq
#' @description S3 print method for object of class\code{"Pfreq"}
#' @param x object of class \code{"Pfreq"}
#' @param ... further arguments passed to or from other methods.
#' @param digits minimal number of significant digits, see \code{\link{print.default}}.
#' @param quote logical, indicating whether or not strings should be printed with surrounding quotes.
#' @param right logical, indicating whether or not strings should be right aligned.
#' @param row.names logical (or character vector), indicating whether (or what) row names should be printed.
#' @param max numeric or \code{NULL}, specifying the maximal number of entries to be printed. By default, when \code{NULL}, \code{\link{getOption}("max.print")} used.
#' @return output printed to the console


########################### hier die print method #class Pfreq #################
print.Pfreq <- function (x, ..., digits = NULL, quote = FALSE, right = TRUE, row.names = TRUE, max = NULL) 
{
  x$wgt <- attributes(x)$WGT
  x$wgtFreq <- attributes(x)$wgtfreq
  
  ### rest wie data.frame method
  n <- length(row.names(x))
  if (length(x) == 0L) {
    cat(sprintf(ngettext(n, "data frame with 0 columns and %d row", 
                         "data frame with 0 columns and %d rows"), n), "\n", 
        sep = "")
  }
  else if (n == 0L) {
    print.default(names(x), quote = FALSE)
    cat(gettext("<0 rows> (or 0-length row.names)\n"))
  }
  else {
    if (is.null(max)) 
      max <- getOption("max.print", 99999L)
    if (!is.finite(max)) 
      stop("invalid 'max' / getOption(\"max.print\"): ", 
           max)
    omit <- (n0 <- max%/%length(x)) < n
    m <- as.matrix(format.data.frame(if (omit) 
      x[seq_len(n0), , drop = FALSE]
      else x, digits = digits, na.encode = FALSE))
    if (!isTRUE(row.names)) 
      dimnames(m)[[1L]] <- if (isFALSE(row.names)) 
        rep.int("", if (omit) 
          n0
          else n)
    else row.names
    print(m, ..., quote = quote, right = right, max = max)
    
    if (omit) 
      cat(" [ reached 'max' / getOption(\"max.print\") -- omitted", 
          n - n0, "rows ]\n")
  }
  invisible(x)
}