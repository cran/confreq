\name{design_cfg_cfa}
\alias{design_cfg_cfa}
\title{Designmatrix for log linear CFA models}
\usage{
  design_cfg_cfa(kat,
    form = paste("~", paste(paste("V", 1:length(kat), sep = ""), collapse = " + ")),
    ...)
}
\arguments{
  \item{kat}{a numerical vector containing kardinal
  numbers, giving the number of categories for each
  variable of a dataset (in the respective order of the
  variables in such a dataset) which corresponds to the
  requested designmatrix. So the length of this numerical
  vector represents the number of variables.}

  \item{form}{a character string which can be coerced into
  a model formulae with the function \code{as.formula} in
  the package \code{stats}. If this argument is left empty
  the function \code{design_cfg_cfa()} will return a
  designmatrix coding only main effects and no interactions
  -- for a designmatrix refering to three variables for
  example, leaving the argument \code{form} empty will be
  equivalent to assigning the character \code{"~ V1 + V2 +
  V3"} to the argument (\code{form="~ V1 + V2 + V3"}).

  A special Case is to define a null-model or rather a cfa
  model of order zero. In such a model no (main) effects
  are considered. This can be achieved bei passing the
  character expression \code{"null"} to the argument
  \code{form} -- so: \code{form = "null"}}

  \item{...}{additional parameters passed through to
  function \code{model.matrix} in package \code{stats}.}
}
\value{
  A designmatrix (an object of class "matrix") for the
  formula therm given in argument\code{form}.
}
\description{
  Calculates the designmatrix corresponding to a dataset
  with \code{length(kat)} columns (variables).
}
\details{
  This function internaly calls the function
  \code{pos_cfg_cfa}.

  For further information on designmatrices see decription
  on function \code{model.matrix} in the package
  \code{stats}.
}
\examples{
#######################################
# designmatrix with three main effects.
# three variables with two categories each.
design_cfg_cfa(kat=c(2,2,2))
# two variables with two categories each and one variable
# with 7 categories (Linert LSD example).
design_cfg_cfa(kat=c(2,2,7))
###########
# designmatrix with three main effects an three interactions.
# three variables with two categories each.
design_cfg_cfa(kat=c(2,2,2),form="~ V1 + V2 + V3 + V1:V2 + V1:V3 + V2:V3")
# two variables with two categories each and one variable
# with 7 categories (Linert LSD example).
design_cfg_cfa(kat=c(2,2,7),form="~ V1 + V2 + V3 + V1:V2 + V1:V3 + V2:V3")
#######################################
}
\references{
  No references in the moment
}

