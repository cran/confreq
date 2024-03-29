% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/df_des_cfa.R
\name{df_des_cfa}
\alias{df_des_cfa}
\title{Degrees of freedom}
\usage{
df_des_cfa(des)
}
\arguments{
\item{des}{a designmatrix (object of class "matrix") as returned by function \code{design_cfg_cfa}.}
}
\value{
An object of class "integer" giving the degrees of freedom for the designmatrix defined in argument \code{des}.
}
\description{
Calculates the degrees of freedom based on an designmatrix for a (log liniear) CFA model.
}
\details{
No details
}
\examples{
#######################################
# degrees of freedom for designmatrix with three main effects.
# three variables with two categories each.
df_des_cfa(design_cfg_cfa(kat=c(2,2,2)))
# two variables with two categories each and one variable
# with 7 categories (e.g. The Krauth & Lienert suicide Data).
df_des_cfa(design_cfg_cfa(kat=c(2,2,7)))
###########
# degrees of freedom for designmatrix with three main effects
# and three 'two by two'interactions.
# and tripple interaction --> saturated model --> df=0
# three variables with two categories each.
df_des_cfa(design_cfg_cfa(kat=c(2,2,2),form="~ V1 + V2 + V3 + V1:V2 + V1:V3 + V2:V3 + V1:V2:V3"))
####################################### 
}
\references{
No references in the moment
}
\keyword{misc}
