#' @title Configural Frequencies Analysis for two Samples.
#' @keywords mainfunction 
#' @export S2CFA
#' @description Calculates coefficients for the two-sample CFA. Instead of differentiating between 'Types' and 'Antitypes', two-sample CFA looks for discrimination types, that is configurations with significant differences in frequencies between two sub samples.
#' @details With survey weights, observed and expected frequencies use sums
#' of weights, without rounding. Inference compares the two weighted proportions
#' by an asymptotic Wald chi-squared test (one degree of freedom) using a
#' linearized variance for independent sampling units, with replacement.
#' No clustering, stratification, calibration or finite-population correction
#' can be incorporated from the weight vector alone.
#'
#' The influence for record i in group g is
#' \code{w_i * (I_i - p_g) / sum(w_g)}. The two group contributions are
#' disjoint. Their squared influences are summed and multiplied by
#' \code{n / (n - 1)}, where n counts all positive-weight records in both
#' groups. The variance is obtained from \code{wgtfreq2} and
#' \code{wgtpositive}, stored by the revised \code{dat2fre()}.
#' Zero-variance configurations are not tested and have an NA p-value.
#' For old weighted objects without these attributes, only weighted descriptive
#' results are returned with a warning; recreate them with \code{dat2fre()}.
#'
#' The weighted p-value is \code{p.wald}, not a Fisher exact test and not a
#' Pearson test on population expansion totals. Classical test columns remain
#' NA in weighted results. Multiplying all weights by a common positive
#' factor leaves the Wald statistic and p-value unchanged.
#' The unweighted statistical calculation retains its historical behavior.
#' Both groups are paired by categories, even after a row permutation.
#' Incomplete or duplicated configuration grids are rejected explicitly.
#' 
#' @param patternfreq an object of class \code{"Pfreq"}, which is data in pattern frequencies representation - see function \code{\link{dat2fre}}. The variable defining the two sub samples (a variable with max. two categories) must be located in the last but one column of the object of class \code{"Pfreq"} 
#' 
#' @param alpha a numeric giving the alpha level for testing (default set to \code{alpha=.05})
#' 
#' @param ccor For unweighted analysis only; TRUE is rejected for survey weights. A logical (TRUE / FALSE) determining whether to apply a continuity correction or not. When set to \code{ccor=TRUE} continuity correction is applied. For \code{ccor=FALSE} no continuity correction is applied.
#'   
#' @param ... Reserved for compatibility; not used by this implementation.
#' @param weighted Logical or NULL. By default, use survey weights if the
#'   \code{wgtfreq} attribute is present. TRUE explicitly requests weighted
#'   analysis; FALSE explicitly uses raw frequencies. This argument follows
#'   \code{...} and must be named in full.
#' @return An object of class \code{S2CFA}. Weighted results additionally
#'   contain \code{weighted = TRUE} and a \code{weighting} description.
#'   The first ten local result columns retain the historical layout; survey
#'   columns are \code{Wald}, \code{df.wald}, \code{p.wald}, \code{prop.diff}
#'   and \code{SE.diff}. Frequencies in weighted results are sums of weights.
#' @references Stemmler, M. (2020). \emph{Person-Centered Methods – Configural Frequency Analysis (CFA) and Other Methods for the Analysis of Contingency Tables.} Cham Heidelberg New York Dordrecht London: Springer.
#' @references Stemmler, M., & Hammond, S. (1997). Configural frequency analysis of dependent samples for intra-patient treatment comparisons. \emph{Studia Psychologica, 39}, 167–175.
#' 
#' @examples #######################################
#' ############### some examples #########
#' ######### example from Marks Textbook
#' data(Lienert1978)
#' res1 <- S2CFA(Lienert1978)
#' summary(res1)
#' res2 <- S2CFA(Lienert1978, ccor=TRUE) # with continuity correction
#' summary(res2)
#' ######### example with biger numbers
#' data(suicide)
#' ftab(suicide) # 'Epoche' may divide the sample into 2 subsamples  
#' suicide_2s <- suicide[, c(1,3,2) ] # reorder data that 'Epoche' is the last column
#' ftab(suicide_2s) # check reordering
#' suicide_2s_fre <- dat2fre(suicide_2s)
#' res3 <- S2CFA(suicide_2s_fre)
#' summary(res3)
#' res4 <- S2CFA(suicide_2s_fre, ccor=TRUE) # with continuity correction
#' summary(res4)
###############################################################
############### start of function definition ##################
S2CFA<-function(patternfreq, alpha=.05, ccor=FALSE, ..., weighted=NULL){
if(any(class(patternfreq)=="Pfreq") != TRUE){stop("patternfreq must be an object of class 'Pfreq'","\n","see func. dat2fre()", call. = TRUE) }

  # Eine Teilmenge bleibt eventuell eine Pfreq-Tabelle, aber die Analyse
  # braucht genau eine Zeile fuer jede Kombination der Kategorien.
  # Faktorlevels zaehlen auch dann mit, wenn ihre Nullzellen entfernt wurden.
  configurations <- as.data.frame(patternfreq)[seq_len(ncol(patternfreq) - 1L)]
  if (!ncol(configurations) || !nrow(configurations))
    stop("The analysis requires a nonempty configuration table.")
  categories <- lapply(configurations, function(z)
    if (is.factor(z)) levels(z) else unique(z))
  category_sizes <- lengths(categories) + as.integer(
    vapply(configurations, anyNA, logical(1)) & !vapply(categories, anyNA, logical(1)))
  if (anyDuplicated(configurations) ||
      prod(category_sizes) != nrow(configurations))
    stop("The analysis requires one row per complete configuration; recreate missing zero cells with dat2fre().")

  # Survey-Gewichte sind keine Replikationshaeufigkeiten. Der folgende Zweig
  # verwendet gewichtete Anteile und eine linearisierte Varianz unter der
  # Annahme unabhaengiger Beobachtungen (ohne Cluster/Schichten/FPC).
  # Der ungewichtete Originalzweig darunter bleibt unveraendert.
  if (!is.null(weighted) && (!is.logical(weighted) || length(weighted) != 1L || is.na(weighted)))
    stop("weighted must be NULL, TRUE or FALSE.")
  has_weights <- !is.null(attr(patternfreq, "wgtfreq", exact = TRUE))
  if (is.null(weighted)) weighted <- has_weights
  if (weighted) {
    z <- .cf_frequency(patternfreq, weighted = TRUE)
    configuration <- z$data
    if (ncol(configuration) < 2L)
      stop("Need at least one configuration variable followed by the grouping variable.")
    if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1)
      stop("alpha must be between zero and one.")
    if (!is.logical(ccor) || length(ccor) != 1L || is.na(ccor))
      stop("ccor must be TRUE or FALSE.")
    if (ccor) stop("ccor is a count-table correction, not used for the survey Wald test.")
    group_name <- names(configuration)[ncol(configuration)]
    predictors <- names(configuration)[-ncol(configuration)]
    group_values <- configuration[[group_name]]
    if (anyNA(group_values)) stop("Missing group membership is not supported.")
    groups <- names(table(group_values))
    if (length(groups) != 2L) stop("The grouping variable must have two categories.")
    i1 <- which(as.character(group_values) == groups[1L])
    i2 <- which(as.character(group_values) == groups[2L])
    if (!length(i1) || length(i1) != length(i2))
      stop("Both groups need the same complete configuration table.")

    # Paarung ueber die Kategorien, nicht ueber die zufaellige Zeilenposition.
    # Mehrstellige Codes oder Labels mit Leerzeichen bleiben eindeutig.
    paired <- integer(length(i1))
    for (j in seq_along(i1)) {
      hit <- rep(TRUE, length(i2))
      for (nm in predictors) {
        a <- as.character(configuration[[nm]][i1[j]])
        b <- as.character(configuration[[nm]][i2])
        same <- (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
        hit <- hit & same
      }
      if (sum(hit) != 1L) stop("Configurations must occur once in each group.")
      paired[j] <- i2[which(hit)]
    }
    if (anyDuplicated(paired)) stop("Duplicate configurations in the first group.")
    i2 <- paired
    patterns <- function(rows) do.call(paste, c(lapply(configuration[rows, predictors, drop = FALSE], as.character), sep = " "))
    observed_1 <- z$freq[i1]
    observed_2 <- z$freq[i2]
    total_1 <- sum(observed_1)
    total_2 <- sum(observed_2)
    if (!is.finite(total_1 + total_2) || total_1 <= 0 || total_2 <= 0)
      stop("Both groups need a positive finite total weight.")
    p1 <- observed_1 / total_1
    p2 <- observed_2 / total_2
    pooled <- (observed_1 + observed_2) / (total_1 + total_2)
    expected_1 <- pooled * total_1
    expected_2 <- pooled * total_2

    # Fuer die Varianz werden die zweiten Gewichtsmomente benoetigt.
    # Sie lassen sich NICHT aus Zellgewichtssumme und Fallzahl rekonstruieren.
    square_sums <- attr(patternfreq, "wgtfreq2", exact = TRUE)
    positive_n <- attr(patternfreq, "wgtpositive", exact = TRUE)
    moments_available <- !is.null(square_sums) && !is.null(positive_n)
    k <- length(p1)
    variance <- rep(NA_real_, k)
    statistic <- rep(NA_real_, k)
    p_wald <- rep(NA_real_, k)
    if (moments_available) {
      if (!is.numeric(square_sums) || length(square_sums) != nrow(configuration) ||
          anyNA(square_sums) || any(!is.finite(square_sums)) || any(square_sums < 0) ||
          !is.numeric(positive_n) || length(positive_n) != nrow(configuration) ||
          anyNA(positive_n) || any(!is.finite(positive_n)) || any(positive_n < 0 | positive_n != floor(positive_n)))
        stop("Invalid survey weight moments; recreate the table with dat2fre().")
      n_units <- sum(positive_n[c(i1, i2)])
      if (n_units > 1L) {
        # Einflussfunktion des Gruppenanteils: w_i * (I_i - p_g) / W_g.
        # Die Gruppen sind disjunkt, daher summieren sich ihre Varianzbeitraege.
        # n/(n-1): Stichprobenkorrektur fuer unabhaengige positive Gewichtstraeger.
        q1 <- square_sums[i1]
        q2 <- square_sums[i2]
        v1 <- (q1 * (1 - p1)^2 + pmax(0, sum(q1) - q1) * p1^2) / total_1^2
        v2 <- (q2 * (1 - p2)^2 + pmax(0, sum(q2) - q2) * p2^2) / total_2^2
        variance <- n_units / (n_units - 1) * (v1 + v2)
        usable <- is.finite(variance) & variance > 0
        statistic[usable] <- (p1[usable] - p2[usable])^2 / variance[usable]
        p_wald[usable] <- stats::pchisq(statistic[usable], df = 1, lower.tail = FALSE)
        # Eine Konfiguration ohne Varianz ist nicht testbar. Es wird weder
        # ein kuenstlicher p-Wert noch eine Signifikanz aus Nullvarianz erzeugt.
      }
    } else {
      warning("Survey moments are missing: weighted summaries only. Recreate with dat2fre(..., wgt=...).", call. = FALSE)
    }

    # Die ersten zehn Spalten behalten das historische Schema. Die klassischen
    # Tests bleiben bei Surveygewichten bewusst NA; kein Runden/Pseudo-Fisher.
    erg <- data.frame(pattern_1 = patterns(i1), expected_1, observed_1,
                      pattern_2 = patterns(i2), expected_2, observed_2,
                      ex.fisher.test = rep(NA_real_, k), Chi = rep(NA_real_, k),
                      df = rep(NA_real_, k), pChi = rep(NA_real_, k),
                      Wald = statistic, df.wald = ifelse(is.na(p_wald), NA_real_, 1),
                      p.wald = p_wald, prop.diff = p1 - p2, SE.diff = sqrt(variance))
    names(erg)[1:6] <- paste(rep(groups, each = 3L), rep(c("pat.", "exp.", "obs."), 2L), sep = ".")
    result <- list(local.test = erg, bonferroni.alpha = alpha / k,
                   global.test = NULL, inputdata = patternfreq, alpha = alpha,
                   ccor = ccor, grouping = list(var = group_name, cats = groups),
                   pat = predictors, weighted = TRUE,
                   weighting = list(type = "survey", method = "independent-unit linearized Wald",
                     assumptions = "Independent sampling units; no clustering, stratification, calibration or finite-population correction.",
                     totals = stats::setNames(c(total_1, total_2), groups),
                     variance.available = moments_available))
    class(result) <- c("S2CFA", "list")
    comment(result) <- "using survey weighted data; independent-unit linearized Wald test"
    return(result)
  }

#---------- Abfragen von steuerargumenten -----------------
# bisher keine
#---------- Aufbereiten von patternfreq -----------------
subsamples <- names(table(patternfreq[,(ncol(patternfreq)-1)]))
if(length(subsamples)!=2){stop("wrong definition of subsamples! - check argument 'patternfreq'")}
# Konfigurationen der beiden Gruppen unabhaengig von ihrer Zeilenfolge paaren.
if (anyNA(configurations)) stop("Missing configuration categories are not supported by the unweighted count-table tests.")
raw_i1 <- which(as.character(patternfreq[[ncol(patternfreq)-1L]]) == subsamples[1L])
raw_i2 <- which(as.character(patternfreq[[ncol(patternfreq)-1L]]) == subsamples[2L])
raw_predictors <- seq_len(ncol(patternfreq)-2L)
raw_match <- integer(length(raw_i1))
for (j in seq_along(raw_i1)) {
  hit <- rep(TRUE, length(raw_i2))
  for (v in raw_predictors)
    hit <- hit & as.character(patternfreq[[v]][raw_i2]) == as.character(patternfreq[[v]][raw_i1[j]])
  if (sum(hit) != 1L) stop("Configurations cannot be paired uniquely across groups.")
  raw_match[j] <- raw_i2[which(hit)]
}
raw_i2 <- raw_match
pattern_1 <- do.call(paste, c(lapply(as.data.frame(patternfreq)[raw_i1,raw_predictors,drop=FALSE],as.character),sep=" "))
pattern_2 <- do.call(paste, c(lapply(as.data.frame(patternfreq)[raw_i2,raw_predictors,drop=FALSE],as.character),sep=" "))
names(pattern_1) <- row.names(patternfreq)[raw_i1]
names(pattern_2) <- row.names(patternfreq)[raw_i2]
observed_1 <- patternfreq[[ncol(patternfreq)]][raw_i1]
observed_2 <- patternfreq[[ncol(patternfreq)]][raw_i2]
oplist <- mapply(FUN=function(o1,o2){a=o1;b=o2; c=sum(observed_1)-o1; d=sum(observed_2)-o2; A=a+b; B=c+d;C=a+c; D=b+d; n=A+B; data.frame(a,b,c,d,A,B,C,D,n)  },o1=observed_1 ,o2=observed_2, SIMPLIFY = FALSE)
names(oplist) <- pattern_1
# oplist
n <- sum(c(observed_1,observed_2))
with.out<-function(x, y) {x[!x %in% y]} # hilfsfunktion
#----------------- calculate expected counts ------------------
expected_1 <- unlist(lapply(oplist, function(x){(x["A"]*x["C"])/x["n"]}))
expected_2 <- unlist(lapply(oplist, function(x){(x["A"]*x["D"])/x["n"]}))
names(expected_1) <- pattern_1
names(expected_2) <- pattern_2
#----------------- calculate exact fisher Test ------------------
#### factorial mit kürzen und package{gmp}  weniger limitiert
oplist_1 <- lapply(oplist, function(x){ list(oben= tabulate(sort(with.out(c((1:x[1,"A"]),(1:x[1,"B"]),(1:x[1,"C"]),(1:x[1,"D"])),0)),n)      , unten= tabulate(sort(with.out(c((1:x[1,"n"]),(1:x[1,"a"]),(1:x[1,"b"]),(1:x[1,"c"]),(1:x[1,"d"])),0)),n) ) } )
# oplist_1
#### kürzen
oplist_2 <- lapply(oplist_1, function(x){x$oben - x$unten}) #  alle mit "-" stehen dann unten
# oplist_2
oplist_3_ind <- lapply(oplist_2, function(x){list(ind_oben = x>0, ind_unten = x<0)}) # index was bleibt oben und unten
# oplist_3_ind
oplist_3 <- lapply(oplist_2, function(x){xv <- abs(x); names(xv) <- 1:n; xv}) # absolut werte haufigkeiten der produkttherme oben und unten
# oplist_3
oplist_4 <- mapply(function(ind, wert){ list(oben=wert[ind$ind_oben],unten=wert[ind$ind_unten]) }, ind=oplist_3_ind, wert=oplist_3, SIMPLIFY = FALSE) # names sind die Produkterme; vector ist die anzahl des jew. Produkterms 
# oplist_4
###### berechnungen mit gekürzter Produkterm liste 'oplist_4'
# library(gmp) # für funkktion 'as.bigz'
oplist_5 <- lapply(oplist_4, function(x){ list(oben=as.bigz(rep(x=as.numeric(names(x$oben)), times=x$oben)), unten=as.bigz(rep(x=as.numeric(names(x$unten)), times=x$unten))      )} )
# oplist_5
oplist_6 <- lapply(oplist_5, function(x){ list(oben= prod(x$oben), unten= prod(x$unten)) } )
# oplist_6
fisher_res <- unlist(lapply(oplist_6, function(x){as.double(div.bigq (x$oben , x$unten)) }))
# fisher_res

#----------------- calculate traditional Chi-square Test ------------------
# oplist
if(ccor==FALSE){
  # x <- oplist[[5]]
oplist_Chi_1 <- lapply(oplist, function(x){ list(oben=(prod(as.bigz(c(x[1,"n"], ((x[1,"a"]*x[1,"d"]) - (x[1,"b"]*x[1,"c"])) , ((x[1,"a"]*x[1,"d"]) - (x[1,"b"]*x[1,"c"])))))), unten=(prod(as.bigz(c(x[1,"A"], x[1,"B"], x[1,"C"], x[1,"D"]))))  ) } )
# x <- oplist_Chi_1[[5]]
  Chi <- unlist(lapply(oplist_Chi_1, function(x){ if(as.double(x$unten)!=0){as.double(div.bigq( x$oben, x$unten))}else{NaN} })) 
}
if(ccor==TRUE){
  # x <- oplist[[2]]
  # 36* (abs(3*4-14*15)  -   36 * 0.5)    ^2  / (17*19*18*18)
  # Chi <- unlist(lapply(oplist, function(x){ ( x[1,"n"] * (abs((x[1,"a"]*x[1,"d"]) - (x[1,"b"]*x[1,"c"])) - (x[1,"n"]) * 0.5)^2 )  / (x[1,"A"]*x[1,"B"]*x[1,"C"]*x[1,"D"])   } )) #  as.bigz fehlt noch
  oplist_Chi_1_cont <- lapply(oplist, function(x){ list(oben=(prod(as.bigz(c(x[1,"n"],(abs((x[1,"a"]*x[1,"d"])-(x[1,"b"]*x[1,"c"]))-x[1,"n"]*0.5),(abs((x[1,"a"]*x[1,"d"])-(x[1,"b"]*x[1,"c"]))-x[1,"n"]*0.5))))) , unten=(prod(as.bigz(c(x[1,"A"], x[1,"B"], x[1,"C"], x[1,"D"]))))  ) } )
  Chi <- unlist(lapply(oplist_Chi_1_cont, function(x){ if(as.double(x$unten)!=0){as.double(div.bigq( x$oben, x$unten))}else{NaN} }))  
}

df <- rep(1, times=length(Chi))

pChi <- (1-pchisq(Chi,df))

bonferroni <- alpha/length(expected_1)

erg <- data.frame(pattern_1, expected_1, observed_1, pattern_2, expected_2, observed_2, ex.fisher.test=fisher_res, Chi=Chi, df=df, pChi=pChi)
names(erg)[1:6] <- paste(c( rep(subsamples[1],3),rep(subsamples[2],3) ), rep(c("pat.","exp.","obs."),2)    ,sep=".")

result <- list( local.test = erg, bonferroni.alpha=bonferroni, global.test = NULL, inputdata=patternfreq, alpha=alpha, ccor=ccor, grouping= list(var=names(patternfreq)[ncol(patternfreq)-1], cats=subsamples), pat=names(patternfreq)[-c(ncol(patternfreq), ncol(patternfreq)-1)] ) 

class(result)<-c("S2CFA","list")

return(result)
}
# End of function ------ 
