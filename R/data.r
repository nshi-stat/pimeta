#' Set-shifting data
#'
#' Higgins et al. (2009) re-analyzed data (Roberts et al., 2007)
#' that included 14 studies evaluating the set-shifting ability
#' in people with eating disorders by using a prediction interval.
#' Standardized mean differences in the time taken to complete
#' Trail Making Test between subjects with eating disorders and
#' healthy controls were collected.
#' Positive estimates indicate impairment in set shifting ability
#' in people with eating disorders.
#'
#' \itemize{
#'   \item \code{y}: Standardized mean difference
#'   \item \code{sigmak}: Standard error
#' }
#'
#' @docType data
#' @keywords datasets
#' @name setshift
#' @usage data(setshift)
#' @format A data frame with 14 rows and 2 variables
#' @references
#' Roberts, M. E., Tchanturia, K., Stahl, D.,
#' Southgate, L., and Treasure, J. (2007).
#' A systematic review and meta-analysis of set-shifting
#' ability in eating disorders.
#' \emph{Psychol Med.}
#' \strong{37}(8): 1075-1084.
#' \url{https://doi.org/10.1017/S0033291707009877}
#' 
#' Higgins, J. P. T, Thompson, S. G., Spiegelhalter, D. J. (2009).
#' A re-evaluation of random-effects meta-analysis.
#' \emph{J R Stat Soc Ser A Stat Soc.}
#' \strong{172}(1): 137-159.
#' \url{https://doi.org/10.1111/j.1467-985X.2008.00552.x}
NULL


#' Pain data
#'
#' The pain data (Riley et al., 2011; Hauser et al., 2009) included
#' 22 studies comparing the treatment effect of antidepressants on
#' reducing pain in patients with fibromyalgia syndrome.
#' The treatment effects were summarized using standardized mean
#' differences on a visual analog scale for pain between the
#' antidepressant group and control group.
#' Negative estimates indicate the reduction of pain in the
#' antidepressant group.
#'
#' \itemize{
#'   \item \code{y}: Standardized mean difference
#'   \item \code{sigmak}: Standard error
#' }
#'
#' @docType data
#' @keywords datasets
#' @name pain
#' @usage data(pain)
#' @format A data frame with 22 rows and 2 variables
#' @references
#' Hauser, W., Bernardy, K, Uceyler, N., and Sommer, C. (2009).
#' Treatment of fibromyalgia syndrome with antidepressants: a meta-analysis.
#' \emph{JAMA.}
#' \strong{301}(2): 198-209.
#' \url{https://doi.org/10.1001/jama.2008.944}
#' 
#' Riley, R. D., Higgins, J. P. T, and Deeks, J. J. (2011).
#' Interpretation of random effects meta-analyses.
#' \emph{BMJ.}
#' \strong{342}: d549.
#' \url{https://doi.org/10.1136/bmj.d549}
NULL


#' Systolic blood pressure (SBP) data
#'
#' Riley et al. (2011) analyzed a hypothetical meta-analysis.
#' They generated a data set of 10 studies examining the same
#' antihypertensive drug.
#' Negative estimates suggested reduced blood pressure in the
#' treatment group.
#'
#' \itemize{
#'   \item \code{y}: Standardized mean difference
#'   \item \code{sigmak}: Standard error
#' }
#'
#' @docType data
#' @keywords datasets
#' @name sbp
#' @usage data(sbp)
#' @format A data frame with 10 rows and 2 variables
#' @references
#' Riley, R. D., Higgins, J. P. T, and Deeks, J. J. (2011).
#' Interpretation of random effects meta-analyses.
#' \emph{BMJ.}
#' \strong{342}: d549.
#' \url{https://doi.org/10.1136/bmj.d549}
NULL


#' Hypertension data
#'
#' The hypertension data (Wang et al., 2005) included
#' 7 studies comparing the treatment effect of
#' anti-hypertensive treatment versus control on reducing
#' diastolic blood pressure (DBP) in patients with hypertension.
#' Negative estimates indicate the reduction of DBP in the
#' anti-hypertensive treatment group.
#' 
#' \itemize{
#'   \item \code{y}: Standardized mean difference
#'   \item \code{se}: Standard error
#'   \item \code{label}: Labels for each study
#' }
#'
#' @docType data
#' @keywords datasets
#' @name hyp
#' @usage data(hyp)
#' @format A data frame with 10 rows and 2 variables
#' @references
#' Wang, J. G., Staessen, J. A., Franklin, S. S., Fagard, R.,
#' and Gueyffier, F. (2005).
#' Systolic and diastolic blood pressure lowering as determinants
#' of cardiovascular outcome.
#' \emph{Hypertension.}
#' \strong{45}(5): 907-913.
#' \url{https://doi.org/10.1161/01.HYP.0000165020.14745.79}
NULL
