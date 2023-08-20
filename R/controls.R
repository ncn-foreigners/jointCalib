#' @title control parameters
#' @author Maciej Beręsewicz
#'
#' @description
#' \code{control_calib} is function that contains control parameters for \code{joint_calib_create_matrix}
#'
#' @param interpolation type of interpolation: \code{logit} or \code{linear},
#' @param logit_const constant for \code{logit} interpolation,
#' @param survey_sparse whether to use sparse matrices via \code{Matrix} package in [survey::grake()] (currently not supported),
#' @param ebal_constraint_tolerance This is the tolerance level used by ebalance to decide if the moments in the reweighted data are equal to the target moments (see [ebal::ebalance()]),
#' @param ebal_print_level Controls the level of printing: 0 (normal printing), 2 (detailed), and 3 (very detailed) (see [ebal::ebalance()]).
#'
#' @return a list with parameters
#'
#' @export
control_calib <- function(interpolation = c("logit", "linear"),
                          logit_const = -1000,
                          survey_sparse = FALSE,
                          ebal_constraint_tolerance=1,
                          ebal_print_level=0) {

  if (missing(interpolation)) interpolation <- "logit"

  return(list(interpolation = interpolation,
              logit_const = logit_const,
              survey_sparse = survey_sparse,
              ebal_constraint_tolerance = ebal_constraint_tolerance,
              ebal_print_level = ebal_print_level))
}
