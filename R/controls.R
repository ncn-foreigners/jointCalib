#' @title control parameters
#' @author Maciej Beręsewicz
#'
#' @description
#' \code{control_calib} is function that contains control parameters for \code{joint_calib_create_matrix}
#'
#' @param interpolation type of interpolation: \code{logit} or \code{linear}
#' @param logit_const constant for \code{logit} interpolation
#'
#' @return a list with parameters
#'
#' @export
control_calib <- function(interpolation=c("logit", "linear"),
                          logit_const=-1000) {

  if (missing(interpolation)) interpolation <- "logit"

  return(list(interpolation=interpolation,
              logit_const=logit_const))
}
