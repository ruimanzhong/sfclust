
#' Calculate Log Marginal Likelihood for All Clusters
#'
#' This function computes the log marginal likelihood for each cluster specified in the membership vector.
#' It applies a specified model to each cluster to calculate the likelihood.
#'
#' @param stars_obj A stars object containing spatial data, including covariates and response variables.
#' @param membership Integer vector indicating the cluster membership for each observation.
#' @param formula An object of class \code{\link[stats]{formula}} specifying the model to be used in INLA.
#' @param family Character string specifying the family of distributions to use for the model. Defaults to "normal".
#' @param correction Logical indicating whether a correction for dispersion or other factors should be applied. Defaults to FALSE.
#' @param time_var The name of the time variable in the stars object.
#' @param N_var The name of the variable used as exposure or trials (for Poisson, binomial, etc.).
#' @param ... Additional arguments passed to the underlying `log_mlik_each` function.
#'
#' @return A numeric vector containing the log marginal likelihood for each cluster.
#'
#' @examples
#' \dontrun{
#' log_mlikelihoods <- log_mlik_all(stars_obj, membership, family = "poisson", formula = Y ~ temperature + precipitation, time_var = "time", N_var = "expected_cases")
#' }
#' @export
log_mlik_all <- function(stars_obj, membership, formula, family = "normal", correction = FALSE, time_var, N_var, ...) {
  k <- max(membership)
  sapply(1:k, log_mlik_each, stars_obj, membership, formula, family, correction, time_var, N_var, ...)
}
