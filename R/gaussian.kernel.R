#' A function to calculate Gaussian kernel-based weights for geographically weighted regression
#'
#' This function lets you calculate distance-based observation weghts using a Gaussian kernel
#' @param distance matrix of distances between locations in a defined spatial area
#' @param bandwidth bandwidth value used in calculating Gaussian weights; bandwidth selected using standard AIC or CV error minimization
#' @keywords geographically weighted regression
#' @export
#' @examples
#' gaussian.kernel()

gaussian.kernel = function(distance, bandwidth) exp(-0.5*(distance/bandwidth)^2)
