#' A GWR model-fitting function adjusting for measurement error
#'
#' This function lets you fit a geographically weighted regression model incorporating measurement error on independent or dependent variables that may vary spatially
#' @param X design matrix of independent variable values (variables contained in columns)
#' @param Y vector of dependent variable values
#' @param u design matrix of measurement error standard deviations on X (first column, corresponding to the intercept term, is usually a vector of 0s). This should have the same dimensions as X with corresponding elements
#' @param q vector of measurement error standard deviations on Y
#' @param DM n x n distance matrix (matrix of distances between each n locations)
#' @param p.ci percent confidence intervals
#' @param error.type classical or non-classical measurement error model
#' @keywords measurement error, sampling error, geographically weighted regression
#' @export
#' @examples
#' spatialME()

spatialME=function(X,Y,u,q,DM,bw,p.ci=0.95,error.type="classical") {
  n=nrow(DM)
  est_mat=se_mat=matrix(NA,nrow=n,ncol=ncol(X))
  for(i in 1:n) {
    dist_weights=gaussian.kernel(DM[i,], bw)
    me_fit=fitME(X=X,Y=Y,u=u,q=q,weights=dist_weights,error.type=error.type)
    est_mat[i,]=me_fit$Est
    se_mat[i,]=me_fit$SE
  }
  colnames(est_mat)=colnames(se_mat)=colnames(X)
  Lo=est_mat-qnorm(p.ci+(1-p.ci)/2)*se_mat; Up=est_mat+qnorm(p.ci+(1-p.ci)/2)*se_mat
  return(list("Est"=est_mat, "SE"=se_mat, "CI_Lo"=Lo, "CI_Up"=Up))
}
