#' An ME model-fitting function
#'
#' This function lets you fit a regression model incorporating measurement error on independent or dependent variables
#' @param X design matrix of independent variable values (variables contained in columns)
#' @param Y vector of dependent variable values
#' @param u design matrix of measurement error standard deviations on X (first column, corresponding to the intercept term, is usually a vector of 0s). This should have the same dimensions as X with corresponding elements
#' @param q vector of measurement error standard deviations on Y
#' @param weights optional weights for each observation in (X,Y)
#' @param error.type classical or non-classical measurement error model
#' @keywords measurement error, sampling error, regression
#' @export
#' @examples
#' fitME()

fitME=function(X, Y, u, q, weights=rep(1,nrow(X)),error.type='clasical') {
  k=ncol(X)
  n=nrow(X)
  ### jackknife standard errors, so looping
  jack_B=matrix(NA,nrow(X),k)
  for(i in 1:n) {
    # weighted covariance matrix of u (kxk dimensions); have to do this because var(ui)=mean(ui) when
    # each ui is a standard error value
    .X=X[-i,]; .Y=Y[-i]; .u=u[-i,]; .q=q[-i]; .weights=weights[-i]
    wcor_u=cov.wt(.u,wt=.weights,method="unbiased",cor=T)$cor
    wcor_u[wcor_u=="NaN"]=0; diag(wcor_u)=1
    wcmeans_u=apply(.u,2,function(x) weighted.mean(x^2,w=.weights)) # this is the variance vector
    wcov_u=wcor_u*(sqrt(wcmeans_u) %*% t(sqrt(wcmeans_u))) # to make Corrs Covs
    # covariance matrix of u and q (kx1 dimensions)
    wcor_uq=t(apply(.u,2,function(x) cov.wt(cbind(x,.q),wt=.weights,method="unbiased",cor=T)$cor[2,1]))
    wcor_uq[wcor_uq=="NaN"]=0
    wmean_q=weighted.mean(.q^2,w=.weights)
    wcov_uq=wcor_uq*(wcmeans_u*wmean_q)

    # if X and u or y and q are correlated
    if(error.type=="non-classical") {
      # covariance matrix of X and u
      wcors_xu=mapply(function(x, y){cor(x, y)}, as.data.frame(.X), as.data.frame(.u))
      wcors_xu[is.na(wcors_xu)]=0
      wcovs_xu=wcors_xu*sqrt(diag(var(.X))*wcmeans_u)
      betas_xu=wcovs_xu/diag(var(.X)); betas_xu[is.na(betas_xu)]=0
      # covariance of Y and q
      cov_yq=cor(.Y,.q)*sqrt(var(.Y)*mean(.q)); cov_yq[is.na(cov_yq)]=0
      beta_yq=cov_yq/var(.Y)
      # betas from regressing each X on its corresponding error (b1s calculating by covariance and variance)
      .A=.X %*% solve(diag(ncol(.X))+diag(betas_xu)) # new X
      .Q=.Y %*% solve(1+beta_yq)
      # .e=.A-.X # new measurement error on X ... this is a model helpful for reference when doing things below
      # .v=.Q-.Y # new measurement error on Y ... this is a model helpful for reference when doing things below

      # need to calculate the variance-covariance matrix of the new errors (epsilon)
      wcors_au=mapply(function(x, y){cor(x, y)}, as.data.frame(.A), as.data.frame(.u))
      wcors_au[is.na(wcors_au)]=0
      wcovs_au=wcors_au*sqrt(diag(var(.A))*wcmeans_u)
      wcov_e=var(.A)+var(.X)+wcov_u-2*cov(.A,.X)-2*diag(wcovs_xu)-2*diag(wcovs_au)

      # also need the vector of covariances between each u vector and q
      wcov_ev=wcov_uq # derive Cov(v,e) to see Cov(v,e)=Cov(u,q)

      ### make error covariance matrix estimation
      mod.wcov_u=wcov_e
      mod.wcov_uq=wcov_ev
    } else {
      .A=.X; .Q=.Y
      mod.wcov_u=wcov_u
      mod.wcov_uq=wcov_uq
    }
    # for beta calculations, dont need Corr(x,u) or Corr(Y,q)
    Mxx=1/nrow(.A)*(t(.weights*.A) %*% (.A) - mod.wcov_u)
    Mxy=1/nrow(.A)*(t(.weights*.A) %*% .Q - t(mod.wcov_uq))
    jack_B[i,]=solve(Mxx) %*% Mxy
  }
  se=apply(jack_B,2,function(x) sqrt((n-1)*var(x)))
  return(list("Est"=colMeans(jack_B), "SE"=se, "jack_est"=jack_B))
}
