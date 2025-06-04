############################ Linear Regression ##################################

##' @param data_int: 3 columns, Y, X1 and X2 where X1 = T and X2 = NT
##' @param Z: a matrix with other covariates, NULL if no other covariates
##' @param rho_shared = n_shared/N
##' @param alpha_ext = summary statistic gwas
##' @param sd_alpha_ext = sd of alpha from gwas summary statistic
##' @export
calibrated_est_trio_1 <- function(data_int,Z = NULL,rho_shared,alpha_ext,sd_alpha_ext){
  # Creating the internal dataset
  Y_int <- data_int[,1]
  X_1_int <- data_int[,2]
  X_2_int <- data_int[,3]
  #internal model
  mod_int <- lm(Y_int ~ cbind(X_1_int,X_2_int,Z))
  #external model
  mod_unadjusted <- lm(Y_int ~ cbind(X_1_int,Z))
  n <- nrow(data_int)
  # creating the X matrix for internal model
  X <- cbind(rep(1,n),X_1_int,X_2_int,Z)
  # creating the X_tilde matrix for external data
  X_tilde <- cbind(rep(1,n),X_1_int,Z)
  p <- ncol(X)-3   ## no. of columns in Z
  sigma_hat_sq <- sum(mod_int$residuals^2)/(n-(p+3))
  sigma_star_hat_sq <- sum(mod_unadjusted$residuals^2)/(n-(p+2))
  # this returns (X_i*X_i^T) which is a (p+3)*(p+3) matrix as the i-th column of out.prod1
  out.prod1 <- apply(X,1,function(x){outer(x,x)})
  # want  (\tilde{X}_i*\tilde{X}_i^T) which is a 2*2 matrix
  # as the i-th column of out.prod2
  # like out.prod2 <- apply(X_tilde,1,function(x){outer(x,x)})
  ## but for faster computation, remove some indexes from out.prod1 to get out.prod2
  index_rem <- c(seq(3,(p+3)^2,(p+3)),seq(2*(p+3)+1,3*(p+3),1))
  out.prod2 <- out.prod1[-index_rem,]
  D_1_hat <- - (1/(n*sigma_hat_sq))*matrix(rowSums(out.prod1),(p+3),(p+3))
  D_2_hat <- - (1/(n*sigma_star_hat_sq))*matrix(rowSums(out.prod2),(p+2),(p+2))
  C_11_hat <- (1/n)*(1/sigma_hat_sq^2)*
    matrix(rowSums(out.prod1%*%diag((mod_int$residuals)^2)),(p+3),(p+3))
  C_22_hat <- (1/n)*(1/sigma_star_hat_sq^2)*
    matrix(rowSums(out.prod2%*%diag((mod_unadjusted$residuals)^2)),(p+2),(p+2))
  out.prod3 <- out.prod1[-c((2*(p+3)+1):(3*(p+3))),]
  C_12_hat <- (1/n)*(1/(sigma_hat_sq*sigma_star_hat_sq))*
    matrix(rowSums(out.prod3%*%diag((mod_unadjusted$residuals)*(mod_int$residuals))),(p+3),(p+2),byrow = "F")
  D_1_inv <- solve(D_1_hat)
  D_2_inv <- solve(D_2_hat)
  a <- D_1_inv %*% C_11_hat %*% D_1_inv
  b <- D_1_inv %*% C_12_hat %*% D_2_inv
  c <- D_2_inv %*% C_22_hat %*% D_2_inv
  V_1 <- cbind(a[1:3,1:3],b[1:3,1:2],rho_shared*b[1:3,2])
  V_2 <- cbind(t(b[1:3,1:2]),c[1:2,1:2],rho_shared*c[1:2,2])
  V_3 <- c(as.vector(rho_shared*b[1:3,2]),as.vector(rho_shared*c[1:2,2]),n*sd_alpha_ext^2)
  V <- rbind(V_1,V_2,V_3)
  A <- rbind(c(0,1,-1,0,0,0),c(0,0,0,0,1,-1))
  V_hat <- A%*%V%*%t(A)
  tau_raw <- as.vector(mod_int$coef)[2] - as.vector(mod_int$coef)[3]
  tau_cal <-  tau_raw -
    V_hat[1,2]/V_hat[2,2]*(as.vector(mod_unadjusted$coef)[2] - alpha_ext)
  var_tau_cal <- (1/n)*(V_hat[1,1] - V_hat[1,2]^2/V_hat[2,2])
  return (list("Raw Estimator" = tau_raw,"raw_variance" = V_hat[1,1]/n,"calibrated_est" = tau_cal,"variance" = var_tau_cal))
}


##' @param data_int: 3 columns, Y, X1 and X2 where X1 = G and X2 = G_PA
##" @param Z: a matrix with other covariates, NULL if no other covariates
##' @param data_ext: 2 columns, Y, X1
##' @param rho_shared = n_shared/N
##' @param alpha_ext = summary statistic gwas
##' @param sd_alpha_ext
##' @export
calibrated_est_trio_2 <- function(data_int,Z = NULL,rho_shared,alpha_ext,sd_alpha_ext){
  # Creating the internal dataset
  Y_int <- data_int[,1]
  X_1_int <- data_int[,2]
  X_2_int <- data_int[,3]
  #internal model
  mod_int <- lm(Y_int ~ cbind(X_1_int,X_2_int,Z))
  #external model
  mod_unadjusted <- lm(Y_int ~ cbind(X_1_int,Z))
  n <- nrow(data_int)
  # creating the X matrix for internal model
  X <- cbind(rep(1,n),X_1_int,X_2_int,Z)
  # creating the X_tilde matrix for external data
  X_tilde <- cbind(rep(1,n),X_1_int,Z)
  p <- ncol(X)-3
  sigma_hat_sq <- sum(mod_int$residuals^2)/(n-(p+3))
  sigma_star_hat_sq <- sum(mod_unadjusted$residuals^2)/(n-(p+2))
  # this returns (X_i*X_i^T) which is a (p+3)*(p+3) matrix as the i-th column of out.prod1
  out.prod1 <- apply(X,1,function(x){outer(x,x)})
  # want  (\tilde{X}_i*\tilde{X}_i^T) which is a 2*2 matrix
  # as the i-th column of out.prod2
  # like out.prod2 <- apply(X_tilde,1,function(x){outer(x,x)})
  ## but for faster computation, remove some indices from out.prod1 to get out.prod2
  index_rem <- c(seq(3,(p+3)^2,(p+3)),seq(2*(p+3)+1,3*(p+3),1))
  out.prod2 <- out.prod1[-index_rem,]
  D_1_hat <- - (1/(n*sigma_hat_sq))*matrix(rowSums(out.prod1),(p+3),(p+3))
  D_2_hat <- - (1/(n*sigma_star_hat_sq))*matrix(rowSums(out.prod2),(p+2),(p+2))
  C_11_hat <- (1/n)*(1/sigma_hat_sq^2)*
    matrix(rowSums(out.prod1%*%diag((mod_int$residuals)^2)),(p+3),(p+3))
  C_22_hat <- (1/n)*(1/sigma_star_hat_sq^2)*
    matrix(rowSums(out.prod2%*%diag((mod_unadjusted$residuals)^2)),(p+2),(p+2))
  out.prod3 <- out.prod1[-c((2*(p+3)+1):(3*(p+3))),]
  C_12_hat <- (1/n)*(1/(sigma_hat_sq*sigma_star_hat_sq))*
    matrix(rowSums(out.prod3%*%diag((mod_unadjusted$residuals)*(mod_int$residuals))),(p+3),(p+2),byrow = "F")
  D_1_inv <- solve(D_1_hat)
  D_2_inv <- solve(D_2_hat)
  a <- D_1_inv %*% C_11_hat %*% D_1_inv
  b <- D_1_inv %*% C_12_hat %*% D_2_inv
  c <- D_2_inv %*% C_22_hat %*% D_2_inv
  V_1 <- cbind(a[1:3,1:3],b[1:3,1:2],rho_shared*b[1:3,2])
  V_2 <- cbind(t(b[1:3,1:2]),c[1:2,1:2],rho_shared*c[1:2,2])
  V_3 <- c(as.vector(rho_shared*b[1:3,2]),as.vector(rho_shared*c[1:2,2]),n*sd_alpha_ext^2)
  V <- rbind(V_1,V_2,V_3)
  A <- rbind(c(0,1,0,0,0,0),c(0,0,0,0,1,-1))
  V_hat <- A%*%V%*%t(A)
  tau_raw <- as.vector(mod_int$coef)[2]
  tau_cal <-  tau_raw -
    V_hat[1,2]/V_hat[2,2]*(as.vector(mod_unadjusted$coef)[2] - alpha_ext)
  var_tau_cal <- (1/n)*(V_hat[1,1] - V_hat[1,2]^2/V_hat[2,2])
  return (list("Raw Estimator" = tau_raw, "raw_variance" = V_hat[1,1]/n, "calibrated_est" = tau_cal,"variance" = var_tau_cal))
}


##' @param data_int: n by 3 matrix of internal full data, (Y1, X1, X2) with X1 = T, X2 = NT
##' @param Z: a matrix with other covariates, NULL if no other covariates
##' @param alpha_ext: coefficient from unadjusted model on external data
##' @param sd_alpha_ext: sd of alpha_ext from GWAS summary statistic
##' @param rho_shared: proportion of shared samples between internal and external data
##' @returns calibrated estimator for logistic regression in trio data
##' @export


calibrated_est_trio_logistic_1 <- function(data_int, Z = NULL, rho_shared, alpha_ext, sd_alpha_ext){
  Y_int <- data_int[,1]
  X_1_int <- data_int[,2]
  X_2_int <- data_int[,3]
  n <- nrow(data_int)
  ##internal model
  mod_int <- glm(Y_int ~ cbind(X_1_int,X_2_int,Z), family = binomial)
  ##external model
  mod_unadjusted <- glm(Y_int ~ cbind(X_1_int,Z), family = binomial)
  #coef_adjusted <- c(mod_int$coef[2], mod_int$coef[3])
  #coef <- mod_ext$coef[2]

  X <- cbind(rep(1,n), X_1_int, X_2_int,Z)
  X_tilde <- cbind(rep(1,n), X_1_int,Z)
  p <- ncol(X)-3
  # this returns (X_i*X_i^T) which is a 3*3 matrix as the i-th column of out.prod1
  out.prod1 <- apply(X,1,function(x){outer(x,x)})
  # this returns (\tilde{X}_i*\tilde{X}_i^T) which is a 2*2 matrix as the i-th column of out.prod2
  index_rem <- c(seq(3,(p+3)^2,(p+3)),seq(2*(p+3)+1,3*(p+3),1))
  out.prod2 <- out.prod1[-index_rem,]
  D_1_hat <- - (1/n)*matrix(rowSums(out.prod1%*%diag(mod_int$fitted.values*(1 - mod_int$fitted.values))),(p+3),(p+3))
  D_2_hat <- - (1/n)*matrix(rowSums(out.prod2%*%diag(mod_unadjusted$fitted.values*(1 - mod_unadjusted$fitted.values))),(p+2),(p+2))
  C_11_hat <- (1/n)*matrix(rowSums(out.prod1%*%diag((data_int[,1] - mod_int$fitted.values)^2)),(p+3),(p+3))
  C_22_hat <- (1/n)*matrix(rowSums(out.prod2%*%diag((data_int[,1] - mod_unadjusted$fitted.values)^2)),(p+2),(p+2))
  out.prod3 <- out.prod1[-c((2*(p+3)+1):(3*(p+3))),]
  C_12_hat <- (1/n)*matrix(rowSums(out.prod3%*%diag((data_int[,1] - mod_int$fitted.values)*
                                                      (data_int[,1] - mod_unadjusted$fitted.values))),(p+3),(p+2),byrow = "F")
  D_1_inv <- solve(D_1_hat)
  D_2_inv <- solve(D_2_hat)
  a <-  D_1_inv %*% C_11_hat %*% D_1_inv
  b <- D_1_inv %*% C_12_hat %*% D_2_inv
  c <- D_2_inv %*% C_22_hat %*% D_2_inv

  V_1 <- cbind(a[1:3,1:3],b[1:3,1:2],rho_shared*b[1:3,2])
  V_2 <- cbind(t(b[1:3,1:2]),c[1:2,1:2],rho_shared*c[1:2,2])
  V_3 <- c(as.vector(rho_shared*b[1:3,2]),as.vector(rho_shared*c[1:2,2]),n*sd_alpha_ext^2)
  V <- rbind(V_1,V_2,V_3)
  A <- rbind(c(0,1,-1,0,0,0),c(0,0,0,0,1,-1))
  V_hat <- A%*%V%*%t(A)
  # calculate calibrated estimator
  tau_raw <- as.vector(mod_int$coef)[2] - as.vector(mod_int$coef)[3]
  tau_cal <-  tau_raw -
    V_hat[1,2]/V_hat[2,2]*(as.vector(mod_unadjusted$coef)[2] - alpha_ext)
  var_tau_cal <- (1/n)*(V_hat[1,1] - V_hat[1,2]^2/V_hat[2,2])
  return (list("Raw Estimator" = tau_raw,"raw_variance" = V_hat[1,1]/n,"calibrated_est" = tau_cal,"variance" = var_tau_cal))
}


##' @param data_int: n by 3 matrix of internal full data, (Y1, X1, X2) with X1 = G, X2 = G_pa
##' @param Z: a matrix with other covariates, NULL if no other covariates
##' @param alpha_ext: coefficient from unadjusted model on external data
##' @param sd_alpha_ext: sd of alpha_ext from GWAS summary statistic
##' @param rho_shared: proportion of shared samples between internal and external data
##' @returns calibrated estimator for logistic regression in trio data
##' @export
calibrated_est_trio_logistic_2 <- function(data_int,Z,rho_shared, alpha_ext, sd_alpha_ext){
  Y_int <- data_int[,1]
  X_1_int <- data_int[,2]
  X_2_int <- data_int[,3]
  n <- nrow(data_int)
  ##internal model
  mod_int <- glm(Y_int ~ cbind(X_1_int,X_2_int,Z), family = binomial)
  ##external model
  mod_unadjusted <- glm(Y_int ~ cbind(X_1_int,Z), family = binomial)
  #coef_adjusted <- c(mod_int$coef[2], mod_int$coef[3])
  #coef <- mod_ext$coef[2]

  X <- cbind(rep(1,n), X_1_int, X_2_int,Z)
  X_tilde <- cbind(rep(1,n), X_1_int,Z)
  p <- ncol(X)-3
  # this returns (X_i*X_i^T) which is a 3*3 matrix as the i-th column of out.prod1
  out.prod1 <- apply(X,1,function(x){outer(x,x)})
  # this returns (\tilde{X}_i*\tilde{X}_i^T) which is a 2*2 matrix as the i-th column of out.prod2
  index_rem <- c(seq(3,(p+3)^2,(p+3)),seq(2*(p+3)+1,3*(p+3),1))
  out.prod2 <- out.prod1[-index_rem,]
  D_1_hat <- - (1/n)*matrix(rowSums(out.prod1%*%diag(mod_int$fitted.values*(1 - mod_int$fitted.values))),(p+3),(p+3))
  D_2_hat <- - (1/n)*matrix(rowSums(out.prod2%*%diag(mod_unadjusted$fitted.values*(1 - mod_unadjusted$fitted.values))),(p+2),(p+2))
  C_11_hat <- (1/n)*matrix(rowSums(out.prod1%*%diag((data_int[,1] - mod_int$fitted.values)^2)),(p+3),(p+3))
  C_22_hat <- (1/n)*matrix(rowSums(out.prod2%*%diag((data_int[,1] - mod_unadjusted$fitted.values)^2)),(p+2),(p+2))
  out.prod3 <- out.prod1[-c((2*(p+3)+1):(3*(p+3))),]
  C_12_hat <- (1/n)*matrix(rowSums(out.prod3%*%diag((data_int[,1] - mod_int$fitted.values)*
                                                      (data_int[,1] - mod_unadjusted$fitted.values))),(p+3),(p+2),byrow = "F")
  D_1_inv <- solve(D_1_hat)
  D_2_inv <- solve(D_2_hat)
  a <-  D_1_inv %*% C_11_hat %*% D_1_inv
  b <- D_1_inv %*% C_12_hat %*% D_2_inv
  c <- D_2_inv %*% C_22_hat %*% D_2_inv

  V_1 <- cbind(a[1:3,1:3],b[1:3,1:2],rho_shared*b[1:3,2])
  V_2 <- cbind(t(b[1:3,1:2]),c[1:2,1:2],rho_shared*c[1:2,2])
  V_3 <- c(as.vector(rho_shared*b[1:3,2]),as.vector(rho_shared*c[1:2,2]),n*sd_alpha_ext^2)
  V <- rbind(V_1,V_2,V_3)
  A <- rbind(c(0,1,0,0,0,0),c(0,0,0,0,1,-1))
  V_hat <- A%*%V%*%t(A)
  tau_raw <- as.vector(mod_int$coef)[2]
  tau_cal <-  tau_raw -
    V_hat[1,2]/V_hat[2,2]*(as.vector(mod_unadjusted$coef)[2] - alpha_ext)
  var_tau_cal <- (1/n)*(V_hat[1,1] - V_hat[1,2]^2/V_hat[2,2])
  return (list("Raw Estimator" = tau_raw,"raw_variance" = V_hat[1,1]/n,"calibrated_est" = tau_cal,"variance" = var_tau_cal))
}
