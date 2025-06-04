##' Format the data for calibration (shared across SNPs)
##'
##' @param pheno: N-vector of individual phenotype
##' @param F_ind: N-vector of invidual family index
##' @param Z: N by p matrix of covariates
##'
##' @return a list of following quantities
##' (1) N: number of individuals
##' (2) K: number of families
##' (3) Y: N-vector of phenotype
##' (4) Y_tilde: N vector of normalized phenotype
##' (5) Z: N by p matrix of covariates
##' (6) Z_tilde: N by p matrix of normalized covariates
##' (7) F_ind: N vector of family index
##' (8) size_dic: K vector of family size of each family
##' (9) F_size: N vector of family size of each individual
##' (10) dim_Z: dimension of Z
##' (11) family2ind: a K-list such that the kth element contains indices of
##' indicies from family k
##' @export
format_data <- function(pheno, F_ind, Z = NULL){
  Z_tilde = NULL

  F_ind = match(F_ind, unique(F_ind))
  # Number of families
  K = max(F_ind)
  # Number of individuals
  N = length(F_ind)

  size_dic <- as.integer(table(F_ind))
  F_size <- size_dic[F_ind]

  family2ind = list()
  for(i in 1:K){
    family2ind[[i]] = which(F_ind == i)
  }

  if(1 %in% table(F_ind) ){
    warning("Some families only have 1 individual")
  }
  pheno_tilde <- pheno
  family_means <- tapply(pheno_tilde, F_ind, mean)
  pheno_tilde <- pheno_tilde - family_means[F_ind]

  if(! is.null(Z)){
    Z = as.matrix(model.matrix( ~ -1 + ., data = Z))
    Z_tilde = Z
    for(k in 1:ncol(Z)){
      kth_mean = tapply(Z[, k], F_ind, mean)
      Z_tilde[, k] = Z[, k] - kth_mean[F_ind]
    }
  }

  dim_Z = 0
  if(! is.null(Z)){
    dim_Z = length(Z) / N
  }

  return(list(N = length(F_ind), K = max(F_ind),
              Y = pheno, Y_tilde = as.vector(pheno_tilde),  Z = Z, Z_tilde = Z_tilde,
              F_ind = F_ind,  size_dic = size_dic, F_size = F_size, dim_Z = dim_Z,
              family2ind = family2ind))
}

############################ Linear Regression ##################################


##' Compute the calibrated estimator
##' @param X: a vector that contains transmitted allele values
##' @param data: output from format_data
##' @param alpha_ext: estimate of alpha from external data
##' @param alpha_ext_var: variance of alpha from external data
##' @param N_ext: number of samples in external data
##' @param overlap_ratio: proportion of internal data from external data
##'
##' @return a list of following elements
##' (1) beta_cal: calibrated estimate
##' (2) beta_cal_var: calibrated variance
##' (3) beta_int: within-family estimate
##' (4) beta_int_var: within-family variance
##' (5) alpha_int: marginal estimate
##' (6) alpha_int_var: marginal variance
##' @export
calibrated_estimator <- function(X, data, alpha_ext, alpha_ext_var, N_ext,
                                 overlap_ratio = 0){
  with(data, {
    # Adjust external variance
    C11 = (alpha_ext_var / N_ext) * N

    # Compute correct model

    X = matrix(X, N, 1)
    X_tilde = X - matrix(tapply(X, F_ind, mean)[F_ind], N, 1)

    XZ = cbind(X, Z)
    XZ_tilde = cbind(X_tilde, Z_tilde)


    correct_int = lm(Y_tilde ~ -1 +  XZ_tilde)
    beta_int = summary(correct_int)$coefficients[1, 1]

    # Compute the incorrect model for internal data

    incorrect_int = lm( Y ~  XZ)
    alpha_int = summary(incorrect_int)$coefficients[2, 1]
    # Compute the C_{22}, C_{33}, C_{23}

    M = max(F_size)
    Cs = list()
    # Consider each family size separately
    if(M == 2){
      resid = resid(summary(incorrect_int))

      resid_tilde = resid(summary(correct_int))


      C = cbind(XZ_tilde * resid_tilde, resid, XZ * resid)

      C = t(sapply(1:K,
                   function(x) colSums(C[family2ind[[x]], ])))


      Cs[[2]] = cov(C)

    }
    else {
      for(m in 2:M){

        # Only look at family of size m
        XZ_sub = XZ[F_size == m, ]
        resid_sub = resid(summary(incorrect_int))[F_size == m]

        XZ_tilde_sub = XZ_tilde[F_size == m, ]
        resid_tilde_sub = resid(summary(correct_int))[F_size == m]


        ind_sub = F_ind[F_size == m]

        C = cbind(XZ_tilde_sub * resid_tilde_sub, resid_sub, XZ_sub * resid_sub)

        C = t(sapply(unique(ind_sub),
                     function(x) colSums(C[family2ind[[x]], ])))


        Cs[[m]] = cov(C)

      }}


    final_C <- Reduce("+", Cs[size_dic[1:K]])


    ## Compute covariance components





    if(! is.null(Z)){
      qr_decomp <- correct_int$qr
      R <- qr.R(qr_decomp)
      R_inv <- backsolve(R, diag(ncol(R)))
      const1 <-  R_inv %*% t(R_inv)
      C22 = N* (const1 %*% final_C[(1:(dim_Z + 1)), 1:(dim_Z + 1)] %*% const1)[1, 1]

      qr_decomp <- incorrect_int$qr
      R <- qr.R(qr_decomp)
      R_inv <- backsolve(R, diag(ncol(R)))
      const2 <- R_inv %*% t(R_inv)
      C33 = final_C[(dim_Z + 2):dim(final_C)[2],  (dim_Z + 2):dim(final_C)[2]]
      C33 = (const2 %*% C33 %*% const2 )[2, 2] * N


      C23 = final_C[1:(dim_Z + 1),  (dim_Z + 2):dim(final_C)[2]]
      C23 = ( (const1) %*% C23 %*% const2 )[1, 2] * N

    } else{
      const1 = solve(t(XZ_tilde) %*% XZ_tilde)
      C22 = N*  final_C[1, 1] * const1^2

      C33 = final_C[(dim_Z + 2):dim(final_C)[2],  (dim_Z + 2):dim(final_C)[2]]
      design = t(cbind(rep(1, N), XZ)) %*% cbind(rep(1, N), XZ)
      C33 = (solve(design) %*% C33 %*% solve(design) )[2, 2] * N


      C23 = final_C[1:(dim_Z + 1),  (dim_Z + 2):dim(final_C)[2]]
      C23 = ((const1) %*% matrix(C23, 1, 2) %*% solve(design) )[1, 2] * N


    }



    # Compute C12, C13

    C12 = overlap_ratio * C23 * N / N_ext
    C13 = overlap_ratio * C33 * N / N_ext

    # Covariance matrix of (alpha_ext - alpha_int, beta_int - beta)
    result_cov = rbind(c(C11 - 2 * C13+ C33, C12 - C23), c(C12 - C23, C22))
    # Compute the final calibrated estimator
    beta_cal = beta_int +  (C23 - C12) / ( C11 + C33 - 2 * C13) * (alpha_ext - alpha_int)
    beta_cal_var = (C22 - (C23 - C12)^2 / (C11 + C33 - 2 * C13) ) / N

    # Compute internal false model by sampling one sibling from each family
    if (is.null(sib_ind)){

      grouped_indices <- split(seq_along(F_ind), F_ind)

      # Sample one index from each group
      sib_ind <- as.numeric( sapply(grouped_indices, sample, size = 1) )
    }


    return(list(beta_cal  = beta_cal, beta_cal_var = (beta_cal_var), beta_int = beta_int,
                beta_int_var = (C22/N), alpha_int = alpha_int, alpha_int_var = (C33/N)) )
  }


  )
}


######################### Logistic Regression ###################################


##' Compute the calibrated estimator
##' @param X: a vector that contains transmitted allele values
##' @param data: output from format_data
##' @param alpha_ext: estimate of alpha from external data
##' @param alpha_ext_var: variance of alpha from external data
##' @param N_ext: number of samples in external data
##' @param overlap_ratio: proportion of internal data from external data
##'
##' @return a list of following elements
##' (1) beta_cal: calibrated estimate
##' (2) beta_cal_var: calibrated variance
##' (3) beta_int: within-family estimate
##' (4) beta_int_var: within-family variance
##' (5) alpha_int: marginal estimate
##' (6) alpha_int_var: marginal variance
##' @export
calibrated_logistic_estimator <- function(X, data, alpha_ext, alpha_ext_var, N_ext,
                                          overlap_ratio = 0){
  with(data, {
    # Adjust external variance
    C11 = (alpha_ext_var / N_ext) * N

    # Process data
    X = matrix(X, N, 1)
    family_effect = tapply(X, F_ind, sum)
    f = family_effect[F_ind]
    XZ = cbind(X, Z)

    # Compute mis-specified model
    int_mis = glm(Y ~ XZ, family = binomial(link = "logit"))
    alpha_int = as.numeric(int_mis$coefficients[2])
    predict_mis = logit(predict(int_mis))


    # Compute the correct model for internal data
    int_cor = glm(Y ~ XZ + f, family = binomial(link = "logit"))
    beta_int = as.numeric(int_cor$coefficients[2])
    predict_cor = logit(predict(int_cor))


    # Compute the variance

    XZ = cbind(1, XZ)
    XZ_tilde =cbind(XZ, f)

    # Efficiently create diagonal matrices using Matrix package
    coef_cor <- Diagonal(x = predict_cor * (1 - predict_cor))
    coef_mis <- Diagonal(x = predict_mis * (1 - predict_mis))

    # Efficient matrix multiplications and scaling
    D2 <- crossprod(XZ_tilde, coef_cor %*% XZ_tilde) / N
    D3 <- crossprod(XZ, coef_mis %*% XZ) / N



    # Compute the C_{22}, C_{33}, C_{23}

    M = max(F_size)

    correct_var = list()
    incorrect_var = list()
    correct_incorrect_cov = list()

    dim_XZ = dim(XZ)[2]
    # Consider each family size separately
    if(M == 2){
      prod = (Y - predict_mis) * XZ
      prod_tilde =  (Y - predict_cor) * XZ_tilde


      num_F_ind_levels <- length(unique(F_ind))
      est <- matrix(0, nrow = num_F_ind_levels, ncol = (dim_XZ + 1) + dim_XZ)

      # Compute the first part of est
      for(i in 1:(dim_XZ + 1)){
        est[, i] <- tapply(prod_tilde[, i], F_ind, sum)
      }

      # Compute the second part of est
      for(i in 1:dim_XZ){
        est[, (dim_XZ + 1 + i)] <- tapply(prod[, i], F_ind, sum)
      }

      # Calculate the covariance of est
      temp <- cov(est)


      correct_var = temp[1:(dim_XZ + 1), 1:(dim_XZ + 1)]

      # This is for C33
      incorrect_var = temp[(dim_XZ+2): (2 * dim_XZ + 1), (dim_XZ+2): (2 * dim_XZ + 1)]

      # This is for C23
      correct_incorrect_cov = temp[1:(dim_XZ + 1), (dim_XZ+2):(2 * dim_XZ + 1)]

      C22 =  correct_var * K / N
      C33 = incorrect_var * K / N
      C23 =  correct_incorrect_cov  * K / N


    } else{
      for(m in 2:M){

        # Only look at family of size m
        XZ_sub = XZ[F_size == m, ]
        predict_mis_sub = predict_mis[F_size == m]

        XZ_tilde_sub = XZ_tilde[F_size == m, ]
        predict_cor_sub = predict_cor[F_size == m]

        Y_sub = Y[F_size == m]

        prod = (Y_sub - predict_mis_sub) * XZ_sub
        prod_tilde =  (Y_sub - predict_cor_sub) * XZ_tilde_sub


        # This is for C22
        est = NULL

        for(i in 1: (dim_XZ + 1)){
          est = cbind(est, tapply(prod_tilde[,i], F_ind[F_size == m], sum))
        }
        for(i in 1:dim_XZ){
          est = cbind(est,  tapply(prod[,i], F_ind[F_size == m], sum))
        }

        temp = cov(est)


        correct_var[[m]] = temp[1:(dim_XZ + 1), 1:(dim_XZ + 1)]

        # This is for C33
        incorrect_var[[m]] = temp[(dim_XZ+2): (2 * dim_XZ + 1), (dim_XZ+2): (2 * dim_XZ + 1)]

        # This is for C23
        correct_incorrect_cov[[m]] = temp[1:(dim_XZ + 1), (dim_XZ+2):(2 * dim_XZ + 1)]
      }
      C22 =  Reduce("+", correct_var[size_dic[1:K]]) / N
      C33 = Reduce("+", incorrect_var[size_dic[1:K]]) / N
      C23 = Reduce("+", correct_incorrect_cov[size_dic[1:K]]) / N
    }





    C22 = ( solve(D2) %*% C22 %*% solve(D2) )[2,2]
    C33 = (solve(D3) %*% C33 %*% solve(D3) )[2,2]
    C23 = (solve(D2) %*% C23 %*% solve(D3))[2,2]

    # Compute C12, C13

    C12 = overlap_ratio * C23 * N / N_ext
    C13 = overlap_ratio * C33 * N / N_ext

    # Covariance matrix of (alpha_ext - alpha_int, beta_int - beta)
    result_cov = rbind(c(C11 - 2 * C13+ C33, C12 - C23), c(C12 - C23, C22))
    # Compute the final calibrated estimator
    beta_cal = beta_int +  (C23 - C12) / ( C11 + C33 - 2 * C13) * (alpha_ext - alpha_int)
    beta_cal_var = (C22 - (C23 - C12)^2 / (C11 + C33 - 2 * C13) ) / N

    sub_ind = sample_indices(F_ind, max(F_ind))

    int_mis = glm(Y[sub_ind] ~ XZ[sub_ind, ], family = binomial(link = "logit"))
    alpha_int_single = summary(int_mis)$coefficients[2,1]
    alpha_int_var_single = summary(int_mis)$coefficients[2,2]^2

    return(list(beta_cal  = beta_cal, beta_cal_var = (beta_cal_var), beta_int = beta_int,
                beta_int_var = (C22/N), alpha_int = alpha_int, alpha_int_var = (C33)/N))
  } )
}
