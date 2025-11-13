##' Estimate the correlation between external and internal gwas (i.e. sample overlap)
##'
##' @param df: a data frame with these columns: "SNP", "beta_int", "beta_int_var", "alpha_int", "alpha_int_var",
##' "alpha_ext", "alpha_ext_var", "N_int", "N_ext"
##' @param snp_list: list of SNPs that has zero direct effects (can be selected from an independent selection file)
##'
##' @return r: noise correlation between external and internal unadjusted model
noise_correlation <- function(df, snp_list = NA){

  # Select snps with zero effects
  if(! all(is.na(snp_list))){
    sub_df = df[df$SNP %in% snp_list, ]
  } else {
    sub_df = df
  }

  # Standardize estimates
  alpha_int_std = sub_df$alpha_int / sqrt(sub_df$alpha_int_var)
  alpha_ext_std = sub_df$alpha_ext / sqrt(sub_df$alpha_ext_var)

  # Compute correlation
  noise_cor = median(alpha_int_std * alpha_ext_std)

  return(noise_cor)


}


##' Estimate the correlation between external and internal trio gwas (i.e. sample overlap)
##'
##' @param df: a data frame with these columns: "SNP", "beta_int", "beta_int_var", "alpha_int", "alpha_int_var",
##' "alpha_ext", "alpha_ext_var", "N_int", "N_ext"
##' @param rho: sample overlap ratio
##' @param family: type of family structure. Either trio or sibling
##' @param r: assumption on correlation between T and NT. Can be zero, shared, or a provided vector.
##'
##' @return a data frame of snp id, calibrated beta and its standard error
calibration_summary_stats_raw <- function(df, rho = NA, family = "trio", r = "shared"){

  n = df$N_int
  N = df$N_ext

  if (is.na(rho)){
    noise_cor = noise_correlation(df)}

  # Compute value and variance of alpha_int - alpha_ext
  alpha_diff = df$alpha_ext - df$alpha_int
  alpha_diff_var = df$alpha_int_var + df$alpha_ext_var - 2 * noise_cor * sqrt(df$alpha_int_var * df$alpha_ext_var )

  # Standardize estimates
  beta_int_std = df$beta_int / sqrt(df$beta_int_var)
  alpha_diff_std = alpha_diff / sqrt(alpha_diff_var)

  if(family == "trio"){
    rho = noise_cor * sqrt(n/N)

    if(all(r == "zero")){
      beta_alpha_cor = (1 - rho) / sqrt( 2 * (1 + n/N - 2 * rho)  )
    } else if(all(r != "shared")){
      beta_alpha_cor = ( (1 - rho) * sqrt(1 - r) )/ sqrt( 2 * (1 + n/N - 2 * rho)  )
    } else {
      beta_alpha_cor = cov(beta_int_std,  alpha_diff_std)
    }
  } else{
    beta_alpha_cor = cov(beta_int_std,  alpha_diff_std)
  }

  beta_alpha_cov = beta_alpha_cor * sqrt(df$beta_int_var) * sqrt(alpha_diff_var)

  lambda = - beta_alpha_cov  / alpha_diff_var
  beta_cal = df$beta_int + lambda * alpha_diff

  beta_cal_var = df$beta_int_var + lambda^2 * alpha_diff_var + 2 * lambda * beta_alpha_cov

  return(data.frame(SNP = df$SNP, beta_cal = beta_cal, beta_cal_var = beta_cal_var))

}




##' Estimate the correlation between external and internal trio gwas (i.e. sample overlap)
##'
##' @param df: a data frame with these columns: "SNP", "beta_int", "beta_int_var", "alpha_int", "alpha_int_var",
##' "alpha_ext", "alpha_ext_var", "N_int", "N_ext" (N_int, N_ext won't be used if assumed shared r)
##' @param rho: sample overlap ratio
##' @param family: type of family structure. Either trio or sibling
##' @param method: theoretical or empirical
##' @param r: vector of correlation between T and NT for each SNP. Default is 0 (only used ifmethod = "theory" and family == "trio").
##' @param pheno_cor: scalar of phenotype correlation between siblings. (only used if method == "theory" and family == "sibling").
##' @param family_size: scalar of family size (only used if method == "theory" and family == "sibling").
##' @param kappa: vector of correlation between G and Gsib for each SNP. Default is 0 (only used if method = "theory" and family == "sibling").
##'
##' @return a data frame of snp id, calibrated beta and its standard error
calibration_summary_stats <- function(df, rho = NA, family = "trio", method =  "empirical", r = NA,  
                                      pheno_cor = NA, family_size = 2, kappa = .5){
  
  n = df$N_int
  N = df$N_ext
  
  if (is.na(rho)){
    noise_cor = noise_correlation(df)
  } else{
    noise_cor = rho
  }
  
  # Compute value and variance of alpha_int - alpha_ext
  alpha_diff = df$alpha_ext - df$alpha_int
  alpha_diff_var = df$alpha_int_var + df$alpha_ext_var - 2 * noise_cor * sqrt(df$alpha_int_var * df$alpha_ext_var )
  
  # Standardize estimates
  beta_int_std = df$beta_int / sqrt(df$beta_int_var)
  alpha_diff_std = alpha_diff / sqrt(alpha_diff_var)
  
  if(family == "trio"){
    rho = noise_cor * sqrt(n/N)
    
    if(method == "theory"){{
      if (is.na(r)){
        beta_alpha_cor = -(1 - rho) / sqrt( 2 * (1 + n/N - 2 * rho)  )
      } else{
        beta_alpha_cor = -( (1 - rho) * sqrt(1 - r) )/ sqrt( 2 * (1 + n/N - 2 * rho)  )
      }
    }} else {
      beta_alpha_cor = mean(beta_int_std* alpha_diff_std)
    }
  } else{
    
    if(method == "theory"){
      x = 1 + (family_size - 1) * pheno_cor * kappa
      beta_alpha_cor = -(N / sqrt( N * ( (N * x + n))) ) * (1 - pheno_cor) * sqrt( (family_size - 1) * (1 - kappa) / family_size) / sqrt(1 - pheno_cor/(family_size - 1))
    } else{
      beta_alpha_cor = mean(beta_int_std* alpha_diff_std)
    }
  }
  beta_alpha_cov = beta_alpha_cor * sqrt(df$beta_int_var) * sqrt(alpha_diff_var)
  
  lambda = - beta_alpha_cov  / alpha_diff_var
  beta_cal = df$beta_int + lambda * alpha_diff
  
  beta_cal_var = df$beta_int_var + lambda^2 * alpha_diff_var + 2 * lambda * beta_alpha_cov
  return(data.frame(SNP = df$SNP, beta_cal = beta_cal, beta_cal_var = beta_cal_var, lambda = lambda))
  
}
