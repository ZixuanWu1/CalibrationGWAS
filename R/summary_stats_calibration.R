##' Estimate the correlation between external and internal gwas (i.e. sample overlap)
##'
##' @param df: a data frame with these columns: "SNP", "beta_int", "beta_int_var", "alpha_int", "alpha_int_var",
##' "alpha_ext", "alpha_ext_var", "N_int", "N_ext"
##' @param snp_list: list of SNPs that has zero direct effects (can be selected from an independent selection file)
##'
##' @return rho: sample overlap ratio
sample_overlap <- function(df, snp_list = NA){

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
calibration_summary_stats <- function(df, rho, family = "trio", r = "shared"){

  n = df$N_int
  N = df$N_ext


  # Compute value and variance of alpha_int - alpha_ext
  alpha_diff = df$alpha_ext - df$alpha_int
  alpha_diff_var = df$alpha_int_var + df$alpha_ext_var - 2 * rho * sqrt(df$alpha_int_var * df$alpha_ext_var )

  # Standardize estimates
  beta_int_std = df$beta_int / sqrt(df$beta_int_var)
  alpha_diff_std = alpha_diff / sqrt(alpha_diff_var)

  if(family == "trio"){
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




# Example:
#
# Create summary stats
# library(mvtnorm)
# x = rmvnorm(100, mean = c(2, 0), sigma = rbind(c(1, .5), c(.5, 1)))
#
# df = cbind(1:100, x[,1], 1, x[,2], 1, rnorm(100, 0, 0.0001), 0.01, 100, 10000 )
#
#
# df = data.frame(df)
#
# colnames(df) = c("SNP", "beta_int", "beta_int_var", "alpha_int", "alpha_int_var",
#                  "alpha_ext", "alpha_ext_var", "N_int", "N_ext")
#
# Compute calibrated estimator
# rho = sample_overlap(df, 1:100)
# calibrated_estimator = calibration(df, rho)
