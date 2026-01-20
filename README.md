# Installation

To download, use the following command

```
devtools::install_github("ZixuanWu1/CalibrationGWAS")
```

# Calibration with summary-level data

To perform calibration using summary-level data, prepare a data frame `df` containing the following columns:

- `"SNP"`: SNP identifier  
- `"beta_int"` and `"beta_int_var"`: family-adjusted estimate and its variance  
- `"alpha_int"` and `"alpha_int_var"`: internal population estimate and its variance  
- `"alpha_ext"` and `"alpha_ext_var"`: external population estimate and its variance  

Optionally, one may also provide:

- `"N_int"` and `"N_ext"`: internal and external sample sizes  

If these sample sizes are not provided, the summary-level calibration assumes shared T and NT correlations for trio data, or shared G and G_sib correlations for sibling data.

The calibrated estimator

$$
\hat{\beta}_{\text{int}} + \lambda \bigl(\hat{\alpha}_{\text{ext}} - \hat{\alpha}_{\text{int}}\bigr)
$$

can be computed using:

```r
calibrated_df = calibration_summary_stats(
  df,
  rho = NA,
  family = "trio",
  method = "empirical",
  r = NA,
  pheno_cor = NA,
  family_size = 2,
  kappa = 0.5
)
```

Here, `rho` represents the correlation between alpha_int and alpha_ext, and will be automatically estimated if not provided.  
The argument `family` can be either `"trio"` or `"sibling"`.

By default, `method = "empirical"` estimates the optimal lambda from the data, assuming shared T and NT correlations in trio designs or shared G and G_sib correlations in sibling designs.

If `method = "theoretical"`, lambda is chosen based on its theoretical expression.

- For trio designs, one may provide `r` to allow different T and NT correlations, where `r` is either 0 or a vector specifying the correlation between T_j and NT_j.  
- For sibling designs, one must provide `pheno_cor` (phenotypic correlation between siblings), `family_size` (number of siblings per family), and the G and G_sib correlations.

The function returns a data frame with the following columns:

1. `SNP`: SNP identifier  
2. `beta_cal`: calibrated estimate  
3. `beta_cal_var`: estimated variance of the calibrated estimator  
4. `lambda`: the selected lambda for each SNP  

# Calibration with individual-level data


## Trio data

One needs to prepare a data frame `data_int` with the following columns:

- `Y`: phenotype  
- `X1`: T  
- `X2`: NT or `Gpa = (T + NT)`  

Optionally, one may provide a matrix of additional covariates `Z`, such as age and gender.

The family-adjusted and calibrated estimators can then be computed using:

```r
res = calibrated_est_trio(
  data_int,
  alpha_ext,
  alpha_ext_sd,
  Z = NULL,
  rho_shared = 0,
  family_info = "NT",
  model = "linear"
)
```

Here:

- `alpha_ext` and `alpha_ext_sd` are the external population estimate and its standard error.  
- `rho_shared` represents $n_{\text{shared}} / N_{\text{ext}}$, where $n_{\text{shared}}$ is the number of overlapping samples between the internal and external datasets, and $N_{\text{ext}}$ is the external sample size.  
- `family_info` can be either `"NT"` or `"Gpa"`.  
- `model` can be either `"linear"` or `"logistic"`.

The function returns a list containing:

1. `"raw_estimator"`: uncalibrated family-adjusted estimator  
2. `"raw_variance"`: variance of the uncalibrated estimator  
3. `"calibrated_est"`: calibrated estimator  
4. `"variance"`: variance of the calibrated estimator  

## Sibling data

One needs to prepare:

- a vector of phenotypes `pheno`  
- a vector of family indices `F_ind`, indicating the family membership of each individual  
- a matrix of covariates `Z` (optional)  
- a vector of genotypes `X`  
- `alpha_ext` and `alpha_ext_var`, the external population estimate and its variance  
- the external sample size `N_ext`  

Each element of these vectors corresponds to one individual.

First, format the internal data:

```r
data = format_data(pheno, F_ind, Z)
```

Then, for **linear regression**, run:

```r
res = calibrated_est_sib(X, data, alpha_ext, alpha_ext_var, N_ext, overlap_ratio = 0)
```

For **logistic regression**, run:

```r
res = calibrated_est_sib_logistic(X, data, alpha_ext, alpha_ext_var, N_ext, overlap_ratio = 0)
```

Here, `overlap_ratio` represents $n_{\text{shared}} / N_{\text{int}}$, where  
$n_{\text{shared}}$ is the number of overlapping samples between the internal and external datasets, and  
$N_{\text{int}}$ is the internal sample size.

Each function returns a list containing:

1. `beta_cal`: calibrated estimator  
2. `beta_cal_var`: variance of the calibrated estimator  
3. `beta_int`: internal family-adjusted estimator  
4. `beta_int_var`: variance of the internal estimator  
5. `alpha_int`: internal population estimator  
6. `alpha_int_var`: variance of the internal population estimator  
```
