#' SMR Test (Wald Test)
#'
#' @param b_eqtl Effect size of instrument on gene expression (eQTL)
#' @param se_eqtl Standard error of b_eqtl
#' @param b_gwas Effect size of instrument on trait (GWAS)
#' @param se_gwas Standard error of b_gwas
#' @return A list containing b_smr, se_smr, and p_smr
#' @export
smr_test <- function(b_eqtl, se_eqtl, b_gwas, se_gwas) {
    check_auth()
    # Estimate SMR effect (Alpha = Beta_GWAS / Beta_eQTL)
    b_smr <- b_gwas / b_eqtl

    # Estimate SE of SMR effect using Delta Method
    se_smr <- sqrt((se_gwas^2 * b_eqtl^2 + se_eqtl^2 * b_gwas^2) / b_eqtl^4)

    z_stat <- b_smr / se_smr
    p_smr <- 2 * pnorm(abs(z_stat), lower.tail = FALSE)

    return(list(b_SMR = b_smr, se_SMR = se_smr, p_SMR = p_smr))
}

#' HEIDI Test (Heterogeneity in Dependent Instruments)
#'
#' This implementation perfectly mirrors the C++ SMR v1.4.0 and Python smrpy logic.
#'
#' @param b_eqtl Vector of eQTL betas
#' @param se_eqtl Vector of eQTL SEs
#' @param b_gwas Vector of GWAS betas
#' @param se_gwas Vector of GWAS SEs
#' @param ld_matrix LD correlation matrix (r) between SNPs
#' @return A list containing p_heidi and nsnp_heidi
#' @export
heidi_test <- function(b_eqtl, se_eqtl, b_gwas, se_gwas, ld_matrix) {
    check_auth()
    m <- length(b_eqtl)
    if (m < 2) {
        return(list(p_heidi = NA, nsnp_heidi = m))
    }

    # 1. Identify Top SNP (based on eQTL significance Z-score)
    # Note: using abs(z) to find max
    z_eqtl <- b_eqtl / se_eqtl
    top_idx <- which.max(abs(z_eqtl))

    # SMR estimates
    b_smr <- b_gwas / b_eqtl

    # 2. Estimate Covariance Matrix of b_smr (Cov(b_xy))
    # Formula from SMR C++ (SMR_data.cpp: est_cov_bxy)
    # Cov_ij = R_ij * [ (se_Gi * se_Gj)/(b_ei * b_ej) + (b_SMRi * b_SMRj)/(Z_ei * Z_ej) ]
    #          - (b_SMRi * b_SMRj)/(Z_ei^2 * Z_ej^2)

    # Pre-calculate outer products for vectorized calculation
    # term1_part = (se_gwas * se_gwas') / (b_eqtl * b_eqtl')
    term1_num <- outer(se_gwas, se_gwas)
    term1_den <- outer(b_eqtl, b_eqtl)
    term1 <- term1_num / term1_den

    # term2_part = (b_smr * b_smr') / (z_eqtl * z_eqtl')
    term2_num <- outer(b_smr, b_smr)
    term2_den <- outer(z_eqtl, z_eqtl)
    term2 <- term2_num / term2_den

    # term3 = (b_smr * b_smr') / (z_eqtl^2 * z_eqtl'^2)
    z_eqtl_sq <- z_eqtl^2
    term3_den <- outer(z_eqtl_sq, z_eqtl_sq)
    term3 <- term2_num / term3_den

    # Combine terms with LD matrix
    # Note: ld_matrix is R_ij
    cov_bxy <- ld_matrix * (term1 + term2) - term3

    # 3. Calculate Deviations (d)
    # d_i = b_smr_top - b_smr_i
    # We work with the subset of SNPs excluding the top SNP
    keep_idx <- setdiff(seq_len(m), top_idx)

    if (length(keep_idx) == 0) {
        return(list(p_heidi = NA, nsnp_heidi = 1))
    }

    d <- b_smr[top_idx] - b_smr[keep_idx]

    # 4. Construct Variance Matrix of d (V_dev)
    # V_ij = Cov(top - i, top - j)
    #      = Var(top) + Cov(i, j) - Cov(top, i) - Cov(top, j)

    var_top <- cov_bxy[top_idx, top_idx]
    cov_ij <- cov_bxy[keep_idx, keep_idx]
    cov_top_i <- cov_bxy[top_idx, keep_idx] # vector of length (m-1)

    # Create matrices for top covariances to vectorize subtraction
    n_sub <- length(keep_idx)
    mat_cov_top_i <- matrix(rep(cov_top_i, n_sub), nrow = n_sub, ncol = n_sub, byrow = TRUE) # Rows are same
    mat_cov_top_j <- matrix(rep(cov_top_i, n_sub), nrow = n_sub, ncol = n_sub, byrow = FALSE) # Cols are same (transpose)

    V_dev <- var_top + cov_ij - mat_cov_top_i - mat_cov_top_j

    # Add small ridge for stability (matches C++ and Python implementation)
    diag(V_dev) <- diag(V_dev) + 1e-8

    # 5. Calculate Chi-square Statistic (Weighted Sum of Chi-squares)
    # Statistic is sum(d_i^2 / V_ii) -> Note: Only diagonal is used for the statistic sum!
    # This comes from SMR_data.cpp logic where dev[i] is divided by tmp3[i] (variance)

    q_heidi <- sum(d^2 / diag(V_dev))

    # 6. Calculate P-value using Satterthwaite Approximation
    # We need eigenvalues of the *correlation* matrix of V_dev

    # Convert V_dev to Correlation Matrix
    sds <- sqrt(diag(V_dev))
    R_dev <- V_dev / outer(sds, sds)

    # Get eigenvalues of the correlation matrix
    # Use 'eigen' with symmetric=TRUE
    ev <- eigen(R_dev, symmetric = TRUE, only.values = TRUE)$values

    # Filter out small/negative eigenvalues due to precision
    ev <- ev[ev > 1e-8]

    # Satterthwaite approximation
    # Distribution of Q approx a * chi^2(df)
    # a = sum(lambda^2) / sum(lambda)
    # df = (sum(lambda))^2 / sum(lambda^2)
    # P = P(chi^2(df) > q_heidi / a)

    sum_ev <- sum(ev)
    sum_ev_sq <- sum(ev^2)

    if (sum_ev < 1e-10 || sum_ev_sq < 1e-10) {
        # Fallback if eigenvalues imply zero variance
        return(list(p_heidi = 1.0, nsnp_heidi = m))
    }

    a <- sum_ev_sq / sum_ev
    df <- (sum_ev^2) / sum_ev_sq

    p_val <- pchisq(q_heidi / a, df = df, lower.tail = FALSE)

    return(list(p_heidi = p_val, nsnp_heidi = m))
}
