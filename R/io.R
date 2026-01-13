#' Load GWAS Summary Statistics
#' @param path Path to file
#' @param verbose Print messages
#' @return data.table
#' @import data.table
#' @export
load_gwas <- function(path, verbose = TRUE) {
    check_auth()
    if (verbose) message(paste("Reading GWAS:", path))
    df <- data.table::fread(path)

    # Standardize
    cols <- colnames(df)
    map <- c(
        "p_val" = "P", "p" = "P", "beta" = "B", "b" = "B", "se" = "SE",
        "bp" = "POS", "chrom" = "CHR", "ref" = "A2", "alt" = "A1",
        "rsid" = "SNP", "freq" = "FREQ"
    )

    for (orig in names(map)) {
        if (orig %in% cols && !map[orig] %in% cols) setnames(df, orig, map[orig])
    }
    return(df)
}

#' Load eQTL Summary Statistics
#' @param path Path to file
#' @param verbose Print messages
#' @return data.table
#' @import data.table
#' @export
load_eqtl <- function(path, verbose = TRUE) {
    check_auth()
    if (verbose) message(paste("Reading eQTL:", path))
    df <- data.table::fread(path)

    cols <- colnames(df)
    map <- c(
        "gene_id" = "GENE", "pval_nominal" = "P", "slope" = "B", "slope_se" = "SE",
        "rsid" = "SNP", "chromosome" = "CHR", "position" = "POS"
    )

    for (orig in names(map)) {
        if (orig %in% cols && !map[orig] %in% cols) setnames(df, orig, map[orig])
    }

    if ("tss_distance" %in% cols && "POS" %in% cols) {
        df[, TSS := POS - tss_distance]
    }
    return(df)
}

#' Get LD Matrix from PLINK (Robust genio Implementation)
#'
#' Uses the C++ optimized 'genio' package to read binary PLINK files efficiently.
#'
#' @param bfile_prefix Prefix of .bed/.bim/.fam
#' @param snp_list Vector of SNPs to extract
#' @return Correlation matrix (r)
#' @import genio
#' @export
get_ld_matrix <- function(bfile_prefix, snp_list) {
    check_auth()
    bim_path <- paste0(bfile_prefix, ".bim")
    bed_path <- paste0(bfile_prefix, ".bed")
    fam_path <- paste0(bfile_prefix, ".fam")

    if (!file.exists(bim_path) || !file.exists(bed_path) || !file.exists(fam_path)) {
        warning(paste("LD reference files not found:", bfile_prefix))
        return(NULL)
    }

    # 1. Read BIM to map SNPs to Indices
    # genio::read_bim is fast
    bim <- genio::read_bim(bim_path, verbose = FALSE)

    # 2. Read FAM to get number of individuals
    fam <- genio::read_fam(fam_path, verbose = FALSE)
    n_ind <- nrow(fam)
    m_loci <- nrow(bim)

    # Match SNPs
    # bim$id contains the rsIDs
    target_indices <- which(bim$id %in% snp_list)

    if (length(target_indices) == 0) {
        return(NULL)
    }

    # Sort indices ensures sequential seeking optimization in some readers
    target_indices <- sort(target_indices)
    found_snps <- bim$id[target_indices]

    # 3. Read BED subset
    # genio::read_bed does not support subsetting in some versions, so we read full
    # For very large files this is inefficient, but for typical 1000G ref split by Chr it's acceptable (10-50MB)
    X_full <- genio::read_bed(bed_path, n_ind = n_ind, m_loci = m_loci, verbose = FALSE)

    # Subset rows (variants)
    # genio returns variants x samples
    X <- X_full[target_indices, , drop = FALSE]

    # Check dimensions
    if (ncol(X) == 0 || nrow(X) == 0) {
        return(NULL)
    }

    # 4. Calculate Correlation (LD)
    # X is (Variants x Samples)

    # cor() works on columns (variables). We want correlation between variants.
    # So we must transpose X -> (Samples x Variants)
    X_t <- t(X)

    # Calculate correlation dealing with missing data
    ld_r <- cor(X_t, use = "pairwise.complete.obs")

    # Assign names
    rownames(ld_r) <- found_snps
    colnames(ld_r) <- found_snps

    # Reorder to match input snp_list order
    # This is crucial for correct matrix algebra in HEIDI
    final_snps <- snp_list[snp_list %in% found_snps]
    ld_r <- ld_r[final_snps, final_snps, drop = FALSE]

    return(ld_r)
}
