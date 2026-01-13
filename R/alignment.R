#' Align GWAS and eQTL Statistics
#'
#' Aligns GWAS alleles to eQTL alleles, handling strand flips and frequency checks.
#'
#' @param gwas GWAS data.table (Must contain SNP, A1, A2, B, SE, P, FREQ)
#' @param eqtl eQTL data.table (Must contain SNP, GENE, A1, A2, B, SE, P, FREQ)
#' @param freq_thresh Threshold for allele frequency difference check (default 0.2)
#' @return Merged and aligned data.table
#' @import data.table
#' @export
align_stats <- function(gwas, eqtl, freq_thresh = 0.2) {
    check_auth()

    # Merge on SNP
    # Utilizing data.table merge for speed
    # Inner join: we only care about intersection
    merged <- merge(eqtl, gwas, by = "SNP", suffixes = c("_EQTL", "_GWAS"))

    if (nrow(merged) == 0) {
        return(merged)
    }

    # Helper: Complement base
    complement <- function(base) {
        if (length(base) == 0) {
            return(character(0))
        }
        # Vectorized replacement
        map <- c(
            "A" = "T", "T" = "A", "C" = "G", "G" = "C",
            "a" = "t", "t" = "a", "c" = "g", "g" = "c"
        )
        return(map[base])
    }

    # We need to determine alignment status for each row
    # Python smrpy logic:
    # 1. Exact Match: A1_e == A1_g & A2_e == A2_g -> Keep, B_aligned = B_eqtl
    # 2. Flip Match: A1_e == A2_g & A2_e == A1_g -> Keep, B_aligned = -B_eqtl
    # 3. Complement Match (if not palindromic): e.g. A/C vs T/G

    # Initialize outcome columns
    merged[, `:=`(
        B_EQTL_ALIGNED = NA_real_,
        ALIGN_STATUS = "FAIL"
    )]

    # 1. Direct Match
    # A1_E == A1_G & A2_E == A2_G
    merged[A1_EQTL == A1_GWAS & A2_EQTL == A2_GWAS, `:=`(
        B_EQTL_ALIGNED = B_EQTL,
        ALIGN_STATUS = "OK"
    )]

    # 2. Reverse Match (Flip)
    # A1_E == A2_G & A2_E == A1_G
    merged[A1_EQTL == A2_GWAS & A2_EQTL == A1_GWAS, `:=`(
        B_EQTL_ALIGNED = -B_EQTL,
        ALIGN_STATUS = "FLIP"
    )]

    # 3. Complement (Strand Flip)
    # Only consider rows not yet aligned
    # Calculate complements for GWAS alleles
    # Note: This is computationally expensive if loop, use vectorized logic

    # Create temp columns for check
    # We only check if ALIGN_STATUS is FAIL

    # Define Palindromic SNPs (A/T or C/G) - ambiguous if freq is near 0.5
    merged[, is_palindromic := (A1_GWAS == "A" & A2_GWAS == "T") |
        (A1_GWAS == "T" & A2_GWAS == "A") |
        (A1_GWAS == "C" & A2_GWAS == "G") |
        (A1_GWAS == "G" & A2_GWAS == "C")]

    # Attempt Complement Match logic for FAIL rows
    # C_A1_G = complement(A1_GWAS)
    # C_A2_G = complement(A2_GWAS)

    # Since R doesn't have easy vector mapping for strings like Python's dict.get on vector easily,
    # We use chartr for complement
    merged[, A1_GWAS_COMP := chartr("ATCGatcg", "TAGCtagc", A1_GWAS)]
    merged[, A2_GWAS_COMP := chartr("ATCGatcg", "TAGCtagc", A2_GWAS)]

    # Check 3a: Complement Exact
    # A1_E == C_A1_G & A2_E == C_A2_G
    idx_comp_ok <- merged$ALIGN_STATUS == "FAIL" &
        merged$A1_EQTL == merged$A1_GWAS_COMP &
        merged$A2_EQTL == merged$A2_GWAS_COMP

    merged[idx_comp_ok, `:=`(
        B_EQTL_ALIGNED = B_EQTL,
        ALIGN_STATUS = "COMP_OK"
    )]

    # Check 3b: Complement Flip
    # A1_E == C_A2_G & A2_E == C_A1_G
    idx_comp_flip <- merged$ALIGN_STATUS == "FAIL" &
        merged$A1_EQTL == merged$A2_GWAS_COMP &
        merged$A2_EQTL == merged$A1_GWAS_COMP

    merged[idx_comp_flip, `:=`(
        B_EQTL_ALIGNED = -B_EQTL,
        ALIGN_STATUS = "COMP_FLIP"
    )]

    # Frequency Filtering for Palindromic SNPs
    # If is_palindromic, we must check frequency difference to resolve ambiguity
    # or drop if freq ~ 0.5

    # If status matches but is palindromic, verify FREQ
    # Let's assume input has FREQ. If not, can't filter, warning?
    if ("FREQ_GWAS" %in% names(merged) && "FREQ_EQTL" %in% names(merged)) {
        # Ambiguous freq range usually 0.4 - 0.6? Or using diff check.
        # Python code used freq_thresh (0.2) difference check.

        # Calculate freq difference
        # Note: If FLIP/COMP_FLIP, the frequencies are for reference allele.
        # We need to ensure we compare the SAME allele's frequency.
        # A1_GWAS vs A1_EQTL

        # If OK or COMP_OK: A1 matches A1. Diff = |F_g - F_e|
        # If FLIP or COMP_FLIP: A1 matches A2. Diff = |F_g - (1 - F_e)|

        merged[, F_diff := abs(FREQ_GWAS - FREQ_EQTL)]
        merged[ALIGN_STATUS %in% c("FLIP", "COMP_FLIP"), F_diff := abs(FREQ_GWAS - (1 - FREQ_EQTL))]

        # Filter out palindromic SNPs with ambiguous freq or high diff
        # Usually drop palindromic if |Freq - 0.5| < threshold? Or strictly diff?
        # SMR implementation uses diff threshold Check.

        # Mark FAIL if is_palindromic AND difference is too large > freq_thresh
        # Also if simple palindromic and no freq info, we should strictly drop?
        # Here we drop if diff > thresh
        merged[is_palindromic == TRUE & F_diff > freq_thresh, ALIGN_STATUS := "FAIL_FREQ"]
    }

    # Remove Failures
    final <- merged[ALIGN_STATUS != "FAIL" & ALIGN_STATUS != "FAIL_FREQ"]

    # Cleanup temp cols
    cols_to_remove <- c("A1_GWAS_COMP", "A2_GWAS_COMP", "is_palindromic", "F_diff", "ALIGN_STATUS")
    cols_to_remove <- intersect(cols_to_remove, names(final))
    if (length(cols_to_remove) > 0) {
        final[, (cols_to_remove) := NULL]
    }

    return(final)
}
