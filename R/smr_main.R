#' Run Complete SMR Analysis Pipeline
#'
#' @param gwas_path Path to GWAS summary file
#' @param eqtl_path Path to eQTL summary file
#' @param ld_reference Path to LD reference prefix (PLINK binary)
#' @param output_prefix Output filename prefix
#' @param peqtl_smr P-value threshold to select instruments (default 5e-8)
#' @param peqtl_heidi P-value threshold for HEIDI instruments (default 0.05)
#' @param cis_window_kb Cis window size in Kb (default 2000)
#' @param threads Number of threads for parallel processing (default 1)
#' @param plot Logical, whether to generate plots (default TRUE)
#' @param activation_code Optional activation code for one-time use
#' @import data.table
#' @import parallel
#' @import pbapply
#' @importFrom utils write.table
#' @export
SMR_analysis <- function(gwas_path, eqtl_path, ld_reference, output_prefix,
                         peqtl_smr = 5e-8, peqtl_heidi = 0.0157,
                         cis_window_kb = 2000, threads = 1, plot = TRUE,
                         activation_code = NULL) {
    # Check activation
    check_auth(activation_code)

    message("--- SMRr Pipeline Started ---")

    # Force progress bar visibility for Rscript
    pbapply::pboptions(type = "txt", style = 3, char = "=")

    # 1. Load Data
    gwas <- load_gwas(gwas_path)
    eqtl <- load_eqtl(eqtl_path)

    # 3. Alignment
    message("Aligning datasets...")
    data <- align_stats(gwas, eqtl)

    if (nrow(data) == 0) stop("No overlapping SNPs found.")

    message(paste("Total analyzed SNP-Gene pairs:", nrow(data)))

    # 4. Processing Loop (Gene-wise)
    genes <- unique(data$GENE)
    message(paste("Processing", length(genes), "genes..."))

    # Plot directory
    if (plot) {
        plot_dir <- paste0(output_prefix, "_plots")
        if (!dir.exists(plot_dir)) dir.create(plot_dir)
    }

    # Worker Function
    process_one_gene <- function(gene) {
        # Subset data
        gene_df <- data[GENE == gene]

        # Identify Candidates (P < peqtl_smr)
        candidates <- gene_df[P_EQTL <= peqtl_smr]

        if (nrow(candidates) == 0) {
            return(NULL)
        }

        # Top SNP Selection
        # Select min P
        top_idx <- which.min(gene_df$P_EQTL)
        if (gene_df$P_EQTL[top_idx] > peqtl_smr) {
            return(NULL)
        }

        top_snp_row <- gene_df[top_idx, ]

        # SMR Test
        res_smr <- smr_test(
            top_snp_row$B_EQTL_ALIGNED, top_snp_row$SE_EQTL,
            top_snp_row$B_GWAS, top_snp_row$SE_GWAS
        )

        # HEIDI Test
        heidi_candidates <- gene_df[P_EQTL < peqtl_heidi & SNP != top_snp_row$SNP]

        p_heidi <- NA
        nsnp_heidi <- 0

        if (nrow(heidi_candidates) >= 2) {
            # Get LD
            snp_list <- c(top_snp_row$SNP, heidi_candidates$SNP)
            ld_r <- get_ld_matrix(ld_reference, snp_list)

            if (!is.null(ld_r)) {
                # Synchronize data frames with returned LD matrix
                if (top_snp_row$SNP %in% colnames(ld_r)) {
                    valid_snps <- colnames(ld_r)
                    sub_df <- gene_df[SNP %in% valid_snps]
                    sub_df <- sub_df[match(valid_snps, sub_df$SNP)]

                    if (nrow(sub_df) >= 2) {
                        res_heidi <- heidi_test(
                            sub_df$B_EQTL_ALIGNED, sub_df$SE_EQTL,
                            sub_df$B_GWAS, sub_df$SE_GWAS, ld_r
                        )
                        p_heidi <- res_heidi$p_heidi
                        nsnp_heidi <- res_heidi$nsnp_heidi
                    }
                }
            }
        }

        # Result Row
        res_entry <- list(
            probeID = gene,
            ProbeChr = top_snp_row$CHR,
            Probe_bp = top_snp_row$POS,
            topSNP = top_snp_row$SNP,
            b_GWAS = top_snp_row$B_GWAS,
            p_GWAS = top_snp_row$P_GWAS,
            b_eQTL = top_snp_row$B_EQTL_ALIGNED,
            p_eQTL = top_snp_row$P_EQTL,
            b_SMR = res_smr$b_SMR,
            p_SMR = res_smr$p_SMR,
            p_HEIDI = p_heidi,
            nsnp_HEIDI = nsnp_heidi
        )

        # Plotting
        if (plot && res_smr$p_SMR < 0.05) {
            # 1. Effect Size Plot
            p_effect <- plot_effect_size(gene_df, res_entry)
            ggplot2::ggsave(
                filename = paste0(plot_dir, "/", gene, "_effect.png"),
                plot = p_effect, width = 7, height = 7, dpi = 300
            )

            # 2. Locus Plot
            p_locus <- plot_locus(gene_df, res_entry)
            ggplot2::ggsave(
                filename = paste0(plot_dir, "/", gene, "_locus.png"),
                plot = p_locus, width = 10, height = 8, dpi = 300
            )
        }

        return(res_entry)
    }

    # Parallel/Serial Execution with Progress Bar
    if (threads > 1) {
        results <- pbapply::pblapply(genes, process_one_gene, cl = threads)
    } else {
        results <- pbapply::pblapply(genes, process_one_gene)
    }

    # Combine Results
    results <- Filter(Negate(is.null), results)
    final_df <- rbindlist(results)

    # Write Output
    write.table(final_df,
        file = paste0(output_prefix, ".smr"),
        sep = "\t", quote = FALSE, row.names = FALSE
    )

    message(paste("Analysis complete. Saved to", output_prefix))
    return(final_df)
}
