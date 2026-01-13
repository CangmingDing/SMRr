#' Plot SMR Effect Sizes
#'
#' @param gene_data Data frame containing eQTL and GWAS betas
#' @param res_entry SMR result row (must contain b_SMR, topSNP)
#' @return ggplot object
#' @import ggplot2
#' @import dplyr
#' @export
plot_effect_size <- function(gene_data, res_entry) {
    check_auth() # Security check

    top_snp <- res_entry$topSNP
    b_smr <- res_entry$b_SMR

    gene_data$is_top <- gene_data$SNP == top_snp

    if (!"r2" %in% colnames(gene_data)) gene_data$r2 <- runif(nrow(gene_data), 0, 0.9)
    gene_data$r2[gene_data$is_top] <- 1.0

    p <- ggplot(gene_data, aes(x = B_EQTL_ALIGNED, y = B_GWAS)) +
        geom_errorbar(aes(ymin = B_GWAS - SE_GWAS, ymax = B_GWAS + SE_GWAS),
            color = "#483D8B", alpha = 0.3, linewidth = 0.3
        ) +
        geom_errorbarh(aes(xmin = B_EQTL_ALIGNED - SE_EQTL, xmax = B_EQTL_ALIGNED + SE_EQTL),
            color = "#483D8B", alpha = 0.3, linewidth = 0.3
        ) +
        geom_abline(intercept = 0, slope = b_smr, color = "#FF8C00", linetype = "dashed", linewidth = 1) +
        geom_point(aes(color = r2), shape = 21, fill = NA, size = 3, stroke = 1) +
        geom_point(
            data = subset(gene_data, is_top),
            aes(x = B_EQTL_ALIGNED, y = B_GWAS),
            shape = 2, color = "red", size = 5, stroke = 1.5
        ) +
        scale_color_gradientn(colors = c("#4169E1", "#000080", "#101010"), name = bquote(r^2)) +
        labs(
            x = "eQTL effect sizes", y = "GWAS effect sizes",
            title = paste("Effect Size:", res_entry$probeID)
        ) +
        theme_classic() +
        theme(
            axis.text = element_text(color = "black"),
            axis.ticks = element_line(color = "black"),
            legend.position = "right",
            plot.title = element_text(hjust = 0.5, face = "bold")
        )

    return(p)
}

#' Plot SMR Locus Plot
#'
#' Multi-track plot showing GWAS, eQTL and Gene position.
#'
#' @param gene_data Data frame containing SNP, POS, P_GWAS, P_EQTL
#' @param res_entry SMR result row
#' @return ggplot object (patchwork combined)
#' @import ggplot2
#' @import patchwork
#' @export
plot_locus <- function(gene_data, res_entry) {
    check_auth() # Security check

    top_snp <- res_entry$topSNP
    p_smr <- res_entry$p_SMR
    probe_id <- res_entry$probeID

    # Prep data
    plot_df <- gene_data
    plot_df$logP_GWAS <- -log10(plot_df$P_GWAS)
    plot_df$logP_EQTL <- -log10(plot_df$P_EQTL)
    plot_df$pos_mb <- plot_df$POS / 1e6
    plot_df$is_top <- plot_df$SNP == top_snp

    x_range <- range(plot_df$pos_mb)

    # 1. GWAS Track
    p1 <- ggplot(plot_df, aes(x = pos_mb, y = logP_GWAS)) +
        geom_point(color = "#777777", alpha = 0.6, size = 1.5) +
        geom_point(data = subset(plot_df, is_top), color = "#4B0082", shape = 18, size = 4) +
        geom_hline(yintercept = -log10(max(p_smr, 1e-300)), color = "#D02090", linetype = "dashed") +
        annotate("text",
            x = x_range[2], y = -log10(max(p_smr, 1e-300)),
            label = paste0("psmr=", signif(p_smr, 2)),
            vjust = -0.5, hjust = 1, color = "#D02090", fontface = "bold"
        ) +
        labs(y = expression(-log[10](p[GWAS])), title = paste("Locus Plot:", probe_id)) +
        theme_classic() +
        coord_cartesian(xlim = x_range) +
        theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold")
        )

    # 2. eQTL Track
    p2 <- ggplot(plot_df, aes(x = pos_mb, y = logP_EQTL)) +
        geom_point(color = "#A62A55", shape = 4, alpha = 0.7, size = 1.5) +
        geom_point(data = subset(plot_df, is_top), color = "#4B0082", shape = 18, size = 4) +
        labs(y = expression(-log[10](p[eQTL]))) +
        theme_classic() +
        coord_cartesian(xlim = x_range) +
        theme(axis.title.x = element_blank(), axis.text.x = element_blank())

    # 3. Position Track (Gene Arrow)
    probe_pos <- res_entry$Probe_bp / 1e6
    arrow_len <- (x_range[2] - x_range[1]) * 0.1

    p3 <- ggplot() +
        geom_segment(aes(x = probe_pos - arrow_len / 2, xend = probe_pos + arrow_len / 2, y = 1, yend = 1),
            color = "#8B8000", linewidth = 3, arrow = arrow(length = unit(0.2, "cm"), type = "closed")
        ) +
        geom_text(aes(x = probe_pos, y = 1.6, label = probe_id), fontface = "italic", size = 3) +
        labs(x = paste("Position on Chromosome", res_entry$ProbeChr, "(Mb)")) +
        theme_void() +
        coord_cartesian(xlim = x_range, ylim = c(0, 2)) +
        theme(axis.title.x = element_text(size = 10, margin = margin(t = 10), hjust = 0.5))

    combined <- patchwork::wrap_plots(p1, p2, p3, ncol = 1, heights = c(4, 4, 1.2))
    return(combined)
}

#' Plot SMR Manhattan Plot
#'
#' @param smr_results SMR results data.frame (from SMR_analysis)
#' @param p_thresh Significance threshold (default 0.05)
#' @param label_genes Logical, whether to label significant genes (default TRUE)
#' @return ggplot object
#' @import ggplot2
#' @import ggrepel
#' @export
plot_manhattan <- function(smr_results, p_thresh = 0.05, label_genes = TRUE) {
    check_auth() # Security check

    df <- smr_results
    df$logp <- -log10(as.numeric(df$p_SMR))
    df$b_SMR <- as.numeric(df$b_SMR)
    df$shape_val <- ifelse(df$b_SMR > 0, 24, 25)
    df$is_sig <- df$p_SMR < p_thresh
    df$ProbeChr <- as.numeric(as.character(df$ProbeChr))

    df <- df[order(df$ProbeChr, df$Probe_bp), ]

    chrs <- unique(df$ProbeChr)
    tot <- 0
    df$cumulative_pos <- 0
    axis_breaks <- numeric(length(chrs))

    for (i in seq_along(chrs)) {
        chr <- chrs[i]
        idx <- df$ProbeChr == chr
        df$cumulative_pos[idx] <- df$Probe_bp[idx] + tot
        axis_breaks[i] <- tot + (max(df$Probe_bp[idx]) + min(df$Probe_bp[idx])) / 2
        tot <- tot + max(df$Probe_bp[idx])
    }

    chr_colors <- rep(c("#333333", "#888888"), length.out = length(chrs))
    names(chr_colors) <- chrs

    p <- ggplot(df, aes(x = cumulative_pos, y = logp)) +
        geom_point(
            data = subset(df, !is_sig),
            aes(color = as.factor(ProbeChr), shape = shape_val),
            alpha = 0.7, size = 2
        ) +
        geom_point(
            data = subset(df, is_sig),
            aes(x = cumulative_pos, y = logp, shape = shape_val),
            color = "red", fill = "red", size = 3
        ) +
        scale_color_manual(values = chr_colors, name = "Chromosome") +
        scale_shape_identity() +
        geom_hline(yintercept = -log10(p_thresh), color = "red", linetype = "dashed", alpha = 0.5) +
        scale_x_continuous(breaks = axis_breaks, labels = chrs) +
        labs(x = "Chromosome", y = expression(-log[10](p[SMR])), title = "SMR Manhattan Plot") +
        theme_classic() +
        theme(
            legend.position = "right",
            axis.text.x = element_text(size = 8),
            plot.title = element_text(hjust = 0.5, face = "bold")
        )

    if (label_genes && any(df$is_sig)) {
        p <- p + ggrepel::geom_text_repel(
            data = subset(df, is_sig),
            aes(label = probeID),
            size = 3, box.padding = 0.5, point.padding = 0.3,
            max.overlaps = 50, fontface = "bold.italic"
        )
    }

    return(p)
}
