\name{SMR_analysis}
\alias{SMR_analysis}
\title{Run Complete SMR Analysis Pipeline}
\usage{
SMR_analysis(
  gwas_path,
  eqtl_path,
  ld_reference,
  output_prefix,
  peqtl_smr = 5e-08,
  peqtl_heidi = 0.0157,
  cis_window_kb = 2000,
  threads = 1,
  plot = TRUE
)
}
\arguments{
\item{gwas_path}{Character string. Path to the GWAS summary statistics file.}

\item{eqtl_path}{Character string. Path to the eQTL summary statistics file.}

\item{ld_reference}{Character string. Path prefix to the PLINK binary reference files (e.g., "data/1000G"). Expects .bed, .bim, .fam extensions.}

\item{output_prefix}{Character string. Prefix for output files.}

\item{peqtl_smr}{Numeric. P-value threshold to select top cis-eQTL instruments for SMR test (default 5e-8).}

\item{peqtl_heidi}{Numeric. P-value threshold to select instrumental SNPs for HEIDI test (default 0.0157).}

\item{cis_window_kb}{Numeric. Size of the cis-window in Kb (default 2000).}

\item{threads}{Integer. Number of CPU threads to use for parallel processing (default 1).}

\item{plot}{Logical. Whether to generate plot files for significant results (default TRUE).}
}
\value{
A data.table containing the results of the SMR and HEIDI tests for all analyzed genes.
}
\description{
Runs the full Summary-data-based Mendelian Randomization (SMR) analysis pipeline. 
It loads GWAS and eQTL data, aligns alleles, performs the SMR test to detect associations, 
conducts the HEIDI test to check for heterogeneity (pleiotropy vs linkage), and generates visualization plots.
}
\examples{
\dontrun{
SMR_analysis(
  gwas_path = "gwas.txt",
  eqtl_path = "eqtl.txt",
  ld_reference = "ref/1000G",
  output_prefix = "results",
  threads = 4
)
}
}
