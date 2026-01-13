# SMRr: Summary-data-based Mendelian Randomization in R

[中文说明 (Chinese Version)](README.md)

**SMRr** is a high-performance R package designed for integrating GWAS and eQTL summary statistics to identify genes whose expression is pleiotropically associated with complex traits.

This package is an R implementation of the original **SMR** software developed by Prof. Jian Yang's group. The algorithm has been strictly validated to match the precision of the original C++ SMR program ([SMR Official Website](https://yanglab.westlake.edu.cn/software/smr/#Overview)).

## Key Features
- **Precise HEIDI Test**: Uses the exact covariance formula and Satterthwaite approximation (`pchisqsum`) for p-values.
- **Automated Alignment**: Robust allele alignment handling flips, strand complements, and palindromic SNP filtering.
- **Rich Visualizations**:
  - **Locus Plot**: Triple-track visualization of GWAS, eQTL, and gene location.
  - **Effect Size Plot**: Highly customizable beta-beta scattering with $r^2$ gradient coloring.
  - **Manhattan Plot**: Direction-aware (triangle markers) Manhattan plots with chromosome coloring and hit highlighting.
- **High Performance**: Multithreading support and real-time progress bars via `pbapply`.

## Installation

### Method 1: Install from GitHub (Source)
You can install the latest development version directly from GitHub using the `remotes` package:
```R
if (!require("remotes")) install.packages("remotes")
remotes::install_github("CangmingDing/SMRr")
```

### Method 2: Install Pre-compiled Binary (Recommended for Windows/macOS)
Go to the [GitHub Releases](https://github.com/CangmingDing/SMRr/releases) page, download the binary file for your operating system, and install it:

- **Windows**: Download `.zip` file.
- **macOS**: Download `.tgz` file.

```R
# Example for Windows:
install.packages("C:/path/to/SMRr_0.1.0.zip", repos = NULL, type = "win.binary")

# Example for macOS:
install.packages("/path/to/SMRr_0.1.0.tgz", repos = NULL, type = "mac.binary")
```

## Activation

The software requires activation to unlock full functionality. Each machine generates a unique ID.

1. **Get Machine ID**:
   ```R
   library(SMRr)
   get_machine_code()
   ```
2. **Contact Author**: Send your Machine ID to the author to receive an activation code.
   - **Email**: 202201230726@163.com
   - **WeChat**: CangMing-03
3. **Activate**:
   ```R
   activate_smrr("YOUR_ACTIVATION_CODE")
   ```

## Quick Start

```R
library(SMRr)

# 1. Run Analysis
results <- SMR_analysis(
  gwas_path = "gwas_summary.txt",
  eqtl_path = "eqtl_summary.csv",
  ld_reference = "ref_panel_prefix",
  output_prefix = "SMR_out"
)

# 2. Manhattan Plot
plot_manhattan(results)
```

## Authors
- **Cangming** (Main Developer)
- Email: 202201230726@bucm.edu.cn
