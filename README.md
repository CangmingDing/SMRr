# SMRr: 完全基于 R 的SMR分析工具

[English Version](README_EN.md)

**SMRr** 是一款完全基于 R 开发的包，可以实现SMR分析，专门用于整合 GWAS 和 eQTL 总结数据，通过 SMR (Wald test) 和 HEIDI (Heterogeneity in Dependent Instruments) 检验，识别与复杂性状具有多效性关联的基因。

本项目是对杨剑教授团队开发的原始 **SMR** 软件的 R 语言高性能实现。保留了基础SMR分析能力，算法经过严格校对，完全匹配 C++ SMR 原版程序（[SMR Official Website](https://yanglab.westlake.edu.cn/software/smr/#Overview)）

## 核心功能
- **精确的 HEIDI 检验**：采用完整的协方差公式。使用 Satterthwaite 近似法计算 p 值，结果远比简化版 R 脚本精确。
- **自动对齐 (Alignment)**：鲁棒的等位基因对齐逻辑，自动处理转链 (Strand flip)、反向 (Flip) 以及回文 SNP (Palindromic SNP) 及其频率校准。
- **丰富的可视化方案**：
  - **Locus Plot (基因座图)**：三层轨道联展 GWAS, eQTL 及基因物理位置。
  - **Effect Size Plot (效应散点图)**：展示 $\beta_{GWAS}$ vs $\beta_{eQTL}$，支持 $r^2$ 渐变着色。
  - **Manhattan Plot (曼哈顿图)**：支持显示效应方向（正负三角），不同染色体分色，显著基因自动标注。
- **高性能运行**：支持多线程计算，并提供 `pbapply` 实时进度条展示。

## 安装

### 方法 1：从 GitHub 源码安装
您可以使用 `remotes` 包直接从 GitHub 安装最新开发版：
```R
if (!require("remotes")) install.packages("remotes")
remotes::install_github("CangmingDing/SMRr")
```

### 方法 2：安装预编译二进制版本（推荐 Windows/macOS 用户）
请前往 [GitHub Releases](https://github.com/CangmingDing/SMRr/releases) 页面，根据您的操作系统下载对应的安装包：

- **Windows**: 下载 `.zip` 文件。
- **macOS**: 下载 `.tgz` 文件。

```R
# Windows 示例:
install.packages("C:/下载路径/SMRr_0.1.0.zip", repos = NULL, type = "win.binary")

# macOS 示例:
install.packages("/用户路径/下载/SMRr_0.1.0.tgz", repos = NULL, type = "mac.binary")
```

## 软件激活

本软件核心功能需要激活。每台电脑具有唯一的机器码。

1. **获取机器码**：
   ```R
   library(SMRr)
   get_machine_code()
   ```
2. **联系作者获取激活码**：将您的机器码发送给作者。
   - **邮箱**: 202201230726@163.com
   - **微信**: CangMing-03
3. **激活软件**：
   ```R
   # 永久激活（激活后后续使用无需再次输入）
   activate_smrr("您的激活码")
   ```

## 快速上手

```R
library(SMRr)

# 1. 执行分析流程
results <- SMR_analysis(
  gwas_path = "gwas.txt",
  eqtl_path = "eqtl.csv",
  ld_reference = "LDREF/1000G.EUR.22",
  output_prefix = "My_SMR"
)

# 2. 绘制曼哈顿图
plot_manhattan(results)
```

## 作者信息
- **Cangming**
- 邮箱: 202201230726@163.com
- 微信: CangMing-03
