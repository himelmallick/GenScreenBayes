# GenScreenBayes  
### Bayesian Screening for High-Dimensional Genomic Signals Using Tweedie Models

<p align="center">
  <a href="https://github.com/himelmallick/GenScreenBayes">
    <img src="inst/sticker/GenScreenBayes.png" width="220" alt="GenScreenBayes Hex Sticker"/>
  </a>
</p>

---
## 🚀 Overview

**GenScreenBayes** implements a **Bayesian genomic screening framework** based on Tweedie exponential dispersion models. It provides a full pipeline for:

- Empirical-Bayes estimation of Tweedie distribution parameters for *control* samples  
- Construction of corresponding *test* distribution parameters under user-specified transformations  
- Efficient Gauss–Hermite quadrature to compute **marginal densities**, **likelihood ratios**, and **posterior quantities**  
- Screening of genomic features through  
  - posterior π₀ distribution  
  - posterior screening probabilities  
  - difference survival analysis  
- Publication-quality summary plots  

The method allows users to compare distributions of gene expression (or other high-dimensional measurements) between control and test regimes and quantify evidence for differential expression in a fully Bayesian manner.

---

## ✨ Key Features

- **Tweedie-based hierarchical Bayesian model**  
- **Empirical Bayes** estimation of control-regime parameters  
- Flexible **test-regime transformations** via `Tmod`  
- **Marginal test densities** via 3D Gauss–Hermite quadrature  
- **Posterior of π₀** (null proportion) with density & CDF  
- **Screening probabilities** (`P(gamma = 0 | data)`)  
- **Difference survival functions** for test–control contrasts  
- Extensive plotting utilities  
- Fully compatible with high-dimensional genomic workflows  

---

## 📥 Installation

### Install from GitHub

```r
# install.packages("devtools")
devtools::install_github("himelmallick/GenScreenBayes")
```

### Install from source

```bash
git clone https://github.com/himelmallick/GenScreenBayes
cd GenScreenBayes
R CMD build .
R CMD INSTALL GenScreenBayes_*.tar.gz
```

---

## 🔧 Example Workflow

A typical analysis pipeline (interactive mode):

```r
library(GenScreenBayes)

# 1. Interactive menu system to build the analysis input
res_input <- GenScreenMenu.fn()

# 2. Run full Bayesian screening engine
res <- GenScreenCalcs.fn(res_input)

# Outputs include:
res$TstMargDens       # marginal test densities
res$SurvDiffs         # survival function of test-control differences
res$PostDenCDFpi0     # posterior density and CDF of pi0
res$Pgam0             # posterior screening probabilities
res$DenCDFplots       # summary posterior plots
res$CombDenCDFplots   # combined density/CDF plots
```

Users may also construct the `GenScrnList` manually without the menu.

### 🔹 Running in Non-Interactive Mode

`GenScreenMenu.fn` also supports a **non-interactive mode**, allowing scripted or batch execution without prompts:

```r
# Skip all interactive prompts
res_input <- GenScreenMenu.fn(interactive = FALSE)

# Then run the engine
res <- GenScreenCalcs.fn(res_input)
```

This allows GenScreenBayes to be used cleanly in automated pipelines, reproducible research scripts, and HPC environments.

---

## 🧠 Method Summary

GenScreenBayes implements a **Tweedie-based hierarchical Bayesian screening model**:

1. **Control regime (`xC`)**
   - Maximum-likelihood estimation of Tweedie parameters
   - Posterior distribution of transformed parameters
   - Empirical Bayes prior for test regime

2. **Test regime (`xT`)**
   - Parameter shifts defined by `Tmod`
   - Calculation of same-process vs different-process marginal densities

3. **Likelihood ratio products**
   - Combined across features using Gauss–Hermite quadrature

4. **Posterior of π₀**
   - Normalized using a Beta(ζ, 1) prior
   - Used to compute posterior screening probabilities

5. **Difference survival distributions**
   - Analytical Tweedie difference survival used to determine effect shifts

---

## 📊 Visualization

GenScreenBayes includes the following plots:

- Posterior density of π₀
- Posterior CDF of π₀
- Combined π₀ density/CDF overlays
- Gene-level screening probabilities
- Difference survival curves

All plots are designed with **ggplot2** themes suitable for publication.


## 📚 Citation

If you use **GenScreenBayes**, please cite:

**Gould, A. L., Paul, E., Basak, P., Bhattacharyya, A., and Mallick, H.** (2025).  
[*High-dimensional Array Bayesian Screening Based on Distributions with Structural Zeroes.*](https://arxiv.org/abs/2503.13605)  
arXiv:2503.13605.

---

## 🤝 Contributing

We are happy to troubleshoot any issues with the package. Please contact the authors via email or open an issue in the GitHub repository.

