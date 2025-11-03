# Graphical Lasso–guided MICE — Simulation Study on High-Dimensional Missing Data Imputation


## 1. Research Overview

This repository presents a simulation study of the  
**Graphical Lasso–guided Multiple Imputation by Chained Equations (Glasso-MICE)** algorithm,  
designed to improve missing data imputation in **high-dimensional settings**.  

The core idea is to estimate the **conditional dependency structure (precision matrix)**  
using **Graphical Lasso** and incorporate this structure into the  
**MICE predictor matrix** to enhance imputation accuracy.  
Rather than imputing values independently, this approach enables  
**structure-guided imputation** that reflects the underlying network relationships among variables.

---

## 2. Research Objectives

- Improve imputation accuracy under **high-dimensional conditions**  
- Integrate **Graphical Lasso–estimated dependency structures** into the MICE process  
- Compare the effects of different **initial imputations** (mean vs. MICE)  
  and **structural constraints** (excluding diagonal entries in voting)  
- Evaluate performance against **Default MICE**, **Lasso-based MICE**, and **Oracle (complete data)** settings  

---

## 3. Simulation Overview

Each simulation examines how accurately different imputation strategies  
recover the true dependency structure (precision matrix)  
when **10% MCAR (Missing Completely At Random)** values are introduced.  

| Method | Description |
|--------|--------------|
| **Proposed (Glasso-guided MICE)** | Iterative MICE guided by precision matrices estimated via Graphical Lasso |
| **Complete Data (Oracle)** | Benchmark estimation from complete data without missingness |
| **MICE (Default)** | Standard MICE assuming variable independence |
| **MICE (Lasso)** | MICE using lasso-based conditional models (`lasso.select.norm`) |

Each method is repeated across 10 random seeds,  
and evaluated using **FPR (False Positive Rate)**, **TPR (True Positive Rate)**,  
and **AUC (Area Under the ROC Curve)**.

---

## 4. Workflow

| Step | Description |
|------|--------------|
| ① | Load required R packages and set simulation parameters |
| ② | Generate multivariate normal data with a known precision matrix |
| ③ | Introduce 10% MCAR missingness |
| ④ | Perform initial imputation (mean or MICE) |
| ⑤ | Estimate network structure via Graphical Lasso and use it in MICE predictor matrix |
| ⑥ | Compare Proposed, MICE, MICE-Lasso, and Complete Data methods |
| ⑦ | Calculate FPR, TPR, and AUC, and visualize ROC curves |

---

## 5. Repository Structure

| File | Description |
|------|--------------|
| `Graphical Lasso–guided MICE Imputation (p = 20, initial imputation: mean).R` | Baseline experiment using mean initialization |
| `Graphical Lasso–guided MICE Imputation (p = 20, initial imputation: mean, diag(wi)=0 in voting).R` | Version excluding diagonal elements during voting |
| `Graphical Lasso–guided MICE Imputation (p = 20, initial imputation: MICE).R` | Medium-dimensional scenario initialized with MICE |
| `Graphical Lasso–guided MICE Imputation (p = 50, initial imputation: MICE).R` | High-dimensional (𝑝 > 𝑛) extension for scalability testing |

---

## 6. Required Packages

| Package | Purpose |
|----------|----------|
| `MASS` | Generation of multivariate normal data |
| `glasso` | Estimation of sparse precision matrices via Graphical Lasso |
| `mice` | Multiple Imputation by Chained Equations |
| `pracma` | Numerical integration (`trapz`) for AUC computation |
| `ROCR` | ROC curve calculation and performance metrics |
| `ggplot2` | Visualization of ROC curves |
| `gridExtra`, `grid` | Plot arrangement and layout management |

---

## 7. Summary of Results

- **Proposed (Glasso-guided MICE)**  
  Incorporating network structures from Graphical Lasso  
  led to higher AUC values and improved structural recovery compared to standard MICE.  

- **MICE (Default / Lasso)**  
  Provided stable imputations but underperformed in network reconstruction accuracy.  

- **Complete Data (Oracle)**  
  Serves as the upper performance bound for model comparison.  

All results are summarized as averaged ROC curves (across 10 runs),  
AUC comparison tables, and timing summaries for each simulation step.

---

## 8. How to Run

```r
# Install required packages
install.packages(c(
  "MASS", "mice", "glasso", "pracma",
  "ROCR", "ggplot2", "gridExtra", "grid"
))


# Example execution
source("Graphical Lasso–guided MICE Imputation (p = 20, initial imputation: MICE).R")

