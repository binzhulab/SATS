# SATS (Signature Analyzer for Targeted Sequencing)
<br/>

### Introduction
SATS is a panel-aware framework for mutational signature analysis in targeted sequencing data. Unlike tools developed primarily for whole-exome or whole-genome sequencing, SATS models panel-specific sequence context and mutation opportunity, enabling de novo signature extraction, mapping to tumor mutational burden (TMB)-normalized reference signatures, individual-tumor signature refitting and calculation of signature-attributed mutation burdens.

The accompanying manuscript applies SATS to 111,711 tumors from AACR Project GENIE to construct a real-world, panel-calibrated pan-cancer catalogue of targeted sequencing-derived mutational signatures. The package and repository are intended to support analysis of targeted-panel cohorts and use of the catalogue in settings where WES/WGS data are unavailable.

For more information please refer to the [user guide](https://github.com/binzhulab/SATS/blob/main/User_Guide_SATS_v1.0.8.md) or the corresponding [PDF](https://github.com/binzhulab/SATS/blob/main/User_Guide_SATS_v1.0.8.pdf). A redesigned static project webpage is available in [`docs/`](https://github.com/binzhulab/SATS/tree/main/docs) and can be served through GitHub Pages.
<br/>

### Installation
The current source version in this repository is **SATS v1.0.8**. To install SATS directly from GitHub:
```r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("binzhulab/SATS", subdir = "source", upgrade = "never")
```

Alternatively, download `SATS_1.0.8.tar.gz` and install the source archive:
```bash
R CMD INSTALL SATS_1.0.8.tar.gz
```
The source archive can also be installed from within R:
```r
install.packages("./SATS_1.0.8.tar.gz", repos = NULL, type = "source")
```

SATS was formerly available from [CRAN](https://CRAN.R-project.org/package=SATS). CRAN currently lists the package as archived, so the GitHub source installation above is the recommended installation route for the current version.

Once the installation is successful, it can be loaded in **R** by calling 
```
library(SATS)
```

### Repository layout and software quality controls
The repository contains the current R package source under [`source/`](https://github.com/binzhulab/SATS/tree/main/source), a source archive for version 1.0.8, older package archives under [`old_versions/`](https://github.com/binzhulab/SATS/tree/main/old_versions), an example panel-context generation workflow under [`Generating_L/`](https://github.com/binzhulab/SATS/tree/main/Generating_L), and a small Nextflow example under [`nextflow/example1/`](https://github.com/binzhulab/SATS/tree/main/nextflow/example1).

Unit tests are provided in `source/tests/testthat/` for the three main user-facing functions: `CalculateSignatureBurdens()`, `EstimateSigActivity()` and `GeneratePanelSize()`. Each test loads simulated data from `SimData` and compares the returned object with a stored expected result. To run the tests locally:
```bash
cd source
Rscript -e 'testthat::test_dir("tests/testthat")'
```

The repository also includes a GitHub Actions workflow, `.github/workflows/R-CMD-check.yaml`, that runs `R CMD check` on Linux, macOS and Windows. For a local package check:
```bash
R CMD check source --no-manual
```
The check requires the R dependencies listed in `source/DESCRIPTION`, including Bioconductor packages used for genomic interval and reference-genome handling.

The example Nextflow workflow in `nextflow/example1/` calls `GeneratePanelSize()` on an input `.rda` file containing `genomic_information` and optional `Class`, `SBS_order` and `ref.genome` objects. This example is intended as a minimal template for integrating SATS panel-context generation into standardized workflow systems rather than as a complete production pipeline.

### Citation and web resources
If you use SATS or the targeted-sequencing mutational-signature catalogue, please cite the associated manuscript:

Lee et al., "A real-world pan-cancer catalogue of mutational signatures from 111,711 tumors" (submitted).

Interactive catalogue plots displaying the frequencies and etiologies of single base substitution (SBS) and double base substitution (DBS) signatures are available at https://analysistools.cancer.gov/mutational-signatures/#/catalog/STS. An online signature refitting tool is available at https://analysistools.cancer.gov/mutational-signatures/#/refitting.

### Input limitations
SATS currently expects users to provide summarized mutation-count matrices, panel-context matrices and panel-coordinate data frames with the required columns. The package does not directly ingest raw VCF or BED files. VCF/BED-derived data should first be converted into mutation catalogue and panel-coordinate tables before running SATS.

### A schematic workflow of SATS
<img width="1134" alt="image" src="https://github.com/binzhulab/SATS/assets/51965629/64b226ef-58c1-4fc5-aca1-2be4c4a7cf6b">

**a**. The workflow starts with summarizing somatic mutations identified through targeted sequencing, including single base substitutions (SBS), into a mutation type matrix $\mathbf{V}$. 
In addition, SATS requires a panel context matrix $\mathbf{L}$ that specifies the number of trinucleotide contexts for individual panels. 
SATS is based on a Poisson Nonnegative-Matrix Factorization (pNMF) model, approximating $\mathbf{V}$ by $\mathbf{L} \circ \mathbf{W} \times \mathbf{H}$ 
(i.e., $\mathbf{V} \approx \mathbf{L} \circ \mathbf{W} \times \mathbf{H}$, where $\circ$ denotes the element-wise product and $\times$ represents the matrix multiplication operator. <br/>

**b**. The analysis procedure of SATS involves signature detection for a patient cohort and signature refitting for individual patients. 
In this illustrative example, SATS initially identifies de novo tumor mutation burden (TMB) signature 1 and 2 for a patient cohort, and then maps them to reference TMB signatures 1, 2/13 and 5. 
Subsequently, SATS carries out signature refitting for 6 patients (e.g., Pt.1, Pt.2, …, Pt.6), estimating activities of the mapped reference TMB signatures and the expected number of mutations attributed to each signature, namely signature burden.  
For instance, the activities of SBS1, SBS2/13 and SBS5 for patient 3 (Pt.3) are 0.27, 0.84 and 0.18. 
Additionally, we estimate 0.67, 1.16 and 3.17 SBS attributed to signature SBS1, SBS2/13 and SBS5, respectively.

### Example Data
The package includes a simulated dataset: <br/>
- A 96 × 10027 mutation catalog matrix $\mathbf{V}$, representing 10027 targeted sequenced tumors across 96 single base substitution (SBS) types with Panel size matrix $\mathbf{L}$.
- These matrices are stored in `SimData` with corresponding names: $\mathbf{V}$ (`SimData$V`), $\mathbf{L}$ (`SimData$L`) as follows. <br/>
```r
data(SimData, package = "SATS")
SimData$V[1:6, 1:6]
SimData$L[1:6, 1:6]
```

### SATS Quick Usage Guide
#### 1. Main Input matrices
- Mutation catalog matrix $\mathbf{V}$: An **R** dataframe or matrix of size $P \times N$ with non-negative counts, columns represent tumors and rows represent mutation types.
- Panel size matrix $\mathbf{L}$: An **R** dataframe or matrix of size $P \times N$ representing the length of trinucleotide contexts per million base pairs for the corresponding sequencing panel.
- Reference TMB signatures $\mathbf{W}_0$: A predefined reference TMB signatures for refitting stage.

#### 2. Generating the panel size matrix $\mathbf{L}$
- The package provides `GeneratePanelSize()` and `GenerateLMatrix()` to construct panel-context and sample-level panel matrices. The [`Generating_L/`](https://github.com/binzhulab/SATS/tree/main/Generating_L) directory also contains legacy helper scripts and example files for generating a panel context matrix from one or more assays.
- Since the $\mathbf{L}$ matrix contains the panel context associated with each patient, first construct panel-level mutation-context opportunity counts using `GeneratePanelSize()`. The $\mathbf{L}$ matrix is then generated from those panel-level counts by matching `SEQ_ASSAY_ID` with `PATIENT_ID` using `GenerateLMatrix()`.
- `GeneratePanelSize()` can use either the HG19 or HG38 reference genome through `ref.genome = "hg19"` or `ref.genome = "hg38"`. The corresponding Bioconductor package, `BSgenome.Hsapiens.UCSC.hg19` or `BSgenome.Hsapiens.UCSC.hg38`, must be installed.
- `GeneratePanelSize(genomic_information, Class = "SBS", SBS_order = "COSMIC", ref.genome = "hg19")` expects `genomic_information` to contain `Chromosome`, `Start_Position`, `End_Position` and `SEQ_ASSAY_ID`, as below:
  ```r
  > Panel_1
    Chromosome Start_Position End_Position SEQ_ASSAY_ID Hugo_Symbol
  1          9      133738302    133738491    UHN-48-V1        ABL1
  2          9      133747476    133747664    UHN-48-V1        ABL1
  3          9      133748157    133748327    UHN-48-V1        ABL1
  ...
  ```
  - The column `Chromosome` contains the chromosome number, and `Start_Position` and `End_Position` are the start and end positions of the targeted panel.
  - The column `SEQ_ASSAY_ID` distinguishes different sequencing panels in the resulting panel-context matrix.
- **Note**: Please use the column names identical to `Chromosome`, `Start_Position`, `End_Position`, `SEQ_ASSAY_ID` as in the above example (`Hugo_Symbol` is optional and not required by `GeneratePanelSize()`).
-  The `SBS_order` argument of `GeneratePanelSize()` specifies mutation type order as either one of `"COSMIC"` or `"signeR"` where
    - `"COSMIC"` corresponds to the 96 SBS mutation-type order used by COSMIC v3.2 and
    - `"signeR"` corresponds to the order from the `signeR` package
    ```r
    > GeneratePanelSize(Panel_2, Class = "SBS", SBS_order = "COSMIC", ref.genome = "hg19")
            GRCC-CP1 UHN-48-V1
    A[C>A]A 0.000883  0.001487
    A[C>A]C 0.000656  0.001120
    A[C>A]G 0.000278  0.000426
    ...
    > GeneratePanelSize(Panel_2, Class = "SBS", SBS_order = "signeR", ref.genome = "hg19")
            GRCC-CP1 UHN-48-V1
    C>A:ACA 0.000883  0.001487
    C>A:ACC 0.000656  0.001120
    C>A:ACG 0.000278  0.000426
    ...
    ```
    - The entries of the resulting panel context matrix denote the number of trinucleotides per million base pairs.
- To create the $\mathbf{L}$ matrix, call `GenerateLMatrix(Panel_context, Patient_Info)` where `Panel_context` is the panel context matrix generated by `GeneratePanelSize()` above and `Patient_Info` contains patient IDs associated with `SEQ_ASSAY_ID`.
  ```r
  > Panel_context <- GeneratePanelSize(Panel_2, Class = "SBS", SBS_order = "signeR", ref.genome = "hg19")
  > Patient_Info[c(1:2, 11:12), ]
     PATIENT_ID SEQ_ASSAY_ID
  1       UHN_1    UHN-48-V1
  2       UHN_2    UHN-48-V1
  11     GRCC_1     GRCC-CP1
  12     GRCC_2     GRCC-CP1
  > L_mat <- GenerateLMatrix(Panel_context, Patient_Info)
  > L_mat[1:3, 1:3]
             UHN_1    UHN_2    UHN_3
  C>A:ACA 0.001487 0.001487 0.001487
  C>A:ACC 0.001120 0.001120 0.001120
  C>A:ACG 0.000426 0.000426 0.000426
  ```
- **Note**: The extracted $\mathbf{L}$ matrix is used as the opportunity matrix for `signeR()` and as the panel-context matrix for SATS refitting. Its mutation type order must match the input mutation catalogue matrix $\mathbf{V}$ (see Section 3). The `"COSMIC"` option here refers to SBS mutation-type ordering, whereas the COSMIC reference-signature version used for mapping is controlled separately by `MappingSignature(COSMICv=...)`.

#### 3. Mapping *de novo* TMB-based Signatures
- Identify *de novo* TMB-based signatures using the signeR algorithm.
  ```r
  library(signeR)
  signeR_re <- signeR(M=V_sum, Opport=L_sum, nlim=c(1,5))
  signeR_re$Phat
  ```
  - We recommend to group samples in $\mathbf{V}$ and $\mathbf{L}$ matrix for computational feasibility and stability when sample size ($N$) is large.
    These pooled matrices, `V_sum` and `L_sum` are used as inputs for `signeR()` function. 
    See the [user guide](https://github.com/binzhulab/SATS/blob/main/User_Guide_SATS_v1.0.8.md) for details.
  - Once `signeR()` is done, the optimal signature profiles are provided in `signeR_re$Phat` which may be used for the next mapping step.
- Map these signatures to reference TMB signatures using `MappingSignature()` function. <!-- using penalized non-negative least squares (pNNLS) -->
  ```r
  data(RefTMB, package = "SATS")
  W_hat <- signeR_re$Phat
  MappedSig <- MappingSignature(W_hat = W_hat, W_ref = RefTMB$TMB_SBS_v3.4)
  MappedSig
  ```
  - `W_hat` is a *de novo* TMB signatures from signeR (`signeR_re$Phat`) or any other signature analysis tool.
  - `W_ref` is the matrix of reference TMB signature profiles to which the de novo signatures will be mapped. In this example, we use the COSMIC v3.4 SBS TMB reference profiles stored in `RefTMB$TMB_SBS_v3.4`. If `W_ref` is not supplied, `MappingSignature()` defaults to `COSMICv = "v3.4"` and uses `RefTMB$TMB_SBS_v3.4`.
  - `MappedSig` contains mapped reference TMB signatures, such as COSMIC SBS1, SBS2/13, SBS4, SBS5 and other signatures present in the selected reference matrix (`MappedSig$Reference`), with frequencies (`MappedSig$freq`) of coefficients greater than 0.1 out of 100 cross-validated repetitions. 
  
#### 4. Estimating Signature Activities and Burdens
- Utilize the expectation-maximization algorithm to estimate signature activities by running `EstimateSigActivity()` function.
  ```r
  data(RefTMB, package = "SATS")
  SBS.list <- c("SBS1", "SBS2_13", "SBS4", "SBS5", "SBS6", "SBS89")
  W_star <- as.matrix(RefTMB$TMB_SBS_v3.4[, SBS.list])
  H_hat <- EstimateSigActivity(V = SimData$V, L = SimData$L, W = W_star)
  H_hat$H
  ```
  - `V` is the mutation type matrix $\mathbf{V}$, `L` is the panel context matrix $\mathbf{L}$ and `W` is the mapped reference TMB signatures $\mathbf{W}^\*$.
  - The resulting `H_hat$H` is the estimated activity matrix of size $K \times N$, where $K$ is the number of signatures given in `W`.
- Calculate the expected number of mutations attributed to a signature with `CalculateSignatureBurdens()` function.
  ```r
  SigBdn <- CalculateSignatureBurdens(L = SimData$L, W = W_star, H = H_hat$H)
  ```
  - `L` is the panel context matrix $\mathbf{L}$, `W` is the mapped reference TMB signatures $\mathbf{W}^\*$ and `H` is the the estimated activity matrix $\widehat{\mathbf{H}}$.
  - The resulting matrix contains the expected number of mutations attributed to each selected reference signature in each sample.
  ```r
  round(SigBdn[, 1:5], 2)
  ```

#### 5. Signature Refitting for Single Tumors
- An additional benefit of the SATS algorithm is its capability to estimate signature activities and burdens, even when working with a single tumor sample, as long as the set of signatures specific to a particular cancer type is known.
- In cases with limited sample sizes, substitute the mapped reference TMB signatures with those provided for the specific cancer type of the tumor sample.
- `RefTMB$TMB_SBS_v3.4` and `RefTMB$SBS_refSigs` contain the COSMIC v3.4 SBS TMB signature profiles and the list of cancer-specific SBS signature names, respectively.
  Similarly, `RefTMB$TMB_DBS_v3.4` and `RefTMB$DBS_refSigs` contain the COSMIC v3.4 DBS TMB signature profiles and the list of cancer-specific DBS signature names.
- As an example, a single simulated tumor derived from a skin Cancer stored in `SimData$SingleTumorEx` (the `singleV` and `singleV` contain mutation counts and sequencing context respectively).
  ```r
  data(RefTMB, package = "SATS")
  SBS.list <- RefTMB$SBS_refSigs[RefTMB$SBS_refSigs$cancerType == "Skin Cancer or Melanoma", "COSMIC"]
  W_star <- as.matrix(RefTMB$TMB_SBS_v3.4[, SBS.list])
  ## Estimate activity
  V1 <- SimData$SingleTumorEx[, 1, drop = FALSE]
  L1 <- SimData$SingleTumorEx[, 2, drop = FALSE]
  H_hat <- EstimateSigActivity(V = V1, L = L1, W = W_star)
  ## Estimate burden
  SigBdn <- CalculateSignatureBurdens(L = L1, W = W_star, H = H_hat$H)
  ```
  - Similar to step 3, `EstimateSigActivity()` and `CalculateSignatureBurdens()` estimate signature activities and signature burdens.

### Conclusion
SATS provides a comprehensive approach for analyzing mutational signatures in targeted sequenced tumors, addressing the limitations of existing tools and providing detailed steps for analysis in various scenarios. See `source/DESCRIPTION` and the repository license file for licensing information.
