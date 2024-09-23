# SATS (Signature Analyzer for Targeted Sequencing)
<br/>

### Introduction
SATS is a novel method developed for the accurate identification of mutational signatures in tumors sequenced using targeted panels. 
Unlike tools developed for whole-exome or whole-genome sequencing, SATS is specifically designed to address the unique challenges of targeted sequenced tumors. 
It encompasses the detection of de novo signatures, mapping these to reference TMB signatures, estimating signature activities, and calculating signature burdens.

For more information please refer to the [user guide](https://github.com/binzhulab/SATS/blob/main/User_Guide_SATS_v7.pdf).
<br/>

### Installation
To install SATS directly from GitHub:
```r
if (!requireNamespace("devtools", quietly = TRUE))  
	install.packages("devtools")
devtools::install_github("binzhulab/SATS/source")
```
Alternatively, download the package and follow the steps below. Download SATS_0.0.8.tar.gz (for Unix) or SATS_0.0.8.zip (for Windows, R version >= 4.1). To install SATS on Unix, enter the command from a Unix prompt:
```
R CMD INSTALL SATS_0.0.8.tar.gz -l path_to_install_package
```
Alternatively, SATS_0.0.8.tar.gz (for Unix) or SATS_0.0.8.zip (for Windows, R version >= 4.1) from the [Github page](https://github.com/binzhulab/SATS) are available and one may use the following commands:
```
install.packages("./SATS_0.0.8.tar.gz", repos = NULL, type = "source")
install.packages("./SATS_0.0.8.zip", repos = NULL, type = "win.binary")
```
Once the installation is successful, it can be loaded in **R** by calling 
```
library(SATS)
```

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
- We have added an **R** script (*L_matrix_Generation_V2.R*) in [here](https://github.com/binzhulab/SATS/tree/main/Generating_L) containing `Panel_Context_generation()` and `L_matrix_generation()` functions, and panel context informations (*Panel_Info_1_assay.txt* for a single panel information and *Panel_Info_2_assays.txt* for multiple panels) and patient information (*Patient_Info.txt* for the patient IDs associated with each SEQ_ASSAY_ID).
- Since the $\mathbf{L}$ matrix contains the panel context associated with each patient, we first construct referecne panel context using the `Panel_Context_generation()` function. The $\mathbf{L}$ matrix is then generated from the reference panel context by matching SEQ_ASSAY_ID with Patient_ID using the `L_matrix_generation()` function.
- Note also that HG19 reference genome is used to generate panel context matrix (An **R** [package](https://www.bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg19.html) `BSgenome.Hsapiens.UCSC.hg19` will be loaded in *L_matrix_Generation_V2.R* script). For the HG38 reference genome, `BSgenome.Hsapiens.UCSC.hg38` [package](https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg38.html) may be installed and loaded.
- It can be called `Panel_Context_generation(genomic_information, Types)` where `genomic_information` contains `Chromosome`, `Start_Position`, `End_Position`, `SEQ_ASSAY_ID` as belows:
  ```r
  > Panel_1
    Chromosome Start_Position End_Position SEQ_ASSAY_ID Hugo_Symbol
  1          9      133738302    133738491    UHN-48-V1        ABL1
  2          9      133747476    133747664    UHN-48-V1        ABL1
  3          9      133748157    133748327    UHN-48-V1        ABL1
  ...
  ```
  - The column `Chromosome` contains chromsome number where `Start_Position` and `End_Position` columns are start and end positions of targeted panel.
  - The last column `SEQ_ASSAY_ID` distinguishes different panels consisting of the resulting $\mathbf{L}$ matrix (column names in the result).
- **Note**: Please use the column names identical to `Chromosome`, `Start_Position`, `End_Position`, `SEQ_ASSAY_ID` as in the above example (`Hugo_Symbol` is optional and not required to use `Panel_Context_generation()` function).
-  The second arguement on `Panel_Context_generation()` specifies mutation type order as either one of `"COSMIC"` or `"signeR"` where
    - `"COSMIC"` corresponds to the order from the COSMIC database v3.2 and
    - `"signeR"` corresponds to the order from the `signeR` package
    ```r
    > Panel_Context_generation(Panel_1) #not working
    > Panel_Context_generation(Panel_2, Types = "COSMIC")
            GRCC-CP1 UHN-48-V1
    A[C>A]A 0.000883  0.001487
    A[C>A]C 0.000656  0.001120
    A[C>A]G 0.000278  0.000426
    ...
    > Panel_Context_generation(Panel_2, Types = "signeR")
            GRCC-CP1 UHN-48-V1
    C>A:ACA 0.000883  0.001487
    C>A:ACC 0.000656  0.001120
    C>A:ACG 0.000278  0.000426
    ...
    ```
    - The entries of the resulting $\mathbf{L}$ matrix denote the number of trinucleotides per million base pairs.
- **Note**: The extracted $\mathbf{L}$ matrix is used as the opportunity matrix for `signeR()` function and its mutation type order should be the same as input matrix (Mutation catalog matrix $\mathbf{V}$; see Section 3). Thus we highly recommend to confirm that both $\mathbf{V}$ and $\mathbf{L}$ matrices have the same mutation type order corresponding to one of COSMIC database v3.2 or `signeR` package (both have the same order but have different expression) to conduct the consistent analysis.

#### 3. Mapping *de novo* TMB-based Signatures
- Identify *de novo* TMB-based signatures using the signeR algorithm.
  ```r
  library(signeR)
  signeR_re <- signeR(M=V_sum, Opport=L_sum, nlim=c(1,5))
  signeR_re$Phat
  ```
  - We recommend to group samples in $\mathbf{V}$ and $\mathbf{L}$ matrix for computational feasibility and stability when sample size ($N$) is large.
    These pooled matrices, `V_sum` and `L_sum` are used as inputs for `signeR()` function. 
    See [user guide](https://github.com/binzhulab/SATS/blob/main/User_Guide_SATS_v7.pdf) for details.
  - Once `signeR()` is done, the optimal signature profiles are provided in `signeR_re$Phat` which may be used for the next mapping step.
- Map these signatures to reference TMB signatures using `MappingSignature()` function. <!-- using penalized non-negative least squares (pNNLS) -->
  ```r
  W_hat <- signeR_re$Phat
  MappedSig <- MappingSignature(W_hat = W_hat, W_ref = RefTMB$SBS_W)
  MappedSig
  ```
  - `W_hat` is a *de novo* TMB signatures from signeR (`signeR_re$Phat`) or any other signature analysis tool.
  - `W_ref` is for the reference TMB signature profiles which will be mapped to. In this example code, we used the COSIMC reference TMB signature profiles (stored in `RefTMB$SBS_W`).
  - `MappedSig` contains mapped reference TMB signatures, e.g., COSMIC SBS1, SBS2/13, SBS4, SBS5, SBS40, and SBS89 signatures (`MappedSig$Reference`), with frequencies (`MappedSig$freq`) of coefficients greater than 0.1 out of 100 cross-validated repetitions. 
  
#### 4. Estimating Signature Activities and Burdens
- Utilize the expectation-maximization algorithm to estimate signature activities by running `EstimateSigActivity()` function.
  ```r
  W_star <- as.matrix(RefTMB$SBS_W[,SBS.list])
  H_hat <- EstimateSigActivity(V = SimData$V, L = SimData$L, W = W_star)
  H_hat$H
  ```
  - `V` is the mutation type matrix $\mathbf{V}$, `L` is the panel context matrix $\mathbf{L}$ and `W` is the mapped reference TMB signatures $\mathbf{W}^\*$.
  - The resulting `H_hat$H` is the estimated activity matrix of size $K \times N$, where $K$ is the number of signatures given in `W`.
- Calculate the expected number of mutations attributed to a signature with `CalculateSigExpectancy()` function.
  ```r
  SigBdn <- CalculateSigExpectancy(L = SimData$L, W = W_star, H = H_hat$H)
  ```
  - `L` is the panel context matrix $\mathbf{L}$, `W` is the mapped reference TMB signatures $\mathbf{W}^\*$ and `H` is the the estimated activity matrix $\widehat{\mathbf{H}}$.
  - The resulting matrix contains the number of mutations due to the identified signatures.
    For example, Sample 1 has 14 mutations, about 4.8 caused by signature 1, 7 by signature 4, about 0.2 by signature 5 and two by signature 40 as described below.
  ```r
  > round(SigBdn[, 1:5], 2)
          Sample1 Sample2 Sample3 Sample4 Sample5
  SBS1       4.79    0.00       0    0.00    0.00
  SBS2_13    0.00    0.00       1    1.08    0.82
  SBS4       7.04    0.00       0    0.00    0.77
  SBS5       0.21    0.00       0    2.93    0.07
  SBS40      1.95    3.87       0    0.99    3.34
  SBS89      0.00    2.13       0    0.00    0.00
  ```

#### 5. Signature Refitting for Single Tumors
- An additional benefit of the SATS algorithm is its capability to estimate signature activities and burdens, even when working with a single tumor sample, as long as the set of signatures specific to a particular cancer type is known.
- In cases with limited sample sizes, substitute the mapped reference TMB signatures with those provided for the specific cancer type of the tumor sample.
- `RefTMB$SBS_W` and `RefTMB$SBS_refSigs` contain the COSMIC SBS TMB signature profiles and the list of cancer specific SBS signature names, respectively.
  Similarly, `RefTMB$DBS_W` and `RefTMB$DBS_refSigs` contain the COSMIC DBS TMB signature profiles and the list of cancer specific DBS signature names.
- As an example, a single simulated tumor derived from a skin Cancer stored in `SimData$SingleTumorEx` (the `singleV` and `singleV` contain mutation counts and sequencing context respectively).
  ```r
  SBS.list <- RefTMB$SBS_refSigs[RefTMB$SBS_refSigs$cancerType == "Skin Cancer or Melanoma", "COSMIC"]
  W_star <- as.matrix(SimData$W_TMB[,SBS.list])
  ## Estimate activity
  V1 <- SimData$SingleTumorEx[, 1, drop = FALSE]
  L1 <- SimData$SingleTumorEx[, 2, drop = FALSE]
  H_hat <- EstimateSigActivity(V = V1, L = L1, W = W_star)
  ## Estimate burden
  SigBdn <- CalculateSigExpectancy(L = L1, W = W_star, H = H_hat$H)
  ```
  - Similar to the step 3, functions `EstimateSigActivity()` and `CalculateSigExpectancy()` work to estimate signature activities and signature burdens.

### Conclusion
SATS provides a comprehensive approach for analyzing mutational signatures in targeted sequenced tumors, addressing the limitations of existing tools and providing detailed steps for analysis in various scenarios. This work is under the license of CC BY-NC 4.0. 

