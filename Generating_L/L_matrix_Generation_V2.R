## This file contains two R functions, Panel_Context_generation() defined in line 22 and L_matrix_generation() defined in line 123
## In order to produce proper L matrix, first construct Panel context matrix and then match patient ID with this panel context
## Examples are provided from line 135

### Required package
library(BSgenome.Hsapiens.UCSC.hg19) ## This package is for the HG19 reference genome
## library(BSgenome.Hsapiens.UCSC.hg38) ## This package is for the HG38 reference genome
library(tidyverse)

## Example data and results
path <- "path\to\data"
Panel_1 <- read.table(file.path(path, "Panel_Info_1_assay.txt"), header = T, quote = "", sep="\t", stringsAsFactors = FALSE)
Panel_2 <- read.table(file.path(path, "Panel_Info_2_assays.txt"), header = T, quote = "", sep="\t", stringsAsFactors = FALSE)
Patient_Info <- read.table(file.path(path, "Patient_Info.txt"), header = T, quote = "", sep="\t", stringsAsFactors = FALSE)

## The following Panel_Context_generation() function can be used to create panel context matrix
## Please use the column names of the argument "genomic_information" identical to Chromosome, Start_Position, End_Position, SEQ_ASSAY_ID as in the above example
## Chromosome: chromsome number
## Start_Position: start position of targeted panel
## End_Position: end position of targeted panel
## SEQ_ASSAY_ID: distinguish different panels
## The unit of the returned L matrix is the number of trinucleotides per million base pairs
Panel_Context_generation <- function(genomic_information, Types = c("COSMIC", "signeR")){
  
  if(!all(c("Chromosome", "Start_Position", "End_Position", "SEQ_ASSAY_ID") %in% colnames(genomic_information))){
    stop("Please specify the column names of the argument genomic_information in (Chromosome, Start_Position, End_Position, SEQ_ASSAY_ID)")
  } 
  if(length(Types) == 2){
    stop("Please specifiy mutation types: avialable mutatioins types are either \"COSMIC\" or \"signeR\"")
  }
  
  # mutation type categories
  # Alexandrov v3.2 mutation type order
  COSMIC <- c("A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T", "G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T", 
              "T[C>A]A", "T[C>A]C", "T[C>A]G", "T[C>A]T", "A[C>G]A", "A[C>G]C", "A[C>G]G", "A[C>G]T", "C[C>G]A", "C[C>G]C", "C[C>G]G", "C[C>G]T", 
              "G[C>G]A", "G[C>G]C", "G[C>G]G", "G[C>G]T", "T[C>G]A", "T[C>G]C", "T[C>G]G", "T[C>G]T", "A[C>T]A", "A[C>T]C", "A[C>T]G", "A[C>T]T", 
              "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T", "G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T", "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T", 
              "A[T>A]A", "A[T>A]C", "A[T>A]G", "A[T>A]T", "C[T>A]A", "C[T>A]C", "C[T>A]G", "C[T>A]T", "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T", 
              "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T", "A[T>C]A", "A[T>C]C", "A[T>C]G", "A[T>C]T", "C[T>C]A", "C[T>C]C", "C[T>C]G", "C[T>C]T", 
              "G[T>C]A", "G[T>C]C", "G[T>C]G", "G[T>C]T", "T[T>C]A", "T[T>C]C", "T[T>C]G", "T[T>C]T", "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T", 
              "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T", "G[T>G]A", "G[T>G]C", "G[T>G]G", "G[T>G]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T")

  signeR <- c("C>A:ACA", "C>A:ACC", "C>A:ACG", "C>A:ACT", "C>A:CCA", "C>A:CCC", "C>A:CCG", "C>A:CCT", "C>A:GCA", "C>A:GCC", "C>A:GCG", "C>A:GCT", 
              "C>A:TCA", "C>A:TCC", "C>A:TCG", "C>A:TCT", "C>G:ACA", "C>G:ACC", "C>G:ACG", "C>G:ACT", "C>G:CCA", "C>G:CCC", "C>G:CCG", "C>G:CCT", 
              "C>G:GCA", "C>G:GCC", "C>G:GCG", "C>G:GCT", "C>G:TCA", "C>G:TCC", "C>G:TCG", "C>G:TCT", "C>T:ACA", "C>T:ACC", "C>T:ACG", "C>T:ACT", 
              "C>T:CCA", "C>T:CCC", "C>T:CCG", "C>T:CCT", "C>T:GCA", "C>T:GCC", "C>T:GCG", "C>T:GCT", "C>T:TCA", "C>T:TCC", "C>T:TCG", "C>T:TCT", 
              "T>A:ATA", "T>A:ATC", "T>A:ATG", "T>A:ATT", "T>A:CTA", "T>A:CTC", "T>A:CTG", "T>A:CTT", "T>A:GTA", "T>A:GTC", "T>A:GTG", "T>A:GTT", 
              "T>A:TTA", "T>A:TTC", "T>A:TTG", "T>A:TTT", "T>C:ATA", "T>C:ATC", "T>C:ATG", "T>C:ATT", "T>C:CTA", "T>C:CTC", "T>C:CTG", "T>C:CTT", 
              "T>C:GTA", "T>C:GTC", "T>C:GTG", "T>C:GTT", "T>C:TTA", "T>C:TTC", "T>C:TTG", "T>C:TTT", "T>G:ATA", "T>G:ATC", "T>G:ATG", "T>G:ATT", 
              "T>G:CTA", "T>G:CTC", "T>G:CTG", "T>G:CTT", "T>G:GTA", "T>G:GTC", "T>G:GTG", "T>G:GTT", "T>G:TTA", "T>G:TTC", "T>G:TTG", "T>G:TTT")
  
  # To handle L matrix
  # define trinucleotide seq from 5' to 3'. They are the same order as the function trinucleotideFrequency()
  nt = c("A","C","G","T")
  tri_nt = paste(rep(nt,each=16,times=1),rep(nt,each=4,times=4),rep(nt,each=1,times=16), sep="")
  tri_nt_idx = setNames(1:64, tri_nt)
  
  #index for where is A, C, G, T in the middle of trinucleotide
  A_idx = c(1:4,17:20,33:36,49:52)
  C_idx = A_idx + 4
  G_idx = A_idx + 4*2
  T_idx = A_idx + 4*3
  
  tri_nt_CT = tri_nt[c(C_idx,T_idx)]
  tri_nt_AG = tri_nt[c(A_idx,G_idx)]
  
  tri_nt_CT_idx = setNames(1:32,tri_nt_CT)
  tri_nt_AG_idx = setNames(1:32,tri_nt_AG)
  
  #The complementary trinucleotide. Note: flip of 1st and 3rd position be consistent with from 5' to 3'
  tri_nt_comp = paste(rep(rev(nt),each=1,times=16),rep(rev(nt),each=4,times=4),rep(rev(nt),each=16,times=1), sep="")
  
  #This index could help extract columns A or G base mutations, and combine with C or T base mutations
  AGtoCT_idx = tri_nt_CT_idx[tri_nt_comp[c(A_idx,G_idx)]]
  
  
  Seq_assay_GRanges <- GRanges(seqnames = paste0("chr",genomic_information$Chromosome),
                               IRanges(start = genomic_information$Start_Position-1, end=genomic_information$End_Position+1), strand = "+")
  
  Seq_assay_n = length(Seq_assay_GRanges)
  
  Seq_assay_64 = matrix(0,Seq_assay_n,64)
  colnames(Seq_assay_64) = tri_nt
  
  myseq = getSeq(Hsapiens,Seq_assay_GRanges) # get sequences
  Seq_assay_64 = trinucleotideFrequency(myseq) # key step: calculate trinucleotide frequency
  
  Seq_assay_64_CT = Seq_assay_64[,c(C_idx,T_idx)]
  Seq_assay_64_AG = Seq_assay_64[,c(A_idx,G_idx)]
  
  Seq_assay_CT_32 = Seq_assay_64_CT+Seq_assay_64_AG[,AGtoCT_idx]
  
  genomic_information_CT_32 = bind_cols(genomic_information, as.data.frame(Seq_assay_CT_32))
  
  genomic_information_CT_32 = as_tibble(genomic_information_CT_32)
  
  #sum over genes for the whole assay
  assay_CT_32 = genomic_information_CT_32 %>% 
    group_by(SEQ_ASSAY_ID) %>% 
    summarise_at(tri_nt_CT,sum)
  assay_CT_32$assaySize = rowSums(assay_CT_32[,-1])
  
  if(Types == "COSMIC"){
    assay_CT_32_2 <- data.frame(t(assay_CT_32[, paste0(substr(COSMIC, 1, 1), substr(COSMIC, 3, 3), substr(COSMIC, 7, 7)), drop = FALSE]))/10^6
    rownames(assay_CT_32_2) <- COSMIC
  } else if(Types == "signeR"){
    assay_CT_32_2 <- data.frame(t(assay_CT_32[, substr(signeR, 5, 7), drop = FALSE]))/10^6
    rownames(assay_CT_32_2) <- signeR  
  } 
  
  colnames(assay_CT_32_2) <- assay_CT_32$SEQ_ASSAY_ID
  
  return(assay_CT_32_2)
  
}

##################################################################################################################################################################
## The following L_matrix_generation() function can be used to generate "L" matrix.
## Please use two inputs, the first argument is the panel context matrix generated by Panel_Context_generation() function, 
## while the second argument is the Patient information with sequence assay ID (Please refer to the example included in the folder)
## For the second argument, please use the column names as below:
## PATIENT_ID: patient ID corresponds to SEQ_ASSAY_ID
## SEQ_ASSAY_ID: SEQ_ASSAY_ID contained in Panel_context
L_matrix_generation <- function(Panel_context, Patient_Info){
  
  idx <- Patient_Info$SEQ_ASSAY_ID %in% colnames(Panel_context)
  L <- Panel_context[, Patient_Info$SEQ_ASSAY_ID[idx]]
  colnames(L) <- Patient_Info$PATIENT_ID[idx]
  if(sum(idx) != nrow(Patient_Info)){
    warning(sprintf("There are patients for whom the panel context has not been provided in Patient_Info. 
    Patient_Info contains %d sequence assays, but Panel_context has only %d sequence assays", length(unique(Patient_Info$SEQ_ASSAY_ID)), ncol(Panel_context)))
  } 
  
  return(L)
}

##################################################################################################################################################################
## Examples:
Panel_Context_generation(Panel_1) #not working
Panel_Context_generation(Panel_1, Types = "COSMIC")
## Part of the results:
#         UHN-48-V1
# A[C>A]A  0.001487
# A[C>A]C  0.001120
# A[C>A]G  0.000426
# A[C>A]T  0.001155
# C[C>A]A  0.001654
# C[C>A]C  0.001300
# C[C>A]G  0.000543
# C[C>A]T  0.001464

Panel_Context_generation(Panel_1, Types = "signeR")
#         UHN-48-V1
# C>A:ACA  0.001487
# C>A:ACC  0.001120
# C>A:ACG  0.000426
# C>A:ACT  0.001155
# C>A:CCA  0.001654
# C>A:CCC  0.001300
# C>A:CCG  0.000543
# C>A:CCT  0.001464

###########################################################
Panel_Context_generation(Panel_2, Types = "COSMIC")
Panel_Context_generation(Panel_2, Types = "signeR")

###########################################################
Panel_context <- Panel_Context_generation(Panel_2, Types = "signeR")
L_mat <- L_matrix_generation(Panel_context, Patient_Info)
head(L_mat)
# UHN_1    UHN_2    UHN_3    UHN_4    UHN_5    UHN_6    UHN_7    UHN_8    UHN_9   UHN_10   GRCC_1   GRCC_2   GRCC_3   GRCC_4   GRCC_5   GRCC_6   GRCC_7   GRCC_8   GRCC_9  GRCC_10
# C>A:ACA 0.001487 0.001487 0.001487 0.001487 0.001487 0.001487 0.001487 0.001487 0.001487 0.001487 0.000883 0.000883 0.000883 0.000883 0.000883 0.000883 0.000883 0.000883 0.000883 0.000883
# C>A:ACC 0.001120 0.001120 0.001120 0.001120 0.001120 0.001120 0.001120 0.001120 0.001120 0.001120 0.000656 0.000656 0.000656 0.000656 0.000656 0.000656 0.000656 0.000656 0.000656 0.000656
# C>A:ACG 0.000426 0.000426 0.000426 0.000426 0.000426 0.000426 0.000426 0.000426 0.000426 0.000426 0.000278 0.000278 0.000278 0.000278 0.000278 0.000278 0.000278 0.000278 0.000278 0.000278
# C>A:ACT 0.001155 0.001155 0.001155 0.001155 0.001155 0.001155 0.001155 0.001155 0.001155 0.001155 0.000687 0.000687 0.000687 0.000687 0.000687 0.000687 0.000687 0.000687 0.000687 0.000687
# C>A:CCA 0.001654 0.001654 0.001654 0.001654 0.001654 0.001654 0.001654 0.001654 0.001654 0.001654 0.001022 0.001022 0.001022 0.001022 0.001022 0.001022 0.001022 0.001022 0.001022 0.001022
# C>A:CCC 0.001300 0.001300 0.001300 0.001300 0.001300 0.001300 0.001300 0.001300 0.001300 0.001300 0.000653 0.000653 0.000653 0.000653 0.000653 0.000653 0.000653 0.000653 0.000653 0.000653