### Required package
library(BSgenome.Hsapiens.UCSC.hg19)
library(tidyverse)

## Example data and results
path <- "path\to\data"
Panel_1 <- read.table(file.path(path, "Panel_Info_1_assay.txt"), header = T, quote = "", sep="\t", stringsAsFactors = FALSE)
Panel_2 <- read.table(file.path(path, "Panel_Info_2_assays.txt"), header = T, quote = "", sep="\t", stringsAsFactors = FALSE)


L_matrix_generation <- function(genomic_information){
  
  if(!all(c("Chromosome", "Start_Position", "End_Position", "SEQ_ASSAY_ID") %in% colnames(genomic_information))){
    stop("Please specify the column names of the argument genomic_information in (Chromosome, Start_Position, End_Position, SEQ_ASSAY_ID)")
  } 
  
  # mutation type categories
  # Alexandrov v3.4 mutation type order
  Alex_type <- c("A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "A[C>G]A", "A[C>G]C", "A[C>G]G", "A[C>G]T", "A[C>T]A", "A[C>T]C", "A[C>T]G", "A[C>T]T", 
                 "A[T>A]A", "A[T>A]C", "A[T>A]G", "A[T>A]T", "A[T>C]A", "A[T>C]C", "A[T>C]G", "A[T>C]T", "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T", 
                 "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T", "C[C>G]A", "C[C>G]C", "C[C>G]G", "C[C>G]T", "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T", 
                 "C[T>A]A", "C[T>A]C", "C[T>A]G", "C[T>A]T", "C[T>C]A", "C[T>C]C", "C[T>C]G", "C[T>C]T", "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T", 
                 "G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T", "G[C>G]A", "G[C>G]C", "G[C>G]G", "G[C>G]T", "G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T", 
                 "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T", "G[T>C]A", "G[T>C]C", "G[T>C]G", "G[T>C]T", "G[T>G]A", "G[T>G]C", "G[T>G]G", "G[T>G]T", 
                 "T[C>A]A", "T[C>A]C", "T[C>A]G", "T[C>A]T", "T[C>G]A", "T[C>G]C", "T[C>G]G", "T[C>G]T", "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T", 
                 "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T", "T[T>C]A", "T[T>C]C", "T[T>C]G", "T[T>C]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T")
  
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
  
  assay_CT_32_2 <- data.frame(t(assay_CT_32[, paste0(substr(Alex_type, 1, 1), substr(Alex_type, 3, 3), substr(Alex_type, 7, 7)), drop = FALSE]))/10^6
  rownames(assay_CT_32_2) <- Alex_type
  colnames(assay_CT_32_2) <- assay_CT_32$SEQ_ASSAY_ID
  
  return(assay_CT_32_2)
  
}



L_matrix_generation(Panel_1)
## Part of the results:
        UHN-48-V1
A[C>A]A  0.001487
A[C>A]C  0.001120
A[C>A]G  0.000426
A[C>A]T  0.001155
A[C>G]A  0.001487
A[C>G]C  0.001120
A[C>G]G  0.000426
A[C>G]T  0.001155
A[C>T]A  0.001487
A[C>T]C  0.001120
A[C>T]G  0.000426

L_matrix_generation(Panel_2)
## Part of the results:
        GRCC-CP1 UHN-48-V1
A[C>A]A 0.000883  0.001487
A[C>A]C 0.000656  0.001120
A[C>A]G 0.000278  0.000426
A[C>A]T 0.000687  0.001155
A[C>G]A 0.000883  0.001487
A[C>G]C 0.000656  0.001120
A[C>G]G 0.000278  0.000426
A[C>G]T 0.000687  0.001155
A[C>T]A 0.000883  0.001487
A[C>T]C 0.000656  0.001120
A[C>T]G 0.000278  0.000426



