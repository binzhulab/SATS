
## The following function is to create "L" matrix
## Please use the column names of the argument "genomic_information" identical to Chromosome, Start_Position, End_Position, SEQ_ASSAY_ID as in the above example
## Chromosome: chromsome number
## Start_Position: start position of targeted panel
## End_Position: end position of targeted panel
## SEQ_ASSAY_ID: distinguish different panels
## The unit of the returned L matrix is the number of trinucleotides per million base pairs
GeneratePanelSize <- function(genomic_information, Types = c("COSMIC", "signeR")){
  
  # Check input arguments
  check_genomic_info(genomic_information)
  check_Types(Types)
  
  SEQ_ASSAY_ID <- NULL

  # Mutation type categories, there are 96 possible single base substitutions.
  # Alexandrov v3.2 mutation type order. .
  COSMIC <- c("A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "C[C>A]A", 
              "C[C>A]C", "C[C>A]G", "C[C>A]T", "G[C>A]A", "G[C>A]C", 
              "G[C>A]G", "G[C>A]T", "T[C>A]A", "T[C>A]C", "T[C>A]G", 
              "T[C>A]T", "A[C>G]A", "A[C>G]C", "A[C>G]G", "A[C>G]T", 
              "C[C>G]A", "C[C>G]C", "C[C>G]G", "C[C>G]T", "G[C>G]A", 
              "G[C>G]C", "G[C>G]G", "G[C>G]T", "T[C>G]A", "T[C>G]C", 
              "T[C>G]G", "T[C>G]T", "A[C>T]A", "A[C>T]C", "A[C>T]G", 
              "A[C>T]T", "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T", 
              "G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T", "T[C>T]A", 
              "T[C>T]C", "T[C>T]G", "T[C>T]T", "A[T>A]A", "A[T>A]C", 
              "A[T>A]G", "A[T>A]T", "C[T>A]A", "C[T>A]C", "C[T>A]G", 
              "C[T>A]T", "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T", 
              "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T", "A[T>C]A", 
              "A[T>C]C", "A[T>C]G", "A[T>C]T", "C[T>C]A", "C[T>C]C", 
              "C[T>C]G", "C[T>C]T", "G[T>C]A", "G[T>C]C", "G[T>C]G", 
              "G[T>C]T", "T[T>C]A", "T[T>C]C", "T[T>C]G", "T[T>C]T", 
              "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T", "C[T>G]A", 
              "C[T>G]C", "C[T>G]G", "C[T>G]T", "G[T>G]A", "G[T>G]C", 
              "G[T>G]G", "G[T>G]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", 
              "T[T>G]T")

  signeR <- c("C>A:ACA", "C>A:ACC", "C>A:ACG", "C>A:ACT", "C>A:CCA", 
              "C>A:CCC", "C>A:CCG", "C>A:CCT", "C>A:GCA", "C>A:GCC", 
              "C>A:GCG", "C>A:GCT", "C>A:TCA", "C>A:TCC", "C>A:TCG", 
              "C>A:TCT", "C>G:ACA", "C>G:ACC", "C>G:ACG", "C>G:ACT", 
              "C>G:CCA", "C>G:CCC", "C>G:CCG", "C>G:CCT", "C>G:GCA", 
              "C>G:GCC", "C>G:GCG", "C>G:GCT", "C>G:TCA", "C>G:TCC", 
              "C>G:TCG", "C>G:TCT", "C>T:ACA", "C>T:ACC", "C>T:ACG", 
              "C>T:ACT", "C>T:CCA", "C>T:CCC", "C>T:CCG", "C>T:CCT", 
              "C>T:GCA", "C>T:GCC", "C>T:GCG", "C>T:GCT", "C>T:TCA", 
              "C>T:TCC", "C>T:TCG", "C>T:TCT", "T>A:ATA", "T>A:ATC", 
              "T>A:ATG", "T>A:ATT", "T>A:CTA", "T>A:CTC", "T>A:CTG", 
              "T>A:CTT", "T>A:GTA", "T>A:GTC", "T>A:GTG", "T>A:GTT", 
              "T>A:TTA", "T>A:TTC", "T>A:TTG", "T>A:TTT", "T>C:ATA", 
              "T>C:ATC", "T>C:ATG", "T>C:ATT", "T>C:CTA", "T>C:CTC", 
              "T>C:CTG", "T>C:CTT", "T>C:GTA", "T>C:GTC", "T>C:GTG", 
              "T>C:GTT", "T>C:TTA", "T>C:TTC", "T>C:TTG", "T>C:TTT", 
              "T>G:ATA", "T>G:ATC", "T>G:ATG", "T>G:ATT", "T>G:CTA", 
              "T>G:CTC", "T>G:CTG", "T>G:CTT", "T>G:GTA", "T>G:GTC", 
              "T>G:GTG", "T>G:GTT", "T>G:TTA", "T>G:TTC", "T>G:TTG", 
              "T>G:TTT")
  
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
                               IRanges(start = genomic_information$Start_Position-1, 
                               end=genomic_information$End_Position+1), strand = "+")
  
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
    assay_CT_32_2 <- data.frame(t(assay_CT_32[, paste0(substr(COSMIC, 1, 1), substr(COSMIC, 3, 3), 
                                  substr(COSMIC, 7, 7)), drop = FALSE]))/10^6
    rownames(assay_CT_32_2) <- COSMIC
  } else if(Types == "signeR"){
    assay_CT_32_2 <- data.frame(t(assay_CT_32[, substr(signeR, 5, 7), drop = FALSE]))/10^6
    rownames(assay_CT_32_2) <- signeR  
  } 
  
  colnames(assay_CT_32_2) <- assay_CT_32$SEQ_ASSAY_ID
  
  return(assay_CT_32_2)
  
}

