
GeneratePanelSize <- function(genomic_information, Class = c("SBS", "DBS"), SBS_order = c("COSMIC", "signeR"), 
                              ref.genome="hg19"){

  # Check input arguments
  check_genomic_info(genomic_information)
  check_ref.genome(ref.genome)
  check_Class(Class)

  if (Class == "SBS") {
    check_Types(SBS_order)
    ret <- GeneratePanelSize_SBS(genomic_information, Types=SBS_order, ref.genome=ref.genome)
  } else {
    ret <- GeneratePanelSize_DBS(genomic_information, ref.genome=ref.genome)
  }
  ret
}

GeneratePanelSize_DBS <- function(genomic_information, ref.genome="hg19") {

    SEQ_ASSAY_ID <- NULL # To remove warning when compiling

    Seq_assay_GRanges <- GRanges(seqnames=sats_seqnames(genomic_information$Chromosome),
                                 IRanges(start = genomic_information$Start_Position, end=genomic_information$End_Position),
                                 strand = "+")
    
    # get sequences
    if (ref.genome == "hg19") {
      myseq = getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens,Seq_assay_GRanges) 
    } else {
      myseq = getSeq(BSgenome.Hsapiens.UCSC.hg38::Hsapiens,Seq_assay_GRanges)
    }
    Seq_assay_16 = dinucleotideFrequency(myseq) # key step: calculate dinucleotide frequency
    
    # Collapse the 16 possible dinucleotides to the 10 canonical DBS contexts.
    # Non-canonical contexts are reverse-complemented before summing so that
    # the returned rows match the COSMIC DBS78 orientation.
    dinucleotideIncluded = SATS_DBS_INCLUDED_DINUCLEOTIDES
    #all possible 16
    dinucleotide_16 = colnames(Seq_assay_16)
    idx = dinucleotide_16 %in% dinucleotideIncluded
    dinucleotide_10 = dinucleotide_16[idx]
    #inclcuded 10
    assay_DBS = Seq_assay_16[,dinucleotide_10]
    #excluded 6
    assay_DBS_6 = Seq_assay_16[,dinucleotide_16[!idx]]
    #rev names
    colnames(assay_DBS_6) = rev_dinucleotide(colnames(assay_DBS_6))
    #sort alphabetically
    assay_DBS_6 = assay_DBS_6[,sort(colnames(assay_DBS_6))]
    
    dinucleotide_10_set = setNames(1:10,dinucleotide_10)
    dinucleotide_10_set_idx = dinucleotide_10_set[colnames(assay_DBS_6)]
    #sum
    assay_DBS[,dinucleotide_10_set_idx] = assay_DBS[,dinucleotide_10_set_idx]+assay_DBS_6
    
    genomic_information_assay_DBS=bind_cols(genomic_information, as.data.frame(assay_DBS))
    
    
    #sum over genes for the whole assay
    assay_DBS_10 = genomic_information_assay_DBS %>% 
      group_by(SEQ_ASSAY_ID) %>% 
      summarise_at(colnames(assay_DBS),sum)
    
    assay_DBS_re <- data.frame(t(assay_DBS_10[, substr(SATS_DBS_ORDER, 1, 2), drop = FALSE]))/SATS_BASES_PER_MB
    colnames(assay_DBS_re) <- assay_DBS_10$SEQ_ASSAY_ID
    rownames(assay_DBS_re) <- SATS_DBS_ORDER

    return(assay_DBS_re)
}


rev_dinucleotide <- function(nt.seq){
    rev.seq = NULL
    for(i in 1:length(nt.seq)){
      x <- DNAString(nt.seq[i])
      rev.seq = c(rev.seq, as.character(reverseComplement(x)))
    }
    return(rev.seq)
}
  

## The following function is to create "L" matrix
## Please use the column names of the argument "genomic_information" identical to Chromosome, Start_Position, End_Position, SEQ_ASSAY_ID as in the above example
## Chromosome: chromosome number
## Start_Position: start position of targeted panel
## End_Position: end position of targeted panel
## SEQ_ASSAY_ID: distinguish different panels
## The unit of the returned L matrix is the number of trinucleotides per million base pairs
GeneratePanelSize_SBS <- function(genomic_information, Types = c("COSMIC", "signeR"), ref.genome="hg19"){
  
  SEQ_ASSAY_ID <- NULL

  # Enumerate trinucleotide opportunities in the order returned by
  # Biostrings::trinucleotideFrequency(). A/G-centered contexts are folded into
  # their C/T-centered reverse complements to match SBS96 convention.
  nt = SATS_BASES
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
  
  
  Seq_assay_GRanges <- GRanges(seqnames = sats_seqnames(genomic_information$Chromosome),
                               IRanges(start = genomic_information$Start_Position-1,
                               end=genomic_information$End_Position+1), strand = "+")
  
  Seq_assay_n = length(Seq_assay_GRanges)
  
  Seq_assay_64 = matrix(0,Seq_assay_n,64)
  colnames(Seq_assay_64) = tri_nt
  
  # Get sequences
  if (ref.genome == "hg19") {
    myseq = getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens,Seq_assay_GRanges) 
  } else {
    myseq = getSeq(BSgenome.Hsapiens.UCSC.hg38::Hsapiens,Seq_assay_GRanges)
  }
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
    assay_CT_32_2 <- data.frame(t(assay_CT_32[, paste0(substr(SATS_SBS_ORDER_COSMIC, 1, 1), substr(SATS_SBS_ORDER_COSMIC, 3, 3), 
                                  substr(SATS_SBS_ORDER_COSMIC, 7, 7)), drop = FALSE]))/SATS_BASES_PER_MB
    rownames(assay_CT_32_2) <- SATS_SBS_ORDER_COSMIC
  } else if(Types == "signeR"){
    assay_CT_32_2 <- data.frame(t(assay_CT_32[, substr(SATS_SBS_ORDER_SIGNER, 5, 7), drop = FALSE]))/SATS_BASES_PER_MB
    rownames(assay_CT_32_2) <- SATS_SBS_ORDER_SIGNER  
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
  if (any(!idx)) {
    warning(sum(!idx), " sample(s) have SEQ_ASSAY_ID values not present in Panel_context and were removed")
  }
  L <- Panel_context[, Patient_Info$SEQ_ASSAY_ID[idx], drop = FALSE]
  colnames(L) <- Patient_Info$PATIENT_ID[idx]
  #if(sum(idx) != nrow(Patient_Info)){
  #  warning(sprintf("There are patients for whom the panel context has not been provided in Patient_Info. 
  #  Patient_Info contains %d sequence assays, but Panel_context has only %d sequence assays", 
  # length(unique(Patient_Info$SEQ_ASSAY_ID)), ncol(Panel_context)))
  #} 
  
  return(L)
}

GenerateLMatrix <- function(Panel_context, Patient_Info, Class = c("SBS", "DBS"),
                            SBS_order = c("COSMIC", "signeR"), ref.genome = "hg19") {

  check_Patient_Info(Patient_Info)
  Patient_Info <- standardize_Patient_Info(Patient_Info)

  if (is_genomic_info_input(Panel_context)) {
    Class <- match.arg(Class)
    SBS_order <- match.arg(SBS_order)
    Panel_context <- GeneratePanelSize(genomic_information = Panel_context,
                                       Class = Class, SBS_order = SBS_order,
                                       ref.genome = ref.genome)
  } else {
    check_Panel_context(Panel_context)
  }

  ret <- L_matrix_generation(Panel_context, Patient_Info)
  ret
}



