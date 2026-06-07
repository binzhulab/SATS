GenerateVMatrix <- function(mutation_record, Class = c("SBS", "DBS"),
                            ref.genome = "hg19", mutation_order = NULL) {

  Class <- match.arg(Class)

  check_mutation_record(mutation_record)
  check_ref.genome(ref.genome)

  if (is.null(mutation_order)) {
    mutation_order <- get_default_mutation_order(Class)
  }
  check_mutation_order(mutation_order)

  if (Class == "SBS") {
    ret <- GenerateVMatrix_SBS(mutation_record, ref.genome, mutation_order)
  } else {
    ret <- GenerateVMatrix_DBS(mutation_record, ref.genome, mutation_order)
  }

  ret
}

GenerateVMatrix_SBS <- function(mutation_record, ref.genome, mutation_order) {

  x <- prepare_mutation_record(mutation_record)
  x <- x[x$Variant_Type == "SNP", , drop = FALSE]
  if (!nrow(x)) stop("ERROR: mutation_record contains no SNP records")

  x$Start_Position <- suppressWarnings(as.integer(x$Start_Position))

  valid <- !is.na(x$Start_Position) &
           x$Start_Position > 1L &
           nchar(x$Reference_Allele) == 1L &
           nchar(x$Tumor_Seq_Allele2) == 1L &
           x$Reference_Allele %in% c("A", "C", "G", "T") &
           x$Tumor_Seq_Allele2 %in% c("A", "C", "G", "T") &
           x$Reference_Allele != x$Tumor_Seq_Allele2

  x <- drop_invalid_mutations(x, valid, "SBS")
  genome <- get_reference_genome(ref.genome)

  flank5 <- GRanges(seqnames = sats_seqnames(x$Chromosome),
                    IRanges(start = x$Start_Position - 1L,
                            end = x$Start_Position - 1L),
                    strand = "+")
  flank3 <- GRanges(seqnames = sats_seqnames(x$Chromosome),
                    IRanges(start = x$Start_Position + 1L,
                            end = x$Start_Position + 1L),
                    strand = "+")

  seq5 <- as.character(getSeq(genome, flank5))
  seq3 <- as.character(getSeq(genome, flank3))

  ref_tri <- paste0(seq5, x$Reference_Allele, seq3)
  mut_tri <- paste0(seq5, x$Tumor_Seq_Allele2, seq3)

  ref_tri_ct <- ref_tri
  mut_tri_ct <- mut_tri

  flip_idx <- x$Reference_Allele %in% c("A", "G")
  if (any(flip_idx)) {
    ref_tri_ct[flip_idx] <- as.character(reverseComplement(DNAStringSet(ref_tri_ct[flip_idx])))
    mut_tri_ct[flip_idx] <- as.character(reverseComplement(DNAStringSet(mut_tri_ct[flip_idx])))
  }

  category <- paste0(substr(ref_tri_ct, 1L, 1L), "[",
                     substr(ref_tri_ct, 2L, 2L), ">",
                     substr(mut_tri_ct, 2L, 2L), "]",
                     substr(ref_tri_ct, 3L, 3L))

  make_count_matrix(category, x$Tumor_Sample_Barcode, mutation_order)
}

GenerateVMatrix_DBS <- function(mutation_record, ref.genome, mutation_order) {

  x <- prepare_mutation_record(mutation_record)
  x <- x[x$Variant_Type == "DNP", , drop = FALSE]
  if (!nrow(x)) stop("ERROR: mutation_record contains no DNP records")

  x$Start_Position <- suppressWarnings(as.integer(x$Start_Position))
  x$End_Position <- suppressWarnings(as.integer(x$End_Position))

  valid <- !is.na(x$Start_Position) &
           !is.na(x$End_Position) &
           x$End_Position == x$Start_Position + 1L &
           nchar(x$Reference_Allele) == 2L &
           nchar(x$Tumor_Seq_Allele2) == 2L &
           grepl("^[ACGT][ACGT]$", x$Reference_Allele) &
           grepl("^[ACGT][ACGT]$", x$Tumor_Seq_Allele2) &
           x$Reference_Allele != x$Tumor_Seq_Allele2

  x <- drop_invalid_mutations(x, valid, "DBS")
  genome <- get_reference_genome(ref.genome)

  base5 <- GRanges(seqnames = sats_seqnames(x$Chromosome),
                   IRanges(start = x$Start_Position, end = x$Start_Position),
                   strand = "+")
  base3 <- GRanges(seqnames = sats_seqnames(x$Chromosome),
                   IRanges(start = x$End_Position, end = x$End_Position),
                   strand = "+")

  ref_di <- paste0(as.character(getSeq(genome, base5)),
                   as.character(getSeq(genome, base3)))
  mut_di <- x$Tumor_Seq_Allele2

  reference_mismatch <- ref_di != x$Reference_Allele
  if (any(reference_mismatch)) {
    warning(sum(reference_mismatch), " DBS record(s) have Reference_Allele values that do not match the reference genome and were removed")
    ref_di <- ref_di[!reference_mismatch]
    mut_di <- mut_di[!reference_mismatch]
    x <- x[!reference_mismatch, , drop = FALSE]
  }

  if (!nrow(x)) stop("ERROR: no usable DBS records remain after reference-genome checking")

  same_count <- mapply(function(s1, s2) {
    a <- strsplit(s1, "", fixed = TRUE)[[1]]
    b <- strsplit(s2, "", fixed = TRUE)[[1]]
    sum(a == b)
  }, ref_di, mut_di)

  misspecified <- same_count > 0L
  if (any(misspecified)) {
    warning(sum(misspecified), " DBS record(s) altered only one base and were removed")
    ref_di <- ref_di[!misspecified]
    mut_di <- mut_di[!misspecified]
    x <- x[!misspecified, , drop = FALSE]
  }

  if (!nrow(x)) stop("ERROR: no usable DBS records remain after DBS mutation checking")

  ref_di_10 <- ref_di
  mut_di_10 <- mut_di
  included <- c("AC", "AT", "CC", "CG", "CT", "GC", "TA", "TC", "TG", "TT")

  flip_idx <- !ref_di %in% included
  if (any(flip_idx)) {
    ref_di_10[flip_idx] <- as.character(reverseComplement(DNAStringSet(ref_di_10[flip_idx])))
    mut_di_10[flip_idx] <- as.character(reverseComplement(DNAStringSet(mut_di_10[flip_idx])))
  }

  flip_idx_at <- ref_di %in% "AT" & !mut_di_10 %in% c("CA", "CC", "CG", "GA", "GC", "TA")
  flip_idx_ta <- ref_di %in% "TA" & !mut_di_10 %in% c("AT", "CG", "CT", "GC", "GG", "GT")
  flip_idx_cg <- ref_di %in% "CG" & !mut_di_10 %in% c("AT", "GC", "GT", "TA", "TC", "TT")
  flip_idx_gc <- ref_di %in% "GC" & !mut_di_10 %in% c("AA", "AG", "AT", "CA", "CG", "TA")

  flip_palindrome <- flip_idx_at | flip_idx_ta | flip_idx_cg | flip_idx_gc
  if (any(flip_palindrome)) {
    mut_di_10[flip_palindrome] <- as.character(reverseComplement(DNAStringSet(mut_di_10[flip_palindrome])))
  }

  category <- paste0(ref_di_10, mut_di_10)

  make_count_matrix(category, x$Tumor_Sample_Barcode, mutation_order)
}

prepare_mutation_record <- function(mutation_record) {

  x <- as.data.frame(mutation_record, stringsAsFactors = FALSE)
  x$Chromosome <- as.character(x$Chromosome)
  x$Variant_Type <- toupper(as.character(x$Variant_Type))
  x$Reference_Allele <- toupper(as.character(x$Reference_Allele))
  x$Tumor_Seq_Allele2 <- toupper(as.character(x$Tumor_Seq_Allele2))
  x$Tumor_Sample_Barcode <- as.character(x$Tumor_Sample_Barcode)

  x
}

drop_invalid_mutations <- function(x, valid, Class) {

  if (any(!valid)) {
    warning(sum(!valid), " invalid ", Class, " mutation record(s) were removed")
  }
  x <- x[valid, , drop = FALSE]
  if (!nrow(x)) {
    stop("ERROR: no usable ", Class, " records remain after input checking")
  }

  x
}

make_count_matrix <- function(category, sample_id, mutation_order) {

  keep <- category %in% mutation_order
  if (any(!keep)) {
    warning(sum(!keep), " mutation record(s) have categories not present in mutation_order and were removed")
  }

  category <- category[keep]
  sample_id <- as.character(sample_id[keep])

  if (!length(category)) {
    stop("ERROR: no mutation records remain after matching categories to mutation_order")
  }

  sample_order <- unique(sample_id)
  counts <- table(factor(category, levels = mutation_order),
                  factor(sample_id, levels = sample_order))
  ret <- as.matrix(counts)
  mode(ret) <- "numeric"
  rownames(ret) <- mutation_order
  colnames(ret) <- sample_order

  ret
}

get_default_mutation_order <- function(Class) {

  ref <- get_ref_tmb()
  if (Class == "SBS") {
    ret <- rownames(ref[["TMB_SBS_v3.4"]])
  } else {
    ret <- rownames(ref[["TMB_DBS_v3.4"]])
  }

  ret
}

get_ref_tmb <- function() {

  env <- new.env(parent = emptyenv())
  utils::data("RefTMB", package = "SATS", envir = env)
  if (!exists("RefTMB", envir = env, inherits = FALSE)) {
    stop("ERROR: RefTMB data object could not be loaded")
  }

  env$RefTMB
}

get_reference_genome <- function(ref.genome) {

  if (ref.genome == "hg19") {
    ret <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens
  } else {
    ret <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
  }

  ret
}
