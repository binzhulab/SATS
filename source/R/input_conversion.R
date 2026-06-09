ReadVCFAsMutationRecord <- function(vcf_file, sample_id = NULL,
                                    keep_filter = c("PASS", ".")) {

  if (!isString(vcf_file)) {
    stop("ERROR: vcf_file must be a single file path")
  }
  if (!file.exists(vcf_file)) {
    stop("ERROR: vcf_file does not exist")
  }
  if (!is.null(sample_id) && !isString(sample_id)) {
    stop("ERROR: sample_id must be NULL or a single character string")
  }

  lines <- read_text_lines(vcf_file)
  header_idx <- grep("^#CHROM", lines)
  if (!length(header_idx)) {
    stop("ERROR: VCF header line beginning with '#CHROM' was not found")
  }
  header_idx <- header_idx[1]
  header <- sub("^#", "", lines[header_idx])
  body <- if (header_idx == length(lines)) {
    character()
  } else {
    lines[(header_idx + 1L):length(lines)]
  }
  body <- body[!is.na(body) & nzchar(body) & !grepl("^#", body)]
  if (!length(body)) {
    stop("ERROR: VCF file contains no variant records")
  }

  vcf <- utils::read.delim(text = paste(c(header, body), collapse = "\n"),
                           header = TRUE, sep = "\t", quote = "",
                           comment.char = "", stringsAsFactors = FALSE,
                           check.names = FALSE)

  required <- c("CHROM", "POS", "REF", "ALT")
  missing <- setdiff(required, colnames(vcf))
  if (length(missing)) {
    stop("ERROR: VCF file is missing required column(s): ",
         paste(missing, collapse = ", "))
  }

  if (!is.null(keep_filter) && "FILTER" %in% colnames(vcf)) {
    keep <- vcf$FILTER %in% keep_filter
    if (any(!keep)) {
      warning(sum(!keep), " VCF record(s) were removed by FILTER")
    }
    vcf <- vcf[keep, , drop = FALSE]
    if (!nrow(vcf)) stop("ERROR: no VCF records remain after FILTER selection")
  }

  sample_cols <- setdiff(colnames(vcf),
                         c("CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                           "FILTER", "INFO", "FORMAT"))
  if (length(sample_cols) > 1L) {
    stop("ERROR: VCF contains multiple sample columns; provide a single-sample VCF")
  }

  if (is.null(sample_id)) {
    if (length(sample_cols) == 1L) {
      sample_id <- sample_cols
    } else {
      sample_id <- guess_sample_id_from_file(vcf_file)
    }
  }

  pos <- suppressWarnings(as.integer(vcf$POS))
  alt_list <- strsplit(as.character(vcf$ALT), ",", fixed = TRUE)
  alt_n <- lengths(alt_list)
  if (any(alt_n > 1L)) {
    warning(sum(alt_n > 1L), " multiallelic VCF record(s) were split by ALT allele")
  }

  idx <- rep(seq_len(nrow(vcf)), alt_n)
  alt <- unlist(alt_list, use.names = FALSE)
  ref <- as.character(vcf$REF[idx])
  start <- pos[idx]
  end <- start + nchar(ref) - 1L

  ret <- data.frame(
    Chromosome = as.character(vcf$CHROM[idx]),
    Start_Position = start,
    End_Position = end,
    Variant_Type = classify_vcf_variant(ref, alt),
    Reference_Allele = toupper(ref),
    Tumor_Seq_Allele2 = toupper(alt),
    Tumor_Sample_Barcode = sample_id,
    stringsAsFactors = FALSE
  )

  valid <- !is.na(ret$Start_Position) &
           !is.na(ret$End_Position) &
           nzchar(ret$Chromosome) &
           nzchar(ret$Reference_Allele) &
           nzchar(ret$Tumor_Seq_Allele2) &
           ret$Tumor_Seq_Allele2 != "."
  if (any(!valid)) {
    warning(sum(!valid), " invalid VCF-derived mutation record(s) were removed")
  }
  ret <- ret[valid, , drop = FALSE]
  if (!nrow(ret)) {
    stop("ERROR: no usable mutation records were produced from the VCF file")
  }

  ret
}

ReadBEDAsPanelInfo <- function(bed_file, seq_assay_id = NULL,
                               has_header = FALSE,
                               chromosome_col = 1,
                               start_col = 2,
                               end_col = 3,
                               seq_assay_col = NULL,
                               name_col = NULL) {

  if (!isString(bed_file)) {
    stop("ERROR: bed_file must be a single file path")
  }
  if (!file.exists(bed_file)) {
    stop("ERROR: bed_file does not exist")
  }
  if (is.null(seq_assay_col)) {
    if (!isString(seq_assay_id)) {
      stop("ERROR: provide seq_assay_id or seq_assay_col")
    }
  }

  bed <- read_tab_file(bed_file, header = has_header, comment.char = "#")
  if (!nrow(bed)) stop("ERROR: BED file contains no target intervals")

  chrom <- as.character(extract_table_column(bed, chromosome_col, "chromosome_col"))
  bed_start <- suppressWarnings(as.integer(extract_table_column(bed, start_col, "start_col")))
  bed_end <- suppressWarnings(as.integer(extract_table_column(bed, end_col, "end_col")))

  if (is.null(seq_assay_col)) {
    assay <- rep(seq_assay_id, nrow(bed))
  } else {
    assay <- as.character(extract_table_column(bed, seq_assay_col, "seq_assay_col"))
  }

  valid <- !is.na(bed_start) &
           !is.na(bed_end) &
           bed_start >= 0L &
           bed_end > bed_start &
           nzchar(chrom) &
           nzchar(assay)
  if (any(!valid)) {
    warning(sum(!valid), " invalid BED interval(s) were removed")
  }

  ret <- data.frame(
    Chromosome = chrom[valid],
    Start_Position = bed_start[valid] + 1L,
    End_Position = bed_end[valid],
    SEQ_ASSAY_ID = assay[valid],
    stringsAsFactors = FALSE
  )

  if (!is.null(name_col)) {
    name <- as.character(extract_table_column(bed, name_col, "name_col"))
    ret$Hugo_Symbol <- name[valid]
  }

  if (!nrow(ret)) {
    stop("ERROR: no usable panel intervals were produced from the BED file")
  }

  ret
}

read_text_lines <- function(path) {

  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    base::file(path, open = "rt")
  }
  on.exit(close(con))
  readLines(con, warn = FALSE)
}

read_tab_file <- function(path, header, comment.char) {

  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    base::file(path, open = "rt")
  }
  on.exit(close(con))
  utils::read.delim(con, header = header, sep = "\t", quote = "",
                    comment.char = comment.char, stringsAsFactors = FALSE,
                    check.names = FALSE)
}

extract_table_column <- function(x, col, arg_name) {

  if (is.numeric(col)) {
    if (length(col) != 1L || is.na(col) || col < 1L || col > ncol(x)) {
      stop("ERROR: ", arg_name, " is outside the available column range")
    }
    return(x[[col]])
  }

  if (isString(col)) {
    if (!(col %in% colnames(x))) {
      stop("ERROR: ", arg_name, " column was not found in the input table")
    }
    return(x[[col]])
  }

  stop("ERROR: ", arg_name, " must be a column index or column name")
}

guess_sample_id_from_file <- function(file) {

  base <- basename(file)
  sub("\\.vcf(\\.gz)?$", "", base, ignore.case = TRUE)
}

classify_vcf_variant <- function(ref, alt) {

  ref <- toupper(ref)
  alt <- toupper(alt)
  ret <- rep("OTHER", length(ref))

  simple <- grepl("^[ACGT]+$", ref) & grepl("^[ACGT]+$", alt)
  ret[simple & nchar(ref) == 1L & nchar(alt) == 1L] <- "SNP"
  ret[simple & nchar(ref) == 2L & nchar(alt) == 2L] <- "DNP"
  ret[simple & nchar(ref) > nchar(alt)] <- "DEL"
  ret[simple & nchar(ref) < nchar(alt)] <- "INS"
  ret[simple & nchar(ref) == nchar(alt) & nchar(ref) > 2L] <- "ONP"

  ret
}
