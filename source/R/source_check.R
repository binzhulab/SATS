# Functions to check for errors in the input arguments

check_L_W_H <- function(L, W, H) {

  check_mat_df(L, "L")
  check_mat_df(W, "W")
  check_mat_df(H, "H")
  if (nrow(H) != ncol(W)) stop("ERROR: nrow(H) != ncol(W)")
  if (ncol(L) != ncol(H)) stop("ERROR: ncol(H) != ncol(L)")
  if (nrow(L) != nrow(W)) stop("ERROR: nrow(L) != nrow(W)")

  NULL
}

align_L_W_H <- function(L, W, H) {

  L <- align_named_axis(H, L, "H", "L", "columns", "columns", "sample IDs")
  W <- align_named_axis(H, W, "H", "W", "rows", "columns", "signature names")
  W <- align_named_axis(L, W, "L", "W", "rows", "rows", "mutation context names")

  list(L=L, W=W, H=H)
}

check_L_W_V <- function(L, W, V) {

  check_mat_df(L, "L")
  check_mat_df(W, "W")
  check_mat_df(V, "V")
  if (any(dim(L) != dim(V))) stop("ERROR: dim(L) != dim(V)")
  if (nrow(L) != nrow(W)) stop("ERROR: nrow(L) != nrow(W)")

  NULL
}

align_L_W_V <- function(L, W, V) {

  L <- align_named_axis(V, L, "V", "L", "rows", "rows", "mutation context names")
  W <- align_named_axis(V, W, "V", "W", "rows", "rows", "mutation context names")
  L <- align_named_axis(V, L, "V", "L", "columns", "columns", "sample IDs")

  list(L=L, W=W, V=V)
}

align_named_axis <- function(reference, target, reference_name, target_name,
                             reference_axis, target_axis, label) {

  reference_names <- get_axis_names(reference, reference_axis)
  target_names    <- get_axis_names(target, target_axis)

  if (is.null(reference_names) || is.null(target_names)) return(target)

  if (anyDuplicated(reference_names)) {
    stop(paste0("ERROR: duplicate ", label, " found in ", reference_name))
  }
  if (anyDuplicated(target_names)) {
    stop(paste0("ERROR: duplicate ", label, " found in ", target_name))
  }

  missing <- setdiff(reference_names, target_names)
  extra   <- setdiff(target_names, reference_names)
  if (length(missing) || length(extra)) {
    msg <- paste0("ERROR: ", label, " in ", reference_name, " and ",
                  target_name, " do not match")
    if (length(missing)) {
      msg <- paste0(msg, "; missing from ", target_name, ": ",
                    format_name_list(missing))
    }
    if (length(extra)) {
      msg <- paste0(msg, "; extra in ", target_name, ": ",
                    format_name_list(extra))
    }
    stop(msg)
  }

  if (identical(reference_names, target_names)) return(target)

  reorder_axis(target, target_axis, reference_names)
}

get_axis_names <- function(x, axis) {

  if (axis == "rows") return(rownames(x))
  if (axis == "columns") return(colnames(x))
  stop("ERROR: axis must be 'rows' or 'columns'")
}

reorder_axis <- function(x, axis, names_order) {

  if (axis == "rows") return(x[names_order, , drop=FALSE])
  if (axis == "columns") return(x[, names_order, drop=FALSE])
  stop("ERROR: axis must be 'rows' or 'columns'")
}

format_name_list <- function(x, max_names=5) {

  x <- as.character(x)
  if (length(x) > max_names) {
    x <- c(x[seq_len(max_names)], paste0("... plus ", length(x) - max_names, " more"))
  }
  paste(x, collapse=", ")
}


check_mat_df <- function(x, nm, check.finite=1) {

  if (!is.matrix(x) && !is.data.frame(x)) {
    stop(paste0("ERROR: ", nm, " must be a matrix or data frame"))
  }
  if (!nrow(x)) stop(paste0("ERROR: ", nm, " has no rows"))
  if (!ncol(x)) stop(paste0("ERROR: ", nm, " has no columns"))
  if (check.finite) {
    if (any(!is.finite(as.matrix(x)))) {
      stop(paste0("ERROR: ", nm, " contains non-finite values"))
    }
  }

  NULL
}

check_number <- function(x, nm, min=NULL, pos=FALSE, max=NULL) {

  if (length(x) != 1) stop(paste0("ERROR: ", nm, " must be a single numeric value"))
  if (!is.null(min) && (x < min)) {
    stop(paste0("ERROR: ", nm, " must be greater than or equal to ", min))
  }
  if (pos && (x <= 0)) {
    stop(paste0("ERROR: ", nm, " must be positive"))
  }
  if (!is.null(max) && (x > max)) {
    stop(paste0("ERROR: ", nm, " must be less than or equal to ", max))
  }

  NULL

}

isString <- function(x) {
  is.character(x) && (length(x) == 1)
}

check_Types <- function(x) {

  valid <- c("COSMIC", "signeR")
  err   <- "ERROR: Types must be 'COSMIC' or 'signeR'"

  if(length(x) == 2){
    stop("Please specify mutation types: available mutations types are either \"COSMIC\" or \"signeR\"")
  }

  if (!isString(x)) stop(err)
  if (!(x %in% valid)) stop(err)  
    
  NULL
}

check_Class <- function(x) {

  valid <- c("SBS", "DBS")
  err   <- "ERROR: Class must be 'SBS' or 'DBS'"

  if(length(x) != 1) stop(err)
  if (!isString(x)) stop(err)
  if (!(x %in% valid)) stop(err)  
    
  NULL
}


check_COSMICv <- function(x) {

  valid <- c("v3.2", "v3.4")
  err   <- "ERROR: COSMICv must be 'v3.2' or 'v3.4'"

  if (!isString(x)) stop(err)
  if (!(x %in% valid)) stop(err)  
    
  NULL
}

check_ref.genome <- function(x) {

  valid <- c("hg19", "hg38")
  err   <- "ERROR: ref.genome must be 'hg19' or 'hg38'"

  if (!isString(x)) stop(err)
  if (!(x %in% valid)) stop(err)  
    
  NULL
}

check_genomic_info <- function(x) {

  if (!is.data.frame(x)) {
    stop("ERROR: genomic_information must be a data frame")
  }
  cols <- c("Chromosome", "Start_Position", "End_Position", "SEQ_ASSAY_ID")
  tmp  <- !(cols %in% colnames(x))
  if (any(tmp)) {
    miss  <- cols[tmp]
    nmiss <- length(miss)
    mstr  <- paste0("'", miss, "'")
    mstr  <- paste0(mstr, collapse=", ")
    if (nmiss == 1) {
      msg <- paste0("ERROR: column ", mstr, " not found in genomic_information")
    } else {
      msg <- paste0("ERROR: columns ", mstr, " not found in genomic_information")
    }
    stop(msg)
  }

  NULL
}

is_genomic_info_input <- function(x) {

  is.data.frame(x) &&
    all(c("Chromosome", "Start_Position", "End_Position", "SEQ_ASSAY_ID") %in% colnames(x))
}

sats_seqnames <- function(x) {

  ret <- as.character(x)
  has_chr <- grepl("^chr", ret, ignore.case = TRUE)
  ret[has_chr] <- paste0("chr", sub("^chr", "", ret[has_chr], ignore.case = TRUE))
  ret[!has_chr] <- paste0("chr", ret[!has_chr])
  ret[ret %in% c("chrMT", "chrMt")] <- "chrM"

  ret
}

check_mutation_record <- function(x) {

  if (!is.data.frame(x)) {
    stop("ERROR: mutation_record must be a data frame")
  }

  cols <- c("Chromosome", "Start_Position", "End_Position", "Variant_Type",
            "Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode")
  tmp <- !(cols %in% colnames(x))
  if (any(tmp)) {
    miss <- cols[tmp]
    mstr <- paste0("'", miss, "'", collapse = ", ")
    if (length(miss) == 1) {
      stop(paste0("ERROR: column ", mstr, " not found in mutation_record"))
    } else {
      stop(paste0("ERROR: columns ", mstr, " not found in mutation_record"))
    }
  }

  NULL
}

check_clinical_sample <- function(x) {

  if (!is.data.frame(x)) {
    stop("ERROR: clinical_sample must be a data frame")
  }

  cols <- c("SAMPLE_ID", "SEQ_ASSAY_ID")
  tmp <- !(cols %in% colnames(x))
  if (any(tmp)) {
    miss <- cols[tmp]
    mstr <- paste0("'", miss, "'", collapse = ", ")
    if (length(miss) == 1) {
      stop(paste0("ERROR: column ", mstr, " not found in clinical_sample"))
    } else {
      stop(paste0("ERROR: columns ", mstr, " not found in clinical_sample"))
    }
  }

  NULL
}

check_mutation_order <- function(x) {

  if (!is.character(x)) {
    stop("ERROR: mutation_order must be a character vector")
  }
  if (!length(x)) {
    stop("ERROR: mutation_order has no entries")
  }
  if (any(is.na(x))) {
    stop("ERROR: mutation_order contains missing values")
  }
  if (anyDuplicated(x)) {
    stop("ERROR: mutation_order contains duplicated entries")
  }

  NULL
}

check_Panel_context <- function(x) {

  if (!is.data.frame(x)) stop("ERROR: Panel_context must be a data frame")
  check_mat_df(x, "Panel_context", check.finite=0)

  NULL
}

check_Patient_Info <- function(x) {

  if (!is.data.frame(x)) stop("ERROR: Patient_Info must be a data frame")
  check_mat_df(x, "Patient_Info", check.finite=0)
  cx <- colnames(x)
  if (!("SEQ_ASSAY_ID" %in% cx)) {
    stop("ERROR: Patient_Info must contain the column 'SEQ_ASSAY_ID'")
  }
  if (!("PATIENT_ID" %in% cx) && !("SAMPLE_ID" %in% cx)) {
    stop("ERROR: Patient_Info must contain either 'PATIENT_ID' or 'SAMPLE_ID'")
  }

  NULL
}

standardize_Patient_Info <- function(x) {

  ret <- as.data.frame(x, stringsAsFactors = FALSE)
  if (!("PATIENT_ID" %in% colnames(ret)) && "SAMPLE_ID" %in% colnames(ret)) {
    ret$PATIENT_ID <- ret$SAMPLE_ID
  }

  ret <- ret[, c("PATIENT_ID", "SEQ_ASSAY_ID"), drop = FALSE]
  ret$PATIENT_ID <- as.character(ret$PATIENT_ID)
  ret$SEQ_ASSAY_ID <- as.character(ret$SEQ_ASSAY_ID)

  ret
}
