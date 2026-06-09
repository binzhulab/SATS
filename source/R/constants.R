# Internal package constants.
#
# These objects collect fixed biological definitions and supported option values
# in one place. They are not user-facing parameters because changing them would
# change the mutation-channel definitions used by SATS.

SATS_SUPPORTED_CLASSES <- c("SBS", "DBS")
SATS_SUPPORTED_SBS_ORDERS <- c("COSMIC", "signeR")
SATS_SUPPORTED_COSMIC_VERSIONS <- c("v3.2", "v3.4")
SATS_SUPPORTED_GENOMES <- c("hg19", "hg38")

SATS_BASES <- c("A", "C", "G", "T")
SATS_BASES_PER_MB <- 1e6

# COSMIC-style SBS96 order. This is the mutation-channel order used for SATS
# COSMIC-compatible SBS matrices; it is independent of the COSMIC signature
# catalogue version selected during reference-signature mapping.
SATS_SBS_ORDER_COSMIC <- c(
  "A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "C[C>A]A",
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
  "T[T>G]T"
)

# SBS96 order used by the signeR package.
SATS_SBS_ORDER_SIGNER <- c(
  "C>A:ACA", "C>A:ACC", "C>A:ACG", "C>A:ACT", "C>A:CCA",
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
  "T>G:TTT"
)

# Canonical DBS78 order and the 10 pyrimidine-normalized dinucleotide contexts.
SATS_DBS_ORDER <- c(
  "ACCA", "ACCG", "ACCT", "ACGA", "ACGG", "ACGT", "ACTA", "ACTG", "ACTT",
  "ATCA", "ATCC", "ATCG", "ATGA", "ATGC", "ATTA",
  "CCAA", "CCAG", "CCAT", "CCGA", "CCGG", "CCGT", "CCTA", "CCTG", "CCTT",
  "CGAT", "CGGC", "CGGT", "CGTA", "CGTC", "CGTT",
  "CTAA", "CTAC", "CTAG", "CTGA", "CTGC", "CTGG", "CTTA", "CTTC", "CTTG",
  "GCAA", "GCAG", "GCAT", "GCCA", "GCCG", "GCTA",
  "TAAT", "TACG", "TACT", "TAGC", "TAGG", "TAGT",
  "TCAA", "TCAG", "TCAT", "TCCA", "TCCG", "TCCT", "TCGA", "TCGG", "TCGT",
  "TGAA", "TGAC", "TGAT", "TGCA", "TGCC", "TGCT", "TGGA", "TGGC", "TGGT",
  "TTAA", "TTAC", "TTAG", "TTCA", "TTCC", "TTCG", "TTGA", "TTGC", "TTGG"
)

SATS_DBS_INCLUDED_DINUCLEOTIDES <- c(
  "AC", "AT", "CC", "CG", "CT", "GC", "TA", "TC", "TG", "TT"
)

# Mutated dinucleotides that define the COSMIC DBS78 orientation for palindromic
# reference dinucleotides. These are used only to canonicalize DBS mutation
# records; they do not filter biological events.
SATS_DBS_PALINDROME_MUTATIONS <- list(
  AT = c("CA", "CC", "CG", "GA", "GC", "TA"),
  TA = c("AT", "CG", "CT", "GC", "GG", "GT"),
  CG = c("AT", "GC", "GT", "TA", "TC", "TT"),
  GC = c("AA", "AG", "AT", "CA", "CG", "TA")
)

# Internal EM wrapper defaults. User-facing optimization controls remain
# n.start, iter.max and eps in EstimateSigActivity().
SATS_EM_PRINT <- 0L
SATS_EM_DEBUG <- 0L
SATS_EM_INIT_LOWER <- 1e-6
SATS_DOUBLE_MISS <- -9999.0e200
SATS_LOGLIKE_MISS <- -9999.0
