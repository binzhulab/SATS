FROM rocker/r-ver:4.4.3

LABEL org.opencontainers.image.title="SATS"
LABEL org.opencontainers.image.description="Signature Analyzer for Targeted Sequencing R environment"
LABEL org.opencontainers.image.source="https://github.com/binzhulab/SATS"
LABEL org.opencontainers.image.licenses="GPL-2"

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gfortran \
    libcurl4-openssl-dev \
    libgit2-dev \
    libssl-dev \
    libxml2-dev \
    zlib1g-dev \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN R -q -e "install.packages(c('remotes', 'BiocManager', 'testthat', 'glmnet', 'dplyr'), repos='https://cloud.r-project.org')"

RUN R -q -e "BiocManager::install(c('GenomicRanges', 'IRanges', 'Biostrings', 'BSgenome.Hsapiens.UCSC.hg19', 'BSgenome.Hsapiens.UCSC.hg38'), ask=FALSE, update=FALSE)"

COPY source /opt/SATS/source

RUN R CMD INSTALL /opt/SATS/source

RUN Rscript -e 'library(SATS); data(SimData, package="SATS"); data(RefTMB, package="SATS"); stopifnot(is.matrix(SimData$V), is.matrix(SimData$L)); stopifnot(identical(rownames(SimData$V), rownames(SimData$L))); stopifnot(nrow(RefTMB$TMB_SBS_v3.4) == nrow(SimData$V))'

CMD ["R"]
