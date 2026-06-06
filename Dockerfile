FROM rocker/r-ver:4.4.3

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

RUN R -e "install.packages(c('remotes', 'BiocManager'), repos='https://cloud.r-project.org')"

RUN R -e "BiocManager::install(c('GenomicRanges', 'IRanges', 'Biostrings', 'BSgenome.Hsapiens.UCSC.hg19', 'BSgenome.Hsapiens.UCSC.hg38'), ask=FALSE, update=FALSE)"

COPY source /opt/SATS/source

RUN R CMD INSTALL /opt/SATS/source

CMD ["R"]
