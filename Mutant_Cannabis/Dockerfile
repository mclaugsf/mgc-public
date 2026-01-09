FROM ubuntu:22.04

# Prevent tzdata from prompting
ENV DEBIAN_FRONTEND=noninteractive

# ----------------------------
# System dependencies
# ----------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    software-properties-common \
    dirmngr \
    gpg-agent \
    wget \
    curl \
    ca-certificates \
    git \
    build-essential \
    gfortran \
    cmake \
    make \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libpcre2-dev \
    libgit2-dev \
    apt-transport-https \
    bcftools \
    samtools \
    tabix \
    && rm -rf /var/lib/apt/lists/*

# ----------------------------
# Install R (CRAN)
# ----------------------------
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | apt-key add - && \
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/"

RUN apt-get update && apt-get install -y \
    r-base \
    r-base-dev \
    && rm -rf /var/lib/apt/lists/*

# ----------------------------
# Install CRAN R packages
# ----------------------------
RUN Rscript -e "install.packages( \
    c( \
      'BiocManager', \
      'argparse', \
      'ape', \
      'ggplot2', \
      'data.table', \
      'Rcpp', \
      'reticulate' \
    ), \
    repos='https://cran.r-project.org' \
)"

# ----------------------------
# Install Bioconductor packages
# ----------------------------
RUN Rscript -e "BiocManager::install( \
    c( \
      'gdsfmt', \
      'SNPRelate', \
      'GENESIS', \
      'GWASTools', \
      'SeqVarTools', \
      'ggtree' \
    ), \
    ask=FALSE, \
    update=FALSE \
)"

RUN apt-get update && apt-get install -y \
    libcairo2-dev \
    pkg-config \
    libxt-dev \
    libglib2.0-dev \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e "install.packages('ggiraph', repos='https://cran.r-project.org')"
RUN Rscript -e "BiocManager::install( \
    c( \
      'ggtree' \
    ), \
    ask=FALSE, \
    update=FALSE \
)"

RUN Rscript -e "BiocManager::install( \
    c( \
      'cowplot' \
    ), \
    ask=FALSE, \
    update=FALSE \
)"

# ----------------------------
# Default workdir
# ----------------------------
WORKDIR /work

CMD ["/bin/bash"]
