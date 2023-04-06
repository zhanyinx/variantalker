FROM r-base

# Install base utilities
RUN apt-get update && \
    apt-get install -y build-essential  && \
    apt-get install -y wget && \
    apt-get install -y tabix && \
    apt-get install -y libreadline-dev && \
    apt-get clean && \
    apt-get install -y procps g++ && \
    rm -rf /var/lib/apt/lists/* 


# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py39_22.11.1-1-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda

# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

RUN conda install -c bioconda --yes bedtools 
RUN conda install -c bioconda --yes biopython
RUN conda install -c bioconda/label/main --yes gatk4
RUN conda install -c conda-forge --yes pyvcf

RUN pip install cnvkit scipy matplotlib reportlab pyfaidx pysam 
RUN pip install numpy
RUN pip install pandas

RUN apt-get update && apt install -y procps g++ && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# nmd score
RUN pip install pyranges

# RNA seq annotation
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    curl \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libssh2-1-dev \
    zlib1g-dev \
    openssl \
    gdebi-core \
    libgsl* \
    libudunits2-dev \
    libgs-dev \
    imagemagick \
    ghostscript \
    qpdf \
    libfontconfig1-dev \
    libfreetype6-dev \
    r-cran-ragg

RUN mkdir /Rscripts
RUN mkdir data software
COPY resources/docker/packages_signatures.R /Rscripts/ 
COPY resources/docker/utils.R /Rscripts
RUN Rscript /Rscripts/packages_signatures.R 



# pyclone
RUN pip install h5py numba click
RUN apt-get update && apt install -y git
RUN pip install git+https://github.com/Roth-Lab/pyclone-vi.git

# nextflow
RUN conda install -c bioconda nextflow -y

# db updates
RUN pip install bs4 lxml
RUN pip install gsutil

# bcftools
RUN conda install -c bioconda bcftools -y

# copying data
COPY resources/docker/ClassifyCNV /software/ClassifyCNV
COPY resources/docker/escat_tiering.csv /data/escat_tiering.csv
COPY resources/docker/nmd_final_grange.tsv /data

