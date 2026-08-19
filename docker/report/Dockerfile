FROM ubuntu:20.04

WORKDIR /reporter

ENV DEBIAN_FRONTEND=noninteractive

#ENV R_BASE_VERSION=4.1.3-1.2004.0

RUN apt update -y && apt -y upgrade && apt install -y --no-install-recommends \
    ca-certificates software-properties-common  gnupg  gnupg2 gnupg1

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' 

RUN apt update -y && apt -y upgrade && apt install -y --no-install-recommends \
   libfftw3-dev libncurses5-dev libncursesw5-dev \
   wget git gcc g++ libfreetype6-dev libcurl4-openssl-dev curl \
   r-base  r-base-dev \
   libbz2-dev \
   libxml2-dev \
   bioperl \
   libapache-dbi-perl


RUN Rscript -e "install.packages(c('rjson', 'DT'))"

RUN apt update -y && apt -y upgrade && apt install -y --no-install-recommends \
    pandoc-citeproc pandoc \
    libfontconfig1-dev \
    texlive-xetex lmodern texlive texlive-fonts-recommended texlive-plain-generic

RUN Rscript -e "install.packages('rmarkdown')"
RUN Rscript -e "install.packages('BiocManager')"
RUN Rscript -e "BiocManager::install('maftools')"
RUN Rscript -e "install.packages('qqman')"
RUN Rscript -e "install.packages('reshape2')"
