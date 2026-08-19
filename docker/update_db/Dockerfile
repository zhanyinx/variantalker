FROM ubuntu:22.04

LABEL maintainer="yinxiu.zhan@ieo.it"
LABEL description="Docker image for VariantAlker database updates with GATK, vt, tabix, and bgzip"

# Avoid interactive prompts during build
ENV DEBIAN_FRONTEND=noninteractive

# Install basic dependencies
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    build-essential \
    autoconf \
    automake \
    libtool \
    pkg-config \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    python3 \
    python3-pip \
    openjdk-17-jdk \
    perl \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
RUN pip3 install --no-cache-dir \
    psycopg[binary] \
    requests \
    beautifulsoup4 \
    numpy \
    pandas \
    lxml

# Install htslib (provides tabix and bgzip)
WORKDIR /opt
RUN wget https://github.com/samtools/htslib/releases/download/1.19/htslib-1.19.tar.bz2 && \
    tar -xjf htslib-1.19.tar.bz2 && \
    cd htslib-1.19 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && \
    rm -rf htslib-1.19.tar.bz2

# Install vt
WORKDIR /opt
RUN git clone https://github.com/atks/vt.git && \
    cd vt && \
#    make && \
    cp vt /usr/local/bin/ && \
    cd .. && \
    rm -rf vt

# Install GATK 4.5.0.0
WORKDIR /opt
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip && \
    unzip gatk-4.5.0.0.zip && \
    rm gatk-4.5.0.0.zip && \
    ln -s /opt/gatk-4.5.0.0/gatk /usr/local/bin/gatk

# Set up working directory
WORKDIR /workspace

# Set environment variables
ENV PATH="/usr/local/bin:${PATH}"
ENV JAVA_HOME="/usr/lib/jvm/java-17-openjdk-amd64"

# Verify installations
RUN echo "=== Verifying installations ===" && \
    bgzip --version && \
    tabix --version && \
    vt --version && \
    gatk --version && \
    python3 --version && \
    perl --version

CMD ["/bin/bash"]
