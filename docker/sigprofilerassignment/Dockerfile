FROM python:3.9.18

# Install base utilities
RUN apt-get update && \
    apt-get install -y build-essential  && \
    apt-get install -y wget && \
    apt-get install -y tabix && \
    apt-get install -y libreadline-dev && \
    apt-get install -y libcairo2-dev && \
    apt-get install -y git && \
    apt-get clean && \
    apt-get install -y procps g++ && \
    rm -rf /var/lib/apt/lists/* 

RUN apt-get update && apt install -y procps g++ && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN apt-get update && apt install -y ssh rsync && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN pip install SigProfilerAssignment
RUN mkdir -p /usr/local/lib/python3.9/site-packages/SigProfilerMatrixGenerator/references/vcf_files/GRCh38_bench/ /usr/local/lib/python3.9/site-packages/SigProfilerMatrixGenerator/references/vcf_files/GRCh37_bench/

COPY ./install_genomes.py /
