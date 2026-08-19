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

RUN pip install pandas click h5py numba numpy scipy
RUN pip install git+https://github.com/Roth-Lab/pyclone-vi.git
RUN pip install pyranges
RUN pip install pysam
